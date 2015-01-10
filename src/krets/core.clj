(ns krets.core
  (:require [clojure.core.matrix :as x]
            [clojure.core.matrix.protocols :as mp]
            [clojure.string :as s]
            [clojure.walk :as w]
            [clojure.java.shell :as sh]
            [clojure.pprint :as pp])
  (:import [java.awt Color]
           [javax.swing JFrame]
           [org.jfree.chart ChartFactory ChartPanel]
           [org.jfree.chart.plot XYPlot]
           [org.jfree.data.xy XYSeries XYSeriesCollection]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(x/set-current-implementation :vectorz)

(def low-key (comp keyword s/lower-case))

;; Netlist parser

(def spice-metric-prefixes
  {:f 1e-15
   :p 1e-12
   :n 1e-9
   :u 1e-6
   :m 1e-3
   :k 1e3
   :meg 1e6
   :g 1e9
   :t 1e12})

(def spice-number-pattern
  (re-pattern (str "(?i)" "^([+-]?\\d+(?:\\.\\d+)?(?:e[+-]?\\d+)?)"
                   "(" (s/join "|" (sort-by count > (map name (keys spice-metric-prefixes)))) ")?")))

(defn spice-number [s]
  (when-let [[_ x p] (and (string? s)
                          (re-find spice-number-pattern s))]
    (* (double (read-string x))
       (double (spice-metric-prefixes (some-> p low-key) 1)))))

(defn re? [re]
  #(re-find re %))

(def ground? zero?)

(defn number-of-voltage-sources [circuit]
  (count (:v circuit)))

(defn number-of-nodes [circuit]
  (->> (dissoc circuit :.)
       vals
       (apply concat)
       (mapcat #(map % [1 2]))
       (remove ground?)
       set
       count))

(defn commands [circuit]
  (group-by (comp low-key first) (:. circuit)))

(defn sub-commands [m]
  (group-by (comp low-key second) m))

(defn models [circuit]
  (->> (for [[_ n t & kvs] (:.model (commands circuit))]
         {n (with-meta
              (->> (for [[k v] (partition 2 kvs)]
                     [(low-key k) v])
                   (into {}))
              {:element-type (low-key t) :name n})})
       (apply merge)))

(defn time-step [circuit]
  (if-let [[[_ dt]] (:.tran (commands circuit))]
    dt
    Double/NaN))

(def non-linear-elements [:d])

(defn non-linear? [circuit]
  (boolean (some circuit non-linear-elements)))

(defn element-type [[[c]]]
  (low-key c))

(defn options [circuit]
  (dissoc (meta circuit) :netlist :title :compiled-conductance-stamp :compiled-source-stamp))

(defn parse-netlist [netlist]
  (let [[title & lines] (-> netlist
                            (s/replace #"\n\+" "")
                            s/split-lines)
        circuit (->> lines
                     (remove (some-fn (re? #"^\*") (re? #"(?i)^.end$")))
                     (map #(s/split % #"[\s=,()]+"))
                     (w/postwalk (some-fn spice-number identity))
                     (group-by element-type))]
    (with-meta circuit
      (->> '[number-of-nodes number-of-voltage-sources time-step models non-linear?]
           (map (juxt keyword #((ns-resolve 'krets.core %) circuit)))
           (apply merge {:netlist netlist :title title})))))

;; MNA

(defn a-matrix [circuit]
  (let [n (-> circuit meta :number-of-nodes double)
        m (-> circuit meta :number-of-voltage-sources double)]
    (x/zero-matrix (+ n m) (+ n m))))

(defn x-or-z-vector [circuit]
  (let [n (-> circuit meta :number-of-nodes double)
        m (-> circuit meta :number-of-voltage-sources double)]
    (x/zero-vector (+ n m))))

(defmulti conductance-element-fn (fn [opts e] (element-type e)))

(defmethod conductance-element-fn :r [_ [_ _ _ ^double r]]
  (fn [x]
    (/ r)))

(defmethod conductance-element-fn :c [{:keys [^double time-step]} [_ _ _ ^double c]]
  (fn [x]
    (/ c time-step)))

(defmethod conductance-element-fn :d [{:keys [models]} [_ ^double n1 _ model]]
  (let [vt 0.025875
        is (double (get-in models [model :is]))
        is-by-vt (/ is vt)]
    (fn [x]
      (let [vd (double (x/mget x (dec n1)))]
        (* is-by-vt (Math/exp (/ vd vt)))))))

;; This fn doesn't stamp the voltage sources in their rows outside the conductance sub matrix.
(defn compiled-conductance-stamp [circuit linearity]
  (let [a (gensym 'a)
        x (gensym 'x)
        g (gensym 'g)
        opts (gensym 'opts)
        ts (case linearity
            :linear [:r :c]
            :transient []
            :non-linear non-linear-elements)
        es (vec (apply concat (map circuit ts)))]
    `(let [~opts ~(options circuit)
           ~(vec (map (comp symbol first) es)) (map (partial conductance-element-fn ~opts) ~es)]
       (fn [~a ~x]
         ~@(for [t ts
                 [id n1 n2 :as e] (t circuit)]
             `(let [~g (double (~(symbol id) ~x))]
                ~@(for [^long row [n1 n2]
                        ^long col [n1 n2]
                        :when (not (or (ground? row) (ground? col)))
                        :let [row (dec row) col (dec col)]]
                    `(x/mset! ~a ~row ~col (+ (double (x/mget ~a ~row ~col))
                                              ~(if (= row col) `~g `(- ~g)))))))
         ~a))))

(defn conductance-stamp [circuit x linearity]
  (let [a (a-matrix circuit)]
    ((-> circuit meta (get-in [:compiled-conductance-stamp linearity])) a x)))

(defmulti source-element-fn (fn [opts e] (element-type e)))

(defmethod source-element-fn :c [{:keys [^double time-step]} [_ _ _ ^double c]]
  (let [g (/ c time-step)]
    (fn [row x]
      (* g (double (x/mget x row))))))

(defmethod source-element-fn :d [{:keys [models]} [_ ^double n1 _ model]]
  (let [vt 0.025875
        is (double (get-in models [model :is]))
        is-by-vt (/ is vt)]
    (fn [_ x]
      (let [vd (double (x/mget x (dec n1)))
            exp-vd-by-vt (Math/exp (/ vd vt))
            geq (* is-by-vt exp-vd-by-vt)
            id (* is (- exp-vd-by-vt 1))]
        (- id (* geq vd))))))

(defmethod source-element-fn :i [_ [_ _ _ _ ^double i]]
  (constantly i))

(defmethod source-element-fn :v [_ [_ _ _ _ ^double v]]
  (constantly v))

(defn compiled-source-stamp [circuit linearity]
  (let [n (-> circuit meta :number-of-nodes long)
        z (gensym 'z)
        x (gensym 'x)
        opts (gensym 'opts)
        ts (case linearity
             :linear [:v :i]
             :transient [:c]
             :non-linear non-linear-elements)
        es (vec (apply concat (map circuit ts)))]
    `(let [~opts ~(options circuit)
           ~(vec (map (comp symbol first) es)) (map (partial source-element-fn ~opts) ~es)]
       (fn [~z ~x]
         ~@(for [t ts
                 [^long idx [id n1 n2 :as e]] (map-indexed vector (t circuit))
                 [^double row sign] [[n1 `-] [n2 `+]]
                 :when (not (ground? row))
                 :let [row (dec row)
                       real-row (case t
                                  (:i, :c, :d) row
                                  :v (+ idx n))]]
             `(x/mset! ~z ~real-row (+ (double (x/mget ~z ~real-row))
                                       (double (~sign (~(symbol id) ~row ~x))))))
         ~z))))

(defn source-stamp [circuit x linearity]
  (let [z (x-or-z-vector circuit)]
    ((-> circuit meta (get-in [:compiled-source-stamp linearity])) z x)))

(defn compile-circuit [circuit]
  (->> (for [s '[compiled-source-stamp compiled-conductance-stamp]
             :let [compiler (ns-resolve 'krets.core s)]]
         {(keyword s)
          (apply merge (for [t [:linear :transient :non-linear]]
                         {t (eval (compiler circuit t))}))})
       (apply merge)
       (vary-meta circuit merge)))

(defn dc-operating-point
  ([circuit]
   (dc-operating-point circuit (x-or-z-vector circuit)))
  ([circuit x]
   (let [a (conductance-stamp circuit x :linear)
         z (source-stamp circuit x :linear)]
     {:a a :z z :x (mp/solve a z)})))

(def ^:dynamic *newton-tolerance* 1e-7)
(def ^:dynamic *newton-iterations* 500)

(defn non-linear-step [circuit a z x _]
  (let [z (x/add z (source-stamp circuit x :transient))
        newton-step (fn [x]
                      (let [a (x/add a (conductance-stamp circuit x :non-linear))
                            z (x/add z (source-stamp circuit x :non-linear))]
                        (mp/solve a z)))
        within-tolerance? (fn [^double x ^double y]
                            (< (/ (Math/abs (- y x)) (Math/abs y)) (double *newton-tolerance*)))
        converged? (fn [[x y]]
                     (every? true? (map within-tolerance? x y)))
        [_ x] (->> x
                   (iterate newton-step)
                   (take *newton-iterations*)
                   (partition 2 1)
                   (drop-while (complement converged?))
                   first)]
    (assert x "Didn't converge.")
    x))

(defn linear-step [circuit a z x _]
  (let [z (x/add z (source-stamp circuit x :transient))]
    (mp/solve a z)))

(defn transient-analysis
  ([circuit time-step simulation-time]
   (transient-analysis circuit time-step simulation-time (dc-operating-point circuit)))
  ([circuit time-step simulation-time {:keys [a z x]}]
    (let [ts (range 0 simulation-time time-step)
          step (if (-> circuit meta :non-linear?)
                 (partial non-linear-step circuit a z)
                 (partial linear-step circuit a z))]
      (->> ts
           (reductions step x)
           (map vector ts)))))

;; Frontend

(defn node-label [[k v]]
  (str k (long v)))

(defn print-result [circuit series]
  (let [n (-> circuit meta :number-of-nodes double)
        time-label "t"
        dt (time-step circuit)
        time-format (str "%." (.scale (bigdec dt)) "f")
        node-format "%.10f"]
    (doseq [[_ _ & nodes] (:tran (sub-commands (:.print (commands circuit))))
            :let [nodes (partition 2 nodes)]
            :when (seq nodes)]
      (pp/print-table
       (concat [time-label] (map node-label nodes))
       (for [[t x] series]
         (apply merge {time-label (format time-format (double t))}
                (for [[k v :as node] nodes
                      :let [v (double v)]]
                  {(node-label node)
                   (format node-format
                           (x/mget x (case k
                                       "V" (dec v)
                                       "I" (+ n v))))})))))))

(defn xy-series
  ([title x y]
   (reduce
    (fn [^XYSeries xy [x y]]
      (doto xy
        (.add (double x) (double y))))
    (XYSeries. title) (map vector x y)))
  ([title y]
   (xy-series title (range) y)))

(defn xy-series-coll [& xys]
  (reduce
   (fn [^XYSeriesCollection xy-s s]
     (doto xy-s
       (.addSeries s)))
   (XYSeriesCollection.) xys))

(defn plot-theme [^XYPlot plot bg grid]
  (doto plot
    (.setBackgroundPaint bg)
    (.setDomainGridlinePaint grid)
    (.setRangeGridlinePaint grid)))

(defn plot [x-title y-title & series]
  (doto (JFrame.)
    (.setContentPane
     (doto (ChartPanel.
            (doto (ChartFactory/createXYLineChart
                   nil x-title y-title
                   (apply xy-series-coll series))
              (-> .getXYPlot (plot-theme Color/WHITE Color/LIGHT_GRAY))))
       (.setInitialDelay 100)))
    .pack
    (.setVisible true)))

(defn plot-result [circuit series]
  (let [n (-> circuit meta :number-of-nodes double)
        ts (map first series)
        xs (map second series)]
    (doseq [[_ _ & nodes] (:tran (sub-commands (:.plot (commands circuit))))
            :let [nodes (partition 2 nodes)]
            :when (seq nodes)]
      (apply plot
             "t" "V"
             (for [[k v :as node] nodes
                   :let [v (double v)]]
               (xy-series (node-label node)
                          ts
                          (map #(x/mget % (case k
                                            "V" (dec v)
                                            "I" (+ n v))) xs) ))))))

(defn batch [circuit]
  (let [circuit (compile-circuit circuit)]
    (println (-> circuit meta :title))
    (pp/print-table [(options circuit)])
    (println "DC Operating Point Analysis")
    (let [dc-result (time (dc-operating-point circuit))]
      (pp/print-table [dc-result])
      (doseq [[_ dt simulation-time] (:.tran (commands circuit))
              :let [series (do (println "Transient Analysis" dt simulation-time)
                               (time (doall (transient-analysis circuit dt simulation-time dc-result))))]]
        (print-result circuit series)
        (plot-result circuit series)))))

(def ^:dynamic *spice-command* "ngspice")

(defn spice [circuit]
  (let [{:keys [out err ^long exit]} (sh/sh *spice-command* "-b" :in (-> circuit meta :netlist))]
    (when out
      (println out))
    (when (not (zero? exit))
      (println err))))

(defn process-file [f]
  (-> f slurp parse-netlist batch))

(defn -main [& [f]]
  (if f
    (process-file f)
    (println "Need to specify a netlist file")))
