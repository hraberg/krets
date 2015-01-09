(ns krets.core
  (:require [clojure.core.matrix :as x]
            [clojure.core.matrix.protocols :as mp]
            [clojure.string :as s]
            [clojure.walk :as w]
            [clojure.pprint :as pp])
  (:import [java.awt Color]
           [javax.swing JFrame]
           [org.jfree.chart ChartFactory ChartPanel]
           [org.jfree.chart.plot XYPlot]
           [org.jfree.data.xy XYSeries XYSeriesCollection]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(x/set-current-implementation :vectorz)

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
    (* ^double (read-string x)
       ^double (spice-metric-prefixes (some-> p s/lower-case keyword) 1))))

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
  (group-by (comp keyword s/lower-case first) (:. circuit)))

(defn sub-commands [m]
  (group-by (comp keyword s/lower-case second) m))

(defn models [circuit]
  (->> (for [[_ n t & kvs] (:.model (commands circuit))]
         {n (with-meta
              (apply hash-map kvs)
              {:element-type (keyword (s/lower-case t)) :name n})})
       (apply merge)))

(defn time-step [circuit]
  (if-let [[[_ dt]] (:.tran (commands circuit))]
    dt
    Double/NaN))

(def non-linear-elements #{:d})

(defn non-linear? [circuit]
  (boolean (some non-linear-elements (keys circuit))))

(defn parse-netlist [netlist]
  (let [[title & lines] (-> netlist
                            (s/replace #"\n\+" "")
                            s/split-lines)
        circuit (->> lines
                     (remove (some-fn (re? #"^\*") (re? #"(?i)^.end$")))
                     (map #(s/split % #"[\s=,()]+"))
                     (w/postwalk (some-fn spice-number identity))
                     (group-by (fn [[[c]]]
                                 (keyword (s/lower-case c)))))]
    (with-meta circuit
      (->> '[number-of-nodes number-of-voltage-sources time-step models non-linear?]
           (map (juxt keyword #((ns-resolve 'krets.core %) circuit)))
           (apply merge {:netlist netlist :title title})))))

;; MNA

(defn a-matrix [circuit]
  (let [n (-> circuit meta :number-of-nodes double)
        m (-> circuit meta :number-of-voltage-sources double)]
    (x/zero-matrix (+ n m) (+ n m))))

(defn xz-vector [circuit]
  (let [n (-> circuit meta :number-of-nodes double)
        m (-> circuit meta :number-of-voltage-sources double)]
    (x/zero-vector (+ n m))))

(defn conductance-stamp [circuit x linearity]
  (let [a (a-matrix circuit)
        dt (-> circuit meta :time-step)]
    (doseq [k (case linearity
                :linear [:r :c]
                :non-linear non-linear-elements)
            [_ ^double n1 n2 r-or-c-or-model] (k circuit)
            :let [g (double (case k
                              :r (/ ^double r-or-c-or-model)
                              :c (/ ^double r-or-c-or-model ^double dt)
                              :d (let [vt 0.025875
                                       is (-> circuit meta :models (get-in [r-or-c-or-model "IS"]) double)
                                       vd (double (x/mget x (dec n1)))]
                                   (* (/ is vt) (Math/exp (/ vd vt))))))]
            ^long xm [n1 n2]
            ^long ym [n1 n2]
            :when (not (or (ground? xm) (ground? ym)))
            :let [xm (dec xm) ym (dec ym)]]
      (x/mset! a xm ym (+ ^double (x/mget a xm ym)
                          (if (= xm ym) g (- g)))))
    a))

(defn source-stamp [circuit x linearity]
  (let [z (xz-vector circuit)
        dt (-> circuit meta :time-step double)
        n (-> circuit meta :number-of-nodes long)]
    (doseq [k (case linearity
                :linear [:v :i]
                :transient [:c]
                :non-linear non-linear-elements)
            [^long idx [_ ^double n1 ^double n2 c-or-model ^double val]] (map-indexed vector (k circuit))
            [^double xm sign] [[n1 -] [n2 +]]
            :when (not (ground? xm))
            :let [xm (dec xm)]]
      (#(x/mset! z (case k
                     (:i, :c, :d) xm
                     :v (+ idx n))
                 (sign (case k
                         :c (* (/ ^double c-or-model dt) (double (x/mget x xm)))
                         :d (let [vt 0.025875
                                  is (-> circuit meta :models (get-in [c-or-model "IS"]) double)
                                  vd (double (x/mget x (dec n1)))
                                  geq (* (/ is vt) (Math/exp (/ vd vt)))
                                  id (* is (- (Math/exp (/ vd vt)) 1))]
                              (- id (* geq vd)))
                         val)))))
    z))

(defn dc-operation-point
  ([circuit]
   (dc-operation-point circuit (xz-vector circuit)))
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
        within? (fn [^double x ^double y]
                  (< (/ (Math/abs (- y x)) (Math/abs (- y))) (double *newton-tolerance*)))
        converged? (fn [[x y]]
                     (every? true? (map within? x y)))
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
   (transient-analysis circuit simulation-time time-step (dc-operation-point circuit)))
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
  (println (-> circuit meta :title))
  (pp/print-table [(dissoc (meta circuit) :netlist :title)])
  (println "DC Operating Point Analysis")
  (let [cs (commands circuit)
        dc-result (time (dc-operation-point circuit))]
    (pp/print-table [dc-result])
    (doseq [[_ dt simulation-time] (:.tran cs)
            :let [series (do (println "Transient Analysis" dt simulation-time)
                             (time (doall (transient-analysis circuit dt simulation-time dc-result))))]]
      (print-result circuit series)
      (plot-result circuit series))))

(defn process-file [f]
  (-> f slurp parse-netlist batch))

(defn -main [& [f]]
  (if f
    (process-file f)
    (println "Need to specify a netlist file")))
