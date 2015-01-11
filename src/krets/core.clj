(ns krets.core
  (:require [clojure.string :as s]
            [clojure.walk :as w]
            [clojure.java.shell :as sh]
            [clojure.pprint :as pp])
  (:import [java.awt Color]
           [javax.swing JFrame]
           [org.jfree.chart ChartFactory ChartPanel]
           [org.jfree.chart.plot XYPlot]
           [org.jfree.data.xy XYSeries XYSeriesCollection]
           [org.ejml.data DenseMatrix64F]
           [org.ejml.ops CommonOps MatrixFeatures]
           [clojure.lang IFn$OD]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; Matrix ops

(defn mtag [m]
  (with-meta m {:tag `DenseMatrix64F}))

(definline zero-matrix [rows cols]
  `(DenseMatrix64F. ~rows ~cols))

(definline row-count [m]
  `(.numRows ~(mtag m)))

(definline solve [a b]
  `(doto (zero-matrix (row-count ~a) 1)
     (->> (CommonOps/solve ~(mtag a) ~b))))

(definline add! [a b]
  `(doto ~(mtag a) (-> (CommonOps/addEquals ~b))))

(definline equals [a b epsilon]
  `(MatrixFeatures/isEquals ~(mtag a) ~b ~epsilon))

(definline mget [m row col]
  `(.unsafe_get ~(mtag m) ~row ~col))

(definline madd! [m row col v]
  `(.add ~(mtag m) ~row ~col ~v))

;; Netlist parser

(def low-key (comp keyword s/lower-case))

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

(defn number-of-rows [circuit]
  (+ (long (number-of-voltage-sources circuit))
     (long (number-of-nodes circuit))))

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
  (dissoc (meta circuit) :netlist :title :compiled? :compiled-conductance-stamp :compiled-source-stamp))

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
      (->> '[number-of-nodes number-of-voltage-sources number-of-rows time-step models non-linear?]
           (map (juxt keyword #((ns-resolve 'krets.core %) circuit)))
           (apply merge {:netlist netlist :title title})))))

;; MNA

(defn voltage-diff-fn ^IFn$OD [n+ n-]
  (let [x (gensym 'x)
        term (fn [^long n]
               (if (ground? n)
                 0.0
                 `(mget ~x ~(dec n) 0)))]
    (eval `(fn ^double [~x]
             (- ~(term n+) ~(term n-))))))

(defmulti conductance-element-fn (fn [opts e] (element-type e)))

(defmethod conductance-element-fn :r [_ [_ _ _ ^double r]]
  (fn ^double [x]
    (/ 1.0 r)))

(defmethod conductance-element-fn :c [{:keys [^double time-step]} [_ _ _ ^double c]]
  (fn ^double [x]
    (/ c time-step)))

(defmethod conductance-element-fn :d [{:keys [models]} [_ anode cathode model]]
  (let [vt 0.025875
        is (-> model models :is double)
        is-by-vt (/ is vt)
        vd-fn (voltage-diff-fn anode cathode)]
    (fn ^double [x]
      (let [vd (.invokePrim vd-fn x)]
        (* is-by-vt (Math/exp (/ vd vt)))))))

;; This fn doesn't stamp the voltage sources in their rows outside the conductance sub matrix.
(defn compiled-conductance-stamp [circuit linearity]
  (let [[a x g opts] (map gensym '[a x g opts])
        {:keys [^long number-of-rows]} (meta circuit)
        ts (case linearity
             :linear [:r :c]
             :transient []
             :non-linear non-linear-elements)
        es (vec (mapcat circuit ts))]
    `(let [~opts ~(options circuit)
           ~(vec (map (comp symbol first) es)) (map (partial conductance-element-fn ~opts) ~es)]
       (fn [~x]
         (let [~a (zero-matrix ~number-of-rows ~number-of-rows)]
           ~@(for [t ts
                   [id n+ n- :as e] (t circuit)]
               `(let [~g (.invokePrim ~(with-meta (symbol id) {:tag `IFn$OD}) ~x)]
                  ~@(for [^long row [n+ n-]
                          ^long col [n+ n-]
                          :when (not (or (ground? row) (ground? col)))
                          :let [row (dec row) col (dec col)]]
                      `(madd! ~a ~row ~col ~(if (= row col) `~g `(- ~g))))))
           ~a)))))

(defmulti source-element-fn (fn [opts e] (element-type e)))

(defmethod source-element-fn :c [{:keys [^double time-step]} [_ n+ n- ^double c]]
  (let [g (/ c time-step)
        ;; TODO: why do we need the reverse order here?
        ;;       this used to be in the reverse in the actual netlist, but think that was wrong.
        vd-fn (voltage-diff-fn n- n+)]
    (fn ^double [x]
      (* g (.invokePrim vd-fn x)))))

(defmethod source-element-fn :d [{:keys [models]} [_ anode cathode model]]
  (let [vt 0.025875
        is (-> model models :is double)
        is-by-vt (/ is vt)
        vd-fn (voltage-diff-fn anode cathode)]
    (fn ^double [x]
      (let [vd (.invokePrim vd-fn x)
            exp-vd-by-vt (Math/exp (/ vd vt))
            geq (* is-by-vt exp-vd-by-vt)
            id (* is (- exp-vd-by-vt 1.0))]
        ;; TODO: should the other row be (+ (- id) (* geq vd)) - that is, they are in different directions?
        ;;       ie. more generic, simply reverse the terms, (- (* geq vd) id)
        ;;       so, we cannot call this fn just once, alternatively, we need to stamp in here.
        ;;       could this be related to the above problem as well?
        (- id (* geq vd))))))

(defmethod source-element-fn :i [_ [_ _ _ _ ^double i]]
  (fn ^double [x] i))

(defmethod source-element-fn :v [_ [_ _ _ _ ^double v]]
  (fn ^double [x] v))

(defn compiled-source-stamp [circuit linearity]
  (let [{:keys [^long number-of-rows ^long number-of-nodes]} (meta circuit)
        [z x i-or-v opts] (map gensym '[z x i-or-v opts])
        ts (case linearity
             :linear [:v :i]
             :transient [:c]
             :non-linear non-linear-elements)
        es (vec (mapcat circuit ts))]
    `(let [~opts ~(options circuit)
           ~(vec (map (comp symbol first) es)) (map (partial source-element-fn ~opts) ~es)]
       (fn [~x]
         (let [~z (zero-matrix ~number-of-rows 1)]
           ~@(for [t ts
                   [^long idx [id n+ n- :as e]] (map-indexed vector (t circuit))]
               `(let [~i-or-v (.invokePrim ~(with-meta (symbol id) {:tag `IFn$OD}) ~x)]
                  ~@(for [[^long row sign] [[n+ `-] [n- `+]]
                          :when (not (ground? row))
                          :let [row (case t
                                      (:i, :c, :d) (dec row)
                                      :v (+ idx number-of-nodes))]]
                      `(madd! ~z ~row 0 (~sign ~i-or-v)))))
           ~z)))))

(defn compile-circuit [circuit]
  (if (-> circuit meta :compiled?)
    circuit
    (->> (for [s '[compiled-source-stamp compiled-conductance-stamp]
               :let [compiler (ns-resolve 'krets.core s)]]
           {:compiled? true
            (keyword s)
            (apply merge (for [t [:linear :transient :non-linear]]
                           {t (eval (compiler circuit t))}))})
         (apply merge)
         (vary-meta circuit merge))))

(defn dc-operating-point
  ([circuit]
   (dc-operating-point circuit (zero-matrix (-> circuit meta :number-of-rows) 1)))
  ([circuit x]
   (let [{:keys [compiled-source-stamp compiled-conductance-stamp]} (meta circuit)
         a ((compiled-conductance-stamp :linear) x)
         z ((compiled-source-stamp :linear) x)]
     {:a a :z z :x (solve a z)})))

(def ^:dynamic *newton-tolerance* 1e-8)
(def ^:dynamic *newton-iterations* 500)

(defn non-linear-step-fn [circuit a z]
  (let [{:keys [compiled-source-stamp compiled-conductance-stamp]} (meta circuit)
        transient-source-stamp (compiled-source-stamp :transient)
        non-linear-conductance-stamp (compiled-conductance-stamp :non-linear)
        non-linear-source-stamp (compiled-source-stamp :non-linear)
        newton-tolerance (double *newton-tolerance*)
        newton-iterations (long *newton-iterations*)]
    (fn [x]
      (let [z (add! (transient-source-stamp x) z)]
        (loop [xn-1 x iters newton-iterations]
          (if (zero? iters)
            (assert xn-1 "Didn't converge.")
            (let [xn (solve (add! (non-linear-conductance-stamp xn-1) a)
                            (add! (non-linear-source-stamp xn-1) z))]
              (if (equals xn xn-1 newton-tolerance)
                xn
                (recur xn (dec iters))))))))))

(defn linear-step-fn [circuit a z]
  (let [{:keys [compiled-source-stamp]} (meta circuit)
        transient-source-stamp (compiled-source-stamp :transient)]
    (fn [x]
      (solve a (add! (transient-source-stamp x) z)))))

(defn transient-step-fn [circuit a z]
  (let [step-fn (if (-> circuit meta :non-linear?)
                  non-linear-step-fn
                  linear-step-fn)]
    (step-fn circuit a z)))

(defn transient-analysis
  ([circuit ^double time-step ^double simulation-time]
   (transient-analysis circuit time-step simulation-time (dc-operating-point circuit)))
  ([circuit ^double time-step ^double simulation-time {:keys [a z x]}]
   (let [step (transient-step-fn circuit a z)]
     (loop [t 0.0 x x acc (transient [])]
       (if (> t simulation-time)
         (persistent! acc)
         (let [x (step x)]
           (recur (+ t time-step) x (conj! acc [t x]))))))))

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
                (for [[k ^double v :as node] nodes]
                  {(node-label node)
                   (format node-format
                           (mget x (case k
                                     "V" (dec v)
                                     "I" (+ n v)) 0))})))))))

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
             (for [[k ^double v :as node] nodes]
               (xy-series (node-label node)
                          ts
                          (map #(mget % (case k
                                          "V" (dec v)
                                          "I" (+ n v)) 0) xs) ))))))

(defn batch [circuit]
  (let [circuit (compile-circuit circuit)]
    (println (-> circuit meta :title))
    (pp/print-table [(options circuit)])
    (println "DC Operating Point Analysis")
    (let [{:keys [a z x] :as dc-result} (time (dc-operating-point circuit))]
      (println "A")
      (println a)
      (println "z")
      (println z)
      (println "x")
      (println x)
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
