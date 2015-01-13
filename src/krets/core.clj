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
           [org.ejml.ops CommonOps MatrixFeatures]))

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
                          (re-find spice-number-pattern
                                   (cond->> s (re-find #"^\.\d" s) (str "0"))))]
    (* (double (read-string x))
       (double (spice-metric-prefixes (some-> p low-key) 1)))))

(defn re? [re]
  #(re-find re (str %)))

(def ground? zero?)

(defn number-of-voltage-sources ^long [netlist]
  (count (mapcat netlist [:v :e :u])))

(defn number-of-nodes ^long [netlist]
  (->> (dissoc netlist :.)
       vals
       (apply concat)
       (mapcat #(map % [1 2]))
       (remove ground?)
       set
       count))

(defn number-of-rows [netlist]
  (+ (number-of-voltage-sources netlist)
     (number-of-nodes netlist)))

(defn voltage-source->index [netlist]
  (let [number-of-nodes (number-of-nodes netlist)]
    (into {} (for [[^long idx [id]] (map-indexed vector (mapcat netlist [:v :e :u]))]
               [id (+ number-of-nodes idx)]))))

(defn commands [netlist]
  (group-by (comp low-key first) (:. netlist)))

(defn sub-commands [m]
  (group-by (comp low-key second) m))

(defn models [netlist]
  (->> (for [[_ n t & kvs] (:.model (commands netlist))]
         {n (with-meta
              (->> (for [[k v] (partition 2 kvs)]
                     [(low-key k) v])
                   (into {}))
              {:element-type (low-key t) :name n})})
       (apply merge)))

(defn time-step [netlist]
  (if-let [[[_ time-step]] (:.tran (commands netlist))]
    time-step
    Double/NaN))

(defn non-linear? [netlist]
  (boolean (some netlist [:d])))

(defn elements [netlist]
  (mapcat val (dissoc netlist :.)))

(defn element-type [[[c]]]
  (low-key c))

(defn circuit-info [circuit]
  (dissoc circuit :netlist :title :compiled? :conductance-stamp :source-stamp :voltage-source->index))

(defn parse-netlist [netlist-source]
  (let [[title & lines] (-> netlist-source
                            (s/replace #"\n\+" "")
                            s/split-lines)
        parsed-netlist (->> lines
                            (remove (some-fn (re? #"^\*") (re? #"(?i)^.end$")))
                            (map #(s/split % #"[\s=,()]+"))
                            (w/postwalk (some-fn spice-number identity))
                            (group-by element-type))]
    (apply merge (with-meta {:netlist parsed-netlist :title title} {:netlist-source netlist-source})
           (for [k '[number-of-nodes number-of-voltage-sources number-of-rows time-step models non-linear? voltage-source->index]]
             {(keyword k) ((ns-resolve 'krets.core k) parsed-netlist)}))))

;; MNA Compiler

(defprotocol MNAStamp
  (linear-stamp [_ x])
  (transient-stamp [_ x])
  (non-linear-stamp [_ x]))

(defn mna-stamp [{:keys [netlist] :as circuit} rows cols stamp-fn]
  `(reify MNAStamp
     ~@(for [f (-> MNAStamp :method-map keys)
             :let [[s x] (map gensym '[stamp x])
                   k (low-key (s/replace (name f) "-stamp" ""))]]
         `(~(symbol (name f)) [_# ~x]
           (let [~s (zero-matrix ~rows ~cols)]
             ~@(for [e (elements netlist)]
                 (stamp-fn circuit k e x s))
             ~s)))))

(defn voltage-diff [x n+ n-]
  (let [term (fn [^long n]
               `(mget ~x ~(dec n) 0))]
    (cond
     (ground? n+) `(- ~(term n-))
     (ground? n-) (term n+)
     :else `(- ~(term n+) ~(term n-)))))

(defn conductance-stamp [a n+ n- g]
  (let [gs (gensym 'g)]
    `(let [~gs ~g]
       ~@(for [^long row [n+ n-]
               ^long col [n+ n-]
               :when (not (or (ground? row) (ground? col)))
               :let [row (dec row) col (dec col)]]
           `(madd! ~a ~row ~col ~(if (= row col) `~gs `(- ~gs)))))))

(defmulti conductance-element (fn [circuit stamp e x a] [(element-type e) stamp]))

(defmethod conductance-element :default [_ _ _ _ _])

(defmethod conductance-element [:r :linear] [_ _ [_ n+ n- ^double r] _ a]
  (conductance-stamp a n+ n- (/ 1.0 r)))

(defmethod conductance-element [:c :linear] [{:keys [^double time-step]} _ [_ n+ n- ^double c] _ a]
  (conductance-stamp a n+ n- (/ c time-step)))

(defmethod conductance-element [:d :non-linear] [{:keys [models]} _ [_ anode cathode model] x a]
  (let [vt 0.025875
        is (-> model models :is double)
        is-by-vt (/ is vt)]
    (conductance-stamp a anode cathode
     `(let [vd# ~(voltage-diff x anode cathode)]
        (* ~is-by-vt (Math/exp (/ vd# ~vt)))))))

(defn conductance-voltage-stamp [a n+ n- idx]
  `(do ~@(for [[^long n v] [[n+ 1] [n- -1]]
               :when (not (ground? n))
               :let [n (dec n)]]
           `(do (madd! ~a ~n ~idx ~v)
                (madd! ~a ~idx ~n ~v)))))

(defmethod conductance-element [:v :linear] [{:keys [voltage-source->index]} _ [id n+ n-] _ a]
  (conductance-voltage-stamp a n+ n- (voltage-source->index id)))

(defmethod conductance-element [:e :linear] [{:keys [voltage-source->index]} _
                                             [id ^long out+ ^long out- ^long in+ ^long in- ^double gain] _ a]
  (let [idx (voltage-source->index id)]
    `(do ~(conductance-voltage-stamp a out- out+ idx)
         ~@(for [[^long n g] [[in+ gain] [in- (- gain)]]
                 :when (not (ground? n))]
             `(madd! ~a ~idx ~(dec n) ~g)))))

;; All About Circuits model ideal op amps as a vcvs.
(defn u->e [[id ^long in+ ^long in- ^long out+]]
  [(str "e" id) out+ 0 in+ in- 999e3])

;; ideal op amp, this is only the linear version. QUCS technical.pdf p 117. QUCS source has a transient part.
(defmethod conductance-element [:u :linear] [{:keys [voltage-source->index]} _
                                             [id ^long in+ ^long in- ^long out+] _ a]
  (let [idx (inc (long (voltage-source->index id)))]
    `(do ~@(for [[^long n+ ^long n- v] [[idx in+ 1] [idx in- -1] [out+ idx 1]]
                 :when (not (or (ground? n+) (ground? n-)))]
             `(madd! ~a ~(dec n+) ~(dec n-) ~v)))))

(defn source-current-stamp [z n+ n- in out]
  `(do ~@(for [[^long row i] [[n+ in] [n- out]]
               :when (not (ground? row))]
           `(madd! ~z ~(dec row) 0 ~i))))

(defmulti source-element (fn [circuit stamp e x z] [(element-type e) stamp]))

(defmethod source-element :default [_ _ _ _ _])

(defmethod source-element [:c :transient] [{:keys [^double time-step]} _ [_ n+ n- ^double c] x z]
  (let [ieq (gensym 'ieq)
        geq (/ c time-step)]
    `(let [~ieq (* ~geq ~(voltage-diff x n+ n-))]
       ~(source-current-stamp z n+ n- ieq `(- ~ieq)))))

(defmethod source-element [:d :non-linear] [{:keys [models]} _ [_ anode cathode model] x z]
  (let [ieq (gensym 'ieq)
        vt 0.025875
        is (-> model models :is double)
        is-by-vt (/ is vt)]
    `(let [vd# ~(voltage-diff x anode cathode)
           exp-vd-by-vt# (Math/exp (/ vd# ~vt))
           geq# (* ~is-by-vt exp-vd-by-vt#)
           id# (* ~is (- exp-vd-by-vt# 1.0))
           ~ieq (+ (- id#) (* geq# vd#))]
       ~(source-current-stamp z anode cathode ieq `(- ~ieq)))))

(defn source-value ^double [[id n+ n- & opts]]
  (or (first (filter number? opts)) 0.0))

(defmethod source-element [:i :linear] [_ _ [_ n+ n- :as e] _ z]
  (let [i (source-value e)]
    (source-current-stamp z n+ n- (- i) i)))

(def ^:dynamic *voltage-sources* {})

(defmethod source-element [:v :linear] [{:keys [voltage-source->index netlist]} _ [id :as e] _ z]
  (let [dc-sweep? ((low-key id) (-> netlist commands :.dc sub-commands))
        v (source-value e)
        idx (voltage-source->index id)]
    `(madd! ~z ~idx 0 ~(if dc-sweep?
                         `(*voltage-sources* ~id ~v)
                         v))))

(defn compile-circuit [{:keys [number-of-rows] :as circuit}]
  (if (:compiled? circuit)
    circuit
    (assoc circuit
      :compiled? true
      :conductance-stamp (eval (mna-stamp circuit number-of-rows number-of-rows conductance-element))
      :source-stamp (eval (mna-stamp circuit number-of-rows 1 source-element)))))

;; MNA Analysis

(defn dc-operating-point
  ([{:keys [^long number-of-rows] :as circuit}]
   (dc-operating-point circuit (zero-matrix number-of-rows 1)))
  ([{:keys [conductance-stamp source-stamp]} x]
   (let [a (linear-stamp conductance-stamp x)
         z (linear-stamp source-stamp x)]
     {:a a :z z :x (solve a z)})))

(defn dc-analysis [circuit source start stop step]
  (for [v (concat (range start stop step) [stop])]
    (binding [*voltage-sources* (assoc *voltage-sources* source v)]
      [v (-> circuit dc-operating-point :x)])))

(def ^:dynamic *newton-tolerance* 1e-8)
(def ^:dynamic *newton-iterations* 500)

(defn non-linear-step-fn [{:keys [conductance-stamp source-stamp]} a z]
  (let [newton-tolerance (double *newton-tolerance*)
        newton-iterations (long *newton-iterations*)]
    (fn [x]
      (let [z (add! (transient-stamp source-stamp x) z)]
        (loop [xn-1 x iters newton-iterations]
          (if (zero? iters)
            (throw (ex-info "Didn't converge." {:x xn-1}))
            (let [xn (solve (add! (non-linear-stamp conductance-stamp xn-1) a)
                            (add! (non-linear-stamp source-stamp xn-1) z))]
              (if (equals xn xn-1 newton-tolerance)
                xn
                (recur xn (dec iters))))))))))

(defn linear-step-fn [{:keys [source-stamp]} a z]
  (fn [x]
    (solve a (add! (transient-stamp source-stamp x) z))))

(defn transient-step-fn [{:keys [non-linear?] :as circuit} a z]
  (let [step-fn (if non-linear?
                  non-linear-step-fn
                  linear-step-fn)]
    (step-fn circuit a z)))

(defn transient-analysis
  ([circuit ^double time-step ^double simulation-time]
   (transient-analysis circuit time-step simulation-time (dc-operating-point circuit)))
  ([circuit ^double time-step ^double simulation-time {:keys [a z x]}]
   (let [step (transient-step-fn circuit a z)
         end (+ simulation-time time-step)]
     (loop [t 0.0 x x acc (transient [])]
       (if (> t end)
         (persistent! acc)
         (let [x (step x)]
           (recur (+ t time-step) x (conj! acc [t x]))))))))

;; Frontend

(defn report-node-label [[k n+ n-]]
  (format "%s(%s)" k (s/join ","  (cond-> [(long n+)] n- (concat [(long n-)])))))

(defn report-node-voltage [[k & ns] ^long number-of-nodes x]
  (->> (for [^long n ns]
         (if (ground? n)
           0.0
           (let [n (dec n)]
             (mget x (case (low-key k)
                       :v n
                       :i (+ number-of-nodes n)) 0))))
       (reduce -)))

(defn report-nodes [nodes]
  (->> nodes
       (partition-by (re? #"(?i)[vi]"))
       (partition 2)
       (map (partial apply concat))))

(defn print-result [{:keys [^long number-of-nodes time-step netlist]} series series-type head-label]
  (let [transient? (= :tran series-type)
        head-format (if transient?
                      (str "%." (.scale (bigdec time-step)) "f")
                      "%f")
        node-format "%.10f"]
    (doseq [[_ _ & nodes] (-> netlist commands :.print sub-commands series-type)
            :let [nodes (report-nodes nodes)]
            :when (seq nodes)]
      (pp/print-table
       (concat [head-label] (map report-node-label nodes))
       (for [[h x] series]
         (apply merge {head-label (format head-format h)}
                (for [node nodes]
                  {(report-node-label node)
                   (format node-format (report-node-voltage node number-of-nodes x))})))))))

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

(defn plot-result [{:keys [^long number-of-nodes netlist]} series series-type head-label]
  (let [xs (map first series)
        ys (map second series)]
    (doseq [[_ _ & nodes] (-> netlist commands :.plot sub-commands series-type)
            :let [nodes (report-nodes nodes)]
            :when (seq nodes)]
      (apply plot head-label "V"
             (for [node nodes]
               (xy-series (report-node-label node)
                          xs
                          (map #(report-node-voltage node number-of-nodes %) ys)))))))

(defn batch [{:keys [title netlist] :as circuit}]
  (let [circuit (compile-circuit circuit)]
    (println title)
    (pp/print-table [(circuit-info circuit)])
    (let [{:keys [a z x] :as dc-result} (do (println "DC Operating Point Analysis")
                                            (time (dc-operating-point circuit)))]
      (println "A")
      (println a)
      (println "z")
      (println z)
      (println "x")
      (println x)
      (doseq [[_ source start stop step] (:.dc (commands netlist))
              :let [sweep (do (println  "DC Analysis")
                              (time (doall (dc-analysis circuit source start stop step))))]]
        (print-result circuit sweep :dc source)
        (plot-result circuit sweep :dc source))
      (doseq [[_ time-step simulation-time] (:.tran (commands netlist))
              :let [series (do (println "Transient Analysis" time-step simulation-time)
                               (time (doall (transient-analysis circuit time-step simulation-time dc-result))))]]
        (print-result circuit series :tran "t")
        (plot-result circuit series :tran "t")))))

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
