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
           [org.ejml.interfaces.linsol LinearSolver]
           [org.ejml.factory LinearSolverFactory]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; Matrix ops

(defn mtag [m]
  (with-meta m {:tag `DenseMatrix64F}))

(defn stag [s]
  (with-meta s {:tag `LinearSolver}))

(definline zero-matrix [rows cols]
  `(DenseMatrix64F. ~rows ~cols))

(definline row-count [m]
  `(.numRows ~(mtag m)))

(definline solve [s a b]
  `(do (.setA ~(stag s) ~a)
       (doto (zero-matrix (row-count ~a) 1)
         (->> (.solve ~(stag s) ~b)))))

(definline add! [a b]
  `(doto ~(mtag a) (-> (CommonOps/addEquals ~b))))

(definline equals [a b epsilon]
  `(MatrixFeatures/isEquals ~(mtag a) ~b ~epsilon))

(definline mget [m row col]
  `(.unsafe_get ~(mtag m) ~row ~col))

(definline madd! [m row col v]
  `(.add ~(mtag m) ~row ~col ~v))

(definline zero! [m]
  `(doto ~(mtag m) .zero))

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

(defn element-type [[[c]]]
  (low-key c))

(defn element-nodes [[_ & nodes :as e]]
  (take ({:e 4 :u 3} (element-type e) 2) nodes))

(defn number-of-nodes ^long [netlist]
  (->> (dissoc netlist :.)
       vals
       (apply concat)
       (mapcat element-nodes)
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

(defn circuit-info [circuit]
  (dissoc circuit :netlist :title :mna-stamp :voltage-source->index :solver))

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
  (linear-stamp! [_ a z x])
  (transient-stamp! [_ z x t])
  (non-linear-stamp! [_ a z x]))

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

(defn conductance-voltage-stamp [a n+ n- idx]
  `(do ~@(for [[^long n v] [[n+ 1.0] [n- -1.0]]
               :when (not (ground? n))
               :let [n (dec n)]]
           `(do (madd! ~a ~n ~idx ~v)
                (madd! ~a ~idx ~n ~v)))))

(defn sine-source [[vo va ^double freq td thet]]
  (fn [t]
    (let [td (or td 0.0)
          thet (double (or thet 0.0))]
      `(if (< ~t ~td)
         ~vo
         (+ ~vo (* ~va ~(if (zero? thet)
                          1.0
                          `(Math/exp (- (/ (- ~t ~td) ~thet))))
                   (Math/sin (* ~(* 2 Math/PI freq) (+ ~t ~td)))))))))

(defn pulse-source [[^double v1 ^double v2 td ^double tr ^double tf ^double pw per]]
  (fn [t]
    (let [td (double (or td 0.0))]
      `(let [t# (double ~t)
             tp# (rem (- t# ~td) ~per)]
         (cond
          (< t# ~td) ~v1
          (< tp# ~tr) (+ ~v1 (* ~(- v2 v1) (/ tp# ~tr)))
          (< tp# ~(+ tr pw)) ~v2
          (< tp# ~(+ tr pw tf)) (- ~v2 (* ~(- v2 v1) (/ (- tp# ~tr ~pw) ~tf)))
          :else ~v1)))))

(defn independent-source [[id n+ n- & [t & opts :as source]]]
  (let [t (if (string? t)
            (low-key t)
            :dc)
        dc (or (first (filter number? source)) 0.0)
        f (case t
            :dc (constantly dc)
            :sin (sine-source opts)
            :pulse (pulse-source opts))]
    {:dc dc :type t :transient f}))

(defn source-current-stamp [z n+ n- in out]
  `(do ~@(for [[^long row i] [[n+ in] [n- out]]
               :when (not (ground? row))]
           `(madd! ~z ~(dec row) 0 ~i))))

(defmulti stamp-element (fn [circuit {:keys [stamp] :as env} e] [(element-type e) stamp]))

(defmethod stamp-element :default [_ _ _])

(defmethod stamp-element [:r :linear] [_ {:keys [a]} [_ n+ n- ^double r]]
  (conductance-stamp a n+ n- (/ 1.0 r)))

(defmethod stamp-element [:c :linear] [{:keys [^double time-step]} {:keys [a]} [_ n+ n- ^double c]]
  (conductance-stamp a n+ n- (/ c time-step)))

(defmethod stamp-element [:c :transient] [{:keys [^double time-step]} {:keys [z x]} [_ n+ n- ^double c]]
  (let [ieq (gensym 'ieq)
        geq (/ c time-step)]
    `(let [~ieq (- (* ~geq ~(voltage-diff x n+ n-)))]
       ~(source-current-stamp z n+ n- `(- ~ieq) ieq))))

(defmethod stamp-element [:d :non-linear] [{:keys [models]} {:keys [x a z]} [_ anode cathode model]]
  (let [[ieq geq] (map gensym '[ieq geq])
        vt 0.025875
        {:keys [^double is]} (merge {:is 1.0e-14} (models model))]
    `(let [vd# ~(voltage-diff x anode cathode)
           exp-vd-by-vt# (Math/exp (/ vd# ~vt))
           id# (* ~is (- exp-vd-by-vt# 1.0))
           ~geq (* ~(/ is vt) exp-vd-by-vt#)
           ~ieq (+ (- id#) (* ~geq vd#))]
       ~(conductance-stamp a anode cathode geq)
       ~(source-current-stamp z anode cathode ieq `(- ~ieq)))))

(def ^:dynamic *voltage-sources* {})

(defmethod stamp-element [:v :linear] [{:keys [voltage-source->index netlist]} {:keys [a z]} [id n+ n- :as e]]
  (let [dc-sweep? ((low-key id) (-> netlist commands :.dc sub-commands))
        {:keys [dc type]} (independent-source e)
        idx (voltage-source->index id)]
    `(do ~(when (= :dc type)
            `(madd! ~z ~idx 0 ~(if dc-sweep?
                                 `(double (*voltage-sources* ~id ~dc))
                                 dc)))
         ~(conductance-voltage-stamp a n+ n- (voltage-source->index id)))))

(defmethod stamp-element [:v :transient] [{:keys [voltage-source->index]} {:keys [z t]} [id :as e]]
  (let [{:keys [transient type]} (independent-source e)
        idx (voltage-source->index id)]
    (when-not (= :dc type)
      `(madd! ~z ~idx 0 ~(transient t)))))

(defmethod stamp-element [:e :linear] [{:keys [voltage-source->index]} {:keys [a]}
                                       [id ^long out+ ^long out- ^long in+ ^long in- ^double gain]]
  (let [idx (voltage-source->index id)]
    `(do ~(conductance-voltage-stamp a out- out+ idx)
         ~@(for [[^long n g] [[in+ gain] [in- (- gain)]]
                 :when (not (ground? n))]
             `(madd! ~a ~idx ~(dec n) ~g)))))

(defmethod stamp-element [:i :linear] [_ {:keys [z]} [_ n+ n- :as e]]
  (let [{:keys [^double dc type]} (independent-source e)]
    (when (= :dc type)
      (source-current-stamp z n+ n- (- dc) dc))))

(defmethod stamp-element [:i :transient] [_ {:keys [z t]} [_ n+ n- :as e]]
  (let [i (gensym 'i)
        {:keys [transient type]} (independent-source e)]
    (when-not (= :dc type)
      `(let [~i ^double ~(transient t)]
         ~(source-current-stamp z n+ n- `(- i) i)))))

;; All About Circuits model ideal op amps as a vcvs.
(defn u->e [[id ^long in+ ^long in- ^long out+]]
  [(str "e" id) out+ 0 in+ in- 999e3])

;; ideal op amp, this is only the linear version. QUCS technical.pdf p 117. QUCS source has a transient part.
(defmethod stamp-element [:u :linear] [{:keys [voltage-source->index]} {:keys [a]}
                                       [id ^long in+ ^long in- ^long out+]]
  (let [idx (inc (long (voltage-source->index id)))]
    `(do ~@(for [[^long n+ ^long n- v] [[idx in+ 1] [idx in- -1] [out+ idx 1]]
                 :when (not (or (ground? n+) (ground? n-)))]
             `(madd! ~a ~(dec n+) ~(dec n-) ~v)))))

(defn compile-mna-stamp [{:keys [netlist] :as circuit}]
  `(reify MNAStamp
     ~@(for [[f {:keys [arglists]}] (:sigs MNAStamp)
             :let [args (first arglists)
                   syms (map gensym args)
                   stamp (low-key (s/replace (name f) "-stamp!" ""))
                   env (assoc (zipmap (map keyword args) syms) :stamp stamp)]]
         `(~(symbol (name f)) [~@syms]
           ~@(for [e (elements netlist)]
               (stamp-element circuit env e))))))

(defn compile-circuit [{:keys [mna-stamp ^long number-of-rows] :as circuit}]
  (if mna-stamp
    circuit
    (assoc circuit
      :mna-stamp (eval (compile-mna-stamp circuit))
      :solver (LinearSolverFactory/linear number-of-rows))))

;; MNA Analysis

(defn dc-operating-point
  ([{:keys [^long number-of-rows] :as circuit}]
   (dc-operating-point circuit (zero-matrix number-of-rows 1)))
  ([{:keys [mna-stamp solver ^long number-of-rows]} x]
   (let [a (zero-matrix number-of-rows number-of-rows)
         z (zero-matrix number-of-rows 1)]
     (linear-stamp! mna-stamp a z x)
     {:a a :z z :x (solve solver a z)})))

(defn dc-analysis [circuit source start stop step]
  (for [v (concat (range start stop step) [stop])]
    (binding [*voltage-sources* (assoc *voltage-sources* source v)]
      [v (-> circuit dc-operating-point :x)])))

(def ^:dynamic *newton-tolerance* 1e-8)
(def ^:dynamic *newton-iterations* 500)

(defn non-linear-step-fn [{:keys [mna-stamp solver ^long number-of-rows]}]
  (let [newton-tolerance (double *newton-tolerance*)
        newton-iterations (long *newton-iterations*)
        anl (zero-matrix number-of-rows number-of-rows)
        znl (zero-matrix number-of-rows 1)]
    (fn [a z x]
      (loop [xn-1 x
             iters newton-iterations
             anl (zero! anl)
             znl (zero! znl)]
        (if (zero? iters)
          (throw (ex-info "Didn't converge." {:x xn-1}))
          (let [anl (add! anl a)
                znl (add! znl z)]
            (non-linear-stamp! mna-stamp anl znl xn-1)
            (let [xn (solve solver anl znl)]
              (if (equals xn xn-1 newton-tolerance)
                xn
                (recur xn (dec iters) (zero! anl) (zero! znl))))))))))

(defn linear-step-fn [{:keys [solver]}]
  (fn [a z _]
    (solve solver a z)))

(defn transient-step-fn [{:keys [non-linear?] :as circuit}]
  ((if non-linear?
     non-linear-step-fn
     linear-step-fn) circuit))

(defn transient-analysis
  ([circuit ^double time-step ^double simulation-time]
   (transient-analysis circuit time-step simulation-time (dc-operating-point circuit)))
  ([{:keys [mna-stamp number-of-rows] :as circuit}
    ^double time-step ^double simulation-time {:keys [a x] z0 :z}]
   (let [step (transient-step-fn circuit)
         end (+ simulation-time time-step)
         number-of-rows (long number-of-rows)]
     (loop [t 0.0 x x acc (transient []) z (zero-matrix number-of-rows 1)]
       (if (> t end)
         (persistent! acc)
         (do (add! z z0)
             (transient-stamp! mna-stamp z x t)
             (let [x (step a z x)]
               (recur (+ t time-step) x (conj! acc [t x]) (zero! z)))))))))

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
    (println)
    (let [{:keys [a z x] :as dc-result} (do (println "DC Operating Point Analysis")
                                            (time (dc-operating-point circuit)))]
      (println)
      (println "A")
      (println a)
      (println "z")
      (println z)
      (println "x")
      (println x)
      (println)
      (doseq [[_ source start stop step] (:.dc (commands netlist))
              :let [sweep (do (println  "DC Analysis")
                              (time (doall (dc-analysis circuit source start stop step))))]]
        (print-result circuit sweep :dc source)
        (plot-result circuit sweep :dc source))
      (doseq [[_ time-step simulation-time start] (:.tran (commands netlist))
              :let [series (do (println "Transient Analysis" time-step simulation-time)
                               (time (doall (transient-analysis circuit time-step simulation-time dc-result))))
                    series (cond->> series
                                    (number? start) (drop-while (fn [[^double t]] (< t (double start)))))]]
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
