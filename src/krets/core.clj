(ns krets.core
  (:require [clojure.string :as s]
            [clojure.walk :as w]
            [clojure.java.shell :as sh]
            [clojure.java.io :as io]
            [clojure.pprint :as pp])
  (:import [java.awt Color]
           [javax.swing JFrame]
           [org.jfree.chart ChartFactory ChartPanel]
           [org.jfree.chart.plot XYPlot]
           [org.jfree.data.xy XYSeries XYSeriesCollection]
           [org.ejml.data DenseMatrix64F]
           [org.ejml.ops CommonOps MatrixFeatures]
           [org.ejml.interfaces.linsol LinearSolver]
           [org.ejml.factory LinearSolverFactory]
           [javax.sound.sampled AudioSystem AudioFormat
            AudioFileFormat$Type AudioInputStream AudioFormat$Encoding]
           [java.io ByteArrayInputStream]
           [java.nio ByteBuffer]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; Matrix ops

(defn mtag [m]
  (with-meta m {:tag `DenseMatrix64F}))

(defn stag [s]
  (with-meta s {:tag `LinearSolver}))

(definline zero-matrix [rows cols]
  `(DenseMatrix64F. (unchecked-int ~rows) (unchecked-int ~cols)))

(definline linear-solver [d]
  `(LinearSolverFactory/linear ~d))

(definline set-a-matrix! [s a]
  `(.setA ~(stag s) ~a))

(definline solve [s b]
  `(doto (zero-matrix (.numRows ~(mtag b)) 1)
     (->> (.solve ~(stag s) ~b))))

(definline equals [a b epsilon]
  `(MatrixFeatures/isEquals ~(mtag a) ~b ~epsilon))

(definline set-matrix! [a b]
  `(doto ~(mtag a) (.set ~(mtag b))))

(definline mget [m row col]
  `(.unsafe_get ~(mtag m) (unchecked-int ~row) (unchecked-int ~col)))

(definline madd! [m row col v]
  `(.add ~(mtag m) (unchecked-int ~row) (unchecked-int ~col) ~v))

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
  (count (mapcat netlist [:v :e :opamp])))

(defn element-type [[[t] & nodes]]
  (let [et (low-key t)]
    (if (= :x et)
      (low-key (first (filter string? nodes)))
      et)))

(defn element-nodes [[_ & nodes :as e]]
  (take ({:e 4 :opamp 3} (element-type e) 2) nodes))

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
    (into {} (for [[^long idx [id]] (map-indexed vector (mapcat netlist [:v :e :opamp]))]
               [id (inc (+ number-of-nodes idx))]))))

(defn commands [netlist]
  (group-by (comp low-key first) (:. netlist)))

(defn sub-commands [m]
  (group-by (comp low-key second) m))

(defn models [netlist]
  (->> (for [[_ n t & kvs] (:.model (commands netlist))]
         {n (->> (for [[k v] (partition 2 kvs)]
                   [(low-key k) v])
                 (into {:model-type (low-key t)}))})
       (apply merge)))

(defn time-step [netlist]
  (if-let [[[_ time-step]] (:.tran (commands netlist))]
    time-step
    Double/NaN))

(defn options [netlist]
  (->> (for [[k v] (partition 2 (mapcat rest (:.options (commands netlist))))]
         {(low-key k) v})
       (apply merge {:tnom 27.0})))

(defn non-linear? [netlist]
  (boolean (some netlist [:d :opamp])))

(defn elements [netlist]
  (mapcat val (dissoc netlist :.)))

(defn node-ids [netlist]
  (let [{ns true ss false} (group-by number? (mapcat element-nodes (elements netlist)))]
    (merge (zipmap ns ns) (zipmap (set ss) (remove (conj (set ns) 0.0) (map double (range)))))))

(defn reassign-element-nodes [id-map [id & data :as e]]
  (let [[nodes data] (split-at (count (element-nodes e)) data)]
    (vec (concat [id] (replace id-map nodes) data))))

(defn reassign-netlist-nodes [id-map netlist]
  (let [es (elements netlist)]
    (w/postwalk-replace (zipmap es (map (partial reassign-element-nodes id-map) es)) netlist)))

(defn circuit-info [{:keys [options] :as circuit}]
  (dissoc circuit :netlist :title :mna-stamp :voltage-source->index :solver :models :options :node-ids))

(defn parse-netlist [netlist-source]
  (let [[title & lines] (-> netlist-source
                            (s/replace #"\n\+" "")
                            s/split-lines)
        netlist (->> lines
                     (remove (some-fn (re? #"^\*") (re? #"(?i)^.end$")))
                     (map #(s/split % #"[\s=,()]+"))
                     (w/postwalk (some-fn spice-number identity))
                     (group-by element-type))
        id-map (node-ids netlist)
        netlist (reassign-netlist-nodes id-map netlist)]
    (apply merge (with-meta {:netlist netlist :title title :node-ids id-map} {:netlist-source netlist-source})
           (for [k '[number-of-nodes number-of-voltage-sources number-of-rows
                     time-step models non-linear? voltage-source->index options]]
             {(keyword k) ((ns-resolve 'krets.core k) netlist)}))))

;; MNA Compiler

(defprotocol MNAStamp
  (linear-stamp! [_ a z x])
  (transient-stamp! [_ z x t])
  (non-linear-stamp! [_ a z x]))

(defmacro voltage-diff [x n+ n-]
  (let [term (fn [^long n]
               `(mget ~x ~(dec n) 0))]
    (cond
     (ground? n+) `(- ~(term n-))
     (ground? n-) (term n+)
     :else `(- ~(term n+) ~(term n-)))))

(defmacro safe-exp [x]
  (let [limit 70.0
        limit-exp (Math/exp limit)]
    `(if (< ~x ~limit)
       (Math/exp ~x)
       (+ ~limit-exp (+ 1.0 (- ~x ~limit))))))

(defn temprature-voltage [^double t]
  (let [zero-kelvin -273.15
        boltzmann 8.6173e-5]
    (* (- t zero-kelvin) boltzmann)))

(defmacro diode-current [is vd vt]
  `(* ~is (- (safe-exp (/ ~vd ~vt)) 1.0)))

(defmacro diode-conductance [is id vt]
  `(/ (+ ~id ~is) ~vt))

;; takes 1 based indexes, 0 is ground.
(defmacro stamp-matrix [m row col v]
  (when-not (or (ground? row) (ground? col))
    `(madd! ~m ~(dec (long row)) ~(dec (long col)) ~v)))

(defmacro conductance-stamp [a n+ n- g]
  (let [gs (gensym 'g)]
    `(let [~gs ~g]
       ~@(for [^long row [n+ n-]
               ^long col [n+ n-]]
           `(stamp-matrix ~a ~row ~col ~(if (= row col) gs `(- ~gs)))))))

(defmacro conductance-voltage-stamp [a n+ n- idx]
  `(do ~@(for [[^long n v] [[n+ 1.0] [n- -1.0]]]
           `(do (stamp-matrix ~a ~n ~idx ~v)
                (stamp-matrix ~a ~idx ~n ~v)))))

(defmacro source-current-stamp [z n+ n- i]
  `(let [i# ~i]
     (stamp-matrix ~z ~n+ 1 (- i#))
     (stamp-matrix ~z ~n- 1 i#)))

(defmacro code [& body]
  (let [ks (vec (keys &env))]
    `(cons 'do (->> '~body (w/postwalk-replace (zipmap '~ks ~ks))))))

(defmulti independent-source (fn [circuit {:keys [type]} [id n+ n- source-type & opts]]
                               [(if (string? source-type)
                                   (low-key source-type)
                                   :dc) type]))

(defmethod independent-source :default [_ _ _])

(defmethod independent-source [:dc :linear] [_ _ [_ _ _ & [t :as source]]]
  (or (first (filter number? source)) 0.0))

(defmethod independent-source [:sin :transient]
  [_  {:keys [t]} [_ _ _ _ vo va ^double freq td thet]]
  (let [td (or td 0.0)
        thet (double (or thet 0.0))
        freq-2-pi (* 2 Math/PI freq)]
    (code
     (let [t (double t)]
       (if (< t td)
         vo
         (+ vo (* va (if (zero? thet)
                       1.0
                       (Math/exp (- (/ (- t td) thet))))
                  (Math/sin (* freq-2-pi (+ t td))))))))))

(defmethod independent-source [:pulse :transient]
  [_  {:keys [t]} [_ _ _ _ ^double v1 ^double v2 td ^double tr ^double tf ^double pw per]]
  (let [td (double (or td 0.0))
        pw-end (+ tr pw)
        f-end (+ tr pw tf)]
    (code
     (let [t (double t)
           tp (rem (- t td) per)]
       (-> (cond
            (< t td) v1
            (< tp tr) (+ v1 (* (- v2 v1) (/ tp tr)))
            (< tp pw-end) v2
            (< tp f-end) (- v2 (* (- v2 v1) (/ (- tp tr pw) tf)))
            :else v1)
           double)))))

(def ^:dynamic *wave-table* (atom {}))

;; wavefile="test.wav" chan=0 ;; 0 is left, 1 is right etc. range is -1 to +1 v.
(defmethod independent-source [:wavefile :transient]
  [{:keys [^double time-step] :as circuit} {:keys [t]} [_ _ _ _ file _ ^double chan]]
  (let [simulation-sample-rate (/ time-step)
        filename (str (read-string file))
        f (io/file filename)
        netlist-file (-> circuit meta :netlist-file)
        f (if (and (not (.exists f)) (not (.isAbsolute f)) netlist-file)
            (io/file (.getParent (io/file netlist-file)) filename)
            f)
        file-size (.length f)
        in (AudioSystem/getAudioInputStream (io/input-stream f))
        number-of-channels (-> in .getFormat .getChannels)
        out-format (AudioFormat. AudioFormat$Encoding/PCM_FLOAT
                                 simulation-sample-rate 64 number-of-channels
                                 (* number-of-channels 8) simulation-sample-rate true)
        out (AudioSystem/getAudioInputStream out-format in)]
    (.mark out file-size)
    (swap! *wave-table* assoc file out)
    (code (let [w (@*wave-table* file)]
            (when (zero? (double t))
              (.reset w)
              (.mark w file-size))
            (let [bs (byte-array (* number-of-channels 8))]
              (.read w bs)
              (.getDouble (ByteBuffer/wrap bs) chan))))))

(defmulti stamp (fn [circuit {:keys [type] :as env} e] [(element-type e) type]))

(defmethod stamp :default [_ _ _])

(defmethod stamp [:r :linear] [_ {:keys [a]} [_ n+ n- ^double r]]
  (code (conductance-stamp a n+ n- (/ 1.0 r))))

(defmethod stamp [:c :linear] [{:keys [^double time-step]} {:keys [a]} [_ n+ n- ^double c]]
  (code (conductance-stamp a n+ n- (/ c time-step))))

(defmethod stamp [:c :transient] [{:keys [^double time-step]} {:keys [z x]} [_ n+ n- ^double c]]
  (let [geq (/ c time-step)]
    (code (let [ieq (- (* geq (voltage-diff x n+ n-)))]
            (source-current-stamp z n+ n- ieq)))))

(def ^:dynamic *voltage-sources* {})

(defmethod stamp [:v :linear] [{:keys [voltage-source->index netlist] :as circuit} {:keys [a z] :as env} [id n+ n- :as e]]
  (let [dc-sweep? ((low-key id) (-> netlist commands :.dc sub-commands))
        dc (independent-source circuit env e)
        idx (voltage-source->index id)]
    (code (cond
           dc-sweep? (stamp-matrix z idx 1 (double (*voltage-sources* id dc)))
           dc (stamp-matrix z idx 1 dc))
          (conductance-voltage-stamp a n+ n- idx))))

(defmethod stamp [:v :transient] [{:keys [voltage-source->index] :as circuit} {:keys [z] :as env} [id :as e]]
  (let [idx (voltage-source->index id)]
    (when-let [source (independent-source circuit env e)]
      (code (stamp-matrix z idx 1 source)))))

(defmethod stamp [:e :linear] [{:keys [voltage-source->index]} {:keys [a]}
                                       [id ^long out+ ^long out- ^long in+ ^long in- ^double gain]]
  (let [idx (voltage-source->index id)]
    (code (conductance-voltage-stamp a out- out+ idx)
          (stamp-matrix a idx in+ gain)
          (stamp-matrix a idx in- (- gain)))))

(defmethod stamp [:i :linear] [circuit {:keys [z] :as env} [_ n+ n- :as e]]
  (when-let [dc (independent-source circuit env e)]
    (code (source-current-stamp z n+ n- dc))))

(defmethod stamp [:i :transient] [circuit {:keys [z] :as env} [_ n+ n- :as e]]
  (when-let [source (independent-source circuit env e)]
    (code (source-current-stamp z n+ n- source))))

(defmethod stamp [:d :non-linear] [{:keys [models options]} {:keys [x a z]} [_ anode cathode model]]
  (let [defaults {:tnom (:tnom options) :is 1.0e-14}
        {:keys [^double is ^double tnom]} (merge defaults (models model))
        vt (temprature-voltage tnom)]
    (code (let [vd (voltage-diff x anode cathode)
                id (diode-current is vd vt)
                geq (diode-conductance is id vt)
                ieq (- id (* geq vd))]
            (conductance-stamp a anode cathode geq)
            (source-current-stamp z anode cathode ieq)))))

(defmethod stamp [:j :non-linear] [{:keys [models options]} {:keys [x a z]} [_ nd ng ns model]]
  (let [defaults {:tnom (:tnom options) :is 1.0e-14 :vto -2.0 :beta 1.0e-4 :lambda 0.0}
        {:keys [^double is ^double tnom model-type ^double vto ^double beta ^double lambda]} (merge defaults (models model))
        vt (temprature-voltage tnom)
        pol (case model-type
              :njf 1.0
              :pjf -1.0)]
    (code (let [vgs (* (voltage-diff x ng ns) pol)
                igs (diode-current is vgs vt)
                ggs (diode-conductance is igs vt)
                vgd (* (voltage-diff x ng nd) pol)
                igd (diode-current is vgd vt)
                ggd (diode-conductance is igd vt)
                vds (- vgs vgd)
                id-gm-gds (double-array 3)]
            (if (pos? vds)
              ;; normal mode
              (let [vgs-vto (- vgs vto)
                    b (* beta (+ 1.0 (* lambda vds)))]
                (cond
                 ;; cutoff
                 (<= vgs-vto 0.0) (do)
                 ;; saturation
                 (<= vgs-vto vds) (do (aset-double id-gm-gds 0 (* b vgs-vto vgs-vto))
                                      (aset-double id-gm-gds 1 (* b 2.0 vgs-vto))
                                      (aset-double id-gm-gds 2 (* lambda beta vgs-vto vgs-vto)))
                 ;; linear
                 :else (do (aset-double id-gm-gds 0 (* b vds (- (* 2.0 vgs-vto) vds)))
                           (aset-double id-gm-gds 1 (* b 2.0 vds))
                           (aset-double id-gm-gds 2 (+ (* b 2.0 (- vgs-vto vds))
                                                       (* lambda beta vds (- (* 2.0 vgs-vto) vds)))))))
              ;; inverse mode
              (let [vgd-vto (- vgd vto)
                    b (* beta (- 1.0 (* lambda vds)))]
                ;; cutoff
                (cond
                 (<= vgd-vto 0.0) (do)
                 ;; saturation
                 (<= vgd-vto (- vds)) (do (aset-double id-gm-gds 0 (- (* b vgd-vto vgd-vto)))
                                          (aset-double id-gm-gds 1 (- (* b 2.0 vgd-vto)))
                                          (aset-double id-gm-gds 2 (+ (* lambda beta vgd-vto vgd-vto)
                                                                      (* b 2.0 vgd-vto))))
                 ;; linear
                 :else (do (aset-double id-gm-gds 0 (* b vds (- (* 2.0 vgd-vto) vds)))
                           (aset-double id-gm-gds 1 (* b 2.0 vds))
                           (aset-double id-gm-gds 2 (- (* b 2.0 vgd-vto)
                                                       (* lambda beta vds (+ (* 2.0 vgd-vto) vds))))))))
            (let [id (aget id-gm-gds 0)
                  gm (aget id-gm-gds 1)
                  gds (aget id-gm-gds 2)
                  igseq (- igs (* ggs vgs))
                  igdeq (- igd (* ggd vgd))
                  idseq (- gm (* gm vgs) (* gds vds))]
              (stamp-matrix z ng 1 (* (- igseq (- igdeq)) pol))
              (stamp-matrix z nd 1 (* (- igdeq idseq) pol))
              (stamp-matrix z ns 1 (* (+ idseq igseq) pol))
              (stamp-matrix a ng ng (+ ggs ggd))
              (stamp-matrix a ng nd (- ggd))
              (stamp-matrix a ng ns (- ggs))
              (stamp-matrix a nd ng (- gm ggd))
              (stamp-matrix a nd nd (+ gds ggd))
              (stamp-matrix a nd ns (- (- gm) gds))
              (stamp-matrix a ns ng (- (- ggs) gm))
              (stamp-matrix a ns nd (- gds))
              (stamp-matrix a ns ns (+ ggs gds gm)))))))

(defmethod stamp [:q :non-linear] [{:keys [models options]} {:keys [x a z]} [_ nc nb ne model]]
  (let [defaults {:tnom (:tnom options) :is 1.0e-16 :bf 100.0 :ne 1.5 :nc 2.0}
        {:keys [^double is ^double tnom model-type]} (merge defaults (models model))]))

;; All About Circuits model ideal op amps as a vcvs.
(defn opamp->e [[id ^long in+ ^long in- ^long out+]]
  [(str "e" id) out+ 0 in+ in- 999e3])

(defmethod stamp [:opamp :linear] [{:keys [voltage-source->index]} {:keys [a]}
                                   [id ^long in+ ^long in- ^long out+]]
  (let [idx (voltage-source->index id)]
    (code (stamp-matrix a idx out+ -1.0)
          (stamp-matrix a out+ idx 1.0))))

(defmethod stamp [:opamp :non-linear] [{:keys [voltage-source->index]} {:keys [a z x]}
                                       [id ^long in+ ^long in- ^long out+]]
  (let [idx (voltage-source->index id)
        gain 1e6
        vmax 15
        pi-by-2-v-max (/ Math/PI (* 2 vmax))
        vmax-2-by-pi (* vmax (/ 2 Math/PI))]
    (code (let [vd (voltage-diff x in+ in-)
                tmp (* pi-by-2-v-max gain vd)
                v (* vmax-2-by-pi (Math/atan tmp))
                g (/ gain (+ 1 (Math/pow tmp 2)))]
            (stamp-matrix z idx 1 (- (* g vd) v))
            (stamp-matrix a idx in+ g)
            (stamp-matrix a idx in- (- g))))))

(defn compile-mna-stamp [{:keys [netlist] :as circuit}]
  `(reify MNAStamp
     ~@(for [[f {:keys [arglists]}] (:sigs MNAStamp)
             :let [args (first arglists)
                   syms (map gensym args)
                   t (low-key (s/replace (name f) "-stamp!" ""))
                   env (assoc (zipmap (map keyword args) syms) :type t)]]
         `(~(symbol (name f)) [~@syms]
           ~@(for [e (elements netlist)]
               (stamp circuit env e))))))

(defn compile-circuit [{:keys [mna-stamp ^long number-of-rows] :as circuit}]
  (if mna-stamp
    circuit
    (assoc circuit
      :mna-stamp (binding [*ns* (find-ns 'krets.core)]
                   (eval (compile-mna-stamp circuit)))
      :solver (linear-solver number-of-rows))))

;; MNA Analysis

(def ^:dynamic *newton-tolerance* 1e-8)
(def ^:dynamic *newton-iterations* 500)

(defn non-linear-step-fn [{:keys [mna-stamp solver ^long number-of-rows]}]
  (let [newton-tolerance (double *newton-tolerance*)
        newton-iterations (long *newton-iterations*)
        anl (zero-matrix number-of-rows number-of-rows)
        znl (zero-matrix number-of-rows 1)]
    (fn [a z x]
      (loop [xn-1 x
             iters newton-iterations]
        (if (zero? iters)
          (throw (ex-info "Didn't converge." {:x xn-1}))
          (do (set-matrix! anl a)
              (set-matrix! znl z)
              (non-linear-stamp! mna-stamp anl znl xn-1)
              (set-a-matrix! solver anl)
              (let [xn (solve solver znl)]
                (if (equals xn xn-1 newton-tolerance)
                  xn
                  (recur xn (dec iters))))))))))

(defn linear-step-fn [{:keys [solver]}]
  (fn [a z _]
    (solve solver z)))

(defn step-fn [{:keys [non-linear?] :as circuit}]
  ((if non-linear?
     non-linear-step-fn
     linear-step-fn) circuit))

(defn dc-operating-point
  ([{:keys [^long number-of-rows] :as circuit}]
   (dc-operating-point circuit (zero-matrix number-of-rows 1)))
  ([{:keys [mna-stamp solver ^long number-of-rows] :as circuit} x]
   (let [step (step-fn circuit)
         a (zero-matrix number-of-rows number-of-rows)
         z (zero-matrix number-of-rows 1)]
     (linear-stamp! mna-stamp a z x)
     (set-a-matrix! solver a)
     {:a a :z z :x (step a z x)})))

(defn dc-analysis [circuit source start stop step]
  (doall
   (for [v (concat (range start stop step) [stop])]
     (binding [*voltage-sources* (assoc *voltage-sources* source v)]
       [v (-> circuit dc-operating-point :x)]))))

(defn transient-analysis
  ([circuit ^double time-step ^double simulation-time]
   (transient-analysis circuit time-step simulation-time (dc-operating-point circuit)))
  ([{:keys [mna-stamp number-of-rows] :as circuit}
    ^double time-step ^double simulation-time {:keys [a x] z0 :z}]
   (let [step (step-fn circuit)
         end (+ simulation-time time-step)
         number-of-rows (long number-of-rows)
         z (zero-matrix number-of-rows 1)]
     (loop [t 0.0
            idx 0
            x x
            acc (object-array (/ end time-step))]
       (if (> t end)
         acc
         (do (set-matrix! z z0)
             (transient-stamp! mna-stamp z x t)
             (let [x (step a z x)]
               (recur (+ t time-step) (inc idx) x (doto acc (aset idx [t x]))))))))))

;; Frontend

(defn report-node-label [[k n+ n-]]
  (let [f #(cond-> % (number? %) long)]
    (format "%s(%s)" k (s/join ","  (cond-> [(f n+)] n- (concat [(f n-)]))))))

(defn report-node-voltage [[k & ns] {:keys [^long number-of-nodes node-ids]} x]
  (->> (for [n ns
             :let [n (long (node-ids n 0))]]
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

(defn print-result [{:keys [^long number-of-nodes time-step netlist] :as circuit} series series-type head-label]
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
                   (format node-format (report-node-voltage node circuit x))})))))))

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

(defn plot-result [{:keys [^long number-of-nodes netlist] :as circuit} series series-type head-label]
  (let [xs (map first series)
        ys (map second series)]
    (doseq [[_ _ & nodes] (-> netlist commands :.plot sub-commands series-type)
            :let [nodes (report-nodes nodes)]
            :when (seq nodes)]
      (apply plot head-label "V"
             (for [node nodes]
               (xy-series (report-node-label node)
                          xs
                          (map #(report-node-voltage node circuit %) ys)))))))

(def ^:dynamic *normalize-wave* true)

;; .wave "test.wav" 16 48000 v(1,0) v(2,3)
(defn write-wave [file sample-rate bits simulation-sample-rate & channels]
  (let [number-of-channels (count channels)
        number-of-samples (count (first channels))
        max-sample (if *normalize-wave*
                     (double (apply max (apply concat channels)))
                     1.0)
        normalize (fn [^double d]
                    (max (min (/ d max-sample) 1.0) -1.0))
        simulation-sample-rate (Math/round (double simulation-sample-rate))
        out-format (AudioFormat. sample-rate bits number-of-channels true false)
        in-format (AudioFormat. AudioFormat$Encoding/PCM_FLOAT
                                simulation-sample-rate 64 number-of-channels
                                (* number-of-channels 8) simulation-sample-rate true)
        data (.array ^ByteBuffer (reduce (fn [^ByteBuffer b ds]
                                           (doseq [^double d ds]
                                             (doto b (.putDouble (normalize d))))
                                           b)
                                         (ByteBuffer/wrap (byte-array (* 8 number-of-samples number-of-channels)))
                                         (apply map vector channels)))
        in (AudioInputStream. (ByteArrayInputStream. data) in-format number-of-samples)]
    (AudioSystem/write (AudioSystem/getAudioInputStream out-format in)
                       AudioFileFormat$Type/WAVE (io/file (str (read-string file))))))

(defn wave-result [{:keys [^long number-of-nodes netlist ^double time-step] :as circuit} series series-type _]
  (let [ys (map second series)]
    (doseq [[_ file bits sample-rate & nodes] (-> netlist commands :.wave)
            :let [nodes (report-nodes nodes)]
            :when (seq nodes)]
      (apply write-wave file sample-rate bits (/ time-step)
             (for [node nodes]
               (map #(report-node-voltage node circuit %) ys))))))

(defn batch [{:keys [title netlist models] :as circuit}]
  (let [circuit (compile-circuit circuit)]
    (println title)
    (pp/print-table [(circuit-info circuit)])
    (when models
      (pp/print-table [models]))
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
                              (time (dc-analysis circuit source start stop step)))]]
        (print-result circuit sweep :dc source)
        (plot-result circuit sweep :dc source))
      (doseq [[_ time-step simulation-time start] (:.tran (commands netlist))
              :let [series (do (println "Transient Analysis" time-step simulation-time)
                               (time (transient-analysis circuit time-step simulation-time dc-result)))
                    series (cond->> series
                                    (number? start) (drop-while (fn [[^double t]] (< t (double start)))))]
              f [print-result plot-result wave-result]]
        (f circuit series :tran "t")))))

(defn spice
  ([circuit]
   (spice ["ngspice" "-b"] circuit))
  ([cmd circuit]
   (let [{:keys [out err ^long exit]} (apply sh/sh (concat cmd  [:in (-> circuit meta :netlist-source)]))]
     (when out
       (println out))
      (when (not (zero? exit))
        (println err)))))

(defn process-file [f]
  (-> f slurp parse-netlist (with-meta {:netlist-file f}) batch))

(defn -main [& [f]]
  (if f
    (process-file f)
    (println "Need to specify a netlist file")))
