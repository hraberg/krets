(ns krets.bench-test
  (:use [krets.core]
        [clojure.test])
  (:require [criterium.core :as cc]
            [clojure.pprint :as pp]))

(defn micro-bench [f]
  (let [c (-> f slurp parse-netlist compile-circuit)
        stub (constantly nil)]
    (println "running" f (:number-of-nodes c) "nodes")
    (with-redefs [plot-result stub
                  print-result stub
                  println stub
                  pp/print-table stub]
      (cc/quick-bench
       (binding [*out* (proxy [java.io.Writer] [] (write [c]) (flush []))]
         (batch c)))
      (is true))))

(deftest micro-bench-non-linear
  (micro-bench "test/krets/Non-Linear_DC_ckt.cir"))

(deftest micro-bench-transient
  (micro-bench "test/krets/Transient_DC_ckt.cir"))

(deftest micro-bench-fullwave-bridge-rectifier
  (micro-bench "test/krets/fullwave-bridge-rectifier.cir"))

(deftest micro-bench-instrumentation-amplifier
  (micro-bench "test/krets/instrumentation-amplifier.cir"))
