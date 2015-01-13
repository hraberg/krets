(ns krets.bench-test
  (:use [krets.core]
        [clojure.test])
  (:require [criterium.core :as cc]))

(defn micro-bench [f]
  (println "running" f)
  (let [c (-> f slurp parse-netlist compile-circuit)
        stub (constantly nil)]
    (with-redefs [plot-result stub
                  print-result stub
                  println stub]
      (cc/quick-bench (batch c))
      (is true))))

(deftest micro-bench-non-linear
  (micro-bench "test/krets/Non-Linear_DC_ckt.cir"))

(deftest micro-bench-transient
  (micro-bench "test/krets/Transient_DC_ckt.cir"))
