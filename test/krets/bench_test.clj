(ns krets.bench-test
  (:use [krets.core]
        [clojure.test])
  (:require [criterium.core :as cc]))

(defn micro-bench [f]
  (println "running" f)
  (with-redefs [plot-result (constantly nil)
                print-result (constantly nil)
                println (constantly nil)]
    (let [c (-> f slurp parse-netlist compile-circuit)]
      (cc/quick-bench (batch c))
      (is true))))

(deftest micro-bench-non-linear
  (micro-bench "test/krets/Non-Linear_DC_ckt.cir"))

(deftest micro-bench-transient
  (micro-bench "test/krets/Transient_DC_ckt.cir"))
