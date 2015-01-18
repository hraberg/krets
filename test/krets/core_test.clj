(ns krets.core-test
  (:use [krets.core]
        [clojure.test])
  (:require [clojure.java.io :as io]))

(defn test-circuits []
  (->> (io/file "test/krets/")
       .listFiles
       (filter (re? #"\.cir$"))))

(deftest sanity-check
  (doseq [f (test-circuits)
          :let [stub (constantly nil)]]
    (with-redefs [plot-result stub
                  print-result stub
                  wave-result stub]
      (let [fail (volatile! nil)
            out (with-out-str
                  (try
                    (process-file f)
                    (is true)
                    (catch Exception e
                      (vreset! fail e))))]
        (when-let [e @fail]
          (println f)
          (println out)
          (throw e))))))
