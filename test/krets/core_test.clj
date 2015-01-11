(ns krets.core-test
  (:use [krets.core]
        [clojure.test])
  (:require [clojure.java.io :as io]))

(deftest sanity-check
  (doseq [f (file-seq (io/file "test/krets/"))
          :when (re-find #"\.cir$" (str f))]
    (with-redefs [plot (constantly nil)
                  print-result (constantly nil)]
      (let [fail (volatile! nil)
            out (with-out-str
                  (try
                    (process-file f)
                    (is true)
                    (catch Exception e
                      (vreset! fail e))))]
        (when-let [e @fail]
          (println out)
          (throw e))))))
