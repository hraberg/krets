(ns krets.core-test
  (:use [krets.core]
        [clojure.test])
  (:require [clojure.java.io :as io]
            [clojure.core.matrix.implementations :refer [*matrix-implementation*]]))

(deftest sanity-check
  (let [matrix-backend *matrix-backend*]
    (try
      (doseq [matrix-backend [[:core.matrix :vectorz] [:core.matrix :clatrix] [:ejml :native-ejml]]
              f (file-seq (io/file "test/krets/"))
              :when (re-find #"\.cir$" (str f))]
        (load-matrix-backend matrix-backend)
        (with-redefs [plot (fn [xt yt & _]
                             (println "< plotting " xt yt ">"))]
          (let [fail (volatile! nil)
                out (with-out-str
                      (try
                        (process-file f)
                        (is true)
                        (catch Exception e
                          (vreset! fail e))))]
            (when-let [e @fail]
              (println matrix-backend)
              (println out)
              (throw e)))))
      (finally
        (load-matrix-backend matrix-backend)))))
