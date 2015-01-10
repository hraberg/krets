(defproject krets "0.1.0-SNAPSHOT"
  :description "Modified Nodal Analysis"
  :url "http://github.com/hraberg/krets"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.7.0-alpha5"]
                 [org.jfree/jfreechart "1.0.19"]
                 [net.mikera/core.matrix "0.32.1"]
                 [net.mikera/vectorz-clj "0.28.0"]
                 [clatrix "0.4.0"]
                 [com.googlecode.efficient-java-matrix-library/ejml "0.25"]]
  :main krets.core)
