(defproject krets "0.1.0-SNAPSHOT"
  :description "Modified Nodal Analysis"
  :url "http://github.com/hraberg/krets"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.7.0-alpha5"]
                 [org.jfree/jfreechart "1.0.19"]
                 [com.googlecode.efficient-java-matrix-library/ejml "0.25"]]
  :profiles {:dev {:dependencies [[criterium "0.4.3"]]}}
  :jvm-opts ^:replace []
  :main krets.core)
