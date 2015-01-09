(defproject krets "0.1.0-SNAPSHOT"
  :description "Modified Nodal Analysis"
  :url "http://github.com/hraberg/krets"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.7.0-alpha4"]
                 [org.jfree/jfreechart "1.0.19"]
                 [net.mikera/core.matrix "0.31.1"]
                 [net.mikera/vectorz-clj "0.26.2" :exclude [core.matrix]]
                 [clatrix "0.4.0" :exclude [core.matrix]]])
