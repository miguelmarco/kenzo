
(defsystem :kenzo-test
  :description "A simple regression test suite for Kenzo"
  :serial t
  :version "0.1.0"
  :author "Francis Sergeraert <Francis.Sergeraert@ujf-grenoble.fr>"
  :license "GPLv3"
  :depends-on (:kenzo :fiveam)
  :pathname #P"test/"
  :components ((:file "package")
	       (:file "common")
	       ;; === test files ===
	       (:file "bar-test")
	       (:file "bicones-test")
	       (:file "cartesian-products-test")
               (:file "ccn")
	       (:file "chain-complexes-test")
	       (:file "chcm-elementary-op-test")
	       (:file "circle")
 	       (:file "cl-space-efhm-test")
 	       (:file "classifying-spaces-test")
	       (:file "cobar-test")
	       (:file "combinations-test")
	       (:file "cones-test")
	       (:file "cs-twisted-products-test")
	       (:file "delta-test")
	       (:file "diabolo")
	       (:file "effective-homology-test")
	       (:file "homology-groups-test")
	       (:file "k-pi-n-test")
	       (:file "kan-test")
	       (:file "searching-homology-test")
	       (:file "simplicial-groups-test")
	       (:file "simplicial-mrphs-test")
	       (:file "simplicial-sets-test")
	       (:file "tensor-products-test")
	       (:file "various-test")
	       ))
