
(defsystem :kenzo-test
  :description "Test suite for kenzo"
  :serial t
  :version "0.1"
  :depends-on (:kenzo :fiveam)
  :pathname #P"test/"
  :components ((:file "package")
	       (:file "combinations-test")
	       (:file "various-test")))
