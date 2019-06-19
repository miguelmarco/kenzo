;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*

(in-package :kenzo-test)

(in-suite :kenzo)


(test homotopy-s2
      (cat:cat-init)
      (let* ((s2 (cat:sphere 2))
             )
        (is (equalp (cat:homotopy-list s2 2 '(1))))
        (is (equalp (cat:homotopy-list s2 3 '(1))))
        (is (equalp (cat:homotopy-list s2 4 '(2))))
        (is (equalp (cat:homotopy-list s2 5 '(2))))
        (is (equalp (cat:homotopy-list s2 6 '(12))))
        ))


(test homotopy-s3
      (cat:cat-init)
      (let* ((s3 (cat:sphere 3))
             )
        (is (equalp (cat:homotopy-list s3 2 nil)))
        (is (equalp (cat:homotopy-list s3 3 '(1))))
        (is (equalp (cat:homotopy-list s3 4 '(2))))
        (is (equalp (cat:homotopy-list s3 5 '(2))))
        (is (equalp (cat:homotopy-list s3 6 '(12))))
        ))


(test homotopy-s4
      (cat:cat-init)
      (let* ((s4 (cat:sphere 4))
             )
        (is (equalp (cat:homotopy-list s4 2 nil)))
        (is (equalp (cat:homotopy-list s4 3 nil)))
        (is (equalp (cat:homotopy-list s4 4 '(1))))
        (is (equalp (cat:homotopy-list s4 5 '(2))))
        (is (equalp (cat:homotopy-list s4 6 '(2))))
        (is (equalp (cat:homotopy-list s4 7 '(12 1))))
        ))           


(test homotopy-s5
      (cat:cat-init)
      (let* ((s5 (cat:sphere 5))
             )
        (is (equalp (cat:homotopy-list s5 2 nil)))
        (is (equalp (cat:homotopy-list s5 3 nil)))
        (is (equalp (cat:homotopy-list s5 4 nil)))
        (is (equalp (cat:homotopy-list s5 5 '(1))))
        (is (equalp (cat:homotopy-list s5 6 '(2))))
        (is (equalp (cat:homotopy-list s5 7 '(2))))
        (is (equalp (cat:homotopy-list s5 8 '(24))))
        ))

(test homotopy-s6
      (cat:cat-init)
      (let* ((s6 (cat:sphere 6))
             )
        (is (equalp (cat:homotopy-list s6 2 nil)))
        (is (equalp (cat:homotopy-list s6 3 nil)))
        (is (equalp (cat:homotopy-list s6 4 nil)))
        (is (equalp (cat:homotopy-list s6 5 nil)))
        (is (equalp (cat:homotopy-list s6 6 '(1))))
        (is (equalp (cat:homotopy-list s6 7 '(2))))
        (is (equalp (cat:homotopy-list s6 8 '(2))))
        (is (equalp (cat:homotopy-list s6 9 '(24))))
        ))
