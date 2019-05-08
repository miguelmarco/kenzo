;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*

;;;  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY
;;;  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY
;;;  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY

;;; Authors: Jonathan Heras and Ana Romero


(IN-PACKAGE #:cat)

(provide "homotopy")


(DEFUN CHCM-HOMOLOGY-FORMAT (cc n)
  (declare (type chain-complex cc) (type fixnum n))
  (let ((rsl (homologie (chcm-mat cc n) (chcm-mat cc (1+ n))))
        (str nil))
    (declare (type list rsl))
    (dolist (item rsl)
      (declare (type list item))
      (setf str (concatenate 'string str (format nil "Z")))
      (unless (zerop (first item)) 
        (setf str (concatenate 'string str (format nil "/~DZ" (first item)))))
      (setf str (concatenate 'string str " ")))
    str))


(DEFUN HOMOLOGY-FORMAT (chcm degr1 &optional (degr2 (1+ degr1)))
  (declare (fixnum degr1 degr2))
  (let ((result_hom nil))
  (do ((degr degr1 (1+ degr)))
      ((>= degr degr2))
    (declare (fixnum degr))
    (setf result_hom (chcm-homology-format (echcm chcm) degr))
    (terpri) (clock) (terpri))
  result_hom))


(DEFUN FIRST-NON-NULL-HOMOLOGY-GROUP-AUX (chcm n limit)
  (if (> n limit)
      nil
    (let ((hf (homology-format chcm n)))
      (if hf
          (1- n)
        (first-non-null-homology-group-aux chcm (1+ n) limit)))))
  
  
(DEFUN FIRST-NON-NULL-HOMOLOGY-GROUP (chcm  limit)
  (first-non-null-homology-group-aux chcm 1 limit))


(DEFUN COMPUTE-HOMOTOPY-Z-XSLT (n-hom obj indx)
  (let* ((ch (chml-clss (eval obj) indx))
         (fib (z-whitehead (eval obj) (eval ch)))
         (ft (fibration-total (eval fib)))
         (result (homology-format  (eval ft) (1+ indx))))
    (if (= (1+ indx) n-hom)
        (homology-format (eval ft) n-hom)
      (if (string= "NIL" result) result
          (if (string= "Z " result)
              (compute-homotopy-z-xslt n-hom ft (1+ indx))
            (if (string= "Z/2Z " result)
                (compute-homotopy-z2-xslt n-hom ft (1+ indx))
              (if (and (string= (subseq result 0 2) "Z/") (string= (subseq result (search "Z" result :start2 2)) "Z"))
                (compute-homotopy-zp-xslt n-hom ft (1+ indx) (read-from-string (subseq result 2 (search "Z" result :start2 2))))
              (compute-homotopy-several-xslt n-hom ft (1+ indx) result))))))))


(DEFUN COMPUTE-HOMOTOPY-Z2-XSLT (n-hom obj indx)
  (let* ((ch (chml-clss (eval obj) indx))
         (fib (z2-whitehead (eval obj) (eval ch)))
         (ft (fibration-total (eval fib)))
         (result (homology-format (eval ft) (1+ indx))))
    (if (= (1+ indx) n-hom)
        (homology-format (eval ft) n-hom)
      (if (string= "NIL" result) result
          (if (string= "Z " result)
              (compute-homotopy-z-xslt n-hom ft (1+ indx))
            (if (string= "Z/2Z " result)
                (compute-homotopy-z2-xslt n-hom ft (1+ indx))
              (if (and (string= (subseq result 0 2) "Z/")  (string= (subseq result (search "Z" result :start2 2)) "Z"))
                (compute-homotopy-zp-xslt n-hom ft (1+ indx) (read-from-string (subseq result 2 (search "Z" result :start2 2))))
              (compute-homotopy-several-xslt n-hom ft (1+ indx) result))))))))


(DEFUN COMPUTE-HOMOTOPY-ZP-XSLT (n-hom obj indx n)
  (let* ((ch (chml-clss (eval obj) indx))
         (fib (zp-whitehead n (eval obj) (eval ch)))
         (ft (fibration-total (eval fib)))
         (result (homology-format  (eval ft) (1+ indx))))
    (if (= (1+ indx) n-hom)
        (homology-format (eval ft) n-hom)
      (if (string= "NIL" result) result
          (if (string= "Z " result)
              (compute-homotopy-z-xslt n-hom ft (1+ indx))
            (if (string= "Z/2Z " result)
                (compute-homotopy-z2-xslt n-hom ft (1+ indx))
              (if (and (string= (subseq result 0 2) "Z/")  (string= (subseq result (search "Z" result :start2 2)) "Z"))
                  (compute-homotopy-zp-xslt n-hom ft (1+ indx) (read-from-string (subseq result 2 (search "Z" result :start2 2))))
                (compute-homotopy-several-xslt n-hom ft (1+ indx) result))))))))


(DEFUN SPLIT-COMPONENTS (string)
  (let ((term (if (string= string "") "" 
                (if (string= string "NIL") NIL
                (subseq string 0 2)))))
    (if (string= term "Z ")
        (cons 1 (split-components (subseq string 2)))
      (if (string= term "Z/")
          (cons (read-from-string (subseq string 2 (search "Z" string :start2 2)))
                (split-components (subseq string (+ 2 (search "Z" string :start2 2)))))
        nil))))


(DEFUN CONSTRUCT-SPACE-ITERATIVE (chcm list indx)
  (if (endp list)
      chcm
    (cond ((equal (car list) 1) (let* ((ch (chml-clss chcm indx))
                                       (fib (z-whitehead chcm ch))
                                       (ft (fibration-total fib)))
                                  (construct-space-iterative ft (cdr list) indx)))
          ((equal (car list) 2) (let* ((ch (chml-clss chcm indx))
                                       (fib (z2-whitehead chcm ch))
                                       (ft (fibration-total fib)))
                                  (construct-space-iterative ft (cdr list) indx)))
          (t (let* ((ch (chml-clss chcm indx))
                    (fib (zp-whitehead (car list) chcm ch))
                    (ft (fibration-total fib)))
               (construct-space-iterative ft (cdr list) indx)))
          )))   


(DEFUN COMPUTE-HOMOTOPY-SEVERAL-XSLT (n-hom obj indx hom)
  (let* ((ft (construct-space-iterative obj (split-components hom) indx))
         (result (homology-format ft (1+ indx))))
    (if (= (1+ indx) n-hom)
        (homology-format (eval ft) n-hom)
      (if (string= "NIL" result) result
          (if (string= "Z " result)
              (compute-homotopy-z-xslt n-hom ft (1+ indx))
            (if (string= "Z/2Z " result)
                (compute-homotopy-z2-xslt n-hom ft (1+ indx))
              (if (and (string= (subseq result 0 2) "Z/")  (string= (subseq result (search "Z" result :start2 2)) "Z"))
                  (compute-homotopy-zp-xslt n-hom ft (1+ indx) (read-from-string (subseq result 2 (search "Z" result :start2 2))))
                (compute-homotopy-several-xslt n-hom ft (1+ indx) result))))))))


(DEFUN COMPUTE-HOMOTOPY2-XSLT (n-hom obj degree hom)
  (cond
   ((= n-hom 0)
   (format nil "Z")
   )
   ((< n-hom (1+ degree) )
    (format nil ""))
   ((= n-hom (1+ degree))
    (homology-format (eval obj) n-hom))
   ((string= hom "Z ")
    (compute-homotopy-z-xslt n-hom obj (1+ degree)))
   ((string= hom "Z/2Z ")
    (compute-homotopy-z2-xslt n-hom obj (1+ degree)))
   ((and (string= (subseq hom 0 2) "Z/") (string= (subseq hom (search "Z" hom :start2 2)) "Z"))
    (compute-homotopy-zp-xslt n-hom obj (1+ degree) (read-from-string (subseq hom 2 (search "Z" hom :start2 2)))))
   (t 
    (compute-homotopy-several-xslt n-hom obj (1+ degree) hom))))


(DEFUN COMPUTE-HOMOTOPY (smst n-hom degree)
  (let ((hom (homology-format smst (1+ degree))))
    (compute-homotopy2-xslt n-hom smst degree hom)))



(DEFUN HOMOTOPY (smst degr)
  (let (;; we obtain the first non null homology group
        (first-non-null (first-non-null-homology-group smst degr)))
    (if first-non-null
        (progn 
          (let ((result (split-components (compute-homotopy smst degr first-non-null))))
            (format t "~3%Homotopy in dimension ~D :~%" degr) 
            (dolist (item result)
              (format t "~2%Component Z")
              (unless (equal item 1) 
                (format t "/~DZ" item)))
            (terpri) (terpri)))
      (progn 
          (format t "~3%Homotopy in dimension ~D :~2%" degr)
        (terpri) (terpri)))))


(DEFUN HOMOTOPY-LIST (smst degr)
  (let (;; we obtain the first non null homology group
        (first-non-null (first-non-null-homology-group smst degr)))
    (if first-non-null
        (split-components (compute-homotopy smst degr first-non-null))
      nil)))

      
    


