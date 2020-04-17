;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*

;;;  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY
;;;  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY
;;;  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY  HOMOTOPY

;;; Authors: Jonathan Heras and Ana Romero


(IN-PACKAGE #:cat)

(provide "homotopy")

(DEFVAR *test-comments* nil)

(DEFUN CHML-CLSS-INTR-NOT-1REDUCED (chcm first)
  (declare
   (type chain-complex chcm)
   (fixnum first))
  (when *test-comments*
    (print "Function CHML-CLSS-INTR-NOT-1REDUCED called"))
  (let* ((echcm (echcm chcm))
         (cmpr (cmpr echcm))
         (basis (basis echcm))
         (f-basis (funcall basis first))
         (mtrx1 (gnrt-list-to-mtrx (kernel (chcm-mtrx echcm first))))
         (mtrx2 (chcm-mtrx echcm (1+ first)))         
         (smith-list-1 (smith mtrx1))
         (p1 (first smith-list-1))
         (p1-1 (second smith-list-1))
         (d1 (third smith-list-1))
         (p1-1Xmtrx2 (mtrx-prdc p1-1 mtrx2))
         (p1Xd1 (copy-mtrx p1))
         (d1-1Xp1-1Xmtrx2 (copy-mtrx p1-1Xmtrx2)))
    (declare
     (list smith-list-1)
     (type matrix p1 p1-1 d1 p1-1Xmtrx2 p1Xd1 d1-1Xp1-1Xmtrx2))
    (progn
      (let* ((line-n-1 (line-number p1))
             (column-n-1 (column-number p1)))
        (declare (type fixnum line-n-1 column-n-1))
        (dotimes (ic column-n-1)
          (declare (type fixnum ic))
          (if (and (< ic (line-number d1)) (< ic (column-number d1)))
              (let ((d (aref d1 ic ic)))
                (declare (type fixnum d))
                (dotimes (il line-n-1)
                  (declare (type fixnum il))
                  (if (and (not (eql 1 d)) (not (eql 0 d)))
                      (setf (aref p1Xd1 il ic) (* d (aref p1Xd1 il ic)))))))))
      (let* ((line-n-2 (line-number p1-1Xmtrx2))
             (column-n-2 (column-number p1-1Xmtrx2)))
        (declare (type fixnum line-n-2 column-n-2))    
        (dotimes (il line-n-2)
          (declare (type fixnum il))
          (if (and (< il (line-number d1)) (< il (column-number d1)))
              (let ((d (aref d1 il il)))
                (if (and (not (eql 1 d)) (not (eql 0 d)))
                    (dotimes (ic column-n-2)
                      (declare (type fixnum ic))
                      (setf (aref d1-1Xp1-1Xmtrx2 il ic) 
                        (floor (aref d1-1Xp1-1Xmtrx2 il ic) d))))))))
      (let* ((mtrx-list (smith d1-1Xp1-1Xmtrx2))
             (smith (third mtrx-list))
             (p-1 (mtrx-prdc  (second mtrx-list) p1-1))
             
             (n (line-number smith))
             (m (column-number smith))
             (diag-indx 
              (dotimes (indx (min n m)
                             (if (> n m)
                                 m
                               (error "In CHML-CLSS, the cohomology-ring ~@
                                      is null.")))
                (declare (fixnum indx))
                (unless (= 1 (aref smith indx indx))
                  (return indx)))))
        (declare
         (type chain-complex echcm)
         (type cmprf cmpr)
         (type basis basis)
         (fixnum n m diag-indx)
         (list f-basis mtrx-list)
         (type matrix p-1 smith))
        (flet ((rslt (cmbn)
                     (declare (type cmbn cmbn))
                     (with-cmbn
                         (degr list) cmbn
                       (unless (= degr first)
                         (return-from rslt (zero-cmbn (- degr first))))
                       (do ((rslt 0)
                            (bmark f-basis)
                            (ic 0)
                            (cmark list (cdr cmark)))
                           ((endp cmark)
                            (if (zerop rslt)
                                (zero-cmbn 0)
                              (term-cmbn 0 rslt :z-gnrt)))
                         (declare
                          (fixnum rslt)
                          (list bmark cmark))
                         (with--term
                             (cffc gnrt) cmark
                           (loop
                             (when (eq :equal (funcall cmpr gnrt (car bmark)))
                               (return))
                             (pop bmark)
                             (incf ic))
                           (incf rslt (* cffc (aref p-1 diag-indx ic)))
                           (pop bmark)
                           (incf ic))))))
          (the intr-mrph #'rslt))))))


(DEFUN CHML-CLSS-NOT-1REDUCED (chcm first)
  (declare
   (type chain-complex chcm)
   (fixnum first))
  (when *test-comments*
    (print "Function CHML-CLSS-NOT-1REDUCED called"))
  (the morphism
    (build-mrph
     :sorc (echcm chcm) :trgt (z-chcm) :degr (- first)
     :intr (chml-clss-intr-not-1reduced chcm first)
     :strt :cmbn
     :orgn `(chml-clss-not-1reduced ,chcm ,first))))


(DEFUN LIST-TO-REDUCEDSMITHLIST (l)
  (declare (type list l))
  (let ((b 't)
        (n (length l)))
    (declare (fixnum n))
    (progn 
      (mapcar #'(lambda (i)
                  (if (and (eq 0  (nth (1- i) l)) (not (eq 0 (nth i l))))
                      (setf b nil)
                    (if (and b (not (eq 0 (nth (1- i) l))))
                        (if (not (eq (mod (nth i l) (nth (1- i) l)) 0)) 
                            (setf b nil)))))
        (<a-b< 1 n))
      (if b l 
        (let* ((n (length l))
               (mtrx   (make-array (list n n)
                                   :element-type 'fixnum
                                   :initial-element 0)))
          (declare (type matrix mtrx)
                   (fixnum n))
          (progn
            (dotimes (i n)
              (declare (fixnum i))
              (setf (aref mtrx i i) (nth i l)))
            (smith mtrx)
            (let ((rslt nil))
              (dotimes (i n)
                (let ((ni (aref mtrx i i)))
                  (if (not (eq 1 ni))
                      (push ni rslt))))
              (nreverse rslt))))))))              


(DEFUN CHCM-HOMOLOGY-FORMAT (cc n)
  (declare (type chain-complex cc) (type fixnum n))
  (when *test-comments*
    (print "Function CHCM-HOMOLOGY-FORMAT called"))
  (let* ((hom (homologie (chcm-mat cc n) (chcm-mat cc (1+ n))))
         (l (mapcar #'(lambda (x)
                        (first x))
              hom))
         (rsl (list-to-reducedsmithlist l))
         (str nil))
    (declare (type list hom l rsl))
    (dolist (item rsl)
      (declare (type fixnum item))
      (setf str (concatenate 'string str (format nil "Z")))
      (unless (zerop item) 
        (setf str (concatenate 'string str (format nil "/~DZ" item))))
      (setf str (concatenate 'string str " ")))
    str))


(DEFUN HOMOLOGY-FORMAT (chcm degr1 &optional (degr2 (1+ degr1)))
  (declare (fixnum degr1 degr2))
  (when *test-comments*
    (print "Function HOMOLOGY-FORMAT called"))
  (let ((result_hom nil))
    (do ((degr degr1 (1+ degr)))
        ((>= degr degr2))
      (declare (fixnum degr))
      (setf result_hom (chcm-homology-format (echcm chcm) degr))
      (when *homology-verbose*
        (terpri) (clock) (terpri)))
    result_hom))


(DEFUN FIRST-NON-NULL-HOMOLOGY-GROUP-AUX (chcm n limit)
  (when *test-comments*
    (print "Function FIRST-NON-NULL-HOMOLOGY-GROUP-AUX called"))
  (if (> n limit)
      nil
    (let ((hf (homology-format chcm n)))
      (if hf
          (1- n)
        (first-non-null-homology-group-aux chcm (1+ n) limit)))))


(DEFUN FIRST-NON-NULL-HOMOLOGY-GROUP (chcm  limit)
  (when *test-comments*
    (print "Function FIRST-NON-NULL-HOMOLOGY-GROUP called"))
  (first-non-null-homology-group-aux chcm 1 limit))


(DEFUN COMPUTE-HOMOTOPY-Z-XSLT (n-hom obj indx)
  (when *test-comments*
    (print "Function COMPUTE-HOMOTOPY-Z-XSLT called"))
  (let* ((b (= 0 (length (basis (echcm obj) 1))))
         (ch (if b (chml-clss (eval obj) indx)
               (chml-clss-not-1reduced (eval obj) indx)))
         (fib (z-whitehead (eval obj) (eval ch)))
         (ft (fibration-total (eval fib)))
         (result (homology-format  (eval ft) (1+ indx))))
    (if (= (1+ indx) n-hom)
        (homology-format (eval ft) n-hom)
      (progn
        (if b
            (kill-epis ft (1- indx) (1+ indx)))
        (if (string= "NIL" result) ;; NIL
            (let ((first-non-null (first-non-null-homology-group ft n-hom)))
              (if first-non-null (compute-homotopy ft n-hom first-non-null) nil))
          (if (string= "Z " result)              
              (compute-homotopy-z-xslt n-hom ft (1+ indx))
            (if (string= "Z/2Z " result)
                
                (compute-homotopy-z2-xslt n-hom ft (1+ indx))
              (if (and (string= (subseq result 0 2) "Z/") (string= (subseq result (search "Z" result :start2 2)) "Z"))
                  (compute-homotopy-zp-xslt n-hom ft (1+ indx) (read-from-string (subseq result 2 (search "Z" result :start2 2))))
                (compute-homotopy-several-xslt n-hom ft (1+ indx) result)))))))))


(DEFUN COMPUTE-HOMOTOPY-Z2-XSLT (n-hom obj indx)
  (when *test-comments*
    (print "Function COMPUTE-HOMOTOPY-Z2-XSLT called"))
  (let* ((b (= 0 (length (basis (echcm obj) 1))))
         (ch (if b (chml-clss (eval obj) indx)
               (chml-clss-not-1reduced (eval obj) indx)))
         (fib (z2-whitehead (eval obj) (eval ch)))
         (ft (fibration-total (eval fib)))
         (result (homology-format (eval ft) (1+ indx))))
    (if (= (1+ indx) n-hom)
        (homology-format (eval ft) n-hom)
      (progn
        (if b
            (kill-epis ft (1- indx) (1+ indx)))
        (if (string= "NIL" result) ;; nil
            (let ((first-non-null (first-non-null-homology-group ft n-hom)))
              (if first-non-null (compute-homotopy ft n-hom first-non-null) nil))
          (if (string= "Z " result)
              (compute-homotopy-z-xslt n-hom ft (1+ indx))
            (if (string= "Z/2Z " result)                
                (compute-homotopy-z2-xslt n-hom ft (1+ indx))
              (if (and (string= (subseq result 0 2) "Z/")  (string= (subseq result (search "Z" result :start2 2)) "Z"))                  
                  (compute-homotopy-zp-xslt n-hom ft (1+ indx) (read-from-string (subseq result 2 (search "Z" result :start2 2))))
                (compute-homotopy-several-xslt n-hom ft (1+ indx) result)))))))))


(DEFUN COMPUTE-HOMOTOPY-ZP-XSLT (n-hom obj indx n)
  (when *test-comments*
    (print "Function COMPUTE-HOMOTOPY-ZP-XSLT called"))
  (let* ((b (= 0 (length (basis (echcm obj) 1))))
         (ch (if b (chml-clss (eval obj) indx)
               (chml-clss-not-1reduced (eval obj) indx)))
         (fib (zp-whitehead n (eval obj) (eval ch)))
         (ft (fibration-total (eval fib)))
         (result (homology-format  (eval ft) (1+ indx))))
    (if (= (1+ indx) n-hom)
        (homology-format (eval ft) n-hom)
      (progn
        (if b
            (kill-epis ft (1- indx) (1+ indx)))
        (if (string= "NIL" result) ;; nil
            (let ((first-non-null (first-non-null-homology-group ft n-hom)))
              (if first-non-null (compute-homotopy ft n-hom first-non-null) nil))
          (if (string= "Z " result)
              (compute-homotopy-z-xslt n-hom ft (1+ indx))
            (if (string= "Z/2Z " result)              
                (compute-homotopy-z2-xslt n-hom ft (1+ indx))
              (if (and (string= (subseq result 0 2) "Z/")  (string= (subseq result (search "Z" result :start2 2)) "Z"))                  
                  (compute-homotopy-zp-xslt n-hom ft (1+ indx) (read-from-string (subseq result 2 (search "Z" result :start2 2))))
                (compute-homotopy-several-xslt n-hom ft (1+ indx) result)))))))))


(DEFUN SPLIT-COMPONENTS (string)
  (when *test-comments*
    (print "Function SPLIT-COMPONENTS called"))
  (let ((term (if (string= string "") "" 
                (if (string= string "NIL") nil
                  (subseq string 0 2)))))
    (if (string= term "Z ")
        (cons 1 (split-components (subseq string 2)))
      (if (string= term "Z/")
          (cons (read-from-string (subseq string 2 (search "Z" string :start2 2)))
                (split-components (subseq string (+ 2 (search "Z" string :start2 2)))))
        nil))))


(DEFUN CONSTRUCT-SPACE-ITERATIVE (chcm list indx)
  (when *test-comments*
    (print "Function CONSTRUCT-SPACE-ITERATIVE called"))
  (if (endp list)
      chcm
    (let ((b (= 0 (length (basis (echcm chcm) 1)))))
      (cond ((equal (car list) 1) (let* ((ch (if b (chml-clss chcm indx)
                                               (chml-clss-not-1reduced chcm indx)))
                                         (fib (z-whitehead chcm ch))
                                         (ft (fibration-total fib)))
                                    (if (endp (cdr list)) 
                                        (progn
                                          (if b
                                              (kill-epis ft (1- indx) (1+ indx)))
                                          ft)
                                      (progn
                                        (if b
                                            (kill-epi ft (1- indx)))
                                        (construct-space-iterative ft (cdr list) indx)))))
            ((equal (car list) 2) (let* ((ch (if (= 0 (length (basis (echcm chcm) 1))) (chml-clss chcm indx)
                                               (chml-clss-not-1reduced chcm indx)))
                                         (fib (z2-whitehead chcm ch))
                                         (ft (fibration-total fib)))
                                    (if (endp (cdr list)) 
                                        (progn
                                          (if b
                                              (kill-epis ft (1- indx) (1+ indx)))
                                          ft)
                                      (progn
                                        (if b
                                            (kill-epi ft (1- indx)))
                                        (construct-space-iterative ft (cdr list) indx)))))
            (t (let* ((ch (if (= 0 (length (basis (echcm chcm) 1))) (chml-clss chcm indx)
                            (chml-clss-not-1reduced chcm indx)))
                      (fib (zp-whitehead (car list) chcm ch))
                      (ft (fibration-total fib)))
                 (if (endp (cdr list)) 
                     (progn
                       (if b
                           (kill-epis ft (1- indx) (1+ indx)))
                       ft)
                   (progn
                     (if b
                         (kill-epi ft (1- indx)))              
                     (construct-space-iterative ft (cdr list) indx)))))))) ) 


(DEFUN COMPUTE-HOMOTOPY-SEVERAL-XSLT (n-hom obj indx hom)
  (when *test-comments*
    (print "Function COMPUTE-HOMOTOPY-SEVERAL-XSLT called"))
  (let* ((ft (construct-space-iterative obj (split-components hom) indx))
         (result (homology-format ft (1+ indx))))
    (if (= (1+ indx) n-hom)
        (homology-format (eval ft) n-hom)
      (if (string= "NIL" result) nil
        (if (string= "Z " result)
            (compute-homotopy-z-xslt n-hom ft (1+ indx))
          (if (string= "Z/2Z " result)
              (compute-homotopy-z2-xslt n-hom ft (1+ indx))
            (if (and (string= (subseq result 0 2) "Z/")  (string= (subseq result (search "Z" result :start2 2)) "Z"))
                (compute-homotopy-zp-xslt n-hom ft (1+ indx) (read-from-string (subseq result 2 (search "Z" result :start2 2))))
              (compute-homotopy-several-xslt n-hom ft (1+ indx) result))))))))


(DEFUN COMPUTE-HOMOTOPY2-XSLT (n-hom obj degree hom)
  (when *test-comments*
    (print "Function COMPUTE-HOMOTOPY2-XSLT called"))
  (cond
   ((= n-hom 0)
    (format nil "Z"))
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
  (when *test-comments*
    (print "Function COMPUTE-HOMOTOPY called"))
  (let ((hom (homology-format smst (1+ degree))))
    (compute-homotopy2-xslt n-hom smst degree hom)))


(DEFUN HOMOTOPY (smst degr)
  (when *test-comments*
    (print "Function HOMOTOPY called"))
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
  (when *test-comments*
    (print "Function HOMOTOPY-LIST called"))
  (let (;; we obtain the first non null homology group
        (first-non-null (first-non-null-homology-group smst degr)))
    (if first-non-null
        (split-components (compute-homotopy smst degr first-non-null))
      nil)))


