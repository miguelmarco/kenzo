;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*



(IN-PACKAGE #:cat)

(PROVIDE "sage-interface")



;;;; Chain complexes



(DEFVAR *ASCII-NUM-START* 48)



(DEFVAR *ASCII-NUM-END* 57)



(DEFUN ITOA (value)
  (if (not (integerp value)) (error "Wrong type"))
  (if (zerop value)
      (return-from itoa "0"))
  (let ((negative? nil))
    (if (< value 0)
	(progn
	  (setf negative? t)
	  (setf value (abs value))))
    (loop
       for num = value then (floor num 10)
       for reminder = (mod num 10)
       while (> num 0)
       collect
	 (code-char (+ reminder *ASCII-NUM-START*)) into result
       finally
	 (if negative? (push #\- (cdr (last result))))
	 (return (format nil "~{~a~}" (nreverse result))))))



(DEFUN STRCAT (&rest strings) (apply 'concatenate 'string strings))



(DEFMACRO GEN (num1 num2)
  `(read-from-string (strcat "G" (itoa ,num1) "G" (itoa ,num2))))



(DEFUN KBASIS (schcm)
  (flet ((rslt (dim)
           (let ((pair (assoc dim schcm)))
             (if (null pair)
                 NIL
               (loop for k from 0 to (1- (array-dimension (cdr pair) 1))
                     collect (gen dim k))))))
    (the basis #'rslt)))



(DEFUN KDFFR (schcm)
  (flet ((frslt (dim gnr)
           (let* ((sep (1+ (position #\G (subseq (write-to-string gnr) 1))))
                  (dim_gnr (read-from-string (subseq (write-to-string gnr) 1 sep)))
                  (num_gnr (read-from-string (subseq (write-to-string gnr) (1+ sep)))))
               (declare (fixnum dim_gnr num_gnr))
             (unless (equal dim dim_gnr)      ;It's not working
               (error "WRONG DIMENSION!!!"))
             (let ((pair (assoc dim schcm))
                   (rslt (cmbn (1- dim))))
               (if (null pair)
                   rslt
                 (let ((dif (cdr pair)))
                   (unless (< -1 num_gnr (array-dimension dif 1))
                     (error "WRONG GENERATOR!!!"))
                   (setf (cmbn-list rslt) (loop for i from 0 to (1- (array-dimension dif 0))
                                                unless (eq (aref dif i num_gnr) 0)
                                                collect (term (aref dif i num_gnr) (gen (1- dim) i))))
                   rslt))))))
    #'frslt))



(DEFUN G-CMPR (GnGm1 GnGm2)
  (let* ((sep (1+ (position #\G (subseq (write-to-string GnGm1) 1))))
         (m1 (read-from-string (subseq (write-to-string GnGm1) (1+ sep))))
         (m2 (read-from-string (subseq (write-to-string GnGm2) (1+ sep)))))
    (declare (fixnum m1 m2))
    (the cmpr
         (if (< m1 m2)
             :less
           (if (= m1 m2)
               :equal
             :greater)))))



(DEFUN CONVERTMATRICE (matrice)
  (let* ((numfil (nlig matrice))
         (numcol (ncol matrice))
         (rslt (make-array (list numfil numcol) :initial-element 0)))
    (dotimes (i numfil)
      (dotimes (j numcol)
        (let ((Mij (entry matrice (1+ i) (1+ j))))
          (unless (zerop Mij)
            (setf (aref rslt i j) Mij)))))
    rslt))



(DEFUN ENTRY (mat i j)
  (let ((p (left (chercher-hor (baselig mat i) j))))
    (if (= j (icol p))
        (val p)
      0)))



(DEFUN CHCM-MAT2 (chcm n)
  (declare
   (type chain-complex chcm)
   (type fixnum n))
  (let ((sorc (basis chcm n))
        (trgt (basis chcm (1- n)))
        (cmpr (cmpr chcm)))
    (declare
     (list sorc trgt)
     (type cmprf cmpr))
    (let ((sorcl (length sorc))
          (mat (creer-matrice (length trgt) (length sorc)))
          (test #'(lambda (gnrt1 gnrt2)
                    (eq (funcall cmpr gnrt1 gnrt2) :equal))))
      (declare
       (type fixnum sorcl)
       (type matrice mat)
       (function test))
      (do ((i 1 (1+ i))
           (mark sorc (cdr mark)))
          ((endp mark))
        (declare
         (type fixnum i)
         (list mark))
          (maj-colonne mat i
             (mapcar #'(lambda (term)
                          (declare (type term term))
                          (list
                             (1+ (position (gnrt term) trgt :test test))
                             (cffc term)))
               (cmbn-list (? chcm n (car mark))))))
       mat)))



(DEFUN MAKE-ARRAY-TO-LISTS (array)
  (loop for i below (array-dimension array 0)
        collect (loop for j below (array-dimension array 1)
                      collect (aref array i j))))



(DEFUN MAKE-ARRAY-FROM-LISTS (nrows ncols list)
    (make-array (list nrows ncols) :initial-contents list))



(DEFUN BASIS_AUX1 (kchcm dim)
    (basis kchcm dim))



(DEFUN ORGN_AUX1 (kchcm)
    (orgn kchcm))



(DEFUN DFFR_AUX1(kchcm dim cmbn_list)
    (let ((comb (cmbn dim)))
        (setf (cmbn-list comb) cmbn_list)
        (dffr kchcm comb)))



;;;; Simplicial sets



(DEFUN BUILD-FINITE-SS2 (list)
  (declare (list list))
  (let ((bspn (first list))
        (table (finite-ss-table list))
        (ind-smst (gensym)))
    (declare
     (symbol bspn ind-smst)
     (simple-vector table))
    ;;  (vector (vector gmsm-faces-info))
    (let ((rslt (build-smst
                 :cmpr #'s-cmpr
                 :basis (finite-ss-basis table)
                 :bspn bspn
                 :face (finite-ss-face ind-smst table)
                 :intr-bndr (finite-ss-intr-bndr ind-smst table)
                 :bndr-strt :gnrt
                 :orgn `(build-finite-ss ,list))))
      (setf (symbol-value ind-smst) rslt)
      rslt)))




