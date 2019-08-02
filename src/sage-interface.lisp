;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*



(IN-PACKAGE #:cat)

(PROVIDE "sage-interface")



;;;;; Chain complexes ;;;;;


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
#| Provide the slot :basis of the KenzoSimplicialSet obtained by applying the function 'KChainComplex' to 'schcm' |#
  (flet ((rslt (dim)
           (let ((pair (assoc dim schcm)))
             (if (null pair)
                 NIL
               (loop for k from 0 to (1- (array-dimension (cdr pair) 1))
                     collect (gen dim k))))))
    (the basis #'rslt)))


(DEFUN KDFFR (schcm)
#| Provide the slot :intr-dffr of the KenzoSimplicialSet obtained by applying the function 'KChainComplex' to 'schcm' |#
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
#| Provide the slot :cmpr of the KenzoSimplicialSet obtained by applying the function 'KChainComplex' to 'schcm' |#
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


(DEFUN ENTRY (mat i j)
#| Provide the entry in row 'i' and column 'j' of the matrice 'mat' |#
  (let ((p (left (chercher-hor (baselig mat i) j))))
    (if (= j (icol p))
        (val p)
      0)))


(DEFUN CONVERTMATRICE (matrice)
#| Convert a matrix of type matrice to an array |#
  (let* ((numfil (nlig matrice))
         (numcol (ncol matrice))
         (rslt (make-array (list numfil numcol) :element-type 'fixnum :initial-element 0)))
    (dotimes (i numfil)
      (dotimes (j numcol)
        (let ((Mij (entry matrice (1+ i) (1+ j))))
          (unless (zerop Mij)
            (setf (aref rslt i j) Mij)))))
    rslt))


(DEFUN CHCM-MAT2 (chcm n)
#| The same 'chcm-mat' function but it does not print on screen intermediate computations |#
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
#| Convert an array to a list formed by the array's rows |#
  (loop for i below (array-dimension array 0)
        collect (loop for j below (array-dimension array 1)
                      collect (aref array i j))))


(DEFUN MAKE-ARRAY-FROM-LISTS (nrows ncols list)
#| Construct an array from a list formed by the array's rows |#
    (make-array (list nrows ncols) :initial-contents list))


(DEFUN BASIS_AUX1 (kchcm dim)
    (basis kchcm dim))


(DEFUN ORGN_AUX1 (kchcm)
    (orgn kchcm))


(DEFUN DFFR_AUX1(kchcm dim cmbn_list)
#| Differential of a comb with dimension 'dim' and 'cmbn_list' as its cmbn-list |#
    (let ((comb (cmbn dim)))
        (setf (cmbn-list comb) cmbn_list)
        (dffr kchcm comb)))


(DEFUN KCHAINCOMPLEX_AUX1 (chcm str_orgn)
#| Construct a chain complex in Kenzo from the information of the assoc list 'chcm' constructed from a dictionary whose values are matrices |#
  (build-chcm :cmpr #'G-CMPR
  :strt :GNRT
  :basis (kbasis chcm)
  :intr-dffr (kdffr chcm)
  :orgn `(KChainComplex ,str_orgn)))


;;;;; Simplicial Sets ;;;;;


(DEFUN KABSTRACTSIMPLEX_AUX1 (degeneracies name)
  (absm (dgop-ext-int degeneracies) name))


(DEFUN BUILD-FINITE-SS2 (list)
#| The same 'build-finite-ss' function but it does not check the face relations of the simplices in 'list' ('check-smst' omitted) |#
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


(DEFUN SFINITESIMPLICIALSET_AUX1 (finitess limit)
#| Construct a list (L0 L1 ... Lm), with m <= 'limit', where each list Lk is formed by lists of the form (S (a0 ... ak)), where each aj is a face of the k-simplex S |#
  (let ((rslt NIL))
    (do ((k 1 (1+ k))) ((> k limit))
      (let ((dimk NIL))
        (dolist (simplex (basis finitess k))
          (let ((faces NIL))
            (do ((i k (1- i))) ((< i 0))
              (push (face finitess i k simplex) faces))
            (push (cons simplex faces) dimk)))
      (push dimk rslt)))
    (reverse rslt)))


