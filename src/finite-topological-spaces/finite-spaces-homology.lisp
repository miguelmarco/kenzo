
;;  H-REGULAR-HOMOLOGY H-REGULAR-HOMOLOGY H-REGULAR-HOMOLOGY H-REGULAR-HOMOLOGY
;;  H-REGULAR-HOMOLOGY H-REGULAR-HOMOLOGY H-REGULAR-HOMOLOGY H-REGULAR-HOMOLOGY
;;  H-REGULAR-HOMOLOGY H-REGULAR-HOMOLOGY H-REGULAR-HOMOLOGY H-REGULAR-HOMOLOGY

(IN-PACKAGE #:cat)

(provide "finite-spaces-homology")

;;
;;  Computing homology of h-regular finite spaces
;;


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#



(DEFUN NEWSMITH-KERNEL (mtrx)
  (declare (type matrice mtrx))
  (the list
       (let* ((newsmith-smith-list (newsmith-smith mtrx))
              (m (third newsmith-smith-list))
              (q (fourth newsmith-smith-list))
              (min (min (nlig mtrx) (ncol mtrx))))
         (declare (type list newsmith-smith-list rslt)
                  (type matrice m q)
                  (type fixnum min))
         (let ((start (or (loop for ic from 1 to min
                                thereis (and (zerop (terme m ic ic)) ic))
                          (1+ min))))
           (declare (type fixnum start))
           (loop for j from start to (ncol mtrx)
                 collect (newsmith-extract-column q j))))))


#|
 (newsmith-kernel (convertarray #2A((1 0 -1)(-1 0 1))))
 (newsmith-kernel (convertarray #2A((1 0 1)(1 0 1))))
 (newsmith-kernel (convertarray #2A((0 -1 -1 0 -1)(0 1 0 -1 1)(0 0 1 1 0))))
 (newsmith-kernel (identite 6))
|#


(DEFMETHOD H-REGULAR-DIF ((stong matrice))
  (let ((differentials '())
        (edges (binarymatrice-to-ubasis (nilpot stong)))
        (heights (heights stong)))
    (declare (type list differentials))
    (push (creer-matrice 0 (length (car heights))) differentials)
    (do ((dim_Cn-2 NIL (length (car heights)))
         (heights heights (cdr heights)))
        ((endp (cdr heights)))
      (let* ((Cn-1 (car heights))
             (Cn (cadr heights))
             (dim_n (length Cn))
             (rslt (creer-matrice (length Cn-1) dim_n)))
        (if dim_Cn-2
            (push (let ((dn-1 (car differentials)))
                    (do ((group Cn (cdr group))
                         (ic 1 (1+ ic)))
                        ((> ic dim_n) rslt)
                      (let* ((Uic (svref edges (1- (car group))))
                             (iligs (posiciones-1 Uic Cn-1)) 
                             (ker (car (newsmith-kernel (submatrix-cols dn-1 iligs))))
                             (liste (mapcar #'(lambda (x) (list (nth (1- (car x)) iligs) (cadr x))) ker)))
                        (maj-colonne rslt ic liste))))
                  differentials)
          (push (do ((group Cn (cdr group))
                     (ic 1 (1+ ic)))
                    ((> ic dim_n) rslt)
                  (let* ((Uic (svref edges (1- (car group))))
                         (iligs (posiciones-1 Uic Cn-1))
                         (liste (list (list (first iligs) -1) (list (second iligs) 1))))
                    (maj-colonne rslt ic liste)))
                differentials))))
    (return-from H-REGULAR-DIF (reverse differentials))))


(DEFMETHOD H-REGULAR-DIF ((finspace finite-space))
  (H-REGULAR-DIF (stong finspace)))


(DEFMETHOD EXPLICIT-H-REGULAR-DIF ((finspace finite-space))
  #| Computes the list of groups and differentials ((C0 C1 ... CN) (D0 D1 ... Dn))
of the explicit-h-regular-chain complex of 'finspace' where n = height of 'finspace' |#
  (let ((groups NIL)
        (differentials NIL)
        (edges (binarymatrice-to-ubasis (nilpot (stong finspace))))
        (heights (heights finspace)))
    (declare (type list differentials groups))
    (dotimes (dimension (length heights))
      (let* ((Cn (nth dimension heights))
             (dim_n (length Cn)))
        (if (zerop dimension)
            (progn
              (push (creer-matrice 0 (length Cn)) differentials)
              (push (loop for x in Cn
                          collect (list (list 1 x))) groups))
          (let* ((Cn-1 (nth (1- dimension) heights))
                 (dim_n-1 (length Cn-1))
                 (rslt (creer-matrice dim_n-1 dim_n)))
            (case dimension          
              (1 (progn
                   (push (do ((group Cn (cdr group))
                              (ic 1 (1+ ic)))
                             ((> ic dim_n) rslt)
                           (let* ((Uic (svref edges (1- (car group))))
                                  (iligs (posiciones-1 Uic Cn-1))
                                  (liste (list (list (first iligs) -1) (list (second iligs) 1))))
                             (maj-colonne rslt ic liste)))
                         differentials)
                   (push (loop for x in Cn
                               collect (let ((edges_x (svref edges (1- x))))
                                         (list (list -1 (car edges_x) x) (list 1 (cadr edges_x) x)))) groups)))
              (otherwise (push (let ((Cn-2 (<a-b> 1 (length (nth (- dimension 2) heights))))
                                     (dn-1 (car differentials))
                                     (group Cn)
                                     (rslt2 NIL)
                                     (old-group (car groups)))
                                 
                                 (do ((ic 1 (1+ ic)))
                                     ((null group))
                                   (let* ((pseudo NIL)
                                          (initial (car group))
                                          (maxUj (svref edges (1- (pop group))))
                                          (iligs (posiciones-1 maxUj Cn-1))
                                          (ker (car (newsmith-kernel (submatrix-cols dn-1 iligs)))))
                                     (do ((ker ker (cdr ker)))
                                         ((endp ker))
                                       (let* ((ilig (caar ker))
                                              (updt (nth (1- (nth (1- ilig) iligs)) old-group)))
                                         (loop for z in (reverse updt)
                                               do (push (reverse (cons initial (reverse (cons (* (cadar ker) (car z)) (cdr z))))) pseudo))
                                       (insert-term rslt (nth (1- ilig) iligs) ic (cadar ker))))
                                     (push pseudo rslt2)))
                                 (push (reverse rslt2) groups)
                                 rslt)
                               differentials)))))))
        (return-from EXPLICIT-H-REGULAR-DIF (list (reverse groups) (reverse differentials)))))


#|
(h-regular-dif (bar-subdivision (random-finite-space 6 .5)))
(explicit-h-regular-dif (bar-subdivision (random-finite-space 6 .5)))

(dotimes (q 30)
  (setf finspace (bar-subdivision (random-finite-space 6 .5)))
  (setf hreg-list (h-regular-dif finspace))
  (setf exp-list (second (explicit-h-regular-dif finspace)))
  (dotimes (k (length hreg-list))
    (print (equalmatrix (nth k exp-list) (nth k hreg-list)))))
|#


(DEFMETHOD H-REGULAR-HOMOLOGY ((finspace finite-space) &rest rest)
  #| Compute the homology of 'finspace' directly from the list (H-REGULAR-DIF 'finspace') |#
  #| rest must be: NIL, an integer number or two integer numbers |#
  (let ((size (1- (length (heights finspace))))
        (range '()))
    (case (length rest)
      (0 (setf range (<a-b> 0 size)))
      (1 (if (> (car rest) size)
             (return-from H-REGULAR-HOMOLOGY (progn (format stream "~3%Homology in dimension ~D :~%" (car rest))
                                               (terpri) (terpri)
                                               (done)))
           (setf range rest)))
      (2 (if (> (car rest) size)
             (return-from H-REGULAR-HOMOLOGY (progn (format stream "~3%Homology in dimensions in ~D :~%" (<a-b> (first rest) (second rest)))
                                               (terpri) (terpri)
                                               (done)))
           (setf range (<a-b> (first rest) (second rest))))))
    (let ((list (matrices (save-info-aux finspace))))
      (dolist (dim range)
        (if (> dim size)
            (progn (format stream "~3%Homology in dimension ~D :~%" dim)
              (terpri) (terpri)
              (done))
          (let ((Mn (copier-matrice (nth dim list)))
                (Mn+1 NIL))
            (declare (type matrice Mn Mn+1))
            (if (= dim size)
                (setf Mn+1 (creer-matrice (ncol (nth dim list)) 0))
              (setf Mn+1 (copier-matrice (nth (1+ dim) list))))
            (let ((rsl (homologie Mn Mn+1)))
              (declare (type list rsl))    
              (format stream "~3%Homology in dimension ~D :~%" dim)
              (dolist (item rsl)
                (declare (type list item))
                (format stream "~2%Component Z")
                (unless (zerop (first item)) 
                  (format stream "/~DZ" (first item))))
              (terpri) (terpri)
              (done))))))))
              

#|
 (setf example (bar-subdivision (random-finite-space 6 .5)))
 (h-regular-homology example 1)
|#
   

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


(DEFUN STRCAT (&rest strings)
  (apply 'concatenate 'string strings))


(DEFMACRO WITH-X (num)
  `(read-from-string (strcat "X" (itoa ,num))))


(DEFVAR *H-REGULAR-DIF-HASH-TABLE*)
(setq *H-REGULAR-DIF-HASH-TABLE* (make-hash-table))


(DEFSTRUCT (SAVE-INFO (:conc-name nil))
  groups matrices)


(DEFUN SAVE-INFO-AUX (finspace)
  #| Assign to (idnm 'finspace') the info of (h-regular-dif 'finspace') |#
  (or (gethash (idnm finspace) *H-REGULAR-DIF-HASH-TABLE*)
      (setf (gethash (idnm finspace) *H-REGULAR-DIF-HASH-TABLE*)
            (make-save-info
             :matrices (h-regular-dif finspace)))))


(DEFUN X-CMPR (Xn1 Xn2)
  (let ((n1 (read-from-string (subseq (write-to-string Xn1) 1)))
        (n2 (read-from-string (subseq (write-to-string Xn2) 1))))
   (declare (fixnum n1 n2))
   (the cmpr
      (if (< n1 n2)
         :less
         (if (= n1 n2)
            :equal
            :greater)))))


(DEFUN TOP-BASIS (finspace)
  #| Provide the slot :basis of (CHCM-H-REGULAR 'finspace') |#
  (flet ((rslt (dim)
           (let ((heights (heights finspace)))
             (cond ((= dim -1) (list '*))
                   ((< dim -1) (error "Wrong dimension!"))
                   ((< dim (length heights)) (let ((refer (nth dim heights)))
                                               (loop for k from 0 to (1- (length refer))
                                                     collect (with-x (nth k refer)))))
                   (T +empty-list+)))))
    (the basis #'rslt)))


(DEFUN TOP-DIFF (finspace)
  #| Provide the slot :intr-dffr of (CHCM-H-REGULAR 'finspace') |#
  (flet ((frslt (dim gnr)
           (let ((heights (heights finspace)))
             (unless (<= 0 dim (length heights))
               (error "Wrong dimension!"))
             (let ((num_gnr 0))
               (setf num_gnr (read-from-string (subseq (write-to-string gnr) 1))) ; num_gnr of 'X13' is 13
               (let ((j (position num_gnr (nth dim heights))))  
                 (if (null j)
                     (error "Wrong generator!")
                   (if (zerop dim)
                       (cmbn -1)
                     (let* ((rslt (cmbn (1- dim)))
                            (info (gethash (idnm finspace) *H-REGULAR-DIF-HASH-TABLE*))
                            (diff (nth dim (matrices info))))
                       (setf (cmbn-list rslt) (loop for i from 1 to (length (nth (1- dim) heights))
                                                    unless (eq (terme diff i (1+ j)) 0)
                                                    collect (term (terme diff i (1+ j)) (with-x (nth (1- i) (nth (1- dim) heights))))))
                       rslt))))))))
    #'frslt))
         

(DEFMETHOD CHCM-H-REGULAR ((finspace finite-space))
  (save-info-aux finspace)
  (build-chcm :cmpr #'x-cmpr
              :strt :GNRT
              :basis (top-basis finspace)
              :intr-dffr (top-diff finspace)
              :orgn `(CHCM-H-REGULAR ,finspace)))


(DEFUN H-REGULAR-HOMOLOGY-GENERATORS (finspace dim)
  (let ((chcm (chcm-h-regular finspace)))
    (chcm-homology-gen chcm dim)))


#| 
  (setf finspace (bar-subdivision (random-finite-space 7 .6)))
  (setf chcm (chcm-h-regular finspace))
  (basis chcm -1)
  (basis chcm 0)
  (basis chcm 1)
  (basis chcm 2)
  (dffr chcm 2 (first (basis chcm 2)))
  (chcm-homology chcm 0)
  (chcm-homology chcm 1)
  (chcm-homology chcm 2)
  (chcm-homology chcm 24)
  (h-regular-homology-generators finspace 1)
|#


(DEFUN SAVE-INFO-EXPLICIT-AUX (finspace)
  #| Assign to (idnm 'finspace') the info of (h-regular-dif-explicit 'finspace') |#
  (let ((info (gethash (idnm finspace) *H-REGULAR-DIF-HASH-TABLE*)))
    (if info
        (unless (groups info)
          (let ((groups-matrices (explicit-h-regular-dif finspace)))
            (setf (groups (gethash (idnm finspace) *H-REGULAR-DIF-HASH-TABLE*))
                  (first groups-matrices))))
      (let ((groups-matrices (explicit-h-regular-dif finspace)))
        (setf (gethash (idnm finspace) *H-REGULAR-DIF-HASH-TABLE*)
              (make-save-info
               :groups (first groups-matrices)
               :matrices (second groups-matrices)))))))


(DEFUN EXPLICIT-CMPR (gnr1 gnr2)
  (let ((n1 (car (last (first (cmbn-list gnr1)))))
        (n2 (car (last (first (cmbn-list gnr2))))))
   (declare (fixnum n1 n2))
   (the cmpr
      (if (< n1 n2)
         :less
         (if (= n1 n2)
            :equal
            :greater)))))


(DEFUN EXPLICIT-TOP-BASIS (finspace)
  #| Provide the slot :basis of (EXPLICIT-CHCM-H-REGULAR 'finspace') |#
  (flet ((rslt (dim)
           (cond ((= dim -1) (list '*))
                 ((< dim -1) (error "Wrong dimension!"))
                 ((< dim (length (heights finspace)))
                  (let ((groups (groups (gethash (idnm finspace) *H-REGULAR-DIF-HASH-TABLE*))))
                    (mapcar #'(lambda (x)
                                (let ((rslt2 (cmbn dim)))
                                  (setf (cmbn-list rslt2) x)
                                  rslt2)) (nth dim groups))))
                  (T +empty-list+))))
    (the basis #'rslt)))


(DEFUN EXPLICIT-TOP-DIFF (finspace)
  #| Provide the slot :intr-dffr of (EXPLICIT-CHCM-H-REGULAR 'finspace') |#
  (flet ((frslt (dim gnr)
           (let ((heights (heights finspace)))
             (unless (<= 0 dim (length heights))
               (error "Wrong dimension!"))
             (let ((j (position (car (last (first (cmbn-list gnr)))) (nth dim heights))))
               (if (null j)
                   (error "Wrong generator!")
                 (if (zerop dim)
                     (cmbn -1)
                   (let* ((rslt (cmbn (1- dim)))
                          (info (gethash (idnm finspace) *H-REGULAR-DIF-HASH-TABLE*))
                          (dif (nth dim (matrices info)))
                          (base (basis (EXPLICIT-CHCM-H-REGULAR finspace) (1- dim))))
                     (setf (cmbn-list rslt) (loop for i from 1 to (length (nth (1- dim) heights))
                                                  unless (eq (terme dif i (1+ j)) 0)
                                                  collect (term (terme dif i (1+ j)) (nth (1- i) base))))
                     rslt)))))))
    #'frslt))


(DEFMETHOD EXPLICIT-CHCM-H-REGULAR ((finspace finite-space))
  (save-info-explicit-aux finspace)
  (build-chcm :cmpr #'explicit-cmpr
              :strt :GNRT
              :basis (explicit-top-basis finspace)
              :intr-dffr (explicit-top-diff finspace)
              :orgn `(EXPLICIT-CHCM-H-REGULAR ,finspace)))


(DEFUN EXPLICIT-H-REGULAR-HOMOLOGY-GENERATORS (finspace dim)
  (let ((chcm (explicit-chcm-h-regular finspace)))
    (chcm-homology-gen chcm dim)))

#|
  (setf finspace (bar-subdivision (random-finite-space 7 .6)))
  (setf chcm (explicit-chcm-h-regular finspace))
  (basis chcm -1)
  (basis chcm 0)
  (basis chcm 1)
  (basis chcm 2)
  (dffr chcm 2 (first (basis chcm 2)))
  (chcm-homology chcm 0)
  (chcm-homology chcm 1)
  (chcm-homology chcm 2)
  (chcm-homology chcm 24)
  (explicit-h-regular-homology-generators finspace 1)
|#
