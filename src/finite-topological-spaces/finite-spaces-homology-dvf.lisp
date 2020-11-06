
;;  H-REGULAR-DVF H-REGULAR-DVF H-REGULAR-DVF H-REGULAR-DVF H-REGULAR-DVF 
;;  H-REGULAR-DVF H-REGULAR-DVF H-REGULAR-DVF H-REGULAR-DVF H-REGULAR-DVF 
;;  H-REGULAR-DVF H-REGULAR-DVF H-REGULAR-DVF H-REGULAR-DVF H-REGULAR-DVF  

(IN-PACKAGE #:cat)

(provide "finite-spaces-homology-dvf")

;;
;;  Computing homology of h-regular finite spaces using discrete vector fields
;;


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFUN POSICIONES (pos list)
  (loop for x in pos
        collect (nth (1- x) list)))


(DEFUN POSICIONES-1 (list1 list2)
 #| 'list1' must be a sublist of 'list2' |#
  (loop for x in list1
        collect (1+ (position x list2))))

#|
 (posiciones '(1 3 5) (<a-b> 4 10))
 (posiciones-1 '(4 7 9) (<a-b> 2 10))
|#


(DEFUN MTRX-SUBSTRACT (mtrx1 mtrx2)
  (let ((rslt (creer-matrice (nlig mtrx1) (ncol mtrx1))))
    (do ((i 1 (1+ i)))
        ((> i (nlig rslt)) rslt)
      (let ((ptl1 (baselig mtrx1 i))
            (ptl2 (baselig mtrx2 i))
            (pl1 (left (baselig mtrx1 i)))
            (pl2 (left (baselig mtrx2 i)))
            (plr (left (baselig rslt i))))
        (do () ((and (eq pl1 ptl1) (eq pl2 ptl2)))
          (let ((ic1 (icol pl1))
                (ic2 (icol pl2)))
            (cond
             ((< ic1 ic2) (progn (majterme (chercher-hor plr ic1) (basecol rslt ic2) (- (val pl2)))
                            (setf pl2 (left pl2))))
             ((< ic2 ic1) (progn (majterme (chercher-hor plr ic2) (basecol rslt ic1) (val pl1))
                            (setf pl1 (left pl1))))
             (T (progn (majterme (chercher-hor plr ic1) (basecol rslt ic1) (- (val pl1) (val pl2)))
                  (setf pl2 (left pl2)
                        pl1 (left pl1)))))))))))

#|
  (setf mat1 (mat-aleat 300 300 .05 8))
  (setf mat2 (mat-aleat 300 300 .05 8))
  (setf res1 (mtrx-substract mat1 mat2))
|#


(DEFUN D33-D31D21D23 (4blocks size-d21)
  #| From '4blocks' computes d21^-1 and returns d33-d31(d21^-1)d23 |#
  (let* ((nrows (nlig 4blocks))
         (ncols (ncol 4blocks))
         (d31-mtrx (creer-matrice (- nrows size-d21) size-d21))
         (d23-mtrx (creer-matrice size-d21 (- ncols size-d21)))
         (d33-mtrx (creer-matrice (- nrows size-d21) (- ncols size-d21))))
    (do ((index 1 (1+ index)))
        ((> index size-d21))
      #| change the -1's in the diagonal of d21 by 1's |#
      (let ((diag-term (terme 4blocks index index)))
        (unless (eq diag-term 1)
          (newsmith-column-minus 4blocks index)))
      (do ((pl01 (baselig 4blocks index))
           (p1 (left (baselig 4blocks index)) (left p1)))
          ((eq p1 pl01))
        (if (> (icol p1) size-d21)
            #| extract d23 block |#
          (insert-term d23-mtrx index (- (icol p1) size-d21) (val p1))
          #| make 0's in index-th row of d21 |#
          (unless (or (zerop (val p1)) (eq (icol p1) index))
            (newsmith-column-op 4blocks (- (val p1)) index (icol p1))))))
    #| extract d31 and d33 blocks |#
    (do ((il (1+ size-d21) (1+ il)))
        ((> il nrows))
      (do ((pl01 (baselig 4blocks il))
           (p1 (left (baselig 4blocks il)) (left p1)))
          ((eq p1 pl01))
        (let ((p31 (left (baselig d31-mtrx (- il size-d21))))
              (p33 (left (baselig d33-mtrx (- il size-d21))))
              (ic (icol p1)))
          (if (> ic size-d21)
              (progn (majterme (chercher-hor p33 (- ic size-d21)) (basecol d33-mtrx (- ic size-d21)) (val p1))
                (setf p33 (left p33)))
            (progn (majterme (chercher-hor p31 ic) (basecol d31-mtrx ic) (val p1))
              (setf p31 (left p31)))))))
    (mtrx-substract d33-mtrx (newsmith-mtrx-prdc d31-mtrx d23-mtrx))))


(DEFUN DECOMPOSITION (heights elements-heights trgts srcs)
  (let* ((len (length heights))
         (targets-level (loop repeat len collect '()))
         (sources-level (loop repeat len collect '()))
         (critical-level (copy-list heights)))
    #| targets, sources and critical by levels |#
    (do ((targets trgts (cdr targets))
         (sources srcs (cdr sources)))
        ((endp targets))
      (let* ((source (car sources))
             (target (car targets))
             (height-source (svref elements-heights (1- source))))
        (setf (nth (1+ height-source) targets-level) (append (nth (1+ height-source) targets-level) (list target))); Assume h(target) = h(source) + 1
        (setf (nth height-source sources-level) (append (nth height-source sources-level) (list source)))
        (setf (nth (1+ height-source) critical-level) (remove target (nth (1+ height-source) critical-level)))
        (setf (nth height-source critical-level) (remove source (nth height-source critical-level)))))
    (return-from DECOMPOSITION (mapcar #'list targets-level sources-level critical-level))))


#|
  (setf ejemplo (random-2space 7))
  (setf finspace (2h-regularization ejemplo))
  (setf dvf (dvfield finspace))
  (setf a (heights finspace))
  (setf b (elements-height finspace))
  (setf c (mapcar #'cadr dvf))
  (setf d (mapcar #'car dvf))
  (decomposition a b c d)
|#


(DEFMETHOD H-REGULAR-DIF-DVF ((finspace finite-space) &key targets sources)
  (let* ((differentials '())
         (old-differential NIL)
         (old-matrix NIL)
         (edges (binarymatrice-to-ubasis (nilpot (stong finspace))))
         (elements-heights (elements-height finspace))
         (TSC-decomp (decomposition (heights finspace) elements-heights targets sources))) (format T "Bien la construccion")
    (push (creer-matrice 0 (length (third (car TSC-decomp)))) differentials)
    (setf old-matrix (creer-matrice 0 (length (car (heights finspace)))))
    (do ((dim_Cn-2 NIL (length (car heights)))
         (heights (heights finspace) (cdr heights))
         (TSC TSC-decomp (cdr TSC)))
        ((endp (cdr heights))) (format T "Bien el paso de dimensiones")
      (let* ((Cn-1 (car heights)) (TSC_n-1 (car TSC))
             (Cn (cadr heights))  (TSC_n (cadr TSC)) (dim_n (length Cn))
             (rslt (creer-matrice (length Cn-1) (length Cn))))
        (if dim_Cn-2
            (setf old-differential (let ((dn-1 old-matrix))
                                     (do ((group Cn (cdr group))
                                          (ic 1 (1+ ic)))
                                         ((> ic dim_n) rslt)
                                       (let* ((Uic (svref edges (1- (car group))))
                                              (iligs (posiciones-1 Uic Cn-1))
                                              (ker (car (newsmith-kernel (submatrix-cols dn-1 iligs))))
                                              (liste (mapcar #'(lambda (x) (list (nth (1- (car x)) iligs) (cadr x))) ker)))
                                         (maj-colonne rslt ic liste)))))
          (setf old-differential (do ((group Cn (cdr group))
                                      (ic 1 (1+ ic)))
                                     ((> ic dim_n) rslt)
                                   (let* ((Uic (svref edges (1- (car group))))
                                          (iligs (posiciones-1 Uic Cn-1))
                                          (liste (list (list (first iligs) -1) (list (second iligs) 1))))
                                     (maj-colonne rslt ic liste))))) (format T "Bien dn")
        (setf old-matrix old-differential)
        (push (d33-d31d21d23
               (submatrix old-differential
                          (posiciones-1 (append (second TSC_n-1) (third TSC_n-1)) Cn-1)
                          (posiciones-1 (append (first TSC_n) (third TSC_n)) Cn))
               (length (first TSC_n))) differentials)(format T "bien pequena")))
    (return-from H-REGULAR-DIF-DVF (reverse differentials))))


(DEFMETHOD H-REGULAR-DIF-DVF ((stong matrice) &key targets sources)
  (let ((finspace (build-finite-space :stong stong)))
        (H-REGULAR-DIF-DVF finspace :targets targets :sources sources)))
        
        
(DEFUN H-REGULAR-DIF-DVF-AUX (object targets sources)
  (H-REGULAR-DIF-DVF object :targets targets :sources sources))

#|
 (setf example (bar-subdivision (random-finite-space 6 .5)))
 (setf dvfield (dvfield-facets example))
 (setf trgts (mapcar #'cadr dvfield))
 (setf srcs (mapcar #'car dvfield))
 (dotimes (k (length (heights example)))
   (print (equalmatrix
           (nth k (H-REGULAR-DIF-DVF example :targets trgts :sources srcs))
           (nth k (H-REGULAR-DIF example)))))
|#


(DEFMETHOD H-REGULAR-HOMOLOGY-DVF-ALL ((finspace finite-space) &key targets sources)
  (let ((list (matrices-dvf (save-info-dvf-aux finspace targets sources)))
        (size (1- (length (heights finspace)))))
    (dotimes (dim (1+ size))
      (if (> dim size)
          (progn (format t "~3%Homology in dimension ~D :~%" dim)
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
            (format t "~3%Homology in dimension ~D :~%" dim)
            (dolist (item rsl)
              (declare (type list item))
              (format t "~2%Component Z")
              (unless (zerop (first item)) 
                (format t "/~DZ" (first item))))
            (terpri) (terpri)
            (done)))))))


(DEFMETHOD H-REGULAR-HOMOLOGY-DVF ((finspace finite-space) &rest rest)
  #| Compute the homology of 'finspace' directly from the list (H-REGULAR-DIF-DVF 'finspace')
The order in 'rest' must be: a dvfield and then the dimensions to compute homology |#
  (let ((size (1- (length (heights finspace))))
        (range '())
        (trgts '())
        (srcs '()))
    (case (length rest)
      (1 (setf range (<a-b> 0 size) ; rest = (dvfield)
               trgts (mapcar #'cadr (car rest))
               srcs (mapcar #'car (car rest))))
      (2 (if (integerp (cadr rest))  ; rest = (dvfield dim1)
             (if (> (cadr rest) size)
                 (return-from H-REGULAR-HOMOLOGY-DVF (progn (format t "~3%Homology in dimension ~D :~%" (cadr rest))
                                                       (terpri) (terpri)
                                                       (done)))
               (setf range (cdr rest)
                     trgts (mapcar #'cadr (car rest))
                     srcs (mapcar #'car (car rest))))
           (setf range (<a-b> 0 size) ; rest = (trgts srcs)
                 trgts (car rest)
                 srcs (cadr rest))))
      (3 (if (integerp (cadr rest)) ; rest = (dvfield dim1 dim2)
             (if (> (cadr rest) size)
                 (return-from H-REGULAR-HOMOLOGY-DVF (progn (format t "~3%Homology in dimension ~D :~%" (<a-b> (cadr rest) (caddr rest)))
                                                       (terpri) (terpri)
                                                       (done)))
               (setf range (<a-b> (cadr rest) (caddr rest)) 
                     trgts (mapcar #'cadr (car rest))
                     srcs (mapcar #'car (car rest))))
           (setf range (cddr rest) ; rest = (trgts srcs dim1)
                 trgts (car rest)
                 srcs (cadr rest))))
      (4 (setf range (<a-b> (caddr rest) (cadddr rest)) ; rest = (trgts srcs dim1 dim2)
               trgts (car rest)
               srcs (cadr rest))))
    
    (let ((list (matrices-dvf (save-info-dvf-aux finspace trgts srcs))))
      (dolist (dim range)
        (if (> dim size)
            (progn (format t "~3%Homology in dimension ~D :~%" dim)
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
              (format t "~3%Homology in dimension ~D :~%" dim)
              (dolist (item rsl)
                (declare (type list item))
                (format t "~2%Component Z")
                (unless (zerop (first item)) 
                  (format t "/~DZ" (first item))))
              (terpri) (terpri)
              (done))))))))


(DEFMETHOD H-REGULAR-HOMOLOGY-DVF ((stong matrice) &rest rest)
  (let ((finspace (build-finite-space :stong stong)))
    (typecase (second rest)
      (list (eval `(H-REGULAR-HOMOLOGY-DVF ,finspace (quote ,(first rest)) (quote ,(second rest)) ,@(cddr rest))))
      (integer (eval `(H-REGULAR-HOMOLOGY-DVF ,finspace (quote ,(first rest)) ,@(rest rest)))))))


(DEFUN H-REGULAR-HOMOLOGY-DVF-GENERATORS (finspace dvfield dim)
  (let ((chcm (chcm-h-regular-dvf finspace dvfield)))
    (chcm-homology-gen chcm dim)))
  

#|
 (setf example (bar-subdivision (random-finite-space 6 .5)))
 (setf dvfield (dvfield-facets example))
 (setf trgts (mapcar #'cadr dvfield))
 (setf srcs (mapcar #'car dvfield))
 (H-REGULAR-HOMOLOGY-DVF example dvfield 1 3)
 (H-REGULAR-HOMOLOGY-DVF example trgts srcs 1 3)
 (H-REGULAR-HOMOLOGY-DVF (stong example) dvfield 1 3)

 (setf dvfield (dvfield (finite-model-sphere 3)))
 (H-REGULAR-HOMOLOGY-DVF (finite-model-sphere 3) dvfield)

((lambda (n)
   (setf example (bar-subdivision (random-finite-space n .4)))
   (setf dvfield (dvfield-facets example))
   (setf tar (mapcar #'second dvfield))
   (setf sou (mapcar #'first dvfield))
   (top example)
   (heights example)
   (list (tiempo (lambda () (H-REGULAR-HOMOLOGY example)))
         ;(tiempo (lambda () (H-REGULAR-HOMOLOGY-DVF example tar sou)))
         (tiempo (lambda () (H-REGULAR-HOMOLOGY-DVF-ALL example :targets tar :sources sou))))) 6)
|#


(DEFVAR *H-REGULAR-DIF-DVF-HASH-TABLE*)
(setq *H-REGULAR-DIF-DVF-HASH-TABLE* (make-hash-table))


(DEFSTRUCT (SAVE-INFO-DVF (:conc-name nil))
  targets-dvf sources-dvf matrices-dvf)


(DEFUN SAVE-INFO-DVF-AUX (finspace trgts srcs)
  #| Assign to (idnm 'finspace') the info of (h-regular-dif-dvf 'finspace') :targets-dvf 'trgts' :sources-dvf 'srcs') |#
  (let ((info (gethash (idnm finspace) *H-REGULAR-DIF-DVF-HASH-TABLE*)))
    (if (and info (equal (targets-dvf info) trgts) (equal (sources-dvf info) srcs))
        info
      (setf (gethash (idnm finspace) *H-REGULAR-DIF-DVF-HASH-TABLE*)
            (make-save-info-dvf
             :targets-dvf trgts
             :sources-dvf srcs
             :matrices-dvf (h-regular-dif-dvf finspace :targets trgts :sources srcs))))))


(DEFUN TOP-BASIS-DVF (finspace trgts srcs)
  #| Provide the slot :basis of (CHCM-H-REGULAR-DVF 'finspace') |#
  (flet ((rslt (dim)
           (let ((heights (heights finspace)))
             (cond ((= dim -1) (list '*))
                   ((< dim -1) (error "Wrong dimension!"))
                   ((< dim (length heights)) (let ((critical_dim (set-difference (set-difference (nth dim heights) srcs) trgts)))
                                               (loop for k from 0 to (1- (length critical_dim))
                                                     collect (with-x (nth k critical_dim)))))
                   (T +empty-list+)))))
    (the basis #'rslt)))


(DEFUN TOP-DIFF-DVF (finspace trgts srcs)
  #| Provide the slot :intr-dffr of (CHCM-H-REGULAR-DVF 'finspace') |#
  (flet ((frslt (dim gnr)
           (let ((heights (heights finspace)))
             (unless (<= 0 dim (length heights))
               (error "Wrong dimension!"))
             (let ((num_gnr 0)
                   (critical_dim (set-difference (set-difference (nth dim heights) srcs) trgts)))
               (setf num_gnr (read-from-string (subseq (write-to-string gnr) 1))) ; num_gnr of X13 is 13
               (let ((j (position num_gnr critical_dim)))  
                 (if (null j)
                     (error "Wrong generator!")
                   (if (zerop dim)
                       (cmbn -1)
                     (let* ((rslt (cmbn (1- dim)))
                            (critical_dim-1 (set-difference (set-difference (nth (1- dim) heights) srcs) trgts))
                            (info (gethash (idnm finspace) *H-REGULAR-DIF-DVF-HASH-TABLE*))
                            (dif (nth dim (matrices-dvf info))))
                       (setf (cmbn-list rslt) (loop for i from 1 to (length critical_dim-1)
                                                    unless (eq (terme dif i (1+ j)) 0)
                                                    collect (term (terme dif i (1+ j)) (with-x (nth (1- i) critical_dim-1)))))
                       rslt))))))))
    #'frslt))


(DEFMETHOD CHCM-H-REGULAR-DVF ((finspace finite-space) &rest rest)
  (let ((trgts '())
        (srcs '()))
    (case (length rest)
      (1 (setf trgts (mapcar #'cadr (car rest))
               srcs (mapcar #'car (car rest))))
      (2 (setf trgts (car rest) srcs (cadr rest))))
    (save-info-dvf-aux finspace trgts srcs)
    (build-chcm :cmpr #'x-cmpr
                :strt :GNRT
                :basis (top-basis-dvf finspace trgts srcs)
                :intr-dffr (top-diff-dvf finspace trgts srcs)
                :orgn `(CHCM-H-REGULAR-DVF ,finspace ,trgts ,srcs))))


(DEFMETHOD CHCM-H-REGULAR-DVF ((stong matrice) &rest rest)
  (let ((finspace (build-finite-space :stong stong)))
    (eval `(CHCM-H-REGULAR-DVF ,finspace ,@rest))))

#| 
  (setf example (bar-subdivision (random-finite-space 6 .5)))
  (setf dvfield (dvfield-facets example))
  (setf chcm (CHCM-H-REGULAR-DVF example dvfield))
  (basis chcm -1)
  (basis chcm 0)
  (basis chcm 1)
  (basis chcm 2)
  (dffr chcm 2 (car (basis chcm 2)))
|#
