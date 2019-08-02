;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*

;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 
;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 
;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 

(IN-PACKAGE #:cat)

(PROVIDE "finite-topological-spaces-homology")

;;
;;  Homology of h-regular spaces
;;


(DEFUN COLUMN-J (mtrx j)
   (declare (type matrix mtrx)
            (type fixnum j))
   (assert (< j (column-number mtrx)))
   (the list
     (let ((line-n (line-number mtrx)))
        (declare (type fixnum line-n))
        (mapcar
         #'(lambda (i)
             (declare (type fixnum i))
             (aref mtrx i j))
         (<a-b> 0 (1- line-n))))))


(DEFUN KERNEL (mtrx)
   (declare (type matrix mtrx))
   (the list
     (let* ((smith-list (smith mtrx))
            (m (third smith-list))
            (q (fourth smith-list))
            (min (min (line-number mtrx) (column-number mtrx)))
            (rslt nil))
        (declare
          (type list smith-list rslt)
          (type matrix m q)
          (type fixnum min))
        (progn
         (do ((j (1- (column-number mtrx)) (1- j)))
             ((< j min))
            (declare (type fixnum j))
            (setq rslt (nconc (list (column-j q j)) rslt)))
         (do ((i (1- min) (1- i)))
             ((or (< i 0) (not (eql  0 (aref m i i)))))
            (declare (type fixnum i))
            (setq rslt (nconc (list (column-j q i)) rslt)))
         rslt))))


(DEFMETHOD H-REGULAR-DIF ((finspace finite-space))
  (let ((diferenciales NIL)
        (aristas (matrice-to-lmtrx (nilpot (stong finspace))))
        (alturas (heights finspace)))
    (declare (type list diferenciales))
    (dotimes (dimension (length alturas))
      (let* ((Cn (nth dimension alturas))
             (dim_n (length Cn)))
        (if (zerop dimension)
            (push (creer-matrice 0 (length Cn)) diferenciales)
          (let* ((Cn-1 (nth (1- dimension) alturas))
                 (dim_n-1 (length Cn-1))
                 (rslt (creer-matrice dim_n-1 dim_n)))
            (case dimension          
              (1 (progn (push (do ((group Cn (cdr group))
                                   (posicion 1 (1+ posicion)))
                                  ((endp group) rslt)
                                (let ((colj (aref aristas (1- (car group)))))
                                  (insert-term rslt (1+ (position (1+ (car (first colj))) Cn-1))  posicion -1)
                                  (insert-term rslt (1+ (position (1+ (car (second colj))) Cn-1)) posicion  1)))
                              diferenciales)))
              
              (otherwise (push (let ((Cn-2 (<a-b> 1 (length (nth (- dimension 2) alturas))))
                                     (dn-1 (first diferenciales))
                                     (group Cn))
                                 (do ((posicion 1 (1+ posicion)))
                                     ((null group) rslt)
                                   (let* ((maxUj (mapcar #'(lambda (v) (1+ (car v))) (aref aristas (1- (pop group)))))
                                          (cols (posiciones-1 maxUj Cn-1))
                                          (ker (car (kernel (make-array (list (length Cn-2) (length cols))
                                                                        :element-type 'fixnum
                                                                        :initial-contents (loop for i in Cn-2
                                                                                                collect (loop for j in cols
                                                                                                              collect (terme dn-1 i j))))))))
                                     (do ((ker ker (cdr ker))
                                          (cols cols (cdr cols)))
                                         ((endp cols))
                                       (insert-term rslt (car cols) posicion (car ker))))))
                               diferenciales)))))))
    (return-from H-REGULAR-DIF (reverse diferenciales))))
    
    
(DEFMETHOD H-REGULAR-DIF ((topogenous matrice))
  (let ((finspace (build-finite-space :top topogenous)))
    (H-REGULAR-DIF finspace)))


(DEFUN POSICIONES (pos list)
  (loop for x in pos
    collect (nth (1- x) list)))


(DEFUN POSICIONES-1 (list1 list2)   ; 'list1' is a sublist of 'list2'
  (loop for x in list1
        collect (1+ (position x list2 :test #'equal))))


(DEFUN D21-1 (blocks size)
  #| Modify 'blocks' in order to obtain an identity matrix in the submatrix 'size'x'size' |#
  (declare (type fixnum size))
  (the matrice
       (let ((rslt (copier-matrice blocks)))
         (unless (or (< (nlig rslt) size) (< (ncol rslt) size) (= size 0))
           (loop for k from 1 to size
                 do (progn
                      ;(if (eq (terme rslt k k) 0)
                          ;(0column-swap rslt k (loop for j from (1+ k) to size
                                   ;                  when (not (eq (terme rslt k j) 0))
                                    ;                 return j)))
                      (unless (eq (terme rslt k k) 1)
                        (0column-op rslt (1- (terme rslt k k)) k k))
                      (loop for j from 1 to size
                            if (not (eq j k))
                            do (unless (eq (terme rslt k j) 0)
                                 (0column-op rslt (- (terme rslt k j)) k j))))))
         rslt)))



(DEFUN RESTA (mtrx1 mtrx2)
  (let ((rslt (creer-matrice (nlig mtrx1) (ncol mtrx2))))
    (loop for i from 1 to (nlig mtrx1)
          do (loop for j from 1 to (ncol mtrx1)
                   do (insert-term rslt i j (- (terme mtrx1 i j) (terme mtrx2 i j)))))
    rslt))


(DEFUN FINAL (matrix corte) ; 'matrix' is the result of the function D21-1 and 'corte' is the 'size' in D21-1 function
  (let ((d31 (sub-matrix matrix
                        (>a-b> corte (nlig matrix))
                        (<a-b> 1 corte)))
        (d23 (sub-matrix matrix
                        (<a-b> 1 corte )
                        (>a-b> corte (ncol matrix))))
        (d33 (sub-matrix matrix
                        (>a-b> corte (nlig matrix))
                        (>a-b> corte (ncol matrix)))))
    (resta d33 (matrice-prdc d31 d23))))


(DEFSTRUCT (TSC (:conc-name nil))
  name diferencial group targets sources critical)


(DEFMACRO DESCOMP (name)
  `(setf ,name (make-TSC)
         (sources ,name) (reverse (intersection sou summands))
         (targets ,name) (reverse (intersection tar summands))
         (critical ,name) (let ((cri summands))
                            (dolist (x (append (targets ,name) (sources ,name)))
                              (setf cri (remove x cri)))
                            cri)))


(DEFMETHOD DVF-H-REGULAR-DIF ((finspace finite-space) &key sou tar) 
  (let* ((groups NIL)
         (diferenciales NIL)
         (dn-3 NIL)
         (dn-2 NIL)
         (dn-1 NIL)
         (dn NIL)
         (alturas (heights finspace))
         (height (length alturas))) 
    (declare (type list groups diferenciales)
             (type fixnum height))
    (setf dn-3 (make-TSC))
    (dotimes (dimension height)
      (let ((summands (nth dimension alturas)))
        (case dimension
          (0  (progn (descomp dn-2)
                ;;; Computing C0':
                (push (mapcar #'(lambda (u) (list (list 1 u))) summands) groups)
                ;;; Computing d0:
                (push (creer-matrice 0 (length (critical dn-2))) diferenciales)))
          
          (1  (progn (descomp dn-1)
                ;;; Computing C1':
                (push (loop for x in summands
                            collect (let* ((primer1 (find-if #'(lambda (y) (= (terme (top finspace) y x) 1))
                                                             (<a-b< 1 (nlig (top finspace)))))
                                           (segundo1 (find-if #'(lambda (y) (= (terme (top finspace) y x) 1))
                                                              (>a-b> primer1 (nlig (top finspace))))))
                                      (list (list 1 segundo1 x) (list -1 primer1 x)))) groups)
                
                ;;; Computing d1'':
                (push (let ((rslt (creer-matrice (length (second groups)) (length summands))))
                        (dolist (x (apply #'append (first groups)))
                          (insert-term rslt (1+ (position (second x) (first alturas) :test #'equal)) 
                                       (1+ (position (third x) (second alturas) :test #'equal))
                                       (first x)))
                        rslt) diferenciales)))
          
          (otherwise (progn (descomp dn)
                       (let ((rslt NIL)
                             (Cn-1 (first groups))
                             (Cn-2 (second groups)))
                         
                         ;;; Cálculo de Cn' (rslt):
                         (dolist (x (reverse summands))
                           (let ((pseudo NIL)
                                 (list_fil NIL)
                                 (list_col NIL)
                                 (kernel NIL))
                             
                             (setf list_fil (loop for k from 0 to (1- (length Cn-2)) 
                                                   when (= (terme (top finspace) (car (last (first (nth k Cn-2)))) x) 1)
                                                   collect k)
                                   list_col (loop for k from 0 to (1- (length Cn-1))
                                                   when (= (terme (top finspace) (car (last (first (nth k Cn-1)))) x) 1)
                                                   collect k)
                                   
                                 kernel (car (kernel (convertmatrice (sub-matrix (first diferenciales)
                                                                                (mapcar #'(lambda (x) (1+ x)) list_fil)
                                                                                (mapcar #'(lambda (x) (1+ x)) list_col))))))
                             
                             (do ((lst list_col (cdr lst))
                                  (ker kernel (cdr ker)))
                                 ((endp lst))
                               (unless (equal (car ker) 0)
                                 (let ((actualizar (nth (car lst) Cn-1)))
                                   (dolist (z actualizar)
                                     (push (reverse (cons x (reverse (cons (* (car ker) (car z)) (cdr z))))) pseudo)))))
                             
                             (push pseudo rslt)))
                         
                         ;;; Computing dn'' (rslt2):
                         (push (let ((rslt2 (creer-matrice (length Cn-1) (length summands))))
                               (dotimes (j (length summands) rslt2)
                                 (let ((next-to-last 0))
                                   (do ((tj (nth j rslt) (remove-if #'(lambda (elto)
                                                                        (equal next-to-last (second (reverse elto))))
                                                                    tj)))
                                       ((null tj))
                                     
                                     (setf next-to-last (second (reverse (first tj))))
                                     (let* ((g (nth (position next-to-last (nth (1- dimension) alturas)) Cn-1))
                                            (ui (car (last g))))
                                       (dolist (l tj)
                                         (if (equal (cdr ui) (reverse (cdr (reverse (cdr l)))))
                                             (progn
                                               (insert-term rslt2 (1+ (position g Cn-1)) (1+ j) (* (first ui) (first l)))
                                               (return))))))))) diferenciales)
                         (setf groups (list (car groups)))
                         (push rslt groups))
                       
                       (unless (= dimension 2)
                         (setf (third diferenciales)
                               (final
                                (d21-1
                                 (sub-matrix (third diferenciales)
                                            (posiciones-1 (append (sources dn-3) (critical dn-3)) (nth (- dimension 3) alturas))
                                            (posiciones-1 (append (targets dn-2) (critical dn-2)) (nth (- dimension 2) alturas)))
                                 (length (targets dn-2)))
                                (length (targets dn-2)))))
                       
                       (setf dn-3 dn-2
                             dn-2 dn-1
                             dn-1 dn))))))
    
    ;;; Modify the last 2 matrices
    (unless (= height 1)
      (unless (= height 2)
        (setf (second diferenciales)
              (final
               (d21-1
                (sub-matrix (second diferenciales)
                           (posiciones-1 (append (sources dn-3) (critical dn-3)) (nth (- height 3) alturas))
                           (posiciones-1 (append (targets dn-2) (critical dn-2)) (nth (- height 2) alturas)))
                (length (targets dn-2)))
               (length (targets dn-2)))))
      
      (setf (first diferenciales)
            (final
             (d21-1
              (sub-matrix (first diferenciales)
                         (posiciones-1 (append (sources dn-2) (critical dn-2)) (nth (- height 2) alturas))
                         (posiciones-1 (append (targets dn-1) (critical dn-1)) (nth (- height 1) alturas)))
              (length (targets dn-1)))
             (length (targets dn-1)))))

    (return-from DVF-H-REGULAR-DIF (reverse diferenciales))))


(DEFMETHOD DVF-H-REGULAR-DIF ((topogenous matrice) &key sou tar)
  (let ((finspace (build-finite-space :top topogenous
                                      :stong (stongmat topogenous))))
        (DVF-H-REGULAR-DIF finspace :sou sou :tar tar)))

