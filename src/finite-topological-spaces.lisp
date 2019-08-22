;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*

;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 
;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 
;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 

(IN-PACKAGE #:cat)

(PROVIDE "finite-topological-spaces")

;;
;;   Auxiliary functions (most of them in 'NewSmith')
;;


(DEFUN MAJTERME (pl pc val)
  #| Same 'maj-terme' but changing '=' by 'eq' |#
  (if (eq val 0)
      (if (eq (left pl) (up pc)) (supprimer-terme pl pc))
    (if (eq (left pl) (up pc))
          (setf (val (left pl)) val)
      (inserer-terme pl pc val))))


(DEFUN EQUALMATRIX (mtrx1 mtrx2)
  #| Same 'equal-matrix' but changing 'line-number' by 'nlig' and 'column-number' by 'ncol' |#
  (declare (type matrice mtrx1 mtrx2))
  (the boolean
    (let ((line-number (nlig mtrx1)))
      (declare (type fixnum line-number))
      (unless (eql line-number (nlig mtrx2))
        (return-from equalmatrix +false+))
      (unless (eql (ncol mtrx1) (ncol mtrx2))
        (return-from equalmatrix +false+))
      (do ((il 1 (1+ il)))
          ((> il line-number))
        (declare (type fixnum il))
        (let ((p10 (baselig mtrx1 il))
              (p20 (baselig mtrx2 il)))
          (declare (type t-mat p10 p20))
          (do ((p1 (left p10) (left p1))
               (p2 (left p20) (left p2)))
              (nil)
            (declare (type t-mat p1 p2))
            (when (eq p1 p10)
              (unless (eq p2 p20)
                (return-from equalmatrix +false+))
              (return))
            (when (eq p2 p20)
              (unless (eq p1 p10)
                (return-from equalmatrix +false+)))
            (unless (eql (icol p1) (icol p2))
              (return-from equalmatrix +false+))
            (unless (eql (val p1) (val p2))
              (return-from equalmatrix +false+)))))
      (return-from equalmatrix +true+))))


(DEFUN MATRICE-PRDC (mtrx1 mtrx2)
  #| Same 'mtrx-prdc' but changing 'line-number' by 'nlig' and 'column-number' by 'ncol' |#
  (declare (type matrice mtrx1 mtrx2))
  (the matrice
    (let ((nl (nlig mtrx1))
          (nc (ncol mtrx2))
          (n (ncol mtrx1)))
      (declare (type fixnum nl nc n))
      (assert (eql n (nlig mtrx2)))
      (let ((rslt (creer-matrice nl nc)))
        (declare (type matrice rslt))
        (do ((il 1 (1+ il)))
            ((> il nl))
          (declare (type fixnum il))
          (do ((pl01 (baselig mtrx1 il))
               (pl0r (baselig rslt il))               
               (ic 1 (1+ ic)))
              ((> ic nc))
            (declare (type t-mat pl01 pl0r) (type fixnum ic))
            (let ((pl (left pl01))
                  (pc (up (basecol mtrx2 ic)))
                  (sum 0))
              (declare (type t-mat pl pc) (type fixnum sum))
              (loop
                (let ((ilic (icol pl))
                      (icil (ilig pc)))
                  (declare (type fixnum ilic icil))
                  (when (zerop ilic) (return))
                  (when (zerop icil) (return))
                  (cond ((eql ilic icil)
                         (incf sum (safe-* (val pl) (val pc)))
                         (setf pl (left pl) pc (up pc)))
                        ((> ilic icil)
                         (setf pl (left pl)))
                        (t
                         (setf pc (up pc))))))
              (unless (zerop sum)
                (inserer-terme pl0r (basecol rslt ic) sum)))))
        rslt))))


(DEFUN MATRICE-TO-LMTRX (mtrx)
  (let* ((column-n (ncol mtrx))
        (rslt (make-array column-n)))
    (dotimes (j column-n)
      (let ((rsltj '())
            (columnj-pair-list 
             (let ((ptc (basecol mtrx (1+ j)))
                   (res '()))
               (do ((pc (up ptc) (up pc)))
                   ((eq pc ptc))
                 (push (list (ilig pc) (val pc)) res))
               res)))
        
        (mapcar #'(lambda (pair)
                    (push  (list (1- (first pair)) (second pair)) rsltj))
                columnj-pair-list)
        
        (setf (svref rslt j) (nreverse rsltj) )))
    rslt))


(DEFUN SAFE-* (arg1 arg2)
  #| This is in 'NewSmith' |#
  (declare (type fixnum arg1 arg2))
  (let ((rslt (* arg1 arg2)))
    (declare (type integer rslt))
    (if (typep rslt 'fixnum)
        rslt
      (error "Fixnum-overflow in (safe-* ~D ~D)." arg1 arg2))))


(DEFUN SAFE-+ (arg1 arg2)
  #| This is in 'NewSmith' |#
  (declare (type fixnum arg1 arg2))
  (let ((rslt (+ arg1 arg2)))
    (declare (type integer rslt))
    (if (typep rslt 'fixnum)
        rslt
      (error "Fixnum-overflow in (safe-+ ~D ~D)." arg1 arg2))))


(DEFUN SAFE-- (arg)
  #| This is in 'NewSmith' |#
  (declare (type fixnum arg))
  (let ((rslt (- arg)))
    (declare (type integer rslt))
    (if (typep rslt 'fixnum)
        rslt
      (error "Fixnum-overflow in (SAFE-- ~D)." arg))))


(DEFUN 0NEW-COLUMN (mtrx icol list)
  #| This is 'new-column' in 'NewSmith' |#
  (declare (type matrice mtrx) (type fixnum icol) (type list list))
  (the matrice
    (let ((peigne (peigne-hor mtrx (baselig mtrx 0) icol))
          (pc0 (basecol mtrx icol)))
      (declare (type list peigne) (type t-mat pc0))
      (map nil #'(lambda (item)
                   (declare (type t-mat item))
                   (when (eq (up pc0) (left item))
                     (supprimer-terme item pc0)))
        peigne)
      (do ((markp (nreverse peigne) (cdr markp))
           (ilig 1 (1+ ilig))
           (markc list))
          ((endp markc))
        (declare (type list markp markc) (type fixnum ilig))
        (when (eql (caar markc) ilig)
          (inserer-terme (car markp) pc0 (cadar markc))
          (setf markc (cdr markc))))      
      mtrx)))


(DEFUN 0COLUMN-OP (mtrx lambda icol1 icol2)
  #| This is 'column-op' in 'NewSmith' |#
  (declare
   (type matrice mtrx)
   (type fixnum lambda icol1 icol2))
  (the matrice
    (let ((mark1 (extract-column mtrx icol1))
          (mark2 (extract-column mtrx icol2))
          (new-column2 +empty-list+))
      (declare (type list mark1 mark2 new-column2))
      (loop
        (when (endp mark1)
          (setf new-column2 (nreconc new-column2 mark2))
          (return))
        (when (endp mark2)
          (setf new-column2
            (nreconc new-column2
                     (mapcar #'(lambda (item)
                                 (declare (type list item))
                                 (list (first item)
                                       (safe-* lambda (second item))))
                       mark1)))
          (return))
        (let ((ilig1 (caar mark1))
              (ilig2 (caar mark2)))
          (declare (type fixnum ilig1 ilig2))
          (cond ((< ilig1 ilig2)
                 (push (list ilig1 (safe-* lambda (second (pop mark1))))
                       new-column2))
                ((> ilig1 ilig2)
                 (push (list ilig2 (second (pop mark2))) new-column2))
                (t (let ((new-val (safe-+ (safe-* lambda (second (pop mark1)))
                                     (second (pop mark2)))))
                     (declare (type fixnum new-val))
                     (unless (zerop new-val)
                       (push (list ilig1 new-val) new-column2)))))))
      (0new-column mtrx icol2 new-column2)
      mtrx)))


(DEFUN 0COLUMN-SWAP (mtrx icol1 icol2)
  #| This is 'column-swap' in 'NewSmith' |#
  (declare
    (type matrice mtrx) (type fixnum icol1 icol2))
  (the matrice
    (let ((new-column1 (extract-column mtrx icol2))
          (new-column2 (extract-column mtrx icol1)))
      (declare (type list new-column1 new-column2))
      (0new-column mtrx icol1 new-column1)
      (0new-column mtrx icol2 new-column2)
      mtrx)))


(DEFUN 0COLUMN-MINUS (mtrx icol)
  #| This is 'column-minus' in 'NewSmith' |#
  (declare (type matrice mtrx) (type fixnum icol))
  (the matrice
    (0new-column mtrx icol
                (mapcar #'(lambda (item)
                            (declare (type list item))
                            (the list
                              (list (first item)
                                    (safe-- (second item)))))
                  (extract-column mtrx icol)))))


(DEFUN EXTRACT-COLUMN (mtrx icol)
  #| This is in 'NewSmith' |#
  (declare (type matrice mtrx) (type fixnum icol))
  (the list
    (let ((pc0 (basecol mtrx icol)))
      (declare (type t-mat pc0))
      (do ((pc (up pc0) (up pc))
           (rslt +empty-list+ (cons (list (ilig pc) (val pc)) rslt)))
          ((eq pc pc0) rslt)
        (declare (type t-mat pc) (type list rslt))))))


(DEFUN EXTRACT-LINE (mtrx ilig)
  #| This is in 'NewSmith' |#
  (declare (type matrice mtrx) (type fixnum ilig))
  (the list
    (let ((pl0 (baselig mtrx ilig)))
      (declare (type t-mat pl0))
      (do ((pl (left pl0) (left pl))
           (rslt +empty-list+ (cons (list (icol pl) (val pl)) rslt)))
          ((eq pl pl0) rslt)
        (declare (type t-mat pl) (type list rslt))))))


;;
;;  Finite topological spaces class
;;


(DEFMACRO INSERT-TERM (mat ilig icol val)
  #| Insert the value 'val' in the ('ilig' 'icol')-entry of 'mat' |#
  `(majterme (chercher-hor (baselig ,mat ,ilig) ,icol)
             (chercher-ver (basecol ,mat ,icol) ,ilig)
             ,val))


(DEFUN CONVERTARRAY (array)
  #| Convert an 'array' to a 'matrice' |#
  (let* ((numfil (array-dimension array 0))
         (numcol (array-dimension array 1))
         (rslt (creer-matrice numfil numcol)))
    (dotimes (j numcol)
      (dotimes (i numfil)
        (insert-term rslt (1+ i) (1+ j) (aref array i j))))
    rslt))


(DEFUN EXTRACT-COLUMN2 (mtrx icol list)
  #| Modification of 'EXTRACT-COLUMN' in order to extract only the elements in the column 'icol' whose row indexes are in the list 'list' |#
  (declare (type matrice mtrx) (type fixnum icol) (type list list))
  (the list
       (let ((pc0 (basecol mtrx icol))
             (lst list))
         (declare (type t-mat pc0))
         (do ((pc (up pc0) (up pc))
              (rslt +empty-list+ (if (find (ilig pc) lst)
                                     (progn
                                       (setf lst (remove (ilig pc) lst))
                                       (cons (list (ilig pc) (val pc)) rslt))
                                   rslt)))
             ((eq pc pc0) rslt)
           (declare (type t-mat pc) (type list rslt))))))


(DEFUN SUB-MATRIX (mtrx rows cols)
  #| Extract the submatrix of 'mtrx' formed by the rows whose indexes are in 'rows' the columns whose indexes are in 'cols' |#
  (declare (type matrice mtrx) (type list rows cols))
  (let ((rslt (creer-matrice (length rows) (length cols))))
    (do ((k 1 (1+ k))
         (columnas cols (cdr columnas)))
        ((endp columnas) rslt)
      (0new-column rslt k (let ((valor NIL))
                            (do ((j (length rows) (1- j))
                                 (rslt2 +empty-list+ (if (setf valor (second (find (nth (1- j) rows) (extract-column2 mtrx (car columnas) rows) :key #'first)))
                                                         (cons (list j valor) rslt2)
                                                       rslt2)))
                                ((< j 1) rslt2)))))))


(DEFUN NILPOT (matrix)  
  (declare (type matrice matrix))
  (let ((rslt (copier-matrice matrix)))
    (loop for k from 1 to (nlig rslt)
          do (insert-term rslt k k 0))
    rslt))


(DEFUN NILPOT-1 (matrix)  
  (declare (type matrice matrix))
  (let ((rslt (copier-matrice matrix)))
    (loop for k from 1 to (nlig rslt)
          do (insert-term rslt k k 1))
    rslt))


(DEFUN ADD-ROWS (matrice pos num)
  #| Add 'num' rows from position 'pos' |#
  (let ((rslt (creer-matrice (+ (n-lig matrice) num) (n-col matrice))))
    (do ((i 1 (1+ i)))
        ((eq i pos))
      (do ((ptl1 (baselig matrice i))
           (pl1 (left (baselig matrice i)) (left pl1))
           (pl2 (baselig rslt i) (left pl2)))
          ((eq pl1 ptl1))
        (inserer-terme pl2 (basecol rslt (icol pl1)) (val pl1))))
    
    (do ((i pos (1+ i)))
        ((> i (n-lig matrice)) rslt)
      (do ((ptl1 (baselig matrice i))
           (pl1 (left (baselig matrice i)) (left pl1))
           (pl2 (baselig rslt (+ i num)) (left pl2)))
          ((eq pl1 ptl1))
        (inserer-terme pl2 (basecol rslt (icol pl1)) (val pl1))))))


(DEFUN ADD-COLS (matrice pos num)
  #| Add 'num' columns from position 'pos' |#
  (let ((rslt (creer-matrice (n-lig matrice) (+ (n-col matrice) num))))
    (do ((i 1 (1+ i)))
        ((eq i pos))
      (do ((ptc1 (basecol matrice i))
           (pc1 (up (basecol matrice i)) (up pc1))
           (pc2 (basecol rslt i) (up pc2)))
          ((eq pc1 ptc1))
        (inserer-terme (baselig rslt (ilig pc1)) pc2 (val pc1))))

      (do ((i pos (1+ i)))
          ((> i (n-col matrice)) rslt)
        (do ((ptc1 (basecol matrice i))
             (pc1 (up (basecol matrice i)) (up pc1))
             (pc2 (basecol rslt (+ i num)) (up pc2)))
            ((eq pc1 ptc1))
          (inserer-terme (baselig rslt (ilig pc1)) pc2 (val pc1))))))


(DEFMACRO ADD-POINTS (matrice pos num)
  #| Add 'num' points from position 'pos' |#
  `(add-cols (add-rows ,matrice ,pos ,num) ,pos ,num))


(DEFUN STONGMAT (topogenous)
  #| Stong matrice from the topogenous matrice 'topogenous' |#
  (declare (type matrice topogenous))
  (let* ((dim (nlig topogenous))
         (rslt (identite dim)))
    (declare (type fixnum dim))
    (if (or (eq dim 1) (eq dim 2))
        (return-from stongmat topogenous)
      
      (do ((i 1 (1+ i))) ((> i dim) rslt)
        (do ((j (1+ i) (1+ j))) ((> j dim))
          (cond ((eq j (1+ i)) (insert-term rslt i j (terme topogenous i j)))
                ((eq (terme topogenous i j) 0))
                (T (let ((parada 0))
                     (do ((k (1+ i) (1+ k))) ((or (> k (1- j)) (eq parada 1)))
                       (cond ((and (eq (terme topogenous i k) 1) (eq (terme topogenous k j) 1)) (setf parada 1))
                             ((and (eq k (1- j)) (eq parada 0)) (insert-term rslt i j 1))))))))))))


(DEFUN TOPMAT (stong)
  #| Topogenous matrice from the Stong matrice 'stong' |#
  (declare (type matrice stong))
  (let* ((dim (nlig stong))
         (rslt (identite dim)))
    (declare (type fixnum dim))
    (if (or (eq dim 1) (eq dim 2))
        (return-from topmat stong)

      (do ((k 2 (1+ k))) ((> k dim) rslt)
        (if (eq k 2)
            (do ((i 1 (1+ i))) ((> i (1- dim)))
              (insert-term rslt i (1+ i) (terme stong i (1+ i))))
          
          (do ((i 1 (1+ i))) ((> i (1+ (- dim k))))
            
            (if (eq (terme stong i (1- (+ i k))) 1)
                (insert-term rslt i (1- (+ i k)) 1)
              
              (if (eq (let ((indicator 0))
                        (do ((r (1+ i) (1+ r))) ((or (> r (- (+ i k) 2)) (eq indicator 1))) 
                          (if (and (eq (terme rslt i r) 1) (eq (terme rslt r (1- (+ i k))) 1))
                              (setf indicator 1))) indicator) 1)
                  (insert-term rslt i (1- (+ i k)) 1)))))))))


(DEFUN RANDOMTOP (dim dens)
  #| A 'dim' x 'dim' random topogenous matrice |#
  (declare (fixnum dim))
  (let ((limite (truncate (/ (* dens (- dim 2) (- dim 1)) 2)))
        (rslt (identite dim))
        (contador 1)
        (base (make-array dim)))
    (unless (= dim 1)   
      (do () ((> contador limite))
        (let* ((j1 (1+ (random (1- dim))))
               (i1 (random j1)))
          (if (find i1 (aref base j1))
              (incf contador) ; Avoiding an infinite loop...
            (progn
              (push i1 (aref base j1))
              (setf (aref base j1) (union (aref base j1) (aref base i1)))
              (loop for k in (>a-b< j1 dim)
                    if (find j1 (aref base k))
                    do (and (pushnew i1 (aref base k))
                            (setf (aref base k) (union (aref base k) (aref base i1))) ; Improve...
                            (incf contador)))))))
      (loop for j in (<a-b< 1 dim)
            do (mapcar #'(lambda (elem) (insert-term rslt (1+ elem) (1+ j) 1)) (aref base j))))
    rslt))


(DEFMACRO RANDOMSTONG (dim dens)
  #| A 'dim' x 'dim' random Stong matrice |#
  `(stongmat (randomtop ,dim ,dens)))


(DEFUN HEIGHTS-AUX (stong lst)
  (if (null lst)
      NIL
    (let ((minimals (loop for n in lst
                          if (null (loop for i from 0 to (1- (position n lst))
                                         thereis (eq (terme stong (nth i lst) n) 1)))
                          collect n)))
      (cons minimals (heights-aux stong (set-difference lst minimals)))))) ; 'set-difference' must be in ascending order

(DEFMETHOD HEIGHTS ((stong matrice))
  #| List of elements separated by heights |#
  (heights-aux stong (<a-b> 1 (ncol stong))))


(DEFCLASS FINITE-SPACE ()
   ;; TOPogenous matrix
  ((top :type matrice :initarg :top :reader top)
   ;; STONG matrix
   (stong :type matrice :initarg :stong :reader stong)
   ;; HEIGHTS
   (heights :type list :initarg :heights :reader heights)
   ;; IDentification NuMber
   (idnm :type fixnum :initform (incf *idnm-counter*) :reader idnm)
   ;; ORiGiN
   (orgn :type list :initarg :orgn :reader orgn)))


(DEFMETHOD SLOT-UNBOUND (class (finspace finite-space) (slot-name (eql 'heights)))
  (declare (ignore class))
  (the list
       (setf (slot-value finspace 'heights) (heights (stong finspace)))))


(DEFMETHOD SLOT-UNBOUND (class (finspace finite-space) (slot-name (eql 'top)))
  (declare (ignore class))
  (the matrice
       (setf (slot-value finspace 'top) (topmat (stong finspace)))))


(DEFMETHOD SLOT-UNBOUND (class (finspace finite-space) (slot-name (eql 'stong)))
  (declare (ignore class))
  (the matrice
       (setf (slot-value finspace 'stong) (stongmat (top finspace)))))


(DEFMETHOD SLOT-UNBOUND (class (finspace finite-space) (slot-name (eql 'orgn)))
  (declare (ignore class))
  (the list
       (setf (slot-value finspace 'orgn) `(FINITE-SPACE ,(idnm finspace))))) ;(1+ *idnm-counter*)


(DEFMETHOD INITIALIZE-INSTANCE :after ((finspace FINITE-SPACE) &rest rest)
  (set (intern (format nil "K~D" (idnm finspace))) finspace))


(DEFVAR *FINITE-SPACE-LIST*)
(SETF *FINITE-SPACE-LIST* +empty-list+)
(PUSHNEW '*FINITE-SPACE-LIST* *list-list*)


(DEFMETHOD PRINT-OBJECT ((finspace finite-space) stream)
 (the finite-space
   (progn
      (format stream "[K~D Finite-Space]" (idnm finspace))
      finspace)))


(DEFUN BUILD-FINITE-SPACE (&key top stong heights orgn)
  (the finite-space
       (progn
         (let ((already (find orgn *finite-space-list* :test #'equal :key #'orgn)))
           (declare (type (or finite-space null) already))
           (when already
             (return-from build-finite-space already)))
         (let ((finspace (make-instance 'finite-space)))
           (declare (type finite-space finspace))
           (if top     (setf (slot-value finspace 'top) top))
           (if stong   (setf (slot-value finspace 'stong) stong))
           (if heights (setf (slot-value finspace 'heights) heights))
           (if orgn    (setf (slot-value finspace 'orgn) orgn))
           (push finspace *finite-space-list*)
           finspace))))


(DEFUN TOPOLOGY (dmns dens)
  #| Random Finite-Space of size 'dmns' and density 'dens' |#
  (declare (fixnum dmns))
  (unless (plusp dmns)
    (error "In TOPOLOGY, the dimension ~D should be positive." dmns))
  (build-finite-space :top (randomtop dmns dens)
                      :orgn `(RANDOM-TOPOLOGY ,(1+ *idnm-counter*))))


(DEFUN ADMISSIBLE (topogenous i j)
  #| 'T' if Hat-Uj-{i} is contractible (the edge (i j) is admissible) |#
  (if (and (> j i) (eq (terme topogenous i j) 1))
      (the boolean
           (eq (length (core_list topogenous
                                  (loop for k from 1 to (1- j)
                                        unless (or (zerop (terme topogenous k j)) (eq k i))
                                        collect k))) 1))))


(DEFMETHOD ADMISSIBLE-P ((topogenous matrice)) 
  (the boolean
       (let ((fin T)
             (stong (stongmat topogenous)))
         (do ((j 2 (1+ j))) ((or (null fin) (> j (ncol topogenous))) fin)
           (do ((i 1 (1+ i))) ((or (null fin) (> i (1- j))))
             (if (= (terme stong i j) 1)
                 (if (null (admissible topogenous i j))
                     (setf fin NIL))))))))


(DEFMETHOD ADMISSIBLE-P ((finspace finite-space))
  (the boolean
       (let ((fin T))
         (do ((j 2 (1+ j))) ((or (null fin) (> j (ncol (top finspace)))) fin)
           (do ((i 1 (1+ i))) ((or (null fin) (> i (1- j))))
             (if (= (terme (stong finspace) i j) 1)
                 (if (null (admissible (top finspace) i j))
                     (setf fin NIL))))))))


(DEFMETHOD DVFIELD ((finspace finite-space))
  #| Discrete vector field on 'finspace' |#
  (declare (finite-space finspace))
  (let* ((m_stong (matrice-to-lmtrx (nilpot (stong finspace))))
         (m_top (matrice-to-lmtrx (nilpot (top finspace))))
         (rown (length m_stong))
         (status (make-array rown :initial-element '()))
         (used nil) 
         (vf '()))
    (dotimes (i rown) ; Here it is working on the columns. 'rown' is the dimension of the space
      (unless (member i used)
        (let ((coli_stong (svref m_stong i))
              (coli_top (svref m_top i)))
          (when coli_stong
            (let ((rowilist (mapcar #'car coli_stong))
                  (rowilist_top (mapcar #'car coli_top)))
              (declare (type list rowilist rowilist_top))
              (dolist (rowi rowilist)    
                (declare (type fixnum rowi))
                (let ((rowistatus (svref status rowi)))
                  (unless (eql 1 rowistatus)
                    (unless (member rowi used)
                      (when (and (admissible (top finspace) (1+ rowi) (1+ i))
                                 (null-intersection-p rowistatus rowilist))
                        (setf (svref status rowi) 1)
                        (setf rowistatus (insert rowi rowistatus))
                        (push (list (1+ rowi) (1+ i)) vf)
                        (push i used)    
                        (push rowi used)    
                        (dolist (item rowilist_top)
                          (unless (= item rowi)
                            (case (svref status item)
                              (1 (dotimes (rowj rown)
                                   (let ((rowjstatus (svref status rowj)))
                                     (unless (eql 1 rowjstatus)
                                       (when (member item rowjstatus)
                                         (setf (svref status rowj)
                                               (union-merge rowistatus rowjstatus)))))))
                              (otherwise
                               (setf (svref status item)
                                     (union-merge (member item rowilist_top) (union-merge rowistatus (svref status item))))))))
                        (return)))))))))))
    (return-from DVFIELD (nreverse vf))))


(DEFMETHOD DVFIELD ((topogenous matrice))
  #| Discrete vector field on a matrice 'topogenous' |#
  (dvfield (make-instance 'finite-space
                          :top topogenous)))


;;
;;  Barycentric subdivision
;;


(DEFUN RELAC (topogenous)
  #| Chains of length 2 |#
  (declare (type matrice topogenous))
  (let* ((dim (nlig topogenous))
         (rslt (creer-matrice dim dim)))
    (declare (type fixnum dim))   
    (do ((i 1 (1+ i))) ((> i dim) rslt)
      (do ((j (1+ i) (1+ j))) ((> j dim))
        (if (eq (terme topogenous i j) 1)
            (insert-term rslt i j (list (list i j))))))))


(DEFUN BOUNDARYOPERATORS (topogenous)
  (let* ((relac (relac topogenous))
         (diferenciales NIL)
         (nuevo relac)
         (dim (nlig topogenous))
         (parada 0)
         (dim_n 0)
         (Cn NIL))

    (do ((i 1 (1+ i))) ((> i (1- dim)))
      (do ((j (1+ i) (1+ j))) ((> j dim))
        (unless (eq (terme topogenous i j) 0)
          (setf Cn (append Cn (list (list i j)))))))

    (setf dim_n (length Cn))
    
    (push (let ((D1 (creer-matrice dim dim_n)))
            (dotimes (i dim)
              (dotimes (j dim_n)
                (if (find (1+ i) (nth j Cn))
                    (insert-term D1 (1+ i)  (1+ j) 1))))
            D1) diferenciales)
            
    (if (> dim 1)
        (do ((long 2 (1+ long))) ((or (> long dim) (eq parada 1)))
          (let ((dim_n-1 dim_n)
                (Cn-1 Cn)
                (anterior nuevo))
            
            (setf dim_n 0 Cn NIL nuevo (creer-matrice dim dim))
            
            (do ((i 1 (1+ i))) ((> i (- (1+ dim) long)))
              (do ((j (1+ i) (1+ j))) ((> j dim))
                (let ((cadenasij NIL))
                  (do ((r i (1+ r))) ((> r j))
                    (let ((termeir (terme anterior i r)))
                      (unless (or (eq termeir 0) (eq (terme relac r j) 0))
                        (setf cadenasij (append cadenasij (mapcar #'(lambda (x) (append x (list j))) termeir))))))
                  (if cadenasij
                      (progn (insert-term nuevo i j cadenasij)
                        (setf Cn (append Cn cadenasij)))))))
            
            (setf dim_n (length Cn))
            
            (if (eq dim_n 0)
                (setf parada 1)
              (push (let ((Dn (creer-matrice dim_n-1 dim_n)))
                      (dotimes (i dim_n-1)
                        (dotimes (j dim_n)
                          (if (subsetp (nth i Cn-1) (nth j Cn))
                              (insert-term Dn (1+ i) (1+ j) 1))))
                      Dn) diferenciales)))))
    (return-from BOUNDARYOPERATORS (reverse diferenciales))))


(DEFUN BLOCKS (list)
  (let ((rslt (identite (let ((suma (nlig (car list))))
                          (dolist (x list)
                            (setf suma (+ suma (ncol x))))
                          suma))))
    (do* ((lst list (cdr lst)) 
          (kfilas 0 (if (null lst) 0 (+ kfilas (nlig dif))))
          (kcolumnas (nlig (car list)) (if (null lst) 0 (+ kcolumnas (ncol dif))))
          (dif (car lst) (car lst))) ((endp lst))
      
      (do ((i 1 (1+ i))) ((> i (nlig dif)))
        (do ((j 1 (1+ j))) ((> j (ncol dif)))
          (unless (eq (terme dif i j) 0)
            (insert-term rslt (+ i kfilas) (+ j kcolumnas) 1)))))
      rslt))


(DEFMETHOD BAR_SUBDIVISION ((topogenous matrice))
  #| Topogenous matrice of the barycentric subdivision of 'topogenous' |#
  (topmat (blocks (boundaryoperators topogenous))))


(DEFMETHOD BAR_SUBDIVISION ((finspace finite-space))
  #| Finite-Space whose 'top' is the barycentric subdivision of (top 'finspace') |#
  (let ((already (find `(BARYCENTRIC-SUBDIVISION ,finspace) *finite-space-list* :test #'equal :key #'orgn)))
    (declare (type (or finite-space null) already))
    (if already
        already
      (build-finite-space :orgn `(BARYCENTRIC-SUBDIVISION ,finspace)
                          :stong (blocks (boundaryoperators (top finspace)))))))


;;
;;  Special examples of Finite topological spaces
;;


(DEFUN FINITE-MODEL-SPHERE (dmns)
  #| Minimal finite model of the 'dmns'-dimensional sphere |#
  (declare (fixnum dmns))
  (if (minusp dmns)
      (error "In FINITE-MODEL-SPHERE, the dimension ~D must be non-negative integer." dmns))
  (build-finite-space
   :stong (let ((rslt (identite (+ (* 2 dmns) 2))))
            (unless (= dmns 0)
              (do ((i 1 (1+ i))) ((> i (* 2 dmns)))
                (insert-term rslt i (+ i 2) 1)
                (if (= (mod i 2) 0)
                    (insert-term rslt i (+ i 1) 1)
                  (insert-term rslt i (+ i 3) 1))))
            rslt)
   :heights (let ((rslt (<a-b> 1 (+ (* 2 dmns) 2))))
              (dotimes (q (1+ dmns)) 
                (setf rslt (append (cddr rslt) (list (list (pop rslt) (pop rslt))))))
              rslt)
   :orgn `(FINITE-MODEL-SPHERE ,dmns)))


(DEFUN FINITE-MODEL-RP2 ()
  #| Minimal finite model of the 2-dimensional real projective space |#
  (build-finite-space
       :top (convertarray
             #2A((1 0 0 1 1 1 1 0 0 1 1 1 1)
                 (0 1 0 1 1 0 0 1 1 1 1 1 1)
                 (0 0 1 0 0 1 1 1 1 1 1 1 1)
                 (0 0 0 1 0 0 0 0 0 1 0 0 1)
                 (0 0 0 0 1 0 0 0 0 0 1 1 0)
                 (0 0 0 0 0 1 0 0 0 1 1 0 0)
                 (0 0 0 0 0 0 1 0 0 0 0 1 1)
                 (0 0 0 0 0 0 0 1 0 1 0 1 0)
                 (0 0 0 0 0 0 0 0 1 0 1 0 1)
                 (0 0 0 0 0 0 0 0 0 1 0 0 0)
                 (0 0 0 0 0 0 0 0 0 0 1 0 0)
                 (0 0 0 0 0 0 0 0 0 0 0 1 0)
                 (0 0 0 0 0 0 0 0 0 0 0 0 1)))
       :orgn '(FINITE-MODEL-RP2)))


(DEFUN FINITE-MODEL-KLEIN-BOTTLE ()
  #| Minimal finite model of the Klein bottle space |#
  (build-finite-space
       :top (convertarray
             #2A((1 0 0 0 1 0 0 1 1 0 1 0 1 1 1 1)
                 (0 1 0 0 0 1 0 1 1 0 0 1 1 1 1 1)
                 (0 0 1 0 1 0 1 0 0 1 1 0 1 1 1 1)
                 (0 0 0 1 0 1 1 0 0 1 0 1 1 1 1 1)
                 (0 0 0 0 1 0 0 0 0 0 0 0 1 1 0 0)
                 (0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0)
                 (0 0 0 0 0 0 1 0 0 0 0 0 1 0 1 0)
                 (0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 1)
                 (0 0 0 0 0 0 0 0 1 0 0 0 0 1 1 0)
                 (0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 1)
                 (0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1)
                 (0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1)
                 (0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0)
                 (0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0)
                 (0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0)
                 (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1)))
       :orgn '(FINITE-MODEL-KLEIN-BOTTLE)))


(DEFMETHOD FINITE-SUSPENSION ((topogenous matrice))
  #| Non-Hausdorff suspension of 'topogenous' |#
  (let* ((nlig+1 (1+ (nlig topogenous)))
         (rslt (add-points topogenous nlig+1 2)))
    (maj-colonne rslt nlig+1 (loop for x in (<a-b> 1 nlig+1)
                                   collect (list x 1)))
    (maj-colonne rslt (1+ nlig+1) (loop for x in (<a-b> 1 (1- nlig+1))
                                        collect (list x 1)))
    (insert-term rslt (1+ nlig+1) (1+ nlig+1) 1)
    rslt))


(DEFMETHOD FINITE-SUSPENSION ((finspace finite-space))
  (build-finite-space :top (finite-suspension (top finspace))
                      :orgn `(FINITE-SUSPENSION ,(idnm finspace))))


(DEFUN FISHER-YATES (list)
  (let ((rslt list))
  (do ((i (1- (length list)) (1- i))) ((< i 0) rslt)
    (let ((j (random (1+ i)))
          (ai (nth i rslt)))
      (setf (nth i rslt) (nth j rslt)
            (nth j rslt) ai)))))


(DEFUN WEDGE-SPHERES (num-minimals num-spheres)
  #| Topogenous matrix for an h-regular finite model of a wedge of 'num-spheres' with 'num-minimals' minimals |#
  (unless (> num-spheres 0)
    (error "Second parameter must be a positive integer"))
  (let ((minimals (fisher-yates (<a-b> 1 num-minimals)))
        (maximals (fisher-yates (<a-b> (1+ num-minimals) (1- (* 2 num-minimals)))))
        (rslt (identite (+ num-minimals num-minimals num-spheres -1))))

    (do ((k 0 (1+ k))) ((> k (- num-minimals 2)))
      (insert-term rslt (nth k minimals) (nth k maximals) 1)
      (insert-term rslt (nth (1+ k) minimals) (nth k maximals) 1)) ; Constructing a fence to ensure connectedness

    (insert-term rslt (car minimals) (* 2 num-minimals) 1) ; Putting the first 'extra' maximal point
    (insert-term rslt (nth (1+ (random (1- num-minimals))) minimals) (* 2 num-minimals) 1)
    
    (unless (eq num-spheres 1)
      (insert-term rslt (car (last minimals)) (+ (* 2 num-minimals) num-spheres -1) 1)
      (insert-term rslt (nth (random (- num-minimals 2)) minimals) (+ (* 2 num-minimals) num-spheres -1) 1) ; Avoiding beat points
      (do ((i (- num-spheres 2) (1- i))) ((zerop i))
        (let ((col-index (+ (* 2 num-minimals) i))
              (u1 (random num-minimals))
              (u2 (random num-minimals)))
          (if (eq u1 u2)
              (if (zerop u1)
                  (setf u1 (1+ u1))
                (setf u1 (1- u1))))
          (insert-term rslt (nth u1 minimals) col-index 1)
          (insert-term rslt (nth u2 minimals) col-index 1))))
    rslt))


(DEFUN WEDGE-FINITE-SPHERES (num-minimals num-spheres)
  (build-finite-space 
       :top (wedge-spheres num-minimals num-spheres)
       :orgn `(WEDGE-FINITE-SPHERES ,num-minimals ,num-spheres)))


;;
;;   Point reductions: beat points and weak beat points
;;


(DEFUN DOWNBP (topogenous n &optional (list (<a-b> 1 (nlig topogenous))))
  #| T if Xn is a down beat point of the submatrix 'list'x'list' |#
  (the boolean
       (let ((maximal (loop for k in (reverse (subseq list 0 (position n list)))
                            thereis (and (eq (terme topogenous k n) 1) k))))
         
         (unless (or (null maximal) 
                     (loop for i in (subseq list 0 (position n list))
                           thereis (not (eq (terme topogenous i maximal) (terme topogenous i n)))))
           +TRUE+))))


(DEFUN UPBP (topogenous n &optional (list (<a-b> 1 (nlig topogenous))))
  #| T if Xn is an up beat point of the submatrix 'list'x'list' |#
  (the boolean
       (let ((minimal (loop for k in (subseq list (1+ (position n list)))
                            thereis (and (eq (terme topogenous n k) 1) k))))
         (unless (or (null minimal) 
                     (loop for i in (subseq list (1+ (position n list)))
                           thereis (not (eq (terme topogenous minimal i) (terme topogenous n i)))))
           +TRUE+))))


(DEFUN BEATPOINT (topogenous n &optional (list (<a-b> 1 (nlig topogenous)))) 
  #| T if Xn is a beat point of the submatrix 'list'x'list' |#
  (the boolean
       (if (< (1+ n) (/ (ncol topogenous) 2))
           (if (downbp topogenous n list)
               +TRUE+
             (upbp topogenous n list))
         (if (upbp topogenous n list)
             +TRUE+
           (downbp topogenous n list)))))


(DEFUN CORE_LIST (topogenous &optional (list (<a-b> 1 (nlig topogenous))))
  #| List of the elements of the core of the submatrix 'list'x'list' |#
  (let ((eliminar NIL))
    (setf eliminar (loop for n in list
                         thereis (and (beatpoint topogenous n list) n)))
    (if (null eliminar)
        (return-from core_list list)
      (core_list topogenous (remove eliminar list)))))


(DEFUN Un-N (topogenous n &optional (list (<a-b> 1 (nlig topogenous))))
  #| List of the elements in Un-{Xn} |#
  (loop for i in (subseq list 0 (position n list))
       when (eq (terme topogenous i n) 1)
       collect i))


(DEFUN Fn-N (topogenous n &optional (list (<a-b> 1 (nlig topogenous))))
  #| List of the elements in Fn-{Xn} |#
  (loop for j in (subseq list (1+ (position n list)))
       when (eq (terme topogenous n j) 1)
       collect j))


(DEFUN LINK_LIST (topogenous n &optional (list (<a-b> 1 (nlig topogenous))))
  #| List of the elements in the link of Xn |#
  (append (Un-N topogenous n list) (Fn-N topogenous n list)))


(DEFUN WEAKPOINT (topogenous n &optional (list (<a-b> 1 (nlig topogenous))))
  #| T if Xn is a weak beat point of the submatrix 'list'x'list' |#
  (the boolean
       (or (eq (length (core_list topogenous (Un-N topogenous n list))) 1)
           (eq (length (core_list topogenous (Fn-N topogenous n list))) 1))))
  
  
(DEFUN WEAKCORE_LIST (topogenous &optional (list (<a-b> 1 (nlig topogenous))))
  #| List of the elements of a "weak core" of the submatrix 'list'x'list' |#
  (let ((eliminar NIL))
    (setf eliminar (loop for n in list
                         thereis (and (or (beatpoint topogenous n list) (weakpoint topogenous n list)) n)))
    (if (null eliminar)
        (return-from weakcore_list list)
      (weakcore_list topogenous (remove eliminar list)))))


(DEFMETHOD CORE ((topogenous matrice) &optional (list (<a-b> 1 (nlig topogenous))))
  #| Matrice of the core of 'topogenous' |#
  (the matrice
       (let ((corelist (core_list topogenous list)))
         (sub-matrix topogenous corelist corelist))))


(DEFMETHOD CORE ((finspace finite-space) &optional list)
  #| Finite-Space whose 'top' is the topogenous matrix of the core of (top 'finspace') |#
  (let ((already (find (if list `(CORE ,finspace ,list) `(CORE ,finspace)) *finite-space-list* :test #'equal :key #'orgn)))
    (declare (type (or finite-space null) already))
    (if already
        already
      (let ((list2 (or list (<a-b> 1 (nlig (top finspace))))))
        (build-finite-space :orgn (if list `(CORE ,finspace ,list) `(CORE ,finspace))
                            :top (core (top finspace) list2))))))


(DEFMETHOD WEAKCORE ((topogenous matrice) &optional (list (<a-b> 1 (nlig topogenous))))
  #| Matrice of a "weak core" of 'topogenous' |#
  (the matrice
       (let ((weakcorelist (weakcore_list topogenous (core_list topogenous list))))
         (sub-matrix topogenous weakcorelist weakcorelist))))


(DEFMETHOD WEAKCORE ((finspace finite-space) &optional list)
  #| Finite-Space whose 'top' is the topogenous matrix of a "weak core" of (top 'finspace') |#
  (let ((already (find (if list `(WEAK-CORE ,finspace ,list) `(WEAK-CORE ,finspace)) *finite-space-list* :test #'equal :key #'orgn)))
    (declare (type (or finite-space null) already))
    (if already
        already
      (let ((list2 (or list (<a-b> 1 (nlig (top finspace))))))
        (build-finite-space :orgn (if list `(WEAK-CORE ,finspace ,list) `(WEAK-CORE ,finspace))
                            :top (weakcore (top finspace) list2))))))


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


