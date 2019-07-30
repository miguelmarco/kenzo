;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*

;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 
;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 
;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 

(IN-PACKAGE #:cat)

(PROVIDE "finite-topological-spaces-class")

;;
;;  Finite topological spaces
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


(DEFUN CONVERTMATRICE (matrice)
  #| Convert a 'matrice' to an 'array' |#
  (let* ((numfil (nlig matrice))
         (numcol (ncol matrice))
         (rslt (make-array (list numfil numcol) :initial-element 0)))
    (dotimes (i numfil)
      (dotimes (j numcol)
        (let ((Mij (terme matrice (1+ i) (1+ j))))
          (unless (zerop Mij)
            (setf (aref rslt i j) Mij)))))
    rslt))


(DEFUN VER (matrix)
  (let ((mat '()))
    (typecase matrix
      (matrice (setf mat matrix))
      (T (setf mat (convertarray matrix))))
    (format t "~%")
    (format t "           ========== MATRIX ~A row(s) + ~A column(s) ==========~%~%" (nlig mat) (ncol mat))
    (do ((i 1 (1+ i))) ((> i (nlig mat)))
      (format t "           ")
      (do ((j 1 (1+ j))) ((> j (ncol mat)))
        (if (< (terme mat i j) 0)
            (format t " ~A " (terme mat i j))
          (format t "  ~A " (terme mat i j))))
      (format t "  ~%"))))


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
