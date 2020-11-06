
;;  FINITE-SPACES FINITE-SPACES FINITE-SPACES FINITE-SPACES FINITE-SPACES 
;;  FINITE-SPACES FINITE-SPACES FINITE-SPACES FINITE-SPACES FINITE-SPACES 
;;  FINITE-SPACES FINITE-SPACES FINITE-SPACES FINITE-SPACES FINITE-SPACES 

(IN-PACKAGE #:cat)

(provide "finite-spaces-class")

;;
;;  Representation of Finite Topological Spaces
;;


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFUN CONVERTARRAY (array)
  #| Convert an 'array' to a 'matrice' |#
  (let* ((dims (array-dimensions array))
         (nligs (first dims))
         (ncols (second dims))
         (rslt (creer-matrice nligs ncols))
         (liste '()))
    (if (zerop ncols)
        (return-from convertarray rslt)
      (do ((i nligs (1- i)))
          ((< i 1))
        (let ((ligi '()))
          (do ((j ncols (1- j)))
              ((< j 1))
            (let ((array-ij (aref array (1- i) (1- j))))
              (declare (type fixnum array-ij))
              (unless (zerop array-ij)
                (push (list j array-ij) ligi))))
          (if ligi
              (push (cons i (list ligi)) liste)))))
      (maj-matrice rslt liste)))


(DEFUN CONVERTMATRICE (matrice)
  #| Convert a 'matrice' to an 'array' |#
  (let* ((nligs (nlig matrice))
         (ncols (ncol matrice))
         (rslt (make-array (list nligs ncols) :initial-element 0)))
    (if (zerop nligs)
        (return-from convertmatrice rslt)
      (dotimes (j ncols rslt)
        (let ((ptc (basecol matrice (1+ j))))
          (do ((pc (up ptc) (up pc)))
              ((eq pc ptc))
            (setf (aref rslt (1- (ilig pc)) j)  (val pc))))))))


#|
  (setf mat (mat-aleat 100 130 .5 10))
  (setf mat2 (convertarray (convertmatrice mat)))
  (equalmatrix mat mat2)
|#


(DEFUN BINARYMATRICE-TO-UBASIS (matrice)
  #| Simplified version of 'matrice-to-lmtrx' function |#
  (let* ((ncols (ncol matrice))
         (rslt (make-array ncols)))
    (dotimes (j ncols)
      (let ((Uj (let ((ptc (basecol matrice (1+ j)))
                      (res '()))
                  (do ((pc (up ptc) (up pc)))
                      ((eq pc ptc))
                    (push (ilig pc) res))
                  res)))
        (setf (svref rslt j) Uj)))
    rslt))


(DEFUN BINARYMATRICE-TO-FBASIS (matrice)
  #| Similar to 'binarymatrice-to-fbasis' but running over rows |#
  (let* ((nligs (nlig matrice))
         (rslt (make-array nligs)))
    (dotimes (i nligs)
      (let ((Fi (let ((ptl (baselig matrice (1+ i)))
                      (res '()))
                  (do ((pl (left ptl) (left pl)))
                      ((eq pl ptl))
                    (push (icol pl) res))
                  res)))
        (setf (svref rslt i) Fi)))
    rslt))


(DEFUN UBASIS-TO-BINARYMATRICE (ubasis)
  #| Create the upper triangular binary matrice associated to 'ubasis' |#
  (let* ((dimension (length ubasis))
         (rslt (creer-matrice dimension dimension)))
    (do ((j 1 (1+ j)))
        ((> j dimension) rslt)
      (dolist (i (svref ubasis (1- j)))
        (inserer-terme (baselig rslt i) (basecol rslt j) 1)))))


(DEFUN FBASIS-TO-BINARYMATRICE (fbasis)
  #| Create the upper triangular binary matrice associated to the opposite 'fbasis' |#
  (let* ((dimension (length fbasis))
         (rslt (creer-matrice dimension dimension)))
    (do ((i 1 (1+ i)))
        ((> i dimension) rslt)
      (dolist (j (svref fbasis (1- i)))
        (inserer-terme (baselig rslt i) (basecol rslt j) 1)))))


(DEFUN EDGES-TO-STONG-MTRX (edges)
  (nilpot-1 (ubasis-to-binarymatrice edges)))


(DEFUN EDGES-TO-STONG (dim edges)
  (let ((rslt (creer-matrice dim dim)))
    (dolist (pair edges (nilpot-1 rslt))
      (inserer-terme (baselig rslt (car pair)) (basecol rslt (cadr pair)) 1))))


(DEFUN EDGES-TO-MATRICE (dim edges)
 (let ((slambda (gensym)))
    (let ((slambda edges))
      (edges-to-stong dim slambda))))


(DEFUN STONG-TO-EDGES (stong)
  (let ((vector (binarymatrice-to-ubasis (nilpot stong)))
        (dim (ncol stong))
        (rslt '()))
    (dotimes (j dim rslt)
      (dolist (i (aref vector j))
        (push (list i (1+ j)) rslt)))))


(DEFUN LIST-TO-VECTOR (list)
  (map 'vector #'identity list))


(DEFUN VECTOR-TO-LIST (vector)
  (loop for elto across vector collect elto))


(DEFUN REMOVE-FIRST-OCCURRENCE (elto lst)
  (cond ((null lst) nil)
        ((equal (car lst) elto) (cdr lst))
        (T (cons (car lst)
                 (remove-first-occurrence elto (cdr lst))))))


(DEFUN FISHER-YATES (list)
  #| Shuffle list |#
  (let ((rslt list))
    (do ((i (1- (length list)) (1- i)))
        ((< i 0) rslt)
      (let ((j (random (1+ i)))
            (ai (nth i rslt)))
        (setf (nth i rslt) (nth j rslt)
              (nth j rslt) ai)))))
              
(DEFUN LIST-DIFFERENCE (l1 l2)
  #| Remove the elements of l1 from l2 |#
  (remove-if #'(lambda (elem) (member elem l2)) l1))


#|
  (setf example (randomtop 8 .4))
  (setf ubasis (binarymatrice-to-ubasis example))
  (setf example2 (ubasis-to-binarymatrice ubasis))
  (setf fbasis (binarymatrice-to-fbasis example2))
  (setf example3 (fbasis-to-binarymatrice fbasis))
  (and (equalmatrix example example2) (equalmatrix example example3)) ; Always T
|#


(DEFUN SHOW (mtrx)
  (let (mat)
    (typecase mtrx
      (matrice (setf mat (convertmatrice mtrx)))
      ((or matrix array) (setf mat mtrx)))
    (let ((dims (array-dimensions mat)))
      (format t "~%")
      (format t "           ========== MATRIX ~A row(s) + ~A column(s) ==========~%~%" (first dims) (second dims))
      (dotimes (i (first dims))
        (format t "           ")
        (dotimes (j (second dims))
          (let ((term-ij (aref mat i j)))
            (if (< term-ij 0)
                (format t " ~A " term-ij)
              (format t "  ~A " term-ij))))
        (format t "  ~%"))))
  (values))


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


(DEFUN CONSECUTIVE-COLS-SUBMATRIX (mtrx rows consecutive-cols)
  #| Extract the submatrix of 'mtrx' formed by the rows whose indexes are in 'rows' and the columns whose indexes are in 'consecutive-cols', a list of consecutive indexes in ascending order |#
  (let ((rslt (creer-matrice (length rows) (length consecutive-cols))))
    (do ((il 0 (1+ il)))
        ((= il (length rows)) rslt)
      (declare (type fixnum il))
      (let ((ilig (nth il rows)))
        (do ((pl01 (baselig mtrx ilig))
             (p1 (left (baselig mtrx ilig)) (left p1))
             (pl0r (baselig rslt (1+ il)) (left pl0r)))
            ((eq p1 pl01))
          (let ((j (position (icol p1) consecutive-cols)))
            (if j
                (inserer-terme pl0r (basecol rslt (1+ j)) (val p1)))))))))


(DEFUN SUBMATRIX-COLS (mtrx cols)
  #| Extract the submatrix of 'mtrx' formed by the columns whose indexes are in 'cols' |#
  (let ((rslt (creer-matrice (nlig mtrx) (length cols))))
    (do ((il 1 (1+ il)))
        ((> il (nlig mtrx)) rslt)
      (declare (type fixnum il))
      (let ((pr (baselig rslt il)))
        (do ((pl01 (baselig mtrx il))
             (p1 (left (baselig mtrx il)) (left p1)))
            ((eq p1 pl01))
          (let ((j (position (icol p1) cols)))
            (if j
                (maj-terme (chercher-hor pr (1+ j))
                           (chercher-ver (basecol rslt (1+ j)) il)
                           (val p1)))))))))


(DEFUN SUBMATRIX (mtrx rows cols)
  #| Extract the submatrix of 'mtrx' formed by the rows whose indexes are in 'rows' and the columns whose indexes are in 'cols' |#
  (let ((rslt (creer-matrice (length rows) (length cols))))
    (do ((il 0 (1+ il)))
        ((= il (length rows)) rslt)
      (declare (type fixnum il))
      (let ((ilig (nth il rows)))
        (do ((pl01 (baselig mtrx ilig))
             (p1 (left (baselig mtrx ilig)) (left p1)))
            ((eq p1 pl01))
          (let ((j (position (icol p1) cols)))
            (if j
                (maj-terme (chercher-hor (baselig rslt (1+ il)) (1+ j))
                           (chercher-ver (basecol rslt (1+ j)) (1+ il))
                           (val p1)))))))))


#|
  (setf example (mat-aleat 300 300 .05 8))
  (setf rows (fisher-yates (<a-b> 23 257)))
  (setf cols1 (<a-b> 23 257))
  (setf cols (fisher-yates (<a-b> 34 298)))
  (setf sub (submatrix example rows cols))
  (setf sub1 (consecutive-cols-submatrix example rows cols1))
  (setf sub2 (submatrix example rows cols1))
  (equalmatrix sub1 sub2)
|#


(DEFUN INSERT-TERM (mat ilig icol val)
  #| Insert the value 'val' in the ('ilig' 'icol')-entry of 'mat' ('majterme' is in 'finite-spaces-changes') |#
  (majterme (chercher-hor (baselig mat ilig) icol)
             (chercher-ver (basecol mat icol) ilig)
             val))


(DEFUN NILPOT (mat)
  #| Put zeros in the diagonal |#
  (declare (type matrice mat))
  (let ((rslt (copier-matrice mat)))
    (loop for k from 1 to (nlig rslt)
          do (insert-term rslt k k 0))
    rslt))


(DEFUN NILPOT-1 (mat)
  #| Put ones in the diagonal |#
  (declare (type matrice mat))
  (let ((rslt (copier-matrice mat)))
    (loop for k from 1 to (nlig rslt)
          do (insert-term rslt k k 1))
    rslt))


(DEFUN ADD-ROWS (matrice pos num)
  #| Add 'num' rows from position 'pos' |#
  (let ((rslt (creer-matrice (+ (nlig matrice) num) (ncol matrice))))
    (do ((i 1 (1+ i)))
        ((eq i pos))
      (do ((ptl1 (baselig matrice i))
           (pl1 (left (baselig matrice i)) (left pl1))
           (pl2 (baselig rslt i) (left pl2)))
          ((eq pl1 ptl1))
        (inserer-terme pl2 (basecol rslt (icol pl1)) (val pl1))))
    
    (do ((i pos (1+ i)))
        ((> i (nlig matrice)) rslt)
      (do ((ptl1 (baselig matrice i))
           (pl1 (left (baselig matrice i)) (left pl1))
           (pl2 (baselig rslt (+ i num)) (left pl2)))
          ((eq pl1 ptl1))
        (inserer-terme pl2 (basecol rslt (icol pl1)) (val pl1))))))


(DEFUN ADD-COLS (matrice pos num)
  #| Add 'num' columns from position 'pos' |#
  (let ((rslt (creer-matrice (nlig matrice) (+ (ncol matrice) num))))
    (do ((i 1 (1+ i)))
        ((eq i pos))
      (do ((ptc1 (basecol matrice i))
           (pc1 (up (basecol matrice i)) (up pc1))
           (pc2 (basecol rslt i) (up pc2)))
          ((eq pc1 ptc1))
        (inserer-terme (baselig rslt (ilig pc1)) pc2 (val pc1))))
    
    (do ((i pos (1+ i)))
        ((> i (ncol matrice)) rslt)
      (do ((ptc1 (basecol matrice i))
           (pc1 (up (basecol matrice i)) (up pc1))
           (pc2 (basecol rslt (+ i num)) (up pc2)))
          ((eq pc1 ptc1))
        (inserer-terme (baselig rslt (ilig pc1)) pc2 (val pc1))))))


(DEFUN ADD-POINTS (matrice pos num)
  (add-cols (add-rows matrice pos num) pos num))


#|
  (setf m (randomtop 15 .5))
  (setf rows '(1 5 7 2 4))
  (setf cols '(10 1 5 2 15 13 8))
  (show m)
  (show (submatrix m rows cols))
|#


(DEFUN STONGMAT (topogenous)
  #| Stong matrice from the topogenous matrice 'topogenous' |#
  (let* ((dim (nlig topogenous))
         (rslt (identite dim)))
    (declare (type fixnum dim))
    (do ((i 1 (1+ i))) ((> i dim) rslt)
      (do ((j (1+ i) (1+ j))) ((> j dim))
        (unless (zerop (terme topogenous i j))
          (do ((k (1+ i) (1+ k))) ((> k (1- j)) (insert-term rslt i j 1))
            (if (and (eq (terme topogenous i k) 1) (eq (terme topogenous k j) 1))
                (return))))))))


(DEFUN TOPMAT (stong)
 #| Topogenous matrice from the Stong matrice 'stong' |#
  (let ((rslt (binarymatrice-to-ubasis (nilpot stong))))
    (dotimes (j (length rslt))
      (let ((colj (svref rslt j)))
        (dolist (i colj)
          (setf (svref rslt j) (unionmerge (svref rslt (1- i)) (svref rslt j))))))
    (nilpot-1 (ubasis-to-binarymatrice rslt))))


(DEFUN RANDOMTOP (dim dens)
  #| A 'dim' x 'dim' random topogenous matrice with density 'dens' of the number of ones|#
  (let ((rslt (identite dim)))
    (do ((i 1 (1+ i)))
        ((equal i dim))
       (do ((j (1+ i) (1+ j)))
           ((> j dim))
           (if (< (random 1.0) dens)
           	(insert-term rslt i j 1))))
  (topmat rslt)))


(DEFUN RANDOMSTONG (dim dens)
  #| A 'dim' x 'dim' random Stong matrice |#
  (stongmat (randomtop dim dens)))


#|
  ((lambda (n dens)
     (let ((M (randomtop n dens)))                         
       (equalmatrix M (topmat (stongmat M))))) 8 .5) ; Always T
  
  ((lambda (n dens)
     (let ((M (randomstong n dens)))                       
       (equalmatrix M (stongmat (topmat M))))) 8 .7) ; Always T
  
  ((lambda (n dens)
     (do ((q 1 (1+ q))) ((> q 30))
       (let ((M (randomtop n dens))
             (P (randomstong n dens)))
         (unless (and (equalmatrix M (topmat (stongmat M))) (equalmatrix P (stongmat (topmat P))))
           (print 'Wrong))))) 10 .5)                        ; Must never print 'Wrong'
|#


(DEFUN HEIGHTS-AUX (ubasis-stong lst) 
  (if (null lst)
      NIL
    (let ((minimals (loop for n in lst
                          if (null (intersection lst (svref ubasis-stong (1- n))))
                          collect n)))
      (cons minimals (HEIGHTS-AUX ubasis-stong (reverse (set-difference lst minimals)))))))


(DEFMETHOD HEIGHTS ((stong matrice))
  #| List of elements separated by heights |#
  (HEIGHTS-AUX (binarymatrice-to-ubasis (nilpot stong)) (<a-b> 1 (ncol stong))))


(DEFUN ELEMENTS-HEIGHT (finspace)
  #| Vector whose k-th entry equals the height of the element k+1 |#
  (let ((rslt (make-array (cardinality finspace)))
        (heights (heights finspace)))
    (dotimes (k (length heights) rslt)
      (let ((heights_k (nth k heights)))
        (dolist (elto heights_k)
          (setf (svref rslt (1- elto)) k))))))


#|  
  (heights (randomstong 1 .5))
  (heights (randomstong 8 .5))
  (heights (randomstong 10 .2))
  (heights (randomstong 15 .9))
|#


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


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
       (setf (slot-value finspace 'orgn) `(FINITE-SPACE ,(idnm finspace))))) ; (1+ *idnm-counter*)


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
  (declare
   (type matrice top)
   (type matrice stong)
   (type list orgn)
   (type list heights))
  (the finite-space
       (progn
         (let ((already (find orgn *finite-space-list* :test #'equal :key #'orgn)))
           (declare (type (or finite-space null) already))
           (when already
             (return-from build-finite-space already)))
         (let ((finspace (make-instance 'finite-space)))
           (declare (type finite-space finspace))
           (if top     (setf (slot-value finspace 'top)     top))
           (if stong   (setf (slot-value finspace 'stong)   stong))
           (if heights (setf (slot-value finspace 'heights) heights))
           (if orgn    (setf (slot-value finspace 'orgn)    orgn))
           (push finspace *finite-space-list*)
           finspace))))


(DEFUN RANDOM-FINITE-SPACE (dmns dens)
  #| Random Finite-Space of size 'dmns' and density 'dens' |#
  (declare (fixnum dmns))
  (unless (plusp dmns)
    (error "In RANDOM-FINITE-SPACE, the dimension ~D should be positive." dmns))
  (build-finite-space :top (randomtop dmns dens)
                      :orgn `(RANDOM-FINITE-SPACE ,(1+ *idnm-counter*))))


(DEFUN CARDINALITY (finspace)
  (ncol (stong finspace)))


#|
  (build-finite-space :stong (randomstong 10 .5))
  (build-finite-space :stong (randomstong 10 .5)
                      :orgn 'Top1)

  (random-finite-space 10 .5)
  (stong (random-finite-space 10 .5))
  (top (random-finite-space 10 .5))
  (heights (random-finite-space 10 .5))
|#


(DEFUN ADMISSIBLE (topogenous i j)
 #| Determine if Hat-Uj-{i} is contractible (in particular, the edge (i j) is homologically admissible) |#
  (if (and (> j i) (eq (terme topogenous i j) 1))
      (the boolean
           (eq (length (core-list topogenous
                                  (loop for k from 1 to (1- j)
                                        unless (or (zerop (terme topogenous k j)) (eq k i))
                                        collect k))) 1))))


(DEFMETHOD ADMISSIBLE-P ((topogenous matrice))
  #| Determine if all the edges of 'topogenous' are homologically admissibles |#
  (the boolean
       (let ((rslt T)
             (stong (stongmat topogenous)))
         (do ((j 2 (1+ j))) ((or (null rslt) (> j (ncol topogenous))) rslt)
           (do ((i 1 (1+ i))) ((or (null rslt) (> i (1- j))))
             (if (= (terme stong i j) 1)
                 (unless (admissible topogenous i j)
                     (setf rslt NIL))))))))


(DEFMETHOD ADMISSIBLE-P ((finspace finite-space))
  #| Determine if all the edges of (top 'finspace') are homologically admissibles |#
  (the boolean
       (let ((rslt T))
         (do ((j 2 (1+ j))) ((or (null rslt) (> j (ncol (top finspace)))) rslt)
           (do ((i 1 (1+ i))) ((or (null rslt) (> i (1- j))))
             (if (= (terme (stong finspace) i j) 1)
                 (unless (admissible (top finspace) i j)
                     (setf rslt NIL))))))))


#|
  (admissible-p (bar-subdivision (random-finite-space 6 .4)))    ; Always T
  (admissible-p (random-finite-space 6 .4))
|#


(DEFMETHOD DVFIELD ((finspace finite-space) &key random h-admissible)
  #| Compute a discrete vector field on 'stong' |#
  (let* ((m_stong (binarymatrice-to-ubasis (nilpot (stong finspace))))
         (coln (length m_stong))
         (status (make-array coln))
         (matched '()) 
         (vf '())
         (columns (<a-b> 1 coln)))
      (when random (setf columns (fisher-yates columns))) ;random cols
      (dolist (j columns) ; Here it is working on the columns. 'coln' is the dimension of the space
        (unless (member j matched)
          (let ((coljlist (svref m_stong (1- j))))
            (declare (type list coljlist))
            (when coljlist
              (when random (setf coljlist (fisher-yates coljlist))) ;random rows
              (dolist (rowi coljlist) ; Run over the edges (i,j)
                (declare (type fixnum rowi))
                (unless (member rowi matched)
                  (let ((rowistatus (svref status (1- rowi))))
                    #| rowistatus = {x : h(x) = h(i), x --> i} [Forman]|#
                    (unless (eql 1 rowistatus)
                      (when (and (or h-admissible (admissible (top finspace) rowi j)) ; here we need the topogenous matrix
                                 (null-intersection-p rowistatus coljlist)) ; avoids a loop  x --> i -> j -> x
                        (setf (svref status (1- rowi)) 1) ; rowi is a source
                        (setf rowistatus (insert rowi rowistatus))
                        (setf matched (insert j matched))
                        (setf matched (insert rowi matched))
                        (push (list rowi j) vf)
                        (dolist (item coljlist) ; here we have i --> item 
                          (unless (= item rowi)
                            (case (svref status (1- item))
                              (1 (dotimes (rowk coln)
                                   (let ((rowkstatus (svref status rowk)))
                                     (unless (eql 1 rowkstatus)
                                       (when (member item rowkstatus)
                                         (setf (svref status rowk)
                                               (union-merge rowistatus rowkstatus)))))))
                              (otherwise
                               #| if x --> i and i --> item then x --> item |#
                               (setf (svref status (1- item)) 
                                     (union-merge rowistatus (svref status (1- item))))))))
                        (return))))))))))
    (return-from DVFIELD (nreverse vf))))

    
(DEFMETHOD DVFIELD ((topogenous matrice) &key random h-admissible)
   (dvfield (make-instance 'finite-space
                           :top topogenous) :random random :h-admissible h-admissible))
                                   
                                   
(DEFUN DVFIELD-AUX (topogenous random h-admissible)
   (dvfield topogenous :random random :h-admissible h-admissible))


#|
  (dvfield (random-finite-space 20 .5))
  (dvfield (randomtop 20 .5))
|#


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFUN F-AUX1 (z num)
  (mapcar #'(lambda (x) (+ x num)) z))


(DEFUN F-AUX2 (x num)
  (mapcar #'(lambda (lista)
              (if lista
                  (if (= (car lista) 1)
                      (cons 1 (f-aux1 (cdr lista) num))
                    (f-aux1 lista num)))) x))


(DEFUN WEDGE-AT-X1-AUX (lcovers-lst cardinalities-lst card)
  (let ((rslt '())
        (num (1- card)))
    (do ((lcoversk lcovers-lst (cdr lcoversk))
         (cardk cardinalities-lst (cdr cardk)))
        ((endp lcoversk) rslt)
      (let ((list-lcovers1 (cdr (vector-to-list (car lcoversk)))))
        (setf rslt (append rslt (f-aux2 list-lcovers1 num))))
      (setf num (1- (+ num (car cardk)))))))


(DEFUN WEDGE-AT-X1 (finspace &rest rest)
  #| Wedge of finite spaces at the common point X1 |#
  (unless rest (return-from WEDGE-AT-X1 finspace))
  (let ((lcovers-lst '())
        (cardinalities-lst '()))
    (dolist (space rest)
      (setf lcovers-lst (cons (binarymatrice-to-ubasis (stong space)) lcovers-lst))
      (setf cardinalities-lst (cons (ncol (stong space)) cardinalities-lst)))
    (build-finite-space :stong (ubasis-to-binarymatrice (list-to-vector (append (vector-to-list (binarymatrice-to-ubasis (stong finspace)))
                                                                                (WEDGE-AT-X1-AUX lcovers-lst cardinalities-lst (ncol (stong finspace))))))
                        :orgn `(WEDGE-AT-X1 ,finspace ,@(sort rest #'< :key #'idnm)))))


(DEFUN WEDGE-SPHERES (num-minimals num-spheres)
  #| A finite model of a wedge of 'num-spheres' 1-spheres. This space has 'num-minimals' minimal elements and height = 1|#
  (unless (> num-spheres 0)
    (error "Second parameter must be a positive integer"))
  (let ((minimals (fisher-yates (<a-b> 1 num-minimals)))
        (maximals (fisher-yates (<a-b> (1+ num-minimals) (1- (* 2 num-minimals)))))
        (rslt (identite (+ num-minimals num-minimals num-spheres -1))))

    (do ((k 0 (1+ k))) ((> k (- num-minimals 2)))
      (insert-term rslt (nth k minimals) (nth k maximals) 1)
      (insert-term rslt (nth (1+ k) minimals) (nth k maximals) 1)) ; This is to ensure connectedness

    (insert-term rslt (car minimals) (* 2 num-minimals) 1) ; Put the first extra maximal
    (insert-term rslt (nth (1+ (random (1- num-minimals))) minimals) (* 2 num-minimals) 1)
    
    (unless (eq num-spheres 1)
      (insert-term rslt (car (last minimals)) (+ (* 2 num-minimals) num-spheres -1) 1)
      (insert-term rslt (nth (random (- num-minimals 2)) minimals) (+ (* 2 num-minimals) num-spheres -1) 1) ; In order to obtain a space with no beat points
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


(DEFUN MOORE-N1-AUX (n)
  #| Minimal finite model of the Moore space M(Zn , 1).
     Refer to Poset splitting... Cianci & Ottina : vertices d = 1, c_pares = 2, c_impares = 3 |#
  (let* ((contents '())
         (start)
         (2n+5 (+ (* 2 n) 5))
         (4n+5 (+ (* 2 n) 2n+5)))
    (declare (type fixnum 2n+5 4n+5))
    #| Max level (2) |#
    (push (list 5 6 2n+5) contents)
    (do ((k (1- 2n+5) (1- k)))
        ((< k 6))
      (if (zerop (mod k 2))
          (setf start 4)
        (setf start 5))
      (push (list start k (1+ k)) contents))
    #| (Max -1) level (1) |#
    (loop repeat n
          do (push (list 1 3) contents)
          do (push (list 1 2) contents))
    (setf contents (append '(NIL NIL NIL (2 3) (2 3)) contents))
    (return-from MOORE-N1-AUX (list
                               (make-array 4n+5 :initial-contents contents)
                               (list '(1 2 3) (<a-b> 4 2n+5) (<a-b> (1+ 2n+5) 4n+5))))))

            
(DEFUN MOORE-NK-AUX (n k)
  #| Minimal finite model of the Moore space M(Zn , k) for  k > 1 |#
  (declare (type fixnum n k))
  (let* ((contents '())
         (heights '())
         (start)
         (2k (* 2 k)) (2k-1 (1- 2k)) (2k+1 (1+ 2k)) (2k+3 (+ 2k+1 2))
         (2k+4 (1+ 2k+3)) (2n+2k+3 (+ (* 2 n) 2k+3)) (2k-3 (- 2k-1 2))
         (4n+2k+3 (+ (* 2 n) 2n+2k+3)))
    (declare (type fixnum 2k 2k-1 2k+1 2k+3 2k-3 2k+4 2n+2k+3 4n+2k+3)) 
    #| Max level |#
    (push (list 2k+3 2k+4 2n+2k+3) contents)
    (do ((elto (1- 2n+2k+3) (1- elto)))
        ((< elto 2k+4))
      (if (zerop (mod elto 2))
          (setf start (1- 2k+3))
        (setf start 2k+3))
      (push (list start elto (1+ elto)) contents))
    #| (Max -1) level |#
    (loop repeat n
          do (push (list 2k-1 2k+1) contents)
          do (push (list 2k-1 2k) contents))
    (push (list 2k 2k+1) contents)
    (push (list 2k 2k+1) contents)
    #| ''Inverse'' suspensions|#
    (push (list 2k-3 (1+ 2k-3)) contents)
    (loop for i downfrom 2k-3 to 1 by 2
          do (push (list i (1+ i)) heights)
          do (push (list i (1+ i)) contents)
          do (push (list i (1+ i)) contents))
    (push nil contents)
    (push nil contents) 
    (setf heights (append heights (list (list 2k-1 2k 2k+1)
                                        (<a-b> (1+ 2k+1) 2n+2k+3)
                                        (<a-b> (1+ 2n+2k+3) 4n+2k+3))))
    (return-from MOORE-NK-AUX (list (make-array 4n+2k+3 :initial-contents contents) heights))))


(DEFUN FINITE-MODEL-MOORE (n k)
  #| Minimal finite model of the Moore space M(Zn , k) |#
  (if (< n 2)
      (error "In FINITE-MODEL-MOORE, the argument n = ~D must be greater than 1." n)
    (if (< k 1)
        (error "In FINITE-MODEL-MOORE, the argument k = ~D must be greater than 0." k)
      (let ((info '()))
        (if (eq k 1)
            (setf info (moore-n1-aux n))
          (setf info (moore-nk-aux n k)))
        (build-finite-space
         :stong (nilpot-1 (ubasis-to-binarymatrice (first info)))
         :heights (second info)
         :orgn `(FINITE-MODEL-MOORE ,n ,k))))))


(DEFUN FINITE-MODEL-SPHERE (dmns)
  #| Minimal finite model of the 'dmns'-dimensional sphere |#
  (declare (fixnum dmns))
  (if (minusp dmns)
      (error "In FINITE-MODEL-SPHERE, the dimension ~D must be non-negative integer." dmns)
    (build-finite-space
     :stong (let ((rslt (identite (+ (* 2 dmns) 2))))
                (do ((i 1 (1+ i))) ((> i (* 2 dmns)) rslt)
                  (insert-term rslt i (+ i 2) 1)
                  (if (zerop (mod i 2))
                      (insert-term rslt i (+ i 1) 1)
                    (insert-term rslt i (+ i 3) 1))))
     :heights (loop for k from 1 to (1+ (* 2 dmns)) by 2
                    collect (list k (1+ k)))
     :orgn `(FINITE-MODEL-SPHERE ,dmns))))


(DEFUN RANDOM-STONG-2SPACE (dimension)
  (let* ((rslt (identite dimension))
         (card-H0 (random 100))
         (card-H1 (random 100))
         (card-H2 (random 100))
         (sum (+ card-H0 card-H1 card-H2)))
    (setf card-H0 (+ 2 (round (* (/ card-H0 sum) (- dimension 6))))
          card-H1 (+ 2 (round (* (/ card-H1 sum) (- dimension card-H0 4))))
          card-H2 (- dimension (+ card-H0 card-H1)))
    (let* ((H0 (<a-b> 1 card-H0))
           (H1 (<a-b> (1+ card-H0) (+ card-H0 card-H1)))
           (H2 (<a-b> (1+ (- dimension card-H2)) dimension)))
      
      (dolist (elmt H1)
        (let ((card-Uelmt (+ (random (1+ (- card-H0 2))) 2))
              (card-Felmt (+ (random (1+ (- card-H2 2))) 2)))
          (setf H0 (fisher-yates H0))
          (dotimes (k card-Uelmt)
            (insert-term rslt (nth k H0) elmt 1))
          (setf H2 (fisher-yates H2))
          (dotimes (k card-Felmt)
            (insert-term rslt elmt (nth k H2) 1))))
      
      (dolist (elmt H0)
        (let ((card-1 (+ (random (1+ (- card-H1 2))) 2)))
          (setf H1 (fisher-yates H1))
          (dotimes (k card-1)
            (insert-term rslt elmt (nth k H1) 1))))

      (dolist (elmt H2 rslt)
        (let ((card-2 (+ (random (1+ (- card-H1 2))) 2)))
          (setf H1 (fisher-yates H1))
          (dotimes (k card-2)
            (insert-term rslt (nth k H1) elmt 1)))))))
            

(DEFUN RANDOM-TOP-2SPACE (dimension)
  (topmat (random-stong-2space dimension)))


(DEFUN RANDOM-2SPACE (dimension)
  #| Return a random minimal finite space of size 'dimension', height 2 |#
  (declare (fixnum dimension))
  (unless (> dimension 5)
    (error "In RANDOM-2SPACE, the dimension ~D must be greater than or equal to 6." dimension))
  (build-finite-space :stong (random-stong-2space dimension)
                      :orgn `(RANDOM-2SPACE ,dimension ,(1+ *idnm-counter*))))


#|
  (show (top (finite-model-sphere 0)))
  (show (stong (finite-model-sphere 4)))
  (heights (finite-model-sphere 17))
  (dvfield (finite-model-sphere 15))
  
  (setf finspace (wedge-at-X1
                  (finite-model-moore 7 3)
                  (finite-model-moore 10 2)
                  (finite-model-sphere 4)
                  (finite-model-moore 5 1)))
  (h-regular-homology finspace)
|#


(DEFMETHOD NON-HAUSDORFF-SUSPENSION ((topogenous matrice))
  #| Non-Hausdorff suspension of 'topogenous' |#
  (let* ((nlig+1 (1+ (nlig topogenous)))
         (rslt (add-points topogenous nlig+1 2)))
    (maj-colonne rslt nlig+1 (loop for x in (<a-b> 1 nlig+1)
                                   collect (list x 1)))
    (maj-colonne rslt (1+ nlig+1) (loop for x in (<a-b> 1 (1- nlig+1))
                                        collect (list x 1)))
    (insert-term rslt (1+ nlig+1) (1+ nlig+1) 1)
    rslt))


(DEFMETHOD NON-HAUSDORFF-SUSPENSION ((finspace finite-space))
   #| Non-Hausdorff suspension of 'finspace' |#
  (let* ((ubasis (binarymatrice-to-ubasis (top finspace)))
         (card (length ubasis))
         (eltos (<a-b> 1 card)))
    (build-finite-space :top (ubasis-to-binarymatrice (list-to-vector 
                                                         (append (vector-to-list ubasis)
                                                                 (list (append eltos (list (1+ card)))
                                                                       (append eltos (list (+ card 2)))))))
                        :orgn `(NON-HAUSDORFF-SUSPENSION ,(idnm finspace)))))
   
#|
  (setf finspace (random-finite-space 9 .5))
  (top finspace)
  (setf mat1 (NON-HAUSDORFF-SUSPENSION (top finspace)))
  (setf mat2 (top (NON-HAUSDORFF-SUSPENSION finspace)))
  (equalmatrix mat1 mat2)
|#


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


#|
(setf TORUS
      (build-finite-space
       :top (convertarray
             #2A((1 0 0 0 1 1 1 0 1 0 0 0 1 1 1 1)
                 (0 1 0 0 1 1 0 1 0 1 0 0 1 1 1 1)
                 (0 0 1 0 0 0 1 0 1 0 1 1 1 1 1 1)
                 (0 0 0 1 0 0 0 1 0 1 1 1 1 1 1 1)
                 (0 0 0 0 1 0 0 0 0 0 0 0 1 0 1 0)
                 (0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 1)
                 (0 0 0 0 0 0 1 0 0 0 0 0 1 1 0 0)
                 (0 0 0 0 0 0 0 1 0 0 0 0 1 1 0 0)
                 (0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 1)
                 (0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 1)
                 (0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0)
                 (0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1)
                 (0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0)
                 (0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0)
                 (0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0)
                 (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1)))
       :orgn 'TORUS))  ; Finite model of the torus


(setf RP2
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
       :orgn 'RP2))  ; Finite model of the projective plane


(setf KLEIN
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
       :orgn 'KLEIN))  ; Finite model of the Klein Bottle


(setf HOMTRIVIAL
      (build-finite-space
       :stong (convertarray
               #2A((1 0 0 1 1 0 0 0 0)
                   (0 1 0 1 1 1 0 0 0)
                   (0 0 1 0 1 1 0 0 0)
                   (0 0 0 1 0 0 1 1 0)
                   (0 0 0 0 1 0 1 0 1)
                   (0 0 0 0 0 1 0 1 1)
                   (0 0 0 0 0 0 1 0 0)
                   (0 0 0 0 0 0 0 1 0)
                   (0 0 0 0 0 0 0 0 1)))
       :orgn 'HOMTRIVIAL))  ; Smallest homotopically trivial non-contractible space
|#
