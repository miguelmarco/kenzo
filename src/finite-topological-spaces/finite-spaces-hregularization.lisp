
;;  H-REGULARIZATION H-REGULARIZATION H-REGULARIZATION H-REGULARIZATION
;;  H-REGULARIZATION H-REGULARIZATION H-REGULARIZATION H-REGULARIZATION
;;  H-REGULARIZATION H-REGULARIZATION H-REGULARIZATION H-REGULARIZATION

(IN-PACKAGE #:cat)

(provide "finite-spaces-hregularization")

;;
;;  h-regularization of finite spaces of height at most 2
;;


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFMETHOD HAT-Un ((topogenous matrice) n)
  #| List of the elements in Un-{Xn} |#
  (loop for i from 1 to (1- n)
        when (eq (terme topogenous i n) 1)
        collect i))


(DEFMETHOD HAT-Fn ((topogenous matrice) n)
  #| List of the elements in Fn-{Xn} |#
  (loop for j from (1+ n) to (ncol topogenous)
        when (eq (terme topogenous n j) 1)
        collect j))


(DEFMETHOD HAT-Un ((finspace finite-space) n)
  (HAT-Un (top finspace) n))


(DEFMETHOD HAT-Fn ((finspace finite-space) n)
  (HAT-Fn (top finspace) n))


(DEFUN ARRAY-TO-LIST (array)
  (map 'list #'identity array))


(DEFUN EDGES-TO-MATRICE (list size)
  (let ((rslt (identite size)))
    (dolist (x list rslt)
      (insert-term rslt (car x) (second x) 1))))


#|
 (setf ejemplo (randomtop 10 .6))
 (loop for k from 1 to (nlig ejemplo)
       collect (HAT-UN ejemplo k))
|#


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFMETHOD SEPARATE-SEGMENTS ((nil-topogenous matrice) n)
  #| 'nil-topogenous' must be a core and nilpot |#
  (let* ((hat-list (hat-Un nil-topogenous n))
         (len (length hat-list)))
    (if (eq len 2)
        nil-topogenous
      (let* ((rslt (add-points nil-topogenous n (- len 2)))
             (hatFn (newsmith-extract-line rslt (+ n len -2))))
        (do ((k n (1+ k))
             (pairs hat-list (cdr pairs)))
            ((eq (length pairs) 2) rslt)
          (maj-colonne rslt k (list (list (car pairs) 1) (list (cadr pairs) 1)))
          (maj-ligne rslt k hatFn)
          (insert-term rslt (car pairs) (+ n len -2) 0))))))
        

(DEFMETHOD 1-H-REGULARIZE ((nil-topogenous matrice) heights1)
  #| nil-topogenous must be a core and nilpot, heights1 in descending order |#
  (let ((rslt nil-topogenous))
    (dolist (n heights1)
      (setf rslt (separate-segments rslt n)))
    rslt))


(DEFMETHOD SEPARATE-SEGMENTS ((minimal-finspace finite-space) point)
  #| Separate the segments in point (point must be of height 1) |#
  (build-finite-space :top (nilpot-1 (separate-segments (nilpot (top minimal-finspace)) point))
                      :orgn `(SEPARATE-SEGMENTS ,minimal-finspace ,point)))


(DEFUN 1-H-REGULARIZATION (minimal-finspace)
  #| Return the 1-h-regularization of the minimal finite space minimal-finspace |#
  (let ((heights1 (reverse (second (heights minimal-finspace))))
        (finspace (build-finite-space :orgn `(1-H-REGULARIZATION ,minimal-finspace))))
    (if (>= (idnm finspace) *idnm-counter*)
        (setf (slot-value finspace 'top) (nilpot-1 (1-h-regularize (nilpot (top minimal-finspace)) heights1))))
    finspace))


#|
  (setf ejemplo (topmat (edges-to-matrice '((1 5) (2 5) (3 5)
                                            (4 5) (1 6) (2 6)
                                            (4 6) (5 7) (6 7)
                                            (7 9)) 10)))
  (setf heights1 (reverse (second (heights (stongmat ejemplo)))))
  (setf 1hreg (nilpot-1 (1-h-regularize (nilpot ejemplo) heights1)))
  |#


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFMETHOD Un-CONNECT ((h-topogenous matrice) n)
 #| h-topogenous must be 1-h-regularized and nilpot |#
  (let* ((Ux-list (hat-Un h-topogenous n))
         (Ux-mat (submatrix h-topogenous Ux-list Ux-list))
         (d0-d1 (h-regular-dif Ux-mat))
         (components (homologie (copier-matrice (first d0-d1)) (copier-matrice (second d0-d1)))))
    (if (eq (length components) 1)
        h-topogenous
      (let* ((base-points (mapcar #'(lambda (x) (cdr (second x))) components))
             (minimals (car (heights Ux-mat)))
             (rslt (add-points h-topogenous n (1- (length components))))
             (new-n (+ n (1- (length components))))
             (hatFn (newsmith-extract-line rslt new-n)))
        (do ((k n (1+ k))
             (pairs base-points (cdr pairs)))
            ((eq (length pairs) 1) rslt)
          (insert-term rslt k new-n 1)
          (maj-colonne rslt k (list (list (nth (1- (nth (1- (car pairs))  minimals)) Ux-list) 1)
                                    (list (nth (1- (nth (1- (cadr pairs)) minimals)) Ux-list) 1)))
          (maj-ligne rslt k hatFn))))))


(DEFMETHOD Un-CONNECT ((1hreg-finspace finite-space) point)
  #| It makes the minimal open set of 'point' connected (1hreg-finspace must be 1-h-regular and minimal) |#
  (build-finite-space :top (nilpot-1 (Un-connect (nilpot (top 1hreg-finspace)) point))
                      :orgn `(UN-CONNECT ,1hreg-finspace ,point)))
  

#|
  (setf ejemplo (topmat (edges-to-matrice '((1 6) (2 6) (3 7) (4 7) (5 8) (6 8) (7 8) (8 9)) 9)))
  (setf heights1 (reverse (second (heights ejemplo))))
  (setf 1hreg (nilpot-1 (1-h-regularize (nilpot ejemplo) heights1)))
  (setf connected8 (nilpot-1 (Un-connect (nilpot 1hreg) 8)))
|#


(DEFUN CONNECT (h-topogenous heights2)
  #| h-topogenous must be h-regularized, nilpot and heights2 in descending order |#
  (let ((rslt h-topogenous))
    (dolist (n heights2)
      (setf rslt (Un-connect rslt n)))
    rslt))


(DEFUN Un-CONNECTED (1hreg-finspace)
  #| It makes the minimal open sets connected of points of height 2  (1hreg-finspace must be 1-h-regular and minimal) |#
  (let ((heights2 (reverse (third (heights 1hreg-finspace)))))
    (build-finite-space :top (nilpot-1 (connect (nilpot (top 1hreg-finspace)) heights2))
                        :orgn `(UN-CONNECTED ,1hreg-finspace))))


(DEFUN FIND-H-SPHERE (core minimals maximals)
  #| minimals and maximals in ascending order |#
  (let* ((down-max (loop for x in maximals
                         collect (HAT-Un core x)))
         (up-min (loop for x in minimals
                       collect (HAT-Fn core x)))
         (min-sphere (car down-max))
         (max-sphere (list (car maximals)))
         (position NIL))
    
    (do ((k 1 (1+ k))) (position)
      (let* ((v_k (car (remove (car max-sphere) (nth (position (car min-sphere) minimals) up-min))))
             (u_k (car (remove (car min-sphere) (nth (position v_k maximals) down-max)))))
        (push v_k max-sphere)
        (unless (setf position (position u_k min-sphere))
          (push u_k min-sphere))))

    (subseq max-sphere 0 (1+ position))))


(DEFUN FIND-CROWN (minimal-finspace)
  #| It finds a crown in minimal-finspace (this space must be minimal and of height 1) |#
  (let ((minimals (first (heights minimal-finspace)))
        (maximals (second (heights minimal-finspace))))
    (find-h-sphere minimal-finspace minimals maximals)))


#|
  (setf ejem (wedge-spheres 6 5))
  (setf alturas (heights ejem))
  (find-h-sphere ejem (first alturas) (second alturas))
|#


(DEFUN H-SPHERES (1h-topogenous elements times)
  #| 1h-topogenous must be h-regularized of height 1 |#
  (if (zerop times)
      '()
    (let* ((core-list (core-list 1h-topogenous elements))
           (core-mat (submatrix 1h-topogenous core-list core-list))
           (core-min-max (heights core-mat))
           (core-min (first core-min-max))
           (core-max (second core-min-max))
           (cycle (mapcar #'(lambda (x) (nth (1- x) core-list))
                          (find-h-sphere core-mat core-min core-max)))
           (rslt NIL))
        (setf rslt (append (list cycle) (h-spheres 1h-topogenous (remove (car cycle) elements) (1- times))))
        rslt)))


#|
  ((lambda (m n)
     (setf 1h-topogenous (wedge-spheres m n))
     (h-spheres 1h-topogenous (<a-b> 1 (+ (* 2 m) (1- n))) n)) 3 3)
|#


(DEFMETHOD SEPARATE-CYCLES ((nilh-topogenous matrice) n)
  #| nilh-topogenous must be nilpot, core and h-regularized |#
  (let* ((Ux-list (hat-Un nilh-topogenous n))
         (Ux-mat (submatrix nilh-topogenous Ux-list Ux-list))
         (min-max (heights Ux-mat))
         (minimals (first min-max))
         (maximals (second min-max))
         (char-euler+1 (- (length maximals) (length minimals))))
    
    (if (zerop char-euler+1) ; case char-euler+1 = 0 : h-sphere
        nilh-topogenous

      (if (< char-euler+1 0) ; case char-euler+1 = -1 : contractible
          (let* ((rslt (add-points nilh-topogenous n 1))
                 (hatFn (newsmith-extract-line rslt (1+ n))))
            (maj-colonne rslt n (list (list (car Ux-list) 1) (list (cadr Ux-list) 1)))
            (maj-ligne rslt n (push (list (1+ n) 1) hatFn))
            rslt)

        (let* ((rslt (add-points nilh-topogenous n char-euler+1)) ; case char-euler+1 > 0 : wedge
               (hatFn (newsmith-extract-line rslt (+ n char-euler+1)))
               (list-separations (h-spheres (nilpot-1 Ux-mat) (<a-b> 1 (length Ux-list)) char-euler+1))) ; nilpot-1 since core-list is used
          
          (do ((k n (1+ k))) ((eq k (+ n char-euler+1)) rslt)
            (maj-colonne rslt k (loop for x in (nth (- k n) list-separations)
                                      collect (list (nth (1- x) Ux-list) 1)))
            (maj-ligne rslt k hatFn)
            (insert-term rslt (nth (1- (car (nth (- k n) list-separations))) Ux-list) (+ n char-euler+1) 0)))))))


(DEFMETHOD SEPARATE-CYCLES ((1hreg-finspace finite-space) point)
  (build-finite-space :top (nilpot-1 (separate-cycles (nilpot (top 1hreg-finspace)) point))
                      :orgn `(SEPARATE-CYCLES ,1hreg-finspace ,point)))


#|
  (setf ejemplo (topmat (edges-to-matrice '((1 6) (2 6) (3 7)
                                            (4 7) (5 8) (6 8)
                                            (7 8) (8 9)) 9)))
  (setf heights1 (reverse (second (heights ejemplo))))
  (setf 1hreg (nilpot-1 (1-h-regularize (nilpot ejemplo) heights1)))
  (setf topogenous (connect 1hreg (reverse (third (heights 1hreg)))))
  (separate-cycles topogenous 10)

  (setf ejemplo (topmat (edges-to-matrice '((1 6) (1 7) (1 8) (2 6)
                                            (2 7) (2 8) (2 9) (3 9)
                                            (3 10) (3 11) (3 13) (4 10)
                                            (4 12) (5 11) (5 12) (5 13)
                                            (6 14) (7 14) (8 14) (9 14)
                                            (10 14) (11 14) (12 14) (13 14)) 14)))
  (nilpot-1 (separate-cycles (nilpot ejemplo) 14))

  (setf ejemplo (topmat (edges-to-matrice '((1 4) (1 5) (2 4) (2 5)
                                            (2 6) (3 6) (2 7) (3 7)
                                            (2 8) (3 8) (6 9) (7 9)
                                            (8 9)) 9)))
  (nilpot-1  (separate-cycles (nilpot ejemplo) 9))
|#


(DEFMETHOD 2-H-REGULARIZE ((h-topogenous matrice) heights2)
  #| h-topogenous must be h-regularized and heights2 in descending order |#
  (let ((rslt h-topogenous))
    (dolist (n heights2)
      (setf rslt (separate-cycles rslt n)))
    rslt))


(DEFMETHOD 2-H-REGULARIZE ((finspace finite-space) heights2)
  #| heights2 in descending order |#
  (build-finite-space :top (topmat (2-h-regularize (top finspace) heights2))
                      :orgn `(2-H-REGULARIZATION ,finspace)))



(DEFUN ORIENTATION (matrix)
  (let ((index 0))
    (loop for j from 1 to (ncol matrix)
          do (and (setf index (loop for i from 1 to (nlig matrix)
                                    thereis (and (eq (terme matrix i j) 1) i)))
                  (insert-term matrix index j -1))))
    matrix)


(DEFMETHOD 2H-REGULARIZATION ((core-topogenous matrice))
  #| core-topogenous must not have beat points |#
  (let* ((rslt (nilpot core-topogenous))
         (level (reverse (second (heights rslt)))))
    #| Separate segments |#
    (dolist (elto-H1 level)
      (let* ((hat-list (hat-Un rslt elto-H1))
             (len-2 (- (length hat-list) 2)))
        (unless (zerop len-2)
          (setf rslt (add-points rslt elto-H1 len-2))
          (let ((hatFn (newsmith-extract-line rslt (+ elto-H1 len-2))))
            (do ((k elto-H1 (1+ k))
                 (pairs hat-list (cdr pairs)))
                ((eq (length pairs) 2))
              (maj-colonne rslt k (list (list (car pairs) 1) (list (cadr pairs) 1)))
              (maj-ligne rslt k hatFn)
              (insert-term rslt (car pairs) (+ elto-H1 len-2) 0))))))

    (setf level (reverse (third (heights rslt))))

    #| Connect Ux |#
    (dolist (elto-H2 level)
      (let* ((Ux-list (hat-Un rslt elto-H2))
             (Ux-mat (submatrix rslt Ux-list Ux-list))
             (heights-U (heights Ux-mat))
             (minimals (car heights-U))
             (lmins (length minimals))
             (d1 (submatrix Ux-mat (<a-b> 1 lmins) (>a-b> lmins (ncol Ux-mat))))
             (components (homologie (creer-matrice 0 lmins) (orientation (copier-matrice d1)))))
        (unless (equal (length components) 1)
          (setf rslt (add-points rslt elto-H2 (1- (length components))))
          (let* ((base-points (mapcar #'(lambda (x) (cdr (second x))) components))
                 (new-n (+ elto-H2 (1- (length components))))
                 (hatFn (newsmith-extract-line rslt new-n)))
            (do ((k elto-H2 (1+ k))
                 (pairs base-points (cdr pairs)))
                ((eq (length pairs) 1))
              (insert-term rslt k new-n 1)
              (maj-colonne rslt k (list (list (nth (1- (nth (1- (car pairs))  minimals)) Ux-list) 1)
                                        (list (nth (1- (nth (1- (cadr pairs)) minimals)) Ux-list) 1)))
              (maj-ligne rslt k hatFn))))))
    
    (setf level (reverse (third (heights rslt)))) 

    #| Separate-cycles |#
    (dolist (elto-H2 level)
      (let* ((Ux-list (hat-Un rslt elto-H2))
             (Ux-mat (submatrix rslt Ux-list Ux-list))
             (min-max (heights Ux-mat))
             (minimals (first min-max))
             (maximals (second min-max))
             (char-euler+1 (- (length maximals) (length minimals))))


        (unless (zerop char-euler+1) ; case char-euler+1 = 0  : h-sphere          
          (if (< char-euler+1 0)     ; case char-euler+1 = -1 : contractible
              (progn 
                (setf rslt (add-points rslt elto-H2 1)) 
                (let ((hatFn (newsmith-extract-line rslt (1+ elto-H2))))
                  (maj-colonne rslt elto-H2 (list (list (car Ux-list) 1) (list (cadr Ux-list) 1)))
                  (maj-ligne rslt elto-H2 (push (list (1+ elto-H2) 1) hatFn))))

            (progn                    ; case char-euler+1 > 0 : wedge
              (setf rslt (add-points rslt elto-H2 char-euler+1))  
              (let ((hatFn (newsmith-extract-line rslt (+ elto-H2 char-euler+1))) 
                    (list-separations (h-spheres (nilpot-1 Ux-mat) (<a-b> 1 (length Ux-list)) char-euler+1))) ; nilpot-1 since core-list is used
                
                (do ((k elto-H2 (1+ k))) ((eq k (+ elto-H2 char-euler+1)))
                  (maj-colonne rslt k (loop for x in (nth (- k elto-H2) list-separations)
                                            collect (list (nth (1- x) Ux-list) 1)))
                  (maj-ligne rslt k hatFn)
                  (insert-term rslt (nth (1- (car (nth (- k elto-H2) list-separations))) Ux-list) (+ elto-H2 char-euler+1) 0))))))))
    (nilpot-1 rslt)))


(DEFMETHOD 2H-REGULARIZATION ((finspace finite-space))
  #| 2h-regularization of a finite-space (must be of height 2) |#
  (build-finite-space :top (2H-REGULARIZATION (top finspace))
                      :orgn `(2H-REGULARIZATION ,finspace)))


(DEFUN 2-H-REGULARIZATION (minimal-finspace)
  #| Return the 2-h-regularization of the minimal finite space minimal-finspace |#
  (let ((finspace (build-finite-space :orgn `(2-H-REGULARIZATION ,minimal-finspace))))
    (if (>= (idnm finspace) *idnm-counter*)
        (setf (slot-value finspace 'top) (2h-regularization (top minimal-finspace))))
    finspace))


(DEFUN 2-SKELETON (topogenous)
  (let ((heights (heights topogenous)))
    (if (< (length heights) 3)
        topogenous)
    (let ((list (sort (append (first heights) (second heights) (third heights)) #'<)))
      (submatrix topogenous list list))))


#|
  (dotimes (q 20)
    (setf ejemplo (nilpot (core (randomtop 100 .3))))
    (setf heights1 (reverse (second (heights ejemplo))))
    (setf 1hreg (1-h-regularize ejemplo heights1))
    (setf heights2 (reverse (third (heights 1hreg))))
    (setf connected (connect 1hreg heights2)) 
    (setf heights3 (reverse (third (heights connected))))
    (setf 2hreg (nilpot-1 (2-h-regularize connected heights3)))
    (format t "~A ------------------------------------------~A~%" (ncol 2hreg) q))
|#


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFMACRO H-REGULAR-HOMOLOGY-PRUEBA (object dim)
  `(let ((size (typecase ,object
                 (finite-space (1- (length (heights ,object))))
                 (matrice (1- (length (heights ,object)))))))
     
     (let ((liss (H-REGULAR-DIF ,object)))
       (if (= ,dim size)
           (homologie (copier-matrice (nth ,dim liss))
                      (creer-matrice (ncol (nth ,dim liss)) 0))
         
         (homologie (copier-matrice (nth ,dim liss))
                    (copier-matrice (nth (1+ ,dim) liss)))))))


(DEFUN H-REGULAR-HOMOLOGY-SIM (finspace)
  (let ((size (1- (length (heights finspace))))
        (list (matrices (save-info-aux finspace))))

      (dotimes (dim size)
        (let ((Mn (copier-matrice (nth dim list)))
              (Mn+1 (copier-matrice (nth (1+ dim) list))))
          (declare (type matrice Mn Mn+1))          
          (let ((rsl (homologie Mn Mn+1)))
            (declare (type list rsl))    
            (format T "~%Homology in dimension ~D :" dim)
            (dolist (item rsl)
              (declare (type list item))
              (format T " Z")
              (unless (zerop (first item)) 
                (format T "/~DZ " (first item)))))))

      (let ((rsl (homologie (copier-matrice (nth size list))
                            (creer-matrice (ncol (nth size list)) 0))))
        (declare (type list rsl))    
        (format T "~%Homology in dimension ~D : " size)
        (format T " Z ^ ~D" (length rsl))))
  (values))

#|
  (dotimes (q 30)
    (setf ejemplo (random-2space 7))
    (setf core (core ejemplo))
    (print (eq (cardinality ejemplo) (cardinality core)))) ; It checks if ejemplo is a minimal space
|#


#|
  (dotimes (q 20)
    (setf ejemplo (random-2space 7))
    (format t "Space ready...~%")
    (setf bar (bar-subdivision ejemplo))
    (format t "Subdivision ready...~%")
    (setf hreg (2-h-regularization ejemplo))
    (format t "h-regularization ready...~%")
    (dotimes (k (length (heights hreg)))
      (unless (equal (mapcar #'car (H-REGULAR-HOMOLOGY-prueba bar k))    
                     (mapcar #'car (H-REGULAR-HOMOLOGY-prueba hreg k)))
        (error "Wrong!!")))
    (H-REGULAR-HOMOLOGY-SIM hreg)
    (H-REGULAR-HOMOLOGY-SIM bar)

    (format t "~%~A------------------------------------------(~A, ~A, ~A)~%" q (nlig (stong ejemplo)) (nlig (stong hreg)) (nlig (stong bar))))
|#
