
;;  POINT-REDUCTIONS POINT-REDUCTIONS POINT-REDUCTIONS POINT-REDUCTIONS POINT-REDUCTIONS
;;  POINT-REDUCTIONS POINT-REDUCTIONS POINT-REDUCTIONS POINT-REDUCTIONS POINT-REDUCTIONS
;;  POINT-REDUCTIONS POINT-REDUCTIONS POINT-REDUCTIONS POINT-REDUCTIONS POINT-REDUCTIONS

(IN-PACKAGE #:cat)

(provide "finite-spaces-point-reductions")

;;
;;  Computing some point reductions in Finite Topological Space
;;


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFUN DOWNBP (topogenous n &optional (list (<a-b> 1 (nlig topogenous))))
  #| Determine if Xn is a down beat point of the submatrix 'list'x'list' |#
  (the boolean
       (let ((maximal (loop for k in (reverse (subseq list 0 (position n list)))
                            thereis (and (eq (terme topogenous k n) 1) k))))
         
         (unless (or (null maximal) 
                     (loop for i in (subseq list 0 (position n list))
                           thereis (not (eq (terme topogenous i maximal) (terme topogenous i n)))))
           +TRUE+))))


(DEFMETHOD DOWN-BEAT-POINT ((ubasis vector) point &optional (list '()))
  #| Determine if 'point' is a down beat point of the submatrix 'list'x'list' |#
  (the boolean 
       (let ((Un (svref ubasis (1- point))))
         (declare (type list Un))
         (if list (setf Un (intersection Un list)))
         (if (eq (length Un) 1)
             (return-from down-beat-point +TRUE+)
           (return-from down-beat-point +FALSE+)))))


(DEFMETHOD DOWN-BEAT-POINT ((topogenous matrice) point &optional (list (<a-b> 1 (nlig topogenous))))
  #| Determine if Xn is a down beat point of the submatrix 'list'x'list' |#
  (downbp topogenous point list))


(DEFMETHOD DOWN-BEAT-POINT ((finspace finite-space) point &optional (list '()))
  #| Determine if 'point' is a down beat point of the subspace of 'finspace' of indexes in 'list' |#
  (let ((ubasis (binarymatrice-to-ubasis (nilpot (stong finspace)))))
    (down-beat-point ubasis point list)))


(DEFUN UPBP (topogenous n &optional (list (<a-b> 1 (nlig topogenous)))) 
  #| Deetrmine if Xn is an up beat point of the submatrix 'list'x'list' |#
  (the boolean
       (let ((minimal (loop for k in (subseq list (1+ (position n list)))
                            thereis (and (eq (terme topogenous n k) 1) k))))
         (unless (or (null minimal) 
                     (loop for i in (subseq list (1+ (position n list)))
                           thereis (not (eq (terme topogenous minimal i) (terme topogenous n i)))))
           +TRUE+))))


(DEFMETHOD UP-BEAT-POINT ((fbasis vector) point &optional (list '()))
  #| Determine if 'point' is an up beat point of the submatrix 'list'x'list' |#
  (the boolean 
       (let ((Fn (svref fbasis (1- point))))
         (declare (type list Fn))
         (if list (setf Fn (intersection Fn list)))
         (if (eq (length Fn) 1)
             (return-from up-beat-point +TRUE+)
           (return-from up-beat-point +FALSE+)))))


(DEFMETHOD UP-BEAT-POINT ((topogenous matrice) point &optional (list (<a-b> 1 (nlig topogenous))))
  #| Determine if Xn is an up beat point of the submatrix 'list'x'list' |#
  (upbp topogenous point list))
             

(DEFMETHOD UP-BEAT-POINT ((finspace finite-space) point &optional (list '()))
  #| Determine if 'point' is a down beat point of the subspace of 'finspace' of indexes in 'list' |#
  (let ((fbasis (binarymatrice-to-fbasis (nilpot (stong finspace)))))
    (up-beat-point fbasis point list)))


(DEFUN BEATPOINT (topogenous n &optional (list (<a-b> 1 (nlig topogenous))))
  #| Determine if Xn is a beat point of the submatrix 'list'x'list' |#
  (the boolean
       (if (< (1+ n) (/ (ncol topogenous) 2))
           (if (downbp topogenous n list)
               +TRUE+
             (upbp topogenous n list))
         (if (upbp topogenous n list)
             +TRUE+
           (downbp topogenous n list)))))


(DEFMETHOD BEAT-POINT ((finspace finite-space) point &optional (list '()))
  #| Determine if Xn is a beat point of the subspace of 'finspace' of indexes in 'list' |#
  (the boolean
       (if (< (1+ point) (/ (cardinality finspace) 2))
           (if (down-beat-point finspace point list)
               +TRUE+
             (up-beat-point finspace point list))
         (if (up-beat-point finspace point list)
             +TRUE+
           (down-beat-point finspace point list)))))

(DEFMETHOD BEAT-POINT ((topogenous matrice) point &optional (list (<a-b> 1 (nlig topogenous))))
  #| Determine if Xn is a beat point of the submatrix 'list'x'list' |#
  (beatpoint topogenous point list))


#|
  (show (setf M (randomtop 10 .5)))
  (downbp M 4 (<a-b> 1 10))
  (upbp M 4 (<a-b> 1 10))
  (beatpoint M 4 (<a-b> 1 10))
  (downbp M 5 '(1 3 5 8))
  (upbp M 5 '(1 3 5 8))
  (beatpoint M 5 '(1 3 5 8))
  
  (setf finspace (random-finite-space 30 .5))
  (dotimes (k (1- (cardinality finspace)))
    (print (eq (beat-point finspace (1+ k)) (beatpoint (top finspace) (1+ k)))))
|#


(DEFMETHOD CORE-LIST ((topogenous matrice) &optional (list (<a-b> 1 (nlig topogenous))))
  #| List of the elements of the core of the submatrix 'list'x'list' |#
  (let ((list-bps NIL))
    (setf list-bps (loop for n in list
                         thereis (and (beatpoint topogenous n list) n)))
    (if (null list-bps)
        (return-from core-list list)
      (core-list topogenous (list-difference list (list list-bps))))))


(DEFMETHOD CORE-LIST ((finspace finite-space) &optional (list (<a-b> 1 (cardinality finspace))))
  #| List of the elements of the core of the subspace of 'finspace' of indexes in 'list' |#
  (core-list (top finspace) list))


(DEFUN Un-N (topogenous point &optional (list (<a-b> 1 (nlig topogenous))))
  #| List of the elements in Un-{Xn} |#
  (loop for i in (subseq list 0 (position point list))
        when (eq (terme topogenous i point) 1)
        collect i))


(DEFUN Fn-N (topogenous point &optional (list (<a-b> 1 (nlig topogenous))))
  #| List of the elements in Fn-{Xn} |#
  (loop for j in (subseq list (1+ (position point list)))
        when (eq (terme topogenous point j) 1)
        collect j))


(DEFUN LINK-LIST (topogenous point &optional (list (<a-b> 1 (nlig topogenous))))
  #| List of the elements in the link of Xn |#
  (append (Un-N topogenous point list) (Fn-N topogenous point list)))


(DEFUN Cn-N (topogenous point &optional (list (<a-b> 1 (nlig topogenous))))
  #| List of the elements in the link of Xn |#
  (link-list topogenous point list))


(DEFUN WEAKPOINT (topogenous n &optional (list (<a-b> 1 (nlig topogenous))))
  #| Determine if Xn is a weak beat point of the submatrix 'list'x'list' |#
  (the boolean
       (or (eq (length (core-list topogenous (Un-N topogenous n list))) 1)
           (eq (length (core-list topogenous (Fn-N topogenous n list))) 1))))


(DEFMETHOD WEAK-POINT ((finspace finite-space) point &optional (list (<a-b> 1 (cardinality finspace))))
  #| Determine if Xn is a weak beat point of the subspace of 'finspace' of indexes in 'list' |#
  (weakpoint (top finspace) point list))


(DEFMETHOD WEAK-POINT ((topogenous matrice) point &optional (list (<a-b> 1 (nlig topogenous))))
  #| Determine if Xn is a weak beat point of the submatrix 'list'x'list' |#
  (weakpoint topogenous point list))

  
(DEFMETHOD WEAKCORE-LIST ((topogenous matrice) &optional (list (<a-b> 1 (nlig topogenous))))
  #| List of the elements of a "weak core" of the submatrix 'list'x'list' |#
  (let ((list-wps NIL))
    (setf list-wps (loop for n in list
                         thereis (and (or (beatpoint topogenous n list) (weakpoint topogenous n list)) n)))
    (if (null list-wps)
        (return-from weakcore-list list)
      (weakcore-list topogenous (list-difference list (list list-wps))))))


(DEFMETHOD WEAKCORE-LIST ((finspace finite-space) &optional (list (<a-b> 1 (cardinality finspace))))
  #| List of the elements of the "weak core" of the subspace of 'finspace' of indexes in 'list' |#
  (weakcore-list (top finspace) list))


(DEFMETHOD CORE ((topogenous matrice) &optional (list (<a-b> 1 (nlig topogenous))))
  #| Matrice of a core of the submatrix 'list'x'list' |#
  (the matrice
       (let ((corelist (core-list topogenous list)))
         (submatrix topogenous corelist corelist))))


(DEFMETHOD CORE ((finspace finite-space) &optional list)
  #| Finite-Space whose 'top' is the topogenous matrix of a core of the submatrix 'list'x'list' of (top 'finspace') |#
  (let ((already (find `(CORE ,finspace ,list) *finite-space-list* :test #'equal :key #'orgn)))
    (declare (type (or finite-space null) already))
    (if already
        already
      (let ((list2 (or list (<a-b> 1 (nlig (top finspace))))))
        (build-finite-space :orgn (if list `(CORE ,finspace ,list) `(CORE ,finspace))
                            :top (core (top finspace) list2))))))


(DEFMETHOD WEAKCORE ((topogenous matrice) &optional (list (<a-b> 1 (nlig topogenous))))
  #| Matrice of a "weak core" of the submatrix 'list'x'list' |#
  (the matrice
       (let ((weakcorelist (weakcore-list topogenous (core-list topogenous list))))
         (submatrix topogenous weakcorelist weakcorelist))))


(DEFMETHOD WEAKCORE ((finspace finite-space) &optional list) 
  #| Finite-Space whose 'top' is the topogenous matrix of a "weak core" of the submatrix 'list'x'list' of (top 'finspace') |#
  (let ((already (find `(WEAK-CORE ,finspace ,list) *finite-space-list* :test #'equal :key #'orgn)))
    (declare (type (or finite-space null) already))
    (if already
        already
      (let ((list2 (or list (<a-b> 1 (nlig (top finspace))))))
        (build-finite-space :orgn (if list `(WEAK-CORE ,finspace ,list) `(WEAK-CORE ,finspace)) 
                            :top (weakcore (top finspace) list2))))))


#|
  (show (core (randomtop 15 .3)))
  (show (weakcore (randomtop 15 .3)))
  (core (random-finite-space 15 .3))
  (top (core (random-finite-space 15 .3)))
  (top (weakcore (random-finite-space 15 .3)))

  (setf id (idnm (random-finite-space 20 .4)))
  (show (core (top (k id))))
  (show (top (core (k id)))) 
  (show (weakcore (top (k id))))
  (show (top (weakcore (k id))))   
|#
