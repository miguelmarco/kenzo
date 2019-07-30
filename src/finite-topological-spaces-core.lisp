;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*

;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 
;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 
;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 

(IN-PACKAGE #:cat)

(PROVIDE "finite-topological-spaces-core")

;;
;;   Point reductions: beat points and weak beat points
;;


(DEFUN DOWNBP (topogenous n list)
  #| T if Xn is a down beat point of the submatrix 'list'x'list' |#
  (the boolean
       (let ((maximal (loop for k in (reverse (subseq list 0 (position n list)))
                            thereis (and (eq (terme topogenous k n) 1) k))))
         
         (unless (or (null maximal) 
                     (loop for i in (subseq list 0 (position n list))
                           thereis (not (eq (terme topogenous i maximal) (terme topogenous i n)))))
           +TRUE+))))


(DEFUN UPBP (topogenous n list)
  #| T if Xn is an up beat point of the submatrix 'list'x'list' |#
  (the boolean
       (let ((minimal (loop for k in (subseq list (1+ (position n list)))
                            thereis (and (eq (terme topogenous n k) 1) k))))
         (unless (or (null minimal) 
                     (loop for i in (subseq list (1+ (position n list)))
                           thereis (not (eq (terme topogenous minimal i) (terme topogenous n i)))))
           +TRUE+))))


(DEFUN BEATPOINT (topogenous n list) 
  #| T if Xn is a beat point of the submatrix 'list'x'list' |#
  (the boolean
       (if (< (1+ n) (/ (ncol topogenous) 2))
           (if (downbp topogenous n list)
               +TRUE+
             (upbp topogenous n list))
         (if (upbp topogenous n list)
             +TRUE+
           (downbp topogenous n list)))))


(DEFUN CORE_LIST (topogenous list)
  #| List of the elements of the core of the submatrix 'list'x'list' |#
  (let ((eliminar NIL))
    (setf eliminar (loop for n in list
                         thereis (and (beatpoint topogenous n list) n)))
    (if (null eliminar)
        (return-from core_list list)
      (core_list topogenous (remove eliminar list)))))


(DEFUN Un-N (topogenous n list)
  #| List of the elements in Un-{Xn} |#
  (loop for i in (subseq list 0 (position n list))
       when (eq (terme topogenous i n) 1)
       collect i))


(DEFUN Fn-N (topogenous n list)
  #| List of the elements in Fn-{Xn} |#
  (loop for j in (subseq list (1+ (position n list)))
       when (eq (terme topogenous n j) 1)
       collect j))


(DEFUN LINK_LIST (topogenous n list)
  #| List of the elements in the link of Xn |#
  (append (Un-N topogenous n list) (Fn-N topogenous n list)))


(DEFUN WEAKPOINT (topogenous n list)
  #| T if Xn is a weak beat point of the submatrix 'list'x'list' |#
  (the boolean
       (or (eq (length (core_list topogenous (Un-N topogenous n list))) 1)
           (eq (length (core_list topogenous (Fn-N topogenous n list))) 1))))
  
  
(DEFUN WEAKCORE_LIST (topogenous list)
  #| List of the elements of a "weak core" of the submatrix 'list'x'list' |#
  (let ((eliminar NIL))
    (setf eliminar (loop for n in list
                         thereis (and (or (beatpoint topogenous n list) (weakpoint topogenous n list)) n)))
    (if (null eliminar)
        (return-from weakcore_list list)
      (weakcore_list topogenous (remove eliminar list)))))


(DEFMETHOD CORE ((topogenous matrice) &optional list)
  #| Matrice of the core of 'topogenous' |#
  (the matrice
       (let* ((list2 (if list
                         list
                       (<a-b> 1 (nlig topogenous))))
              (corelist (core_list topogenous list2)))
         (sub-matrix topogenous corelist corelist))))


(DEFMETHOD CORE ((finspace finite-space) &optional list)
  #| Finite-Space whose 'top' is the topogenous matrix of the core of (top 'finspace') |#
  (let ((already (find `(CORE ,finspace ,list) *finite-space-list* :test #'equal :key #'orgn)))
    (declare (type (or finite-space null) already))
    (if already
        already
      (build-finite-space :orgn `(CORE ,finspace ,list)
                          :top (core (top finspace) list)))))


(DEFMETHOD WEAKCORE ((topogenous matrice) &optional list)
  #| Matrice of a "weak core" of 'topogenous' |#
  (the matrice
       (let* ((list2 (if list
                         list
                       (<a-b> 1 (nlig topogenous))))
              (weakcorelist (weakcore_list topogenous (core_list topogenous list2))))
         (sub-matrix topogenous weakcorelist weakcorelist))))


(DEFMETHOD WEAKCORE ((finspace finite-space) &optional list)
  #| Finite-Space whose 'top' is the topogenous matrix of a "weak core" of (top 'finspace') |#
  (let ((already (find `(WEAK-CORE ,finspace ,list) *finite-space-list* :test #'equal :key #'orgn)))
    (declare (type (or finite-space null) already))
    (if already
        already
      (build-finite-space :orgn `(WEAK-CORE ,finspace ,list)
                          :top (weakcore (top finspace) list)))))


