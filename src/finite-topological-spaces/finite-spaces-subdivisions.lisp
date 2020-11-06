
;;  BARYCENTRIC-SUBDIVISION BARYCENTRIC-SUBDIVISION BARYCENTRIC-SUBDIVISION 
;;  BARYCENTRIC-SUBDIVISION BARYCENTRIC-SUBDIVISION BARYCENTRIC-SUBDIVISION 
;;  BARYCENTRIC-SUBDIVISION BARYCENTRIC-SUBDIVISION BARYCENTRIC-SUBDIVISION 

(IN-PACKAGE #:cat)

(provide "finite-spaces-subdivisions")

;;
;;  Computing the barycentric subdivision of a finite space
;;


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


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
  #| Matrices relating consecutive dimensions in the order complex |#
  (let* ((relac (relac topogenous))
         (diferentials NIL)
         (nuevo relac)
         (dim (nlig topogenous))
         (stop 0)
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
            D1) diferentials)
            
    (if (> dim 1)
        (do ((long 2 (1+ long))) ((or (> long dim) (eq stop 1)))
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
                (setf stop 1)
              (push (let ((Dn (creer-matrice dim_n-1 dim_n)))
                      (dotimes (i dim_n-1)
                        (dotimes (j dim_n)
                          (if (subsetp (nth i Cn-1) (nth j Cn))
                              (insert-term Dn (1+ i) (1+ j) 1))))
                      Dn) diferentials)))))
    (return-from BOUNDARYOPERATORS (reverse diferentials))))


(DEFUN BLOQUES (lista)
  (let ((rslt (identite (let ((suma (nlig (car lista))))
                          (dolist (x lista)
                            (setf suma (+ suma (ncol x))))
                          suma))))
    (do* ((lst lista (cdr lst))
          (kfilas 0 (if (null lst) 0 (+ kfilas (nlig dif))))
          (kcolumnas (nlig (car lista)) (if (null lst) 0 (+ kcolumnas (ncol dif))))
          (dif (car lst) (car lst))) ((endp lst))
      
      (do ((i 1 (1+ i))) ((> i (nlig dif)))
        (do ((j 1 (1+ j))) ((> j (ncol dif)))
          (unless (eq (terme dif i j) 0)
            (insert-term rslt (+ i kfilas) (+ j kcolumnas) 1)))))
      rslt))


#|
  ((lambda (m n p)
     (let ((M (randomtop m .5))
           (N (randomtop n .5))
           (P (randomtop p .5)))
       (show M)
       (show N)
       (show P)
       (show (bloques (list M N P))))) 4 4 4)
|#


(DEFMETHOD BAR-SUBDIVISION ((topogenous matrice))
  #| Topogenous matrice of the barycentric subdivision of 'topogenous' |#
  (topmat (bloques (boundaryoperators topogenous))))


(DEFMETHOD BAR-SUBDIVISION ((finspace finite-space)) 
  #| Finite-Space whose 'top' is the barycentric subdivision of (top 'finspace') |#
  (let ((already (find `(BARYCENTRIC-SUBDIVISION ,finspace) *finite-space-list* :test #'equal :key #'orgn)))
    (declare (type (or finite-space null) already))
    (if already
        already
      (build-finite-space :orgn `(BARYCENTRIC-SUBDIVISION ,finspace)
                          :stong (bloques (boundaryoperators (top finspace)))))))


#|
  (bar-subdivision (randomtop 6 .5))
  (bar-subdivision (random-finite-space 6 .5))
  (bar-subdivision (finite-model-sphere 4))
  (bar-subdivision (finite-model-sphere 3))
|#

