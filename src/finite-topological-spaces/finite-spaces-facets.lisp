
;;  FACETS FACETS FACETS FACETS FACETS FACETS FACETS FACETS FACETS FACETS
;;  FACETS FACETS FACETS FACETS FACETS FACETS FACETS FACETS FACETS FACETS
;;  FACETS FACETS FACETS FACETS FACETS FACETS FACETS FACETS FACETS FACETS 

(IN-PACKAGE #:cat)

(provide "finite-spaces-facets")

;;
;;  Construction of facet posets from facets data
;;


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFUN MAX-DIMENSION (facets)
  #| Dimension of the simplicial complex given by 'facets' |#
  (let ((max_dim 0))
    (dolist (facet facets)
      (if (> (length facet) max_dim)
          (setf max_dim (length facet)))) max_dim))


(DEFUN SUBLISTS (list n)
  #| Sublists of 'list' of length 'n' |#
  (cond ((zerop n) nil)
        ((null list) nil)
        ((= n 1) (mapcar #'list list))
        (T (append (mapcar #'(lambda (subset) (cons (car list) subset))
                         (sublists (cdr list) (1- n)))
                 (sublists (cdr list) n)))))


(DEFUN ALL-FACES (facets max_dim)
  #| All the faces of length up to 'max_dim' of the simplicial complex given by 'facets' |#
  (let ((faces (loop for i from 1 to max_dim collect '())))
    (dotimes (i max_dim)
      (dolist (facet facets)
        (unless (< (length facet) i)
          (setf (nth i faces) (union (nth i faces) (sublists facet (1+ i)) :test #'equal)))))
    faces))


(DEFUN REMOVE-NTH (n list)
  #| Remove the 'n' position of the 'list' |#
  (if (zerop n)
      (cdr list)
    (cons (car list) (remove-nth (1- n) (cdr list)))))


(DEFUN FACES-TO-UBASIS (allfaces)
  #| Construct the ubasis for the Stong matrice associated to 'allfaces' and the heights of the space|#
  (let ((dimension 0)
        (heights NIL))
    (dolist (level allfaces)
      (let ((past dimension))
        (setf dimension (+ dimension (length level))
              heights (push (<a-b> (1+ past) dimension) heights))))
    (let ((rslt (make-array dimension)))
      (do* ((n 1 (1+ n))
            (sum_n-1 0 (+ sum_n-1 (length (nth (- n 2) allfaces))))
            (sum_n (length (car allfaces)) (+ sum_n (length (nth (1- n) allfaces)))))
           ((= n (length allfaces)) (list rslt (nreverse heights)))
        (let ((level_n-1 (nth (1- n) allfaces))
              (level_n (nth n allfaces)))
          (dotimes (k (length level_n))
            (let ((face (nth k level_n)))
              (setf (svref rslt (+ sum_n k))
                    (append (svref rslt (+ sum_n k))
                            (loop for rmv from 0 to n ; An n-simplex has (n+1) co-faces
                                  collect (+ sum_n-1 (1+ (position (remove-nth rmv face) level_n-1 :test #'equal)))))))))))))

#|
  (setf facets (list (<a-b> 1 4)))
  (faces-to-ubasis (all-faces facets 4))
|#

(DEFUN FACETS-TO-FINITE-SPACE (facets)
  (let ((info (faces-to-ubasis (all-faces facets (max-dimension facets)))))
    (build-finite-space :stong (nilpot-1 (ubasis-to-binarymatrice (car info)))
                        :heights (second info)
                        :orgn `(FACETS ,facets))))


(DEFUN IMPORT-FACETS-TO-FINITE-SPACE (string_file &key random)
  (let ((facets (with-open-file (jlcr (concatenate 'string data_folder string_file ".txt"))
                  (read jlcr))))
    (when random (setf facets (fisher-yates facets)))
    (facets-to-finite-space facets)))

    
(DEFMETHOD DVFIELD-FACETS ((stong matrice) &key random)
  #| Compute a discrete vector field on 'stong' |#
  (let* ((m_stong (binarymatrice-to-ubasis (nilpot stong)))
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
                      (when (null-intersection-p rowistatus coljlist) ; avoids a loop  x --> i -> j -> x
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
    (return-from DVFIELD-FACETS (nreverse vf))))

    
(DEFMETHOD DVFIELD-FACETS ((finspace finite-space) &key random)
   (dvfield-facets (stong finspace) :random random))
