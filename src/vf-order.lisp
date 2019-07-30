;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*

;;  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  
;;  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  
;;  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  VF-ORDER  

(IN-PACKAGE #:cat)

(PROVIDE "vf-order")

#|

Principe de l'algorithme.

L'ETAT est un vecteur indexé par les lignes.
Chaque ligne peut avoir pour valeur:
  1 = source (peut avoir des prédécesseurs et des successeurs, non notés)
  (...) = liste des prédécesseurs, éventuellement vide
        => non source utilisABLE à condition que l'ensemble des successeurs attribués
           soit disjoint de l'ensemble des prédécesseurs.

Etat initial : tout à ().

Les colonnes sont traitées successivement.

On cherche sur cette colonne un indice de ligne non source avec inversible.

La liste des successeurs désignés par la colonne doit avoir une intersection
vide avec la liste des précécesseurs.

Si oui:
  Cet indice devient source.
  Pour chaque enfant DE CETTE NOUVELLE SOURCE.
    Si c'est une source, tous les sans-enfants CITANT CETTE SOURCE sont mis à jour:
      ajout de cette nouvelle source et de ses prédécesseurs.
    Sinon cet enfant seul mis à jour avec 
      ajout de cette nouvelle source et de ses prédécesseurs.
  

|#


(DEFUN RANDOM-EXP (d cut)
  (floor (log (- 1 cut)) (log (- 1 d))))

#|
()
(dotimes (i 20)
  (print (random-exp 0.5 (random 1.0))))
|#

(DEFUN RANDOM-M (rown coln d1 d2)
  ; d1 = exponential-probability to get a 0
  ; d2 = probability a non-null is invertible.
  (let ((rslt (make-array coln :element-type 'list)))
    (dotimes (i coln)
      (let ((rslt2 '())
            (mark -1))
        (loop
          (incf mark (1+ (random-exp d1 (random 1.0))))
          (when (>= mark rown)
            (setf (svref rslt i) (nreverse rslt2))
            (return))
          (push (list mark (if (> (random 1.0) d2) 1 2))
                rslt2))))
    rslt))

(DEFUN DISPLAY-M (m)
  (dotimes (i (length m))
    (print (svref m i))))

#|
()
(display-m (random-m 5 10 0.5 0.5))
|#

(DEFUN NULL-INTERSECTION-P (list1 list2)
  (declare (type list list1 list2))
  (the boolean (progn
    (unless (and list1 list2) (return-from null-intersection-p t))
    (let ((mark1 list1)
          (mark2 list2))
      (declare (type list list1 list2))
      (loop
        (let ((delta (- (car mark1) (car mark2))))
          (declare (type fixnum delta))
          (if (plusp delta)
              (or (setf mark2 (cdr mark2))
                  (return-from null-intersection-p t))
            (if (minusp delta)
                (or (setf mark1 (cdr mark1))
                    (return-from null-intersection-p t))
              (return-from null-intersection-p nil)))))))))

#|
()
(null-intersection-p '(1 2 3) '())
(null-intersection-p '() '(4 5 6))
(null-intersection-p '(1 2 3) '(4 5 6))
(null-intersection-p '(4 5 6) '(1 2 3))
(null-intersection-p '(1 2 6) '(3 4 5))
(null-intersection-p '(1 2 3) '(3 4 5))
(null-intersection-p '(3 4 5) '(1 2 3))
|#

(DEFUN INSERT (n list)
  (declare (type fixnum n) (type list list))
  (the list
    (do ((rslt nil (cons (car mark) rslt))
         (mark list (cdr mark)))
        ((endp mark) (nreverse (cons n rslt)))
      (when (< n (car mark))
        (return-from insert
          (nreconc rslt (cons n mark)))))))

#|
()
(insert 4 '())
(insert 5 '(6 7 8))
(insert 5 '(3 7 8))
(insert 5 '(1 2 3))        
|#

(DEFUN UNION-MERGE (list1 list2)
  (declare (type list list1 list2))
  (the list
    (let ((rslt nil))
      (declare (type list rslt))
      (unless list1 (return-from union-merge list2))
      (unless list2 (return-from union-merge list1))
      (loop
        (let ((n1 (car list1))
              (n2 (car list2)))
          (declare (type fixnum n1 n2))
          (cond ((< n1 n2)
                 (push n1 rslt)
                 (or (setf list1 (cdr list1))
                     (return-from union-merge (nreconc rslt list2))))
                ((= n1 n2)
                 (push n1 rslt)
                 (or (setf list1 (cdr list1))
                     (return-from union-merge (nreconc rslt list2)))
                 (or (setf list2 (cdr list2))
                     (return-from union-merge (nreconc rslt list1))))
                (t
                 (push n2 rslt)
                 (or (setf list2 (cdr list2))
                     (return-from union-merge (nreconc rslt list1))))))))))

#|
()
(union-merge '() '(5 6 7))
(union-merge '(5 6 7) '())
(union-merge '(4 6 8) '(5 6 7))
(union-merge '(5 6 7) '(4 6 8))
|#

(DEFUN M-VF-HARD (rown m)
  (let ((status (make-array rown :initial-element '()))
        (coln (length m))
        (vf '()))
    (dotimes (i coln)
      (let ((coli (svref m i)))
        (when coli
          (let ((rowilist (mapcar #'car coli)))
            (declare (type list rowilist))
            (do ((rowilistmark rowilist (cdr rowilistmark))
                 (termlistmark (mapcar #'cadr coli) (cdr termlistmark)))
                ((endp rowilistmark))
              (let ((rowi (car rowilistmark)))
                (declare (type fixnum rowi))
                (unless (eql 2 (car termlistmark))
                  (let ((rowistatus (svref status rowi)))
                    (unless (eql 1 rowistatus)
                      (when (null-intersection-p rowistatus rowilist)
                        (setf (svref status rowi) 1)
                        (setf rowistatus (insert rowi rowistatus))
                        (push (list rowi i) vf)
                        (dolist (item rowilist)
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
                                 (union-merge rowistatus (svref status item)))))))
                        (return)))))))))))
    (return-from m-vf-hard (nreverse vf))))

(DEFUN M-VF (data)
  (destructuring-bind (coln rown m) data
      (declare (type fixnum rown) (ignore coln) (list m))
    (let ((status (make-array rown :initial-element 0))
          ;; 0 = unused
          ;; 1 = source
          ;; 2 = minimal
          (vf '()))
      (declare (type simple-vector status) (type list vf))
      (do ((mark m (cdr mark))
           (i 0 (1+ i)))
          ((endp mark))
          (declare (type list mark) (type fixnum i))
        (let ((coli (car mark)))
          (declare (type list coli))
          (when coli
            (let ((rowilist (mapcar #'car coli))
                  (termlist (mapcar #'cadr coli)))
              (declare (type list rowilist termlist))
              (do ((rowimark rowilist (cdr rowimark))
                   (termmark termlist (cdr termmark)))
                  ((endp rowimark))
                (declare (type list rowimark termmark))
                (let ((rowi (car rowimark)))
                  (declare (type fixnum rowi))
                  (block block1
                         (when (eql 2 (car termmark))
                           (return-from block1))
                         (let ((rowistatus (svref status rowi)))
                           (declare (type fixnum rowistatus))
                           (when (eql 1 rowistatus)
                             (return-from block1))
                           (when (eql 2 rowistatus)
                             (dolist (item rowilist)
                               (declare (type fixnum item))
                               (when (eql 1 (svref status item))
                                 (return-from block1))))
                           (push (list rowi i) vf)
                           (dolist (item rowilist)
                             (declare (type fixnum item))
                             (when (eql 0 (svref status item))
                               (setf (svref status item) 2)))
                           (setf (svref status rowi) 1)
                           (return)))))))))
      (return-from m-vf (nreverse vf)))))

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

#|
()
(setf d8 (delta 8))
(setf m (chcm-mat d8 3))
(m-vf-hard 84 (matrice-to-lmtrx m))

(setf m (list 6 6 '(
 ((0 2) (5 1)) 
 ((1 1) (5 2)) 
 ((0 2)) 
 ((1 1) (2 1) (5 1)) 
 ((3 1)) 
 ((2 2)))))
(m-vf m)

(setf m (make-array 6 :initial-contents '(
 ((0 1) (1 2) (4 2) (5 2)) 
 ((0 1)) 
 NIL 
 ((0 1) (1 1) (3 2)) 
 ((3 1) (5 2)) 
 ((3 2) (4 1)) )))
(m-vf-hard 6 m)

(setf m (make-array 8 :initial-contents '(  ;; bug m-vf
 ((0 1) (1 2) (2 2) (5 1) (6 2) (7 2))      ;; Différence entre m-vf et m-vf-soft
 ((0 2) (3 1) (5 2) (6 2)) 
 ((2 2) (3 2) (5 2) (6 2) (7 2)) 
 ((1 1) (3 2) (5 2) (6 2)) 
 ((0 2) (1 2) (2 1) (4 2) (5 2)) 
 ((1 1) (6 2)) 
 ((3 2) (5 2) (6 1) (7 1)) 
 ((1 1) (2 1) (6 2)))))

(display-m (setf m (random-m 8 8 0.5 0.5)))
(m-vf-hard 8 m)
(m-vf (list 8 8 (convertvector m)))
|#


#|
()
(setf m (make-array 6 :initial-contents '(
 ((0 2) (5 1)) 
 ((1 1) (5 2)) 
 ((0 2)) 
 ((1 1) (2 1) (5 1)) 
 ((3 1)) 
 ((2 2)))))

(setf m (make-array 6 :initial-contents '(
 ((0 1) (1 2) (4 2) (5 2)) 
 ((0 1)) 
 NIL 
 ((0 1) (1 1) (3 2)) 
 ((3 1) (5 2)) 
 ((3 2) (4 1)))))

(display-m (setf m (random-m 6 6 0.4 0.5)))
(m-vf 6 m)
|#

(DEFUN READ-D-2 (file)
  (declare (type string file))
  (the list
    (with-open-file (file file :direction :input)
      (declare (type stream file))
      (let ((rslt nil)
            (next nil))
        (declare (type list rslt next))
        (loop
          (setf next (read file nil nil))
          (if next
              (push next rslt)
            (return)))
        (return-from read-d-2 (nreverse rslt))))))
#|
()
(read-d "Vector-Fields\\D6.txt")
|#

(DEFUN D-M-1 (coln d)
  ;; (col-i row-i entry-i)
  (declare (type list d) (type fixnum coln))
  (the simple-vector
    (let ((rslt (make-array coln :element-type 'list :initial-element nil)))
      (declare (type simple-vector rslt))
      (dolist (item d)
        (declare (type list item))
        (destructuring-bind (coli rowi vali) item
          (declare (type fixnum coli rowi) (ignore vali))
          (push (list rowi 1) (svref rslt coli))))
      (return-from d-m-1 rslt))))

(DEFUN READ-D (file)
  (declare (type string file))
  (the list
    (with-open-file (file file :direction :input)
      (declare (type stream file))
      (let ((rslt nil)
            (subrslt nil)
            (next nil)
            (coli 0)
            (rowmax 0))
        (declare (type list rslt subrslt next) (type fixnum coli))
        (loop
          (setf next (read file nil nil))
          (unless next
            (push (nreverse subrslt) rslt)
            (return-from read-d (list (1+ coli) (1+ rowmax) (nreverse rslt))))
          (destructuring-bind (colj rowi vali) next
            (declare (type fixnum colj rowi) (ignore vali))
            (setf rowmax (max rowmax rowi))
            (cond ((eql coli colj)
                   (push (list rowi 1) subrslt))
                  (t
                   (push (nreverse subrslt) rslt)
                   (incf coli)
                   (setf subrslt (list (list rowi 1)))))))))))

#|
()
(setf d (read-d "Vector-Fields\\D1.txt"))
(time (m-vf d))
(length *)
(setf d (read-d "Vector-Fields\\D2.txt"))
(time (m-vf d))
(length *)
(setf d (read-d "Vector-Fields\\D3.txt"))
(time (m-vf d))
(length *)
(setf d (read-d "Vector-Fields\\D4.txt"))
(time (m-vf d))
(length *)
(setf d (read-d "Vector-Fields\\D5.txt"))
(time (m-vf d))
(length *)
(setf d (read-d "Vector-Fields\\D6.txt"))
(time (m-vf d))
(length *)
|#
