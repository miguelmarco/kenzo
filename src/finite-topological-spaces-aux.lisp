;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*

;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 
;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 
;;  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES  FINITE-TOPOLOGICAL-SPACES 

(IN-PACKAGE #:cat)

(PROVIDE "finite-topological-spaces-aux")

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

