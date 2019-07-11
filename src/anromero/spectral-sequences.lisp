;; SPECTRAL-SEQUENCES   SPECTRAL-SEQUENCES   SPECTRAL-SEQUENCES   SPECTRAL-SEQUENCES
;; SPECTRAL-SEQUENCES   SPECTRAL-SEQUENCES   SPECTRAL-SEQUENCES   SPECTRAL-SEQUENCES
;; SPECTRAL-SEQUENCES   SPECTRAL-SEQUENCES   SPECTRAL-SEQUENCES   SPECTRAL-SEQUENCES


(IN-PACKAGE #:cat)

(provide "spectral-sequences")


;; Auxiliary functions to deal with matrices

(DEFUN SUBMATRIX (mtrx first-line last-line first-col last-col)
  (declare (type matrix mtrx)
           (type fixnum first-line last-line first-col last-col))
  (the matrix
    (let ((line-n (line-number mtrx))
          (column-n (column-number mtrx)))
      (declare (type fixnum line-n column-n))
      (assert (<= first-line last-line (1- line-n)))
      (assert (<= first-col last-col (1- column-n)))
      (let* ((line-num (+ (- last-line first-line) 1))
             (column-num (+ (- last-col first-col) 1))
             (rslt (make-array (list line-num column-num)
                               :element-type 'fixnum)))
        (declare (type fixnum line-num column-num) 
                 (type matrix rslt))
        (dotimes (il line-num)
          (declare (type fixnum il))
          (dotimes (ic column-num)
            (declare (type fixnum ic))
            (setf (aref rslt il ic)
              (aref mtrx (+ il first-line) (+ ic first-col)))))
        rslt))))

#|
(setf m (random-matrix 4 6 10))
(submatrix m 1 3 2 4)
|#


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

#|
(setf m (random-matrix 4 6 10))
(column-j m 2)
|#


;; Function that concatenates two matrices (if their line number is different, the 
;; needed zeros are added. 
(DEFUN MTRX-CONC (mtrx1 mtrx2)
  (declare (type matrix mtrx1 mtrx2))
  (the matrix
    (let* ((line-n-1 (line-number mtrx1))
           (column-n-1 (column-number mtrx1))
           (line-n-2 (line-number mtrx2))
           (column-n-2 (column-number mtrx2))
           (line-n (max line-n-1 line-n-2))
           (column-n (+ column-n-1 column-n-2)))
      (declare (type fixnum line-n-1 column-n-1 line-n-2
                     column-n-2 line-n column-n))
      (let ((rslt 
             #-ACLPC
             (make-array (list line-n column-n)
                         :element-type 'fixnum
                         :initial-element 0)
             #+ACLPC
             (if (or (zerop line-n) (zerop column-n))
                 (make-array (list line-n column-n)
                             :element-type 'fixnum)
               (make-array (list line-n column-n)
                           :element-type 'fixnum
                           :initial-element 0))))
        (declare (type matrix rslt))
        (do ((j 0 (1+ j)))
            ((>= j column-n-1))
          (declare (type fixnum j))
          (do ((i 0 (1+ i)))
              ((>= i line-n-1))
            (declare (type fixnum i))
            (setf (aref rslt i j) (aref mtrx1 i j))))
        (do ((j 0 (1+ j)))
            ((>= j column-n-2))
          (declare (type fixnum j))
          (do ((i 0 (1+ i)))
              ((>= i line-n-2))
            (declare (type fixnum i))
            (setf (aref rslt i (+ column-n-1 j)) (aref mtrx2 i j))))
        rslt))))

#|
(setf m1 (random-matrix 3 4 10))
(setf m2 (random-matrix 4 2 10))
(mtrx-conc m1 m2)
(mtrx-conc m2 m1)
|#

 
;; Functions related with submodules and matrices

;; Function that returns a list of generators of the kernel of a matrix.
;; Each generator is a list of integers which represent its coordinates
;; in the initial basis.
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

#|
(setf m (random-matrix 3 4 10))
(kernel m)
|#


;; Returns a matrix which has as columns the elements of the list.
;; (Each element in the list is a list of integers).
(DEFUN GNRT-LIST-TO-MTRX (gnrt-list)
  (declare (type list gnrt-list))
  (the matrix
    (let* ((column-n (length gnrt-list))
           (line-n (length (car gnrt-list)))
           (rslt 
            #-ACLPC
            (make-array (list line-n column-n)
                        :element-type 'fixnum
                        :initial-element 0)
            #+ACLPC
            (if (or (zerop line-n) (zerop column-n))
                (make-array (list line-n column-n)
                            :element-type 'fixnum)
              (make-array (list line-n column-n)
                          :element-type 'fixnum
                          :initial-element 0))))
      (declare 
       (type fixnum column-n line-n)
       (type matrix rslt))
      (do ((j 0 (1+ j))
           (mark2 gnrt-list (cdr mark2)))
          ((endp mark2))
        (declare
         (type fixnum j)
         (list mark2)) 
        (if (eql line-n (length (car mark2)))
            (do ((i 0 (1+ i))
                 (mark1 (car mark2) (cdr mark1)))
                ((endp mark1))
              (declare 
               (type fixnum i)
               (list mark1))
              (setf (aref rslt i j) (car mark1)))))
      rslt)))

#|
(setq gnrt-list (list '(2 4 6 8) '(0 1 0 0) '(0 3 3 3)))
(gnrt-list-to-mtrx gnrt-list)
|# 

;; Returns the list of columns of the matrix mtrx.
(DEFUN GNRT-MTRX-TO-LIST (mtrx)
  (declare (type matrix mtrx))
  (the list
    (let* ((column-n (column-number mtrx)))
      (declare (type fixnum column-n))
      (mapcar
          #'(lambda (j)
              (declare (type fixnum j))
              (column-j mtrx j))
        (<a-b> 0 (1- column-n))))))

#|
(setf m (random-matrix 3 4 10))
(gnrt-mtrx-to-list m)
|#

;; Returns the representation basis-divisors of a quotient of submodules B/C
;; which are the images of the matrices mtrx1 and mtrx2, where basis=list
;; of combinations (c_1...c_m) that forms a basis for B and divisors=list of
;; integers (n_1....n_m) such that (n_1*c_1,...,n_m*c_m) forms a basis for C. 
(DEFUN MTRX-QUOTIENT (mtrx1 mtrx2)
  (declare (type matrix mtrx1 mtrx2))
  (the list
    (let* ((smith-list-1 (smith mtrx1))
           (p1 (first smith-list-1))
           (p1-1 (second smith-list-1))
           (d1 (third smith-list-1))
           (p1-1Xmtrx2 (mtrx-prdc p1-1 mtrx2))
           (p1Xd1 (copy-mtrx p1))
           (d1-1Xp1-1Xmtrx2 (copy-mtrx p1-1Xmtrx2)))
      (declare
       (list smith-list-1)
       (type matrix p1 p1-1 d1 p1-1Xmtrx2 p1Xd1 d1-1Xp1-1Xmtrx2))
      (progn
        (let* ((line-n-1 (line-number p1))
               (column-n-1 (column-number p1)))
          (declare (type fixnum line-n-1 column-n-1))
          (dotimes (ic column-n-1)
            (declare (type fixnum ic))
            (if (and (< ic (line-number d1)) (< ic (column-number d1)))
                (let ((d (aref d1 ic ic)))
                  (declare (type fixnum d))
                  (dotimes (il line-n-1)
                    (declare (type fixnum il))
                    (if (and (not (eql 1 d)) (not (eql 0 d)))
                        (setf (aref p1Xd1 il ic) (* d (aref p1Xd1 il ic)))))))))
        (let* ((line-n-2 (line-number p1-1Xmtrx2))
               (column-n-2 (column-number p1-1Xmtrx2)))
          (declare (type fixnum line-n-2 column-n-2))    
          (dotimes (il line-n-2)
            (declare (type fixnum il))
            (if (and (< il (line-number d1)) (< il (column-number d1)))
                (let ((d (aref d1 il il)))
                  (if (and (not (eql 1 d)) (not (eql 0 d)))
                      (dotimes (ic column-n-2)
                        (declare (type fixnum ic))
                        (setf (aref d1-1Xp1-1Xmtrx2 il ic) 
                          (floor (aref d1-1Xp1-1Xmtrx2 il ic) d))))))))
        (let* ((smith-list-2 (smith d1-1Xp1-1Xmtrx2))
               (p2 (first smith-list-2))
               (basis-mtrx (mtrx-prdc p1Xd1 p2))
               (div-mtrx (third smith-list-2))
               (cmbns nil)
               (divs nil))
          (declare 
           (list smith-list-2 cmbns divs)
           (type matrix p2 basis-mtrx div-mtrx))
          (do ((i 0 (1+ i)))
              ((or (>= i (line-number d1)) (>= i (column-number d1))
                   (eql 0 (aref d1 i i))))
            (declare (type fixnum i))
            (progn
              (setq cmbns (nconc cmbns (list (column-j basis-mtrx i))))
              (if (and (< i (line-number div-mtrx)) (< i (column-number div-mtrx))) 
                  (setq divs (nconc divs (list (aref div-mtrx i i))))
                (setq divs (nconc divs (list 0))))))
          (list cmbns divs))))))

#|
(setf m1 (gnrt-list-to-mtrx '((1 1 1 1) (0 2 0 3) (1 0 0 0))))
(setf m2 (gnrt-list-to-mtrx '((2 3 1 4) (0 1 1 1))))
(mtrx-quotient m1 m2)
|#



;; SPECIFIC FUNCTIONS FOR THE COMPUTATION OF SPECTRAL SEQUENCES 
;; OF EFFECTIVE FILTERED COMPLEXES.


;; GROUPS

;; Returns the matrix that defines Z^r_{p,q} (a submatrix of the 
;; differential matrix M, such that the elements of Z^r_{p,q} are those that
;; when we apply M to them we obtain zero, i.e. Z^r_{p,q}=Ker M). Only can be applied
;; to a finitely generated complex. 
(DEFUN FLTR-CHCM-Z-MTRX (fltrcm r p q)
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum r p q))
  (the matrix
    (let* ((cmpr (cmpr1 fltrcm))
           (dffr (dffr fltrcm))
           (degr (+ p q))
           (sbasis (fltrd-basis fltrcm degr p))
           (tbasis 
            (mapcan 
                #'(lambda (i)
                    (declare (type fixnum i))
                    (bigrd-basis fltrcm i (- (1- degr) i)))
              (<a-b> (1+ (- p r)) p)))) ;; We suppose that the differential
      ;; respects the filtration, i.e., 
      ;; d(F_pX) \in F_p(dX). If not, we must
      ;; consider (<a-b> (1+ (-p r)) t(degr),
      ;; where t(degr)=max{flin(x), x \in X_degr}.
      (declare
       (type cmprf cmpr)
       (type morphism dffr)
       (type fixnum degr)
       (list sbasis tbasis))
      (let ((srank (length sbasis))
            (trank (length tbasis)))
        (declare (type fixnum srank trank))
        (let ((rslt 
               #-ACLPC
               (make-array (list trank srank)
                           :element-type 'fixnum
                           :initial-element 0)
               #+ACLPC
               (if (or (zerop srank) (zerop trank))
                   (make-array (list trank srank)
                               :element-type 'fixnum)
                 (make-array (list trank srank)
                             :element-type 'fixnum
                             :initial-element 0))))               
          (declare (type matrix rslt))
          (do ((j 0 (1+ j))
               (mark sbasis (cdr mark)))
              ((endp mark))
            (declare
             (type fixnum j)
             (list mark))
            (let ((cmbn (gnrt-? dffr degr (car mark))))
              (declare (type cmbn cmbn))
              (do ((mark1 (cmbn-list cmbn) (cdr mark1)))
                  ((endp mark1))
                (declare (list mark1))
                (with--term (cffc gnrt) mark1
                  (declare 
                   (type fixnum cffc)
                   (type gnrt gnrt))
                  (do ((mark2 tbasis (cdr mark2))
                       (i 0 (1+ i))
                       (found nil))
                      ((or (endp mark2) found))
                    (declare 
                     (list mark2)
                     (type fixnum i ))
                    (if (eq :equal (funcall cmpr gnrt (car mark2)))
                        (progn
                          (setq found 1)
                          (setf (aref rslt i j) cffc))))))))
          rslt)))))



(DEFUN FLTR-CHCM-Z-MTRX2 (fltrcm r p q)
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum r p q))
  (the matrix
    (let* ((degr (+ p q))
           (sbasis (fltrd-basis fltrcm degr p))
           (tbasis 
            (mapcan 
                #'(lambda (i)
                    (declare (type fixnum i))
                    (bigrd-basis fltrcm i (- (1- degr) i)))
              (<a-b> (1+ (- p r)) p)))) ;; We suppose that the differential
      ;; respects the filtration, i.e., 
      ;; d(F_pX) \in F_p(dX). If not, we must
      ;; consider (<a-b> (1+ (-p r)) t(degr),
      ;; where t(degr)=max{flin(x), x \in X_degr}.
      (declare
       (type fixnum degr)
       (list sbasis tbasis))
      (let ((srank (length sbasis))
            (trank (length tbasis)))
        (declare (type fixnum srank trank))
        (let ((rslt 
               #-ACLPC
               (make-array (list trank srank)
                           :element-type 'fixnum
                           :initial-element 0)
               #+ACLPC
               (if (or (zerop srank) (zerop trank))
                   (make-array (list trank srank)
                               :element-type 'fixnum)
                 (make-array (list trank srank)
                             :element-type 'fixnum
                             :initial-element 0))))               
          (declare (type matrix rslt))
          
          rslt)))))


#|
(cat-init)
(setq kz2 (k-z2-1))
(setq fcc (bar kz2))
(Change-Chcm-TO-Flcc fcc abar-flin `(abar-flin))
(fltr-chcm-z-mtrx fcc 2 3 2)
(fltr-chcm-z-mtrx fcc 1 2 4)
|#
 

;; Returns a list of generators (list of coordinates) of Z^r_{p,q}
;; (only for finitely generated complexes)
(DEFUN FLTR-CHCM-Z-GNRT-LIST (FltrCm r p q)
  (declare (type filtered-chain-complex fltrcm)
           (type fixnum r p q))
  (the list
    (let* ((degr (+ p q))
           (all-basis-gnrt-list (fltr-basis-gnrt-list fltrcm degr p)))
      (declare
       (type fixnum degr)
       (list all-basis-gnrt-list))
      (if (= 0 r)
          all-basis-gnrt-list
        (let* ((mat (fltr-chcm-z-mtrx fltrcm r p q))
               (line-n (line-number mat))
               (column-n (column-number mat)))
          (declare 
           (type matrix mat)
           (type fixnum line-n column-n))
          (if (= 0 column-n) nil
            (if (= 0 line-n)
                all-basis-gnrt-list
              (kernel mat))))))))
#|
(fltr-chcm-z-gnrt-list fcc 2 3 2)
(fltr-chcm-z-gnrt-list fcc 1 2 4)
|# 

;; Returns a list of lists with the coordinates of the elements of
;; the basis of F_pX_n (with p=fltr-index, n=degr), that are lists with
;; the form (0...0 i 0...0) (only for finitely generated complexes)
(DEFUN FLTR-BASIS-GNRT-LIST (fltrcm degr fltr-index)
  (declare 
   (type filtered-chain-complex fltrcm)
   (type fixnum degr fltr-index))
  (the list
    (let* ((all-basis (fltrd-basis fltrcm degr fltr-index))
           (basis-l (length all-basis))
           (rslt 
            (the list
              (mapcar
                  #'(lambda (i)
                      (nconc (make-list (1- i) :initial-element 0) 
                             (list 1) (make-list (- basis-l i) :initial-element 0)))
                (<a-b> 1 basis-l)))))
      rslt)))

#|
(fltr-basis-gnrt-list fcc 6 2)
|#
 

;; Returns the matrix of the generators of the numerator in E^r_{p,q}
;; (Only for finitely generated complexes)
(DEFUN SPSQ-NUM-MTRX (fltrcm r p q)
  (the matrix
    (let* ((degr (+ p q))
           (fltr-p-basis (fltrd-basis fltrcm degr p))
           (p-basis-l (length fltr-p-basis))
           (z1-list (fltr-chcm-z-gnrt-list fltrcm r p q))
           (fltr-p-1-list (fltr-basis-gnrt-list fltrcm degr (1- p)))
           (z1-mtrx (gnrt-list-to-mtrx z1-list))
           (fltr-p-1-mtrx (gnrt-list-to-mtrx fltr-p-1-list))
           (num1-mtrx (mtrx-conc z1-mtrx fltr-p-1-mtrx))
           (num1-line-n (line-number num1-mtrx)))
      (declare
       (type fixnum degr p-basis-l num1-line-n)
       (type list fltr-p-basis  z1-list fltr-p-1-list)
       (type matrix z1-mtrx fltr-p-1-mtrx num1-mtrx))
      (if (eql num1-line-n p-basis-l)
          num1-mtrx 
        (let* ((num-line-n p-basis-l)
               (num-column-n (column-number num1-mtrx))
               (rslt  
                #-ACLPC (make-array (list num-line-n num-column-n)
                                    :element-type 'fixnum
                                    :initial-element 0)
                #+ACLPC
                (if (or (zerop num-line-n) (zerop num-column-n))
                    (make-array (list num-line-n num-column-n)
                                :element-type 'fixnum)
                  (make-array (list num-line-n num-column-n)
                              :element-type 'fixnum
                              :initial-element 0))))
          (declare
           (type fixnum num-line-n num-column-n)
           (type matrix rslt))
          (do ((j 0 (1+ j)))
              ((>= j num-column-n))
            (declare (type fixnum j))
            (do ((i 0 (1+ i)))
                ((>= i num1-line-n))
              (declare (type fixnum i))
              (setf (aref rslt i j) (aref num1-mtrx i j))))
          rslt)))))

#|
(spsq-num-mtrx fcc 2 3 2)
(spsq-num-mtrx fcc 1 2 4)
|# 

;; Returns the matrix of the generators of the denominator in E^r_{p,q}
;; (Only for finitely generated complexes) 
(DEFUN SPSQ-DEN-MTRX (fltrCm r p q)
  (the matrix
    (let* ((degr (+ p q))
           (fltr-p-basis (fltrd-basis fltrcm degr p))
           (p-basis-l (length fltr-p-basis))
           (fltr-p-1-list (fltr-basis-gnrt-list fltrcm degr (1- p)))
           (z2-list (fltr-chcm-z-gnrt-list fltrcm (1- r) (1- (+ p r)) (+ 2 (- q r))))
           (z2-mtrx (gnrt-list-to-mtrx z2-list))
           (fltr-p-1-mtrx (gnrt-list-to-mtrx fltr-p-1-list))
           (dffr-mtrx (flcc-dffr-mtrx fltrcm (1+ degr) (1- (+ p r))))
           (nil-mtrx 
            #-ACLPC (make-array (list 0 0)
                                :element-type 'fixnum
                                :initial-element 0)
            #+ACLPC (make-array (list 0 0)
                                :element-type 'fixnum))
           (dffr-z2-mtrx 
            (if (or (= 0 (array-total-size z2-mtrx)) (= 0 (array-total-size dffr-mtrx)))
                nil-mtrx
              (mtrx-prdc dffr-mtrx z2-mtrx)))
           (bnd-z2-line-n p-basis-l)
           (bnd-z2-column-n (column-number dffr-z2-mtrx))
           (bnd-z2-mtrx
            (if (or (= 0 bnd-z2-line-n) (= 0 bnd-z2-column-n))
                nil-mtrx
              (submatrix dffr-z2-mtrx 0 (1- bnd-z2-line-n) 0 (1- bnd-z2-column-n))))
           (den1-mtrx (mtrx-conc bnd-z2-mtrx fltr-p-1-mtrx))
           (den-column-n (column-number den1-mtrx))
           (den-line-n p-basis-l)
           (rslt  
            #-ACLPC (make-array (list den-line-n den-column-n)
                                :element-type 'fixnum
                                :initial-element 0)
            #+ACLPC
            (if (or (zerop den-line-n) (zerop den-column-n))
                (make-array (list den-line-n den-column-n)
                            :element-type 'fixnum)
              (make-array (list den-line-n den-column-n)
                          :element-type 'fixnum
                          :initial-element 0))))
      (declare
       (list fltr-p-basis fltr-p-1-list z2-list)
       (type fixnum degr p-basis-l bnd-z2-line-n bnd-z2-column-n den-column-n den-line-n)
       (type matrix z2-mtrx fltr-p-1-mtrx dffr-mtrx nil-mtrx 
             dffr-z2-mtrx bnd-z2-mtrx den1-mtrx rslt))
      (do ((j 0 (1+ j)))
          ((>= j den-column-n))
        (declare (type fixnum j))
        (do ((i 0 (1+ i)))
            ((>= i (min den-line-n (line-number den1-mtrx))))
          (declare (type fixnum i))
          (Setf (aref rslt i j) (aref den1-mtrx i j))))
      rslt)))

#|
(spsq-den-mtrx fcc 2 3 2)
(spsq-den-mtrx fcc 1 2 4)
|#   

;; Function that returns the representation basis-divisors of E^r_{p,q} of a
;; finitely generated complex. 
(DEFUN EFF-SPSQ-BASIS-DVS (FltrCm r p q)
  (declare 
   (type filtered-chain-complex fltrcm)
   (type fixnum r p q))
  (the list
    (let* ((degr (+ p q))
           (fltr-p-basis (fltrd-basis fltrcm degr p))
           (num-mtrx (spsq-num-mtrx fltrcm r p q))
           (den-mtrx (spsq-den-mtrx fltrcm r p q)) 
           (cmbn-list nil)
           (divs nil))
      (declare
       (list fltr-p-basis cmbn-list divs)
       (type fixnum degr)
       (type matrix num-mtrx den-mtrx))
      (progn
        (if (= 0 (array-total-size num-mtrx))
            (progn
              (setq cmbn-list nil)
              (setq divs nil))
          (if (= 0 (array-total-size den-mtrx))
              (progn
                (setq cmbn-list (gnrt-mtrx-to-list num-mtrx))
                (setq divs (make-list (length cmbn-list) :initial-element 0)))
            (let ((quotient-list (mtrx-quotient num-mtrx den-mtrx)))
              (declare (list quotient-list))
              (setq cmbn-list (first quotient-list)
                  divs (second quotient-list)))))
        (list
         (mapcar 
             #'(lambda (int-list)
                 (declare (list int-list))
                 (the cmbn
                   (let* ((cmbn (cmbn degr))
                          (cmpr (cmpr fltrcm)))
                     (declare 
                      (type cmbn cmbn)
                      (type cmprf cmpr))
                     (do ((mark1 int-list (cdr mark1))
                          (mark2 fltr-p-basis (cdr mark2)))
                         ((or (endp mark1) (endp mark2)))
                       (declare (type list mark1 mark2))
                       (if (not (= 0 (car mark1)))
                           (setq cmbn (2cmbn-add cmpr cmbn 
                                                 (cmbn degr (car mark1) (car mark2))))))
                     cmbn))) 
           cmbn-list)
         divs)))))

#|
(do ((r 1 (1+ r)))
    ((> r 3))
  (dotimes (n 8)
    (dotimes (p (+ n 1))
      (let ((q (- n p)))
        (print (list r p q))
        (princ (eff-spsq-basis-dvs fcc r p q))
        ))))
|# 


(DEFUN EFF-SPSQ-GNRTS (FltrCm r p q)
  (declare 
   (type filtered-chain-complex fltrcm)
   (type fixnum r p q))
  (let* ((basis-dvs (eff-spsq-basis-dvs fltrcm r p q))
         (basis (first basis-dvs))
         (dvs (second basis-dvs))
         (rslt nil))
    (declare (type list basis-dvs basis dvs rslt))
    (do ((mark1 basis (cdr mark1))
         (mark2 dvs (cdr mark2)))
        ((endp mark1))
      (if (not (eql 1 (car mark2)))
          (setf rslt (append rslt (list (car mark1))))))
    rslt))


;; Function that presents on the screen the components of the E^r_{p,q}
;; of the spectral sequence of a filtered complex.
(DEFUN EFF-SPSQ-GROUP (fltrcm r p q)
  (declare 
   (type filtered-chain-complex fltrcm)
   (type fixnum r p q))  
  (let* ((basis-dvs (eff-spsq-basis-dvs fltrcm r p q))
         (divs (second basis-dvs)))
    (declare (list basis-dvs divs))
    (format t "Spectral sequence E^~D_{~D,~D}" r p q)
    (dolist (item divs)
      (declare (type fixnum item))
      (if (not (eql 1 item))
          (progn
            (format t "~2%Component Z")
            (unless (zerop item) 
              (format t "/~DZ" item))
            )))
    (terpri)))

#|
(do ((r 1 (1+ r)))
    ((> r 3))
  (dotimes (n 8)
    (dotimes (p (+ n 1))
      (let ((q (- n p)))
        (eff-spsq-group fcc r p q)
        (terpri)
        (terpri)
        ))))
|# 



;; DIFFERENTIAL FUNCTIONS

;; Auxiliary functions to deal with coefficients of elements with regard to 
;; some generators or elements of a basis. 


;; Function that returns the coordinates of a vector vctr with regard to
;; the basis whose elements are the columns of the matrix gnrt-mtrx.
(DEFUN VCTR-COORDINATES (vctr gnrt-mtrx)
  (declare 
   (type list vctr)
   (type matrix gnrt-mtrx))
  (the list
    (let* ((C (gnrt-list-to-mtrx (list vctr)))
           (B (copy-mtrx gnrt-mtrx))
           (smith-list (smith B))
           (p-1 (second smith-list))
           (D (third smith-list))
           (q (fourth smith-list))
           (C2 (mtrx-prdc p-1 C))
           (line-n (column-number B))
           (A2 
            #-ACLPC
            (make-array (list line-n 1)
                        :element-type 'fixnum
                        :initial-element 0)
            #+ACLPC
            (if (zerop line-n)
                (make-array (list line-n 1)
                            :element-type 'fixnum)
              (make-array (list line-n 1)
                          :element-type 'fixnum
                          :initial-element 0)))
           (end nil))
      (declare 
       (type matrix C B p-1 D q  C2 A2)
       (type list smith-list)
       (type fixnum line-n ))
      (do ((i 0 (1+ i)))
          ((or (>= i (column-number D)) (>= i (line-number D)) end))
        (declare (type fixnum i))
        (let* ((di (aref D i i)))
          (declare (type fixnum di))
          (if (not (eql 0 di))
              ;;(setq found 1)
              (progn
                (setf (aref A2 i 0) (floor (aref C2 i 0) di))))))
      (let* ((A (mtrx-prdc Q A2)))
        (declare (type matrix A))
        (column-j A 0)))))

#|
(setf m (gnrt-list-to-mtrx '((1 0 0 0) (2 3 4 5) (0 1 0 3))))
(setq vctr '(5 9 8 19))
(vctr-coordinates vctr m)
|# 

;; Function that returns the list of coefficients of the combination "cmbn"
;; in the basis (list of generators) "basis" (with comparison function "cmpr").
(DEFUN CMBN-CFFC-LIST (cmbn basis cmpr)
  (declare 
   (type cmbn cmbn)
   (type list basis)
   (type cmprf cmpr))
  (the list
    (let* ((rslt nil)
           (mark1 (cmbn-list cmbn))
           (term (car mark1))
           (cffc (cffc term))
           (gnrt (gnrt term)))
      (declare
       (type list rslt mark1)
       (type term term)
       (type cffc cffc)
       (type gnrt gnrt))
      (do ((mark2 basis (cdr mark2)))
          ((endp mark2))
        (declare 
         (type list mark2))
        (if (endp mark1)
            (setq rslt (nconc rslt (list 0)))
          (if (eq :equal (funcall cmpr gnrt (car mark2)))
              (progn
                (setq rslt (nconc rslt (list cffc)))
                (pop mark1)
                (setq term (car mark1))
                (setq cffc (cffc term) gnrt (gnrt term)))
            (setq rslt (nconc rslt (list 0))))))
      rslt)))

#|
(setf basis (list 'a 'b 'c 'd 'e 'f))
(setq cmbn (cmbn 2 3 'a 5 'c 8 'f))
(cmbn-cffc-list cmbn basis #'s-cmpr)
|# 

;; Function that returns the list of coordinates of the combination "cmbn" in the basis
;; formed by the combinations in the list "gnrt-list". All these combinations are combinations
;; of the generators in "basis" ("cmpr" is the comparison function). 
(DEFUN CMBN-COORDINATES (cmbn gnrt-list basis cmpr)
  (declare 
   (type cmbn cmbn)
   (type list gnrt-list basis)
   (type cmprf cmpr))
  (the list
    (let* ((vctr (cmbn-cffc-list cmbn basis cmpr))
           (B (gnrt-list-to-mtrx
               (mapcar
                   #'(lambda (cmbn2)
                       (declare (type cmbn cmbn2))
                       (cmbn-cffc-list cmbn2 basis cmpr))
                 gnrt-list))))
      (declare
       (type list vctr)
       (type matrix B))
      (vctr-coordinates vctr B))))

#|
(setf basis (list 'a 'b 'c 'd 'e 'f))
(setq cmbn (cmbn 2 3 'a 5 'c 6 'f))
(setq gnrt-list (list (cmbn 2 1 'a 2 'd 2 'f) (cmbn 2 5 'c -6 'd)))
(cmbn-coordinates cmbn gnrt-list basis #'s-cmpr)
|# 
    
         
;; Main function to compute the differential of an element in the spectral sequence          


;; Function that returns the list of integers which correspond to the coordinates
;; on each of the components of E^r_{p-r,q+r-1} of the differential 
;; d: E^r_{p,q} --> E^r_{p-r,q+r-1} applied to the element of E^r_{p,q}
;; which has as coordinates the list int-list. (Only for finitely generated complexes)
(DEFUN EFF-SPSQ-DFFR-OF-ONE-ELEMENT (fltrcm r p q int-list)
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum r p q)
   (list int-list))
  (the list
    (let* ((cmpr (cmpr1 fltrcm))
           (degr (+ p q))
           (sbasis (fltrd-basis fltrcm degr p))
           (tbasis (fltrd-basis fltrcm (1- degr) (- p r)))
           (num-mtrx (spsq-num-mtrx fltrcm r p q))
           (den-mtrx (spsq-den-mtrx fltrcm r p q))
           (z-gnrt-list (fltr-chcm-z-gnrt-list fltrcm r p q))
           (s-basis-dvs (mtrx-quotient (copy-mtrx num-mtrx) den-mtrx))
           (s-gnrt-list (first s-basis-dvs))
           (s-dvs (second s-basis-dvs))
           (t-basis-dvs (spsq-basis-dvs fltrcm r (- p r) (1- (+ q r))))
           (t-gnrt-list (first t-basis-dvs))
           (t-dvs (second t-basis-dvs))
           (selement (make-list (length sbasis) :initial-element 0))
           (cmbn (cmbn degr)))
      (declare 
       (type cmprf cmpr)
       (type fixnum degr)
       (type list sbasis tbasis z-gnrt-list s-basis-dvs s-gnrt-list s-dvs
             t-basis-dvs t-gnrt-list t-dvs selement)
       (type matrix num-mtrx den-mtrx)
       (type cmbn cmbn))
      (progn
        (labels ((2list-add (list1 list2)
                            (declare (type list list1 list2))
                            (the list
                              (mapcar #'+ list1 list2)))
                 (n-list (n list)
                         (declare
                          (type list list)
                          (type fixnum n))
                         (the list
                           (mapcar #'(lambda (i)
                                       (* n i))
                             list))))
          (do ((mark1 s-gnrt-list (cdr mark1))
               (mark2 s-dvs (cdr mark2))
               (mark3 int-list))
              ((endp mark1))
            (declare (list mark1 mark2 mark3))
            (if (not (eql 1 (car mark2)))
                (progn
                  (setq selement (2list-add selement (n-list (car mark3) (car mark1))))
                  (pop mark3)))))
        (let* ((selt-coord (vctr-coordinates selement num-mtrx)))
          (declare (list selt-coord))
          (do ((mark1 z-gnrt-list (cdr mark1))
               (mark2 selt-coord (cdr mark2)))
              ((endp mark1))
            (declare (list mark1 mark2))
            (if (not (eql 0 (car mark2)))
                (setq cmbn (2cmbn-add cmpr cmbn (n-cmbn (car mark2)
                                                        (funcall 
                                                         #'(lambda (int-list)
                                                             (declare (list int-list))
                                                             (the cmbn
                                                               (let* ((cmbn2 (cmbn degr)))
                                                                 (declare 
                                                                  (type cmbn cmbn))
                                                                 (do ((mark3 int-list (cdr mark3))
                                                                      (mark4 sbasis (cdr mark4)))
                                                                     ((or (endp mark3) (endp mark4)))
                                                                   (declare (type list mark3 mark4))
                                                                   (if (not (= 0 (car mark3)))
                                                                       (setq cmbn2 (2cmbn-add cmpr cmbn2 
                                                                                              (cmbn degr (car mark3) (car mark4))))))
                                                                 cmbn2)))
                                                         (car mark1))))))))
        (let* ((dffr-cmbn (dffr fltrcm cmbn))
               (crdnt (cmbn-coordinates dffr-cmbn t-gnrt-list tbasis cmpr))
               (rslt nil))
          (declare
           (type cmbn dffr-cmbn)
           (type list crdnt rslt))
          (do ((mark1 t-dvs (cdr mark1))
               (mark2 crdnt (cdr mark2)))
              ((endp mark1))
            (declare (list mark1 mark2))
            (if (not (eql 1 (car mark1)))
                (if (eql 0 (car mark1))
                    (setq rslt (nconc rslt (list (car mark2))))
                  (setq rslt (nconc rslt (list (mod (car mark2) (car mark1))))))))
          rslt)))))          

#|
(eff-spsq-dffr-of-one-element fcc 1 2 3 '(1))
(eff-spsq-dffr-of-one-element fcc 1 2 5 '(0 1))
(eff-spsq-dffr-of-one-element fcc 1 2 5 '(1 0))
(eff-spsq-dffr-of-one-element fcc 1 3 4 '(1 0))
(eff-spsq-dffr-of-one-element fcc 1 3 4 '(0 1))
|# 



;; Function that returns the matrix (as a list of lists) of the differential 
;; d: E^r_{p,q} --> E^r_{p-r,q+r-1} (Only for finitely generated complexes)
(DEFUN EFF-SPSQ-DFFR-MTRX (fltrcm r p q )
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum r p q))
  (the list
    (let* ((cmpr (cmpr1 fltrcm))
           (degr (+ p q))
           (sbasis (fltrd-basis fltrcm degr p))
           (tbasis (fltrd-basis fltrcm (1- degr) (- p r)))
           (num-mtrx (spsq-num-mtrx fltrcm r p q))
           (den-mtrx (spsq-den-mtrx fltrcm r p q))
           (Z-gnrt-list (Fltr-Chcm-Z-Gnrt-List fltrcm r p q))
           (s-basis-dvs (mtrx-quotient (copy-mtrx num-mtrx) den-mtrx))
           (s-gnrt-list (first s-basis-dvs))
           (s-dvs (second s-basis-dvs))
           (t-basis-dvs (eff-spsq-basis-dvs fltrcm r (- p r) (1- (+ q r))))
           (t-gnrt-list (first t-basis-dvs))
           (t-dvs (second t-basis-dvs))
           (s-length (length
                      (mapcan #'(lambda (i)
                                  (if (not (eql 1 i)) (list i) nil))
                        s-dvs)))
           (base-list nil))
      (declare 
       (type cmprf cmpr)
       (type fixnum degr s-length)
       (type list sbasis tbasis z-gnrt-list s-basis-dvs s-gnrt-list s-dvs
             t-basis-dvs t-gnrt-list t-dvs)
       (type matrix num-mtrx den-mtrx))
      (progn
        (dotimes (i s-length)
          (push (nconc (make-list (1- (- s-length i)) :initial-element 0) (list 1)
                       (make-list i :initial-element 0))
                base-list ))        
        (mapcar #'(lambda (int-list)
                    (let* ((selement (make-list (length sbasis) :initial-element 0))
                           (cmbn (cmbn degr))
                           )
                      (declare
                       (type list selement)
                       (type cmbn cmbn))
                      (progn
                        (labels ((2list-add (list1 list2)
                                            (declare (type list list1 list2))
                                            (the list
                                              (mapcar #'+ list1 list2)))
                                 (n-list (n list)
                                         (declare
                                          (type list list)
                                          (type fixnum n))
                                         (the list
                                           (mapcar #'(lambda (i)
                                                       (* n i))
                                             list))))
                          (do ((mark1 s-gnrt-list (cdr mark1))
                               (mark2 s-dvs (cdr mark2))
                               (mark3 int-list))
                              ((endp mark1))
                            (declare (list mark1 mark2 mark3))
                            (if (not (eql 1 (car mark2)))
                                (progn
                                  (setq selement (2list-add selement (n-list (car mark3) (car mark1))))
                                  (pop mark3)))))
                        (let* ((selt-coord (vctr-coordinates selement num-mtrx)))
                          (declare (list selt-coord))
                          (do ((mark1 z-gnrt-list (cdr mark1))
                               (mark2 selt-coord (cdr mark2)))
                              ((endp mark1))
                            (declare (list mark1 mark2))
                            (if (not (eql 0 (car mark2)))
                                (setq cmbn (2cmbn-add cmpr cmbn (n-cmbn (car mark2)
                                                                        (funcall 
                                                                         #'(lambda (int-list)
                                                                             (declare (list int-list))
                                                                             (the cmbn
                                                                               (let* ((cmbn2 (cmbn degr)))
                                                                                 (declare 
                                                                                  (type cmbn cmbn))
                                                                                 (do ((mark3 int-list (cdr mark3))
                                                                                      (mark4 sbasis (cdr mark4)))
                                                                                     ((or (endp mark3) (endp mark4)))
                                                                                   (declare (type list mark3 mark4))
                                                                                   (if (not (= 0 (car mark3)))
                                                                                       (setq cmbn2 (2cmbn-add cmpr cmbn2 
                                                                                                              (cmbn degr (car mark3) (car mark4))))))
                                                                                 cmbn2)))
                                                                         (car mark1))))))))
                        (let* ((dffr-cmbn (dffr fltrcm cmbn))
                               (crdnt (cmbn-coordinates dffr-cmbn t-gnrt-list tbasis cmpr))
                               (rslt nil))
                          (declare
                           (type cmbn dffr-cmbn)
                           (type list crdnt rslt))
                          (do ((mark1 t-dvs (cdr mark1))
                               (mark2 crdnt (cdr mark2)))
                              ((endp mark1))
                            (declare (list mark1 mark2))
                            (if (not (eql 1 (car mark1)))
                                (if (eql 0 (car mark1))
                                    (setq rslt (nconc rslt (list (car mark2))))
                                  (setq rslt (nconc rslt (list (mod (car mark2) (car mark1))))))))
                          rslt))))
          base-list)))))         



;; CONVERGENCE OF THE SPECTRAL SEQUENCE            

;; Function that determines if at level r and for degree n the E^r_{p,q} of the
;; filtered complex fltrcm with p+q=degr are equal to E^\infty_{p,q} (i.e., it 
;; determines if the convergence of the spectral sequence has been reached at 
;; level r for the degree n). The list hmlg-cmpns contains the "components" of
;; the homology group of the complex. It begins with m zeros if the free part of the
;; group is Z^m and the rest of the list contains the torsion coefficients.
;; (only for finitely generated complexes). 
(DEFUN EFF-COMPARE-SPSQ-WITH-HOMOLOGY (fltrcm r degr hmlg-cmpns)
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum r degr)
   (type list hmlg-cmpns))
  (let* ((min-max (flin-min-max fltrcm degr))
         (min (first min-max))
         (max (second min-max))
         (spsq-cmpns
          (mapcan
              #'(lambda (p)
                  (declare (type fixnum p))
                  (let* ((q (- degr p)))
                    (declare (type fixnum q))
                    (second (eff-spsq-basis-dvs fltrcm r p q))))
            (<a-b> min max)))
         (mark2 hmlg-cmpns)
         (hmlg-prdc 1)
         (spsq-prdc 1))
    (declare 
     (type function spsq-cmpns)
     (type list spsq-cmpns mark2)
     (type fixnum hmlg-prdc spsq-prdc))
    (progn
      (do ((mark1 spsq-cmpns (cdr mark1)))
          ((endp mark1))
        (declare (type list mark1))
        (if (eql 0 (car mark1))
            (if (not (eql 0 (car mark2)))
                (return-from eff-compare-spsq-with-homology nil)
              (pop mark2))
          (setq spsq-prdc (* spsq-prdc (car mark1)))))
      (do ((mark3 mark2 (cdr mark3)))
          ((endp mark3))
        (declare (type list mark3))
        (setq hmlg-prdc (* hmlg-prdc (car mark3))))
      (if (eql hmlg-prdc spsq-prdc)
          't
        nil))))

#|
(homology fcc 5) 
(eff-compare-spsq-with-homology fcc  1 5 '(2))
(homology fcc 6)
(eff-compare-spsq-with-homology fcc  1 6 '(2))
|# 

;; Function that determines the level r at which for degree n the E^r_{p,q} of the
;; filtered complex fltrcm with p+q=degr are equal to E^\infty_{p,q} (i.e., it 
;; determines the level r which the convergence of the spectral sequence has been 
;; reached at for the degree n).            
(DEFUN EFF-SPSQ-CNVG (fltrcm degr)
  (declare 
   (type filtered-chain-complex fltrcm)
   (type fixnum degr))
  (the fixnum
    (let* ((hmlg-dvs-gnrs (homologie (chcm-mat fltrcm degr) (chcm-mat fltrcm (1+ degr))))
           (hmlg-cmpns (reverse
                        (mapcar
                            #'(lambda (l)
                                (first l))
                          hmlg-dvs-gnrs)))
           (cnvg-r 0))
      (declare 
       (list hmlg-dvs-gnrs hmlg-cmpns)
       (type fixnum cnvg-r))
      (do ((r 1 (1+ r)))
          ((> cnvg-r 0))
        (declare (type fixnum r))
        (if (eff-compare-spsq-with-homology fltrcm r degr hmlg-cmpns)
            (setq cnvg-r r)))
      cnvg-r)))

#|
(eff-spsq-cnvg fcc 5)
(eff-spsq-cnvg fcc 6)
|# 

;; EFFECTIVE HOMOLOGY


;; Function that returns the order of filtration t (with regard to the filtration in the
;; left complex, which is translated to the top and right ones) of the homotopy in the right reduction,
;; which will determine the level r (=t+1) up to which the spectral sequences E^r_p,q 
;; (with p+q=n) of the left and right complexes are isomorphic. This function can only be used
;; when the top complex is effective (we need to obtain the elements in the basis of this complex). 
(DEFUN HMTP-EQ-FLTR-ORDER (hmtp-eq degr)
  (declare 
   (type homotopy-equivalence hmtp-eq)
   (type fixnum degr))
  (let* ((lcc (lbcc hmtp-eq))
         (tcc (tcc hmtp-eq))
         (max 0))
    (declare
     (type chain-complex lcc tcc)
     (type fixnum max))
    (when (eq (basis tcc) :locally-effective)
      (error "Hmtp-eq-fltr-order cannot work when the top complex of the homotopy equivalence
             is a LOCALLY-EFFECTIVE chain complex."))
    (mapcar 
        #'(lambda (gnrt)
            (declare (type gnrt gnrt))
            (setq max (max max (- (flin lcc (lf hmtp-eq (rh hmtp-eq degr gnrt)))
                                  (flin lcc (lf hmtp-eq degr gnrt))))))
      (basis tcc degr))
    max))

#|
(setf hmtp-eq (efhm fcc)) ;; Trivial homotopy equivalence
(hmtp-eq-fltr-order hmtp-eq 3)

(setq k3 (k-z2 3))
(Setq hmtp-eq3 (efhm k3))
(setq ek3 (rbcc hmtp-eq3))
(change-chcm-to-flcc ek3 :flin abar-flin :orgn `(filtered-chain-complex ,ek3))
(setf k3-flin
  #'(lambda (degr gnrt)
      (flin ek3 (rf hmtp-eq3 (lg hmtp-eq3 degr gnrt)))))
(change-chcm-to-flcc k3 :flin k3-flin :orgn `(filtered-chain-complex ,k3))
(hmtp-eq-fltr-order hmtp-eq3 3)
(hmtp-eq-fltr-order hmtp-eq3 5)
|#

;; Function that changes the right bottom chain complex in a homotopy equivalence into 
;; a filtered complex, with a filtration defined translating the filtration in the left complex 
;; (which must be a filtered complex). 
(DEFUN TRANSLATE-FILTRATION (hmtp-eq)
  (declare
   (type homotopy-equivalence hmtp-eq))
  (let ((fltrcm (lbcc hmtp-eq))
        (eff-chcm (rbcc hmtp-eq)))
    (declare (type chain-complex fltrcm eff-chcm))
    (Change-Chcm-TO-Flcc eff-chcm 
                         #'(lambda (degr gnrt)
                             (flin fltrcm (lf hmtp-eq (rg hmtp-eq degr gnrt))))
                         (orgn eff-chcm))))



;; Function that returns the representation basis-divisors of E^r_{p,q} for a filtered
;; complex fltrcm with effective homology, computing the E^r_{p,q} of the effective complex
;; in the right side of the homotopy equivalence. 
;; If the effective complex is a filtered complex, we supose that both reductions in the equi-
;; valence are compatible with the filtrations (f and g respect the filtration index and h has
;; filtration order < r). If not, we define the filtration in the effective
;; complex translating the filtration in the original filtered complex FltrCm, suposing that the 
;; homotopy in the right side is compatible with the filtration in the left complex for such r (we 
;; supose that flin(rh(lf(y)))-flin(lf(y))<r for every y in the top complex). 
(DEFMETHOD SPSQ-BASIS-DVS ((fltrcm FILTERED-CHAIN-COMPLEX) r p q)
  (declare 
   (type fixnum r p q))
  (let* ((hmtp-eq (efhm fltrcm))
         (eff-chcm (rbcc hmtp-eq))
         )
    (declare 
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm)
     )
    (progn
      (if (not (typep eff-chcm 'filtered-chain-complex))
          (translate-filtration hmtp-eq))
      (let ((spsq-list (eff-spsq-basis-dvs eff-chcm r p q)))
        (declare (list spsq-list))
        (if (eql fltrcm eff-chcm) spsq-list
          (list
           (mapcar
               #'(lambda (cmbn)
                   (lf hmtp-eq (rg hmtp-eq cmbn)))
             (first spsq-list))
           (second spsq-list)))))))

#|
(spsq-basis-dvs k3 2 2 2)
|# 


(DEFUN SPSQ-GNRTS (fltrcm r p q)
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum r p q))
  (let* ((hmtp-eq (efhm fltrcm))
         (eff-chcm (rbcc hmtp-eq))
         )
    (declare
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm)
     )
    (progn
      (if (not (typep eff-chcm 'filtered-chain-complex))
          (translate-filtration hmtp-eq))
      (eff-spsq-gnrts eff-chcm r p q))))

;; Method that presents on the screen the components of the E^r_{p,q}
;; of the spectral sequence of a filtered complex fltrcm with effective homology.
;; If the effective complex is a filtered complex, we supose that both reductions in the equi-
;; valence are compatible with the filtrations (f and g respect the filtration index and h has
;; filtration order < r). If not, we define the filtration in the effective
;; complex translating the filtration in the original filtered complex FltrCm, suposing that the 
;; homotopy in the right side is compatible with the filtration in the left complex for such r (we 
;; supose that flin(rh(lf(y)))-flin(lf(y))<r for every y in the top complex). 
(DEFMETHOD SPSQ-GROUP ((fltrcm FILTERED-CHAIN-COMPLEX) r p q)
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum r p q))
  (let* ((hmtp-eq (efhm fltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm)
     )
    (progn
      (if (not (typep eff-chcm 'filtered-chain-complex))
          (translate-filtration hmtp-eq))
      (eff-spsq-group eff-chcm r p q))))

;; Method that returns the list of integers which correspond to the coordinates
;; on each of the components of E^r_{p-r,q+r-1} of the differential 
;; d: E^r_{p,q} --> E^r_{p-r,q+r-1} applied to the element of E^r_{p,q}
;; which has as coordinates the list int-list, using the effective homology
;; when the complex is not of finite type.
(DEFMETHOD SPSQ-DFFR-OF-ONE-ELEMENT ((fltrcm FILTERED-CHAIN-COMPLEX) r p q int-list)
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum p q r)
   (type list int-list))
  (let* ((hmtp-eq (efhm fltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (progn
      (if (not (typep eff-chcm 'filtered-chain-complex))
          (translate-filtration hmtp-eq))
      (eff-spsq-dffr-of-one-element eff-chcm r p q int-list))))


;; Function that returns the matrix (as a list of lists) of the differential 
;; d: E^r_{p,q} --> E^r_{p-r,q+r-1}, using the effective homology
;; when the complex is not of finite type.
(DEFMETHOD SPSQ-DFFR-MTRX ((fltrcm FILTERED-CHAIN-COMPLEX) r p q)
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum p q r))
  (let* ((hmtp-eq (efhm fltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm))
    (progn
      (if (not (typep eff-chcm 'filtered-chain-complex))
          (translate-filtration hmtp-eq))
      (eff-spsq-dffr-mtrx eff-chcm r p q))))



;; Function that determines the level r at which for degree n the E^r_{p,q} (with p+q=degr)
;; of the filtered complex fltrcm (with effective homology) are equal to E^\infty_{p,q} (i.e., it 
;; determines the level r which the convergence of the spectral sequence has been 
;; reached at for the degree n).
(DEFUN SPSQ-CNVG (fltrcm degr)
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum degr))
  (let* ((hmtp-eq (efhm fltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm)
     (type fixnum degr))
    (progn
      (if (not (typep eff-chcm 'filtered-chain-complex))
          (translate-filtration hmtp-eq))
      (eff-spsq-cnvg eff-chcm degr))))


;; Homology filtration


(DEFUN HOMOLOGY-GNRT-LIST (chcm degr)
  (declare (type chain-complex chcm) 
           (type fixnum degr))
  (let ((src (cons :unused (basis chcm degr)))
        (hom-dvs-gnrts (homologie (chcm-mat chcm degr) (chcm-mat chcm (1+ degr))))
        )
    (declare (type list src hom-dvs-gnrts ))
    (mapcar #'(lambda (item)
                (declare (type list item))
                (let ((cmbn (cmbn degr))
                      (cmpr (cmpr chcm)))
                  (dolist (item2 (rest item))
                    (declare (type cons item2))
                    (setf cmbn (2cmbn-add cmpr cmbn (cmbn degr (car item2) (nth (cdr item2) src)))))
                  cmbn))
      hom-dvs-gnrts)))



#|
(DEFUN HOMOLOGY-FLTR1 (fltrcm degr fltr-index)
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum degr fltr-index))
  (let* ((r (spsq-cnvg fltrcm degr))
         (cmpr (cmpr1 fltrcm))
         (total-basis (fltrd-basis fltrcm degr fltr-index))
         (dffr-mtrx (FLTR-CHCM-Z-MTRX fltrcm r fltr-index (- degr fltr-index)))
         (fltrt-basis 
          (mapcan #'(lambda (p)
                      (let* ((q (- degr p))
                             (num-mtrx (spsq-num-mtrx fltrcm r p q))
                             (den-mtrx (spsq-den-mtrx fltrcm r p q))
                             (Z-gnrt-list (Fltr-Chcm-Z-Gnrt-List fltrcm r p q))
                             (ss-basis-dvs (mtrx-quotient (copy-mtrx num-mtrx) den-mtrx))
                             (ss-gnrt-list (first ss-basis-dvs))
                             (ss-dvs (second ss-basis-dvs))
                             (rslt nil))
                        (declare 
                         (type cmprf cmpr)
                         (type fixnum r q)
                         (type list total-basis z-gnrt-list ss-basis-dvs ss-gnrt-list ss-dvs rslt)
                         (type matrix num-mtrx den-mtrx))
                        (do ((mark1 ss-gnrt-list (cdr mark1))
                             (mark2 ss-dvs (cdr mark2)))
                            ((endp mark1))
                          (declare (type list mark1 mark2))
                          (if (not (eql 1 (car mark2)))
                              (let* ((selt (first mark1))
                                     ;;(coord (cmbn-cffc-list selt total-basis cmpr))
                                     (cmbn (cmbn degr)))
                                (progn
                                  (let* ((selt-coord (vctr-coordinates selt num-mtrx)))
                                    (declare (list selt-coord))
                                    (do ((mark3 z-gnrt-list (cdr mark3))
                                         (mark4 selt-coord (cdr mark4)))
                                        ((endp mark3))
                                      (declare (list mark3 mark4))
                                      (if (not (eql 0 (car mark4)))
                                          (setq cmbn (2cmbn-add cmpr cmbn (n-cmbn (car mark4)
                                                                                  (funcall 
                                                                                   #'(lambda (int-list)
                                                                                       (declare (list int-list))
                                                                                       (the cmbn
                                                                                         (let* ((cmbn2 (cmbn degr)))
                                                                                           (declare 
                                                                                            (type cmbn 
                                                                                                  cmbn2))
                                                                                           (do ((
                                                                                                 mark5 int-list (cdr 
                                                                                                                 mark5))
                                                                                                (mark6 total-basis (cdr mark6)))
                                                                                               ((or (endp mark5) (endp mark6)))
                                                                                             (declare (type list mark5 mark6))
                                                                                             (if (not (= 0 (car mark5)))
                                                                                                 (setq cmbn2 (2cmbn-add cmpr cmbn2 
                                                                                                                        (cmbn degr (car 
                                                                                                                                    mark5) (car 
                                                                                                                                            mark6))))))
                                                                                           cmbn2)))
                                                                                   (car mark3))))))))
                                  (setf rslt (nconc rslt (list cmbn)))
                                  rslt))))))
            (<a-b> 0 fltr-index))))
    (if (not (endp fltrt-basis))
        (let* ((fltrt-mtrx
                (gnrt-list-to-mtrx
                 (mapcar #'(lambda (cmbn)
                             (cmbn-cffc-list cmbn total-basis cmpr))
                   (fltrd-basis))))
               (num-mtrx (mtrx-conc fltrt-mtrx dffr-mtrx)))
          (mtrx-quotient num-mtrx dffr-mtrx))
      nil)))
|#


(DEFUN EFF-HMLG-FLTR (fltrcm degr p)
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum degr p))
  (let* ((min-max-1 (flin-min-max fltrcm (1- degr)))
         (min-1 (first min-max-1))
         (min-max+1 (flin-min-max fltrcm (1+ degr)))
         (min+1 (first min-max+1))
         (r (1+ (- p min-1)))
         (q (- degr p))
         (Z-gnrts (Fltr-Chcm-Z-Gnrt-List FltrCm r p q))
         (dffr-mtrx (FLTR-CHCM-Z-MTRX fltrcm (+ 2 (- degr min+1)) (1+ degr) 0))
         (num-mtrx (mtrx-conc (gnrt-list-to-mtrx z-gnrts) dffr-mtrx))
         (basis-divs (mtrx-quotient num-mtrx dffr-mtrx))
         (divs (second basis-divs)))
    (format t "Filtration F_~D H_~D" p degr)
    (dolist (item divs)
      (declare (type fixnum item))
      (if (not (eql 1 item))
          (progn
            (format t "~2%Component Z")
            (unless (zerop item) 
              (format t "/~DZ" item))
            )))))


(DEFUN HMLG-FLTR (fltrcm degr p)
  (declare
   (type filtered-chain-complex fltrcm)
   (type fixnum degr))
  (let* ((hmtp-eq (efhm fltrcm))
         (eff-chcm (rbcc hmtp-eq)))
    (declare
     (type homotopy-equivalence hmtp-eq)
     (type chain-complex eff-chcm)
     (type fixnum degr))
    (progn
      (if (not (typep eff-chcm 'filtered-chain-complex))
          (translate-filtration hmtp-eq))
      (eff-hmlg-fltr eff-chcm degr p))))


(DEFCLASS SPECTRAL-SEQUENCE ()
  (;; FiLTeRedchainCoMplex
   (fltrcm :type filtered-chain-complex :initarg :fltrcm :reader fltrcm)
   ;; ReSuLTS 
   (group-rslts :type list :initarg :group-rslts :reader group-rslts)
   (dffr-rslts :type list :initarg :dffr-rslts :reader dffr-rslts)
   ;; IDentification NuMber      
   (idnm :type fixnum :initform (incf *idnm-counter*) :reader idnm)      
   ;; ORiGiN      
   (orgn :type list :initarg :orgn :reader orgn)))


(DEFVAR *ss-list*
    "The variable *SS-LIST* is bound to a list of user created spectral sequences")
(SETF *ss-list* +empty-list+)
(PUSHNEW '*ss-list* *list-list*)


#+clisp(eval-when (:compile-toplevel :load-toplevel :execute)
         (setf (ext:package-lock :clos) nil))
(DEFMETHOD PRINT-OBJECT ((ss SPECTRAL-SEQUENCE) stream)
  (the spectral-sequence
    (progn
      (format stream "[K~D Spectral-Sequence]" (idnm ss))
      ss)))
#+clisp(eval-when (:compile-toplevel :load-toplevel :execute)
         (setf (ext:package-lock :clos) t))


(DEFUN SS (idnm)
  (declare (type fixnum idnm))
  (the (or spectral-sequence null)
    (find idnm *ss-list* :key #'idnm)))


;;; Function to build a spectral sequences from a filtered chain complex and an origin
(DEFUN BUILD-SS (fltrcm orgn)
  (declare
   (type filtered-chain-complex fltrcm)
   (type list orgn))
  (the spectral-sequence
    (progn
      (let ((already (find orgn *ss-list* :test #'equal :key #'orgn)))
        (declare (type (or spectral-sequence null) already))
        (when already
          (return-from build-ss already)))
      (let ((ss (make-instance 'spectral-sequence
                  :fltrcm fltrcm
                  :orgn orgn
                  :group-rslts +empty-list+
                  :dffr-rslts +empty-list+)))
        (declare (type spectral-sequence ss))
        (push ss *ss-list*)
        ss))))


;; Method that returns the representation basis-divisors of E^r_{p,q} of the spectral sequence
(DEFMETHOD SPSQ-BASIS-DVS ((ss SPECTRAL-SEQUENCE) r p q)
  (declare
   (type fixnum r p q))
  (let* ((fltrcm (fltrcm ss))
         (rslts (group-rslts ss))
         (pos (position (list r p q nil) rslts :test #'(lambda (l1 l2)
                                                         (and
                                                          (eq (first l1) (first l2))
                                                          (eq (second l1) (second l2))
                                                          (eq (third l1) (third l2)))))))
    (the list
      (if pos (fourth (nth pos rslts))
        (let ((basis-dvs (spsq-basis-dvs fltrcm r p q)))
          (declare (type list basis-dvs))
          (progn
            (push (list r p q basis-dvs) rslts)
            (setf (slot-value ss 'group-rslts) rslts)
            basis-dvs))))))


;; Method that returns the list of abelian invariants of the group E^r_{p,q} of the spectral sequence
(DEFMETHOD SPECTRAL-SEQUENCE-GROUP ((ss SPECTRAL-SEQUENCE) r p q)
  (declare
   (type spectral-sequence ss)
   (type fixnum r p q))  
  (let* ((basis-dvs (SPSQ-BASIS-DVS SS r p q))
         (divs (second basis-dvs)))
    (declare (list basis-dvs divs))
    (mapcan #'(lambda (i)
                (if (eq 1 i) nil
                  (list i)))
      divs)))


;; Method that prints the spectral sequence group
(DEFMETHOD PRINT-SPSQ-GROUP ((ss SPECTRAL-SEQUENCE) r p q)
  (declare
   (type spectral-sequence ss)
   (type fixnum r p q))  
  (let* ((basis-dvs (SPSQ-BASIS-DVS SS r p q))
         (divs (second basis-dvs)))
    (declare (list basis-dvs divs))
    (format t "Spectral sequence E^~D_{~D,~D}" r p q)
    (dolist (item divs)
      (declare (type fixnum item))
      (if (not (eql 1 item))
          (progn
            (format t "~2%Component Z")
            (unless (zerop item) 
              (format t "/~DZ" item))
            )))
    (terpri)))


;; Method that returns the list of integers which correspond to the coordinates
;; on each of the components of E^r_{p-r,q+r-1} of the differential 
;; d: E^r_{p,q} --> E^r_{p-r,q+r-1} applied to the element of E^r_{p,q}
;; which has as coordinates the list int-list.
(DEFMETHOD SPECTRAL-SEQUENCE-DIFFERENTIAL-OF-ONE-ELEMENT ((ss SPECTRAL-SEQUENCE) r p q int-list)
  (declare
   (type spectral-sequence ss)
   (type fixnum p q r)
   (type list int-list))
  (let* ((fltrcm (fltrcm ss))
         (rslts (dffr-rslts ss))
         (pos (position (list r p q int-list nil) rslts :test #'(lambda (l1 l2)
                                                                  (and
                                                                   (eq (first l1) (first l2))
                                                                   (eq (second l1) (second l2))
                                                                   (eq (third l1) (third l2))
                                                                   (eq (fourth l1) (fourth l2)))))))
    (the list
      (if pos (fifth (nth pos rslts))
        (let ((l (spsq-dffr-of-one-element fltrcm r p q int-list)))
          (declare (type list l))
          (progn
            (push (list r p q int-list l) rslts)
            (setf (slot-value ss 'dffr-rslts) rslts)
            l))))))


;; Method that returns the matrix (as a list of lists) of the differential 
;; d: E^r_{p,q} --> E^r_{p-r,q+r-1}, using the effective homology
;; when the complex is not of finite type.
(DEFMETHOD SPECTRAL-SEQUENCE-DIFFERENTIAL-MATRIX ((ss SPECTRAL-SEQUENCE) r p q)
  (declare
   (type spectral-sequence ss)
   (type fixnum p q r))   
  (let* ((fltrcm (fltrcm ss))
         (rslts (dffr-rslts ss))
         (pos (position (list r p q nil) rslts :test #'(lambda (l1 l2)
                                                         (and
                                                          (eq (first l1) (first l2))
                                                          (eq (second l1) (second l2))
                                                          (eq (third l1) (third l2))
                                                          )))))
    (the list
      (if pos (fourth (nth pos rslts))
        (let ((l (spsq-dffr-mtrx fltrcm r p q )))
          (declare (type list l))
          (progn
            (push (list r p q l) rslts)
            (setf (slot-value ss 'dffr-rslts) rslts)
            l))))))


;; Function that constructs the Serre spectral sequence associated with a fibration f
;; Returns an object of the class Spectral-sequence
(DEFUN SERRE-SPECTRAL-SEQUENCE-FIBRATION (f)
  (declare (type fibration f))
  (let* ((fltrcm (fibration-total f))
         (ecc (rbcc (efhm fltrcm))))
    (declare (type simplicial-set fltrcm)
             (type chain-complex ecc))
    (progn
      (change-chcm-to-flcc fltrcm crpr-flin 'crpr-flin)
      (change-chcm-to-flcc ecc tnpr-flin 'tnpr-flin)
      (the spectral-sequence
        (build-ss fltrcm `(Serre-Spectral-Sequence ,f))))))

;; Function that constructs the Serre spectral sequence associated with a simplicial set
;; which is given as a (twisted) cartesian product. Returns an object of the class Spectral-sequence
(DEFUN SERRE-SPECTRAL-SEQUENCE-PRODUCT (x)
  (declare (type simplicial-set x))
  (let ((ecc (rbcc (efhm x))))
    (declare (type chain-complex ecc))
    (progn
      (change-chcm-to-flcc x crpr-flin 'crpr-flin)
      (change-chcm-to-flcc ecc tnpr-flin 'tnpr-flin)
      (the spectral-sequence
        (build-ss x   `(Serre-Spectral-Sequence ,x))))))


;; Function that constructs the Serre spectral sequence of the first fibration of the
;; Whitehead tower of fibrations of a 1-reduced simplicial set
(DEFUN SERRE-WHITEHEAD-SPECTRAL-SEQUENCE (x)
  (let* (;; we obtain the first non null homology group
         (first-non-null (first-non-null-homology-group x 20))
         (degr (1+ first-non-null))
         (hom (homology-format x degr))
         (ft (construct-space-iterative x (split-components hom) degr)))
    (declare (fixnum first-non-null degr)
             (type string hom)
             (type simplicial-set ft))
    (the spectral-sequence
      (serre-spectral-sequence-product ft))))


;; Function that constructs the Eilenberg-Moore spectral sequence associated with a (1-reduced)
;; simplicial set X
;; Returns an object of the class Spectral-sequence
(DEFUN EILENBERG-MOORE-SPECTRAL-SEQUENCE (x)
  (declare (type simplicial-set x))
  (let* ((ox (loop-space x))
         (ecc (rbcc (efhm OX))))
    (declare 
     (type simplicial-set ox)
     (type chain-complex ecc))
    (progn
      (change-chcm-to-flcc ecc cobar-flin 'cobar-flin)
      (the spectral-sequence
        (build-ss ecc `(Eilenberg-Moore-Spectral-Sequence ,x))))))





#|
(cat-init)
(progn
  (setf s3 (sphere 3))
  (setf k3 (chml-clss s3 3))
  (setf F3 (z-whitehead s3 k3)))


(setf ss1 (serre-spectral-sequence-fibration f3))

(spectral-sequence-group ss1 2 0 2)


(setf r 2)

(dotimes (n 8)
  (dotimes (p (1+ n))
    (let ((q (- n p)))
      ;;(print-spsq-group ss1 r p q)
      (format t "Spectral sequence E^~D_{~D,~D}" r p q)
      (print (spectral-sequence-group ss1 r p q))
      (terpri))))

(spectral-sequence-differential-matrix ss1 3 3 0)
(spectral-sequence-differential-matrix ss1 3 3 2)
(spectral-sequence-differential-matrix ss1 3 3 4)


(cat-init)
(setf x (loop-space (sphere 3)))
(setf ss2 (eilenberg-moore-spectral-sequence x))

;;(dotimes (n 6)
;;  (dotimes (p (1+ n))
;;    (let ((q (+ n p)))
;;      (print-spsq-group ss2 1 (- p) q))))

(spectral-sequence-differential-matrix ss2 1 -2 8)
(spectral-sequence-differential-of-one-element ss2 1 -2 8 '(1 0 0))
(spectral-sequence-differential-of-one-element ss2 1 -2 8 '(0 1 0))
(spectral-sequence-differential-of-one-element ss2 1 -2 8 '(0 0 1))


(spectral-sequence-differential-of-one-element ss2 1 -1 6 '(1))
(spectral-sequence-differential-matrix ss2 1 -1 6)
(spectral-sequence-differential-matrix ss2 1 -2 6)


(cat-init)
(setf s3 (sphere 3))
(setf ss3 (serre-whitehead-spectral-sequence s3))

(spectral-sequence-group ss3 2 0 2)


(dotimes (n 8)
  (dotimes (p (1+ n))
    (let ((q (- n p)))
      ;;(print-spsq-group ss3 r p q)
      (format t "Spectral sequence E^~D_{~D,~D}" r p q)
      (print (spectral-sequence-group ss3 r p q))
      (terpri))))


(spectral-sequence-differential-matrix ss3 3 3 0)
(spectral-sequence-differential-matrix ss3 3 3 2)
(spectral-sequence-differential-matrix ss3 3 3 4)



(cat-init)
(setf s3 (sphere 3) s2 (sphere 2))
(setf ss4 (serre-spectral-sequence-product (crts-prdc s3 s2)))

(dotimes (n 8)
  (dotimes (p (1+ n))
    (let ((q (- n p)))
      ;;(print-spsq-group ss1 r p q)
      (format t "Spectral sequence E^~D_{~D,~D}" r p q)
      (print (spectral-sequence-group ss4 r p q))
      (terpri))))




|#
       

