
(IN-PACKAGE #:cat)

(provide "finite-spaces-changes")


(DEFMACRO NEWSMITH-IDNT-MTRX (n)
  `(identite ,n))
  
  
(DEFMACRO NEWSMITH-COPY-MTRX (mtrx)
  `(copier-matrice ,mtrx))
  
  
(DEFUN NEWSMITH-LEFT-SUBMATRIX (mtrx k)
  ;;; keeps the columns of index <= k
  (declare (type matrice mtrx) (type fixnum k))
  (the matrice
    (let ((line-number (line-number mtrx)))
      (declare (type fixnum line-number))
      (let ((rslt (creer-matrice line-number k)))
        (declare (type matrice rslt))
        (do ((basecol (uplig mtrx))
             (il 1 (1+ il)))
            ((> il line-number))
          (declare (type array basecol) (type fixnum il))
          (let ((spl0 (baselig mtrx il)))
            (declare (type t-mat spl0))
            (do ((spl (left (chercher-hor spl0 k)) (left spl))
                 (tpl (baselig rslt il) (left tpl)))
                ((eq spl spl0))
              (declare (type t-mat spl tpl))
              (inserer-terme tpl (aref basecol (icol spl)) (val spl)))))
        rslt))))


(DEFMACRO NEWSMITH-CHCM-MTRX (chcm degr)
  `(chcm-mat ,chcm ,degr))


(DEFUN NEWSMITH-EXTRACT-LINE (mtrx ilig)
  (declare (type matrice mtrx) (type fixnum ilig))
  (the list
    (let ((pl0 (baselig mtrx ilig)))
      (declare (type t-mat pl0))
      (do ((pl (left pl0) (left pl))
           (rslt +empty-list+ (cons (list (icol pl) (val pl)) rslt)))
          ((eq pl pl0) rslt)
        (declare (type t-mat pl) (type list rslt))))))


(DEFUN NEWSMITH-EXTRACT-COLUMN (mtrx icol)
  (declare (type matrice mtrx) (type fixnum icol))
  (the list
    (let ((pc0 (basecol mtrx icol)))
      (declare (type t-mat pc0))
      (do ((pc (up pc0) (up pc))
           (rslt +empty-list+ (cons (list (ilig pc) (val pc)) rslt)))
          ((eq pc pc0) rslt)
        (declare (type t-mat pc) (type list rslt))))))


(DEFUN NEWSMITH-EXTRACT-TERM (matrix il ic)
  (declare (type matrice matrix) (type fixnum ic il))
  (the fixnum
    (do ((p (left (baselig matrix il)) (left p)))
        (nil)
      (declare (type t-mat p))
      (let ((p-ic (icol p)))
        (declare (type fixnum p-ic))
        (when (<= p-ic ic)
          (return-from newsmith-extract-term
            (if (eql p-ic ic)
                (val p)
              0)))))))
              
              
(DEFUN NEWSMITH-NEW-LINE (mtrx ilig list)
  (declare (type matrice mtrx) (type fixnum ilig) (type list list))
  (the matrice
    (let ((peigne (peigne-ver mtrx (basecol mtrx 0) ilig))
          (pl0 (baselig mtrx ilig)))
      (declare (type list peigne) (type t-mat pl0))
      (map nil #'(lambda (item)
                   (declare (type t-mat item))
                   (when (eq (left pl0) (up item))
                     (supprimer-terme pl0 item)))
        peigne)
      (do ((markp (nreverse peigne) (cdr markp))
           (icol 1 (1+ icol))
           (markl list))
          ((endp markl))
        (declare (type list markp markl) (type fixnum icol))
        (when (eql (caar markl) icol)
          (inserer-terme pl0 (car markp) (cadar markl))
          (setf markl (cdr markl))))      
      mtrx)))


(DEFUN NEWSMITH-NEW-COLUMN (mtrx icol list)
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



(DEFUN NEWSMITH-LINE-OP (mtrx lambda iline1 iline2)
  ;; line2 := line2 + lambda * line1
  (declare
   (type matrice mtrx)
   (type fixnum lambda iline1 iline2))
  (the matrice
    (let ((mark1 (newsmith-extract-line mtrx iline1))
          (mark2 (newsmith-extract-line mtrx iline2))
          (new-line2 +empty-list+))
      (declare (type list mark1 mark2 new-line2))
      (loop
        (when (endp mark1)
          (setf new-line2 (nreconc new-line2 mark2))
          (return))
        (when (endp mark2)
          (setf new-line2
            (nreconc new-line2
                     (mapcar #'(lambda (item)
                                 (declare (type list item))
                                 (list (first item)
                                       (safe-* lambda (second item))))
                       mark1)))
          (return))
        (let ((icol1 (caar mark1))
              (icol2 (caar mark2)))
          (declare (type fixnum icol1 icol2))
          (cond ((< icol1 icol2)
                 (push (list icol1 (safe-* lambda (second (pop mark1))))
                       new-line2))
                ((> icol1 icol2)
                 (push (list icol2 (second (pop mark2))) new-line2))
                (t (let ((new-val (safe-+ (safe-* lambda (second (pop mark1)))
                                     (second (pop mark2)))))
                     (declare (type fixnum new-val))
                     (unless (zerop new-val)
                       (push (list icol1 new-val) new-line2)))))))
      (newsmith-new-line mtrx iline2 new-line2)
      mtrx)))


(DEFUN NEWSMITH-COLUMN-OP (mtrx lambda icol1 icol2)
  ;; col2 := col2 + lambda * col1
  (declare
   (type matrice mtrx)
   (type fixnum lambda icol1 icol2))
  (the matrice
    (let ((mark1 (newsmith-extract-column mtrx icol1))
          (mark2 (newsmith-extract-column mtrx icol2))
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
      (newsmith-new-column mtrx icol2 new-column2)
      mtrx)))
      
      
(DEFUN NEWSMITH-EQUAL-MATRIX (mtrx1 mtrx2)
  (declare (type matrice mtrx1 mtrx2))
  (the boolean
    (let ((line-number (nlig mtrx1)))
      (declare (type fixnum line-number))
      (unless (eql line-number (nlig mtrx2))
        (return-from newsmith-equal-matrix +false+))
      (unless (eql (ncol mtrx1) (ncol mtrx2))
        (return-from newsmith-equal-matrix +false+))
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
                (return-from newsmith-equal-matrix +false+))
              (return))
            (when (eq p2 p20)
              (unless (eq p1 p10)
                (return-from newsmith-equal-matrix +false+)))
            (unless (eql (icol p1) (icol p2))
              (return-from newsmith-equal-matrix +false+))
            (unless (eql (val p1) (val p2))
              (return-from newsmith-equal-matrix +false+)))))
      (return-from newsmith-equal-matrix +true+))))


(DEFUN NEWSMITH-LINE-SWAP (mtrx iline1 iline2)
  ;; swaps line1 and line2
  (declare
    (type matrice mtrx) (type fixnum iline1 iline2))
  (the matrice
    (let ((new-line1 (newsmith-extract-line mtrx iline2))
          (new-line2 (newsmith-extract-line mtrx iline1)))
      (declare (type list new-line1 new-line2))
      (newsmith-new-line mtrx iline1 new-line1)
      (newsmith-new-line mtrx iline2 new-line2)
      mtrx)))


(DEFUN NEWSMITH-COLUMN-SWAP (mtrx icol1 icol2)
  ;; swaps column1 and column2
  (declare
    (type matrice mtrx) (type fixnum icol1 icol2))
  (the matrice
    (let ((new-column1 (newsmith-extract-column mtrx icol2))
          (new-column2 (newsmith-extract-column mtrx icol1)))
      (declare (type list new-column1 new-column2))
      (newsmith-new-column mtrx icol1 new-column1)
      (newsmith-new-column mtrx icol2 new-column2)
      mtrx)))
      
      
(DEFUN NEWSMITH-LINE-MINUS (mtrx iline)
  ;; line := - line
  (declare (type matrice mtrx) (type fixnum iline))
  (the matrice
    (newsmith-new-line mtrx iline
              (mapcar #'(lambda (item)
                          (declare (type list item))
                          (the list
                            (list (first item)
                                  (safe-- (second item)))))
                (newsmith-extract-line mtrx iline)))))


(DEFUN NEWSMITH-COLUMN-MINUS (mtrx icol)
  ;; column := - column
  (declare (type matrice mtrx) (type fixnum icol))
  (the matrice
    (newsmith-new-column mtrx icol
                (mapcar #'(lambda (item)
                            (declare (type list item))
                            (the list
                              (list (first item)
                                    (safe-- (second item)))))
                  (newsmith-extract-column mtrx icol)))))
                  

(DEFUN NEWSMITH-MINIMAL-TERM (matrix begin)
  (declare (type matrice matrix) (type fixnum begin))
  (the (values fixnum fixnum fixnum)
    (do ((il (nlig matrix) (1- il))
         (min 0)
         (min-il -1)
         (min-ic -1))
        ((< il begin)
         (return-from newsmith-minimal-term
           (values min min-il min-ic)))
      (declare (type fixnum il min min-il min-ic))
      (do ((p (left (baselig matrix il)) (left p)))
          ((< (icol p) begin))
        (declare (type t-mat p))
        (let ((term (abs (val p))))
          (declare (type fixnum term))
          (when (eql term 1)
            (return-from newsmith-minimal-term (values 1 il (icol p))))
          (when (plusp term)
            (when (or (< term min) (zerop min))
              (setf min term
                min-il il
                min-ic (icol p)))))))))
                

(DEFUN NEWSMITH-MINIMAL-REST-1 (matrix begin)
  ;; Let c (= corner) the term M_{b,b} (b = begin).
  ;; This function looks for the minimal rest of the division
  ;; of M_{il,ic} by c
  ;;              for il = begin and ic > begin
  ;;               or ic = begin and il > begin
  ;; asserts c defined and non-null
  (declare (type matrice matrix) (type fixnum begin))
  (the (values fixnum fixnum fixnum)
    (let ((min 0) (min-il -1) (min-ic -1)
          (corner (do ((p (left (baselig matrix begin)) (left p)))
                      ((eql (icol p) begin) (val p))
                    (declare (type t-mat p))
                    (when (zerop (icol p))
                      (error "Illegal matrix in MINIMAL-REST-1")))))
      (declare (type fixnum min min-il min-ic corner))
      (assert (not (zerop corner)))
      (do ((p (left (baselig matrix begin)) (left p)))
          ((eql (icol p) begin))
        (declare (type t-mat p))
        (let ((term (abs (second (multiple-value-list
                                  (round (val p) corner))))))
          (declare (type fixnum term))
          (when (= term 1)
            (return-from newsmith-minimal-rest-1 (values 1 begin (icol p))))
          (when (plusp term)
            (when (or (< term min) (zerop min))
              (setf min term
                min-il begin
                min-ic (icol p))))))
      (do ((p (up (basecol matrix begin)) (up p)))
          ((eql (ilig p) begin))
        (declare (type t-mat p))
        (let ((term (abs (second (multiple-value-list
                                  (round (val p) corner))))))
          (declare (type fixnum term))
          (when (= term 1)
            (return-from newsmith-minimal-rest-1 (values 1 (ilig p) begin)))
          (when (plusp term)
            (when (or (< term min) (zerop min))
              (setf min term
                min-il (ilig p)
                min-ic begin)))))
      (values min min-il min-ic))))

                
(DEFUN NEWSMITH-MINIMAL-REST-2 (matrix begin)
  ;; Let c (= corner) the term M_{b,b} (b = begin).
  ;; This function looks for the minimal rest of the division
  ;; of M_{il,ic} by c for il > begin and ic > begin.
  (declare (type matrice matrix) (type fixnum begin))
  (the (values fixnum fixnum fixnum)
    (let ((min 0) (min-il -1) (min-ic -1)
          (corner (do ((p (left (baselig matrix begin)) (left p)))
                      ((eql (icol p) begin) (val p))
                    (declare (type t-mat p))
                    (when (zerop (icol p))
                      (error "Illegal matrix in MINIMAL-REST-1")))))
      (declare (type fixnum min min-il min-ic corner))
      (assert (not (zerop corner)))
      (do ((il (nlig matrix) (1- il)))
          ((eql il begin))
        (declare (type fixnum il))
        (do ((p (left (baselig matrix il)) (left p)))
            ((<= (icol p) begin))
          (declare (type t-mat p))
          (let ((term (abs (second (multiple-value-list
                                    (round (val p) corner))))))
            (declare (type fixnum term))
            (when (= 1 term)
              (return-from newsmith-minimal-rest-2 (values 1 il (icol p))))
            (when (plusp term)
              (when (or (< term min) (zerop min))
                (setf min term
                  min-il il
                  min-ic (icol p)))))))
      (values min min-il min-ic))))



(DEFUN NEWSMITH-MINIMAL-TERM-TOP-LEFT (mtrx-list begin il ic)
  (declare (list mtrx-list)  (type fixnum begin il ic))
  (the list
    (progn
      (assert (<= begin il))
      (assert (<= begin ic))
      (when (< begin il)
        (newsmith-line-swap-5 mtrx-list begin il))
      (when (< begin ic)
        (newsmith-column-swap-5 mtrx-list begin ic))
      (let ((corner (newsmith-extract-term (third mtrx-list) begin begin)))
        (declare (type fixnum corner))
        (when (minusp corner)
          (newsmith-line-minus-5 mtrx-list begin)))
      mtrx-list)))


(DEFUN NEWSMITH-PIVOTT (mtrx-list begin &aux (mtrx (third mtrx-list)))
  (declare (type list mtrx-list) (type fixnum begin)
           (type matrice mtrx))
  (the list
    (let ((corner (newsmith-extract-term mtrx begin begin))
          (pc0 (basecol mtrx begin))
          (pl0 (baselig mtrx begin)))
      (declare (type fixnum corner) (type t-mat pc0 pl0))
      (do ((p (up pc0) (up pc0)))
          (nil)
        (declare (type t-mat p))
        (let ((il (ilig p)))
          (declare (type fixnum il))
          (when (eql il begin)
            (return))
          (let ((lambda (safe-- (/ (val p) corner))))
            (declare (type fixnum lambda))
            (newsmith-line-op-5 mtrx-list lambda begin il))))
      (do ((p (left pl0) (left pl0)))
          (nil)
        (declare (type t-mat p))
        (let ((ic (icol p)))
          (declare (type fixnum ic))
          (when (eql ic begin)
            (return))
          (let ((lambda (safe-- (/ (val p) corner))))
            (declare (type fixnum lambda))
            (newsmith-column-op-5 mtrx-list lambda begin ic))))
      mtrx-list)))


(DEFUN NEWSMITH-LIST-SMITH (mtrx-list)
  (declare (list mtrx-list))
  (the list
    (progn
      (let ((matrix (third mtrx-list))
            (begin 1))
        (declare (type matrice matrix) (type fixnum begin))
        (loop
          (multiple-value-bind (term il ic) (newsmith-minimal-term matrix begin)
            (declare (type fixnum term il ic))
            ;; (format t "~%*BEGIN* = ~D ; MIN = ~D." begin term)
            (when (zerop term)
              (return-from newsmith-list-smith mtrx-list))
            (newsmith-minimal-term-top-left mtrx-list begin il ic)
            ; (print (list 1 begin il ic t2))
            ; (break)
            )
          (loop
            (multiple-value-bind (term il ic) (newsmith-minimal-rest-1 matrix begin)
              (declare (type fixnum term il ic))
              (cond ((zerop term)
                     (newsmith-pivott mtrx-list begin)
                     ; (print (list 2 t2))
                     ; (break)
                     (multiple-value-bind (term il ic) (newsmith-minimal-rest-2 matrix begin)
                       (declare (type fixnum term il) (ignore ic))
                       (cond ((zerop term)
                              (return))
                             (t
                              (newsmith-line-op-5 mtrx-list 1 il begin)
                              ; (print (list 3 t2))
                              ; (break)
                              ))))
                    ((= il begin)
                     (newsmith-column-op-5 mtrx-list
                                  (safe-- (round (newsmith-extract-term matrix begin ic)
                                            (newsmith-extract-term matrix begin begin)))
                                  begin ic)
                     ; (print (list 4 t2))
                     ; (break)
                     (newsmith-column-swap-5 mtrx-list begin ic)
                     ; (print (list 5 t2))
                     ; (break)
                     (when (minusp (newsmith-extract-term matrix begin begin))
                       (newsmith-column-minus-5 mtrx-list begin)
                       ; (print (list 6 t2))
                       ; (break)
                       ))
                    (t
                     (newsmith-line-op-5 mtrx-list
                                (safe-- (round (newsmith-extract-term matrix il begin)
                                          (newsmith-extract-term matrix begin begin)))
                                begin il)
                     ; (print (list 7 t2))
                     ; (break)
                     (newsmith-line-swap-5 mtrx-list begin il)
                     ; (print (list 8 t2))
                     ; (break)
                     (when (minusp (newsmith-extract-term matrix begin begin))
                       (newsmith-column-minus-5 mtrx-list begin)
                       ; (print (list 9 t2))
                       ; (break)
                       )))))
          ;; (Format t "~%  Finally the diagonal term is ~D." (aref matrix begin begin))
          (incf begin)))
      mtrx-list)))


(DEFUN NEWSMITH-SMITH (matrix)
  (declare (type matrice matrix))
  (the list
    (let ((line-n (nlig matrix))
          (column-n (ncol matrix)))
      (declare (type fixnum line-n column-n))
      (newsmith-list-smith
       (list (newsmith-idnt-mtrx line-n) (newsmith-idnt-mtrx line-n)
             matrix
             (newsmith-idnt-mtrx column-n) (newsmith-idnt-mtrx column-n))))))


(DEFUN MAJTERME (pl pc val)
  #| Same 'maj-terme' changing "=" by "eq" (NewSmith) |#
  (if (eq val 0)
      (if (eq (left pl) (up pc)) (supprimer-terme pl pc))
    (if (eq (left pl) (up pc))
          (setf (val (left pl)) val)
      (inserer-terme pl pc val))))


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFUN EQUALMATRIX (mtrx1 mtrx2)
  #| Same 'equal-matrix' changing "line-number" by "nlig" and "column-number" by "ncol" (NewSmith) |#
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


#|
  (setf m1 (creer-matrice 2 3))
  (setf m2 (creer-matrice 3 3))
  (equalmatrix m1 m2)
|#


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFUN NEWSMITH-MTRX-PRDC (mtrx1 mtrx2)
   #| Same 'MTRX-PRDC' changing "line-number" by "nlig" and "column-number" by "ncol" (NewSmith) |#
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


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFUN UNIONMERGE (list1 list2)
  #| Same UNION-MERGE but not allowing repetition (list1 and list2 in ascending order) |#
  (declare (type list list1 list2))
  (the list
       (let ((rslt nil))
         (declare (type list rslt))
         (unless list1 (return-from unionmerge list2))
         (unless list2 (return-from unionmerge list1))
         (loop
          (let ((n1 (car list1))
                (n2 (car list2)))
            (declare (type fixnum n1 n2))
            (cond ((< n1 n2)
                   (push n1 rslt)
                   (or (setf list1 (cdr list1))
                       (return-from unionmerge (nreconc rslt list2))))
                  ((= n1 n2)
                   (push n1 rslt)
                   (or (setf list1 (cdr list1))
                       (progn (setf list2 (cdr list2)) ; (1 2) U (2 3) = (1 2 3)
                         (return-from unionmerge (nreconc rslt list2))))
                   (or (setf list2 (cdr list2))
                       (return-from unionmerge (nreconc rslt list1))))
                  (t
                   (push n2 rslt)
                   (or (setf list2 (cdr list2))
                       (return-from unionmerge (nreconc rslt list1))))))))))


#|
  (union-merge '(1 2 3 4) '(4 5 6))
  (unionmerge '(1 2 3 4) '(4 5 6))
  (union-merge '(4 7 9) '(1 2 4))
  (unionmerge '(4 7 9) '(1 2 4))
  (union-merge '(1 2 4) '(1 2 4 5))
  (unionmerge '(1 2 4) '(1 2 4 5))
|#


(DEFUN SAFE-* (arg1 arg2)
   #| This is in 'New-Smith' |#
  (declare (type fixnum arg1 arg2))
  (let ((rslt (* arg1 arg2)))
    (declare (type integer rslt))
    (if (typep rslt 'fixnum)
        rslt
      (error "Fixnum-overflow in (safe-* ~D ~D)." arg1 arg2))))


#|
  (safe-* 23170 23170)
  (safe-* 23170 23171)
|#


(DEFUN SAFE-+ (arg1 arg2)
  #| This is in 'New-Smith' |#
  (declare (type fixnum arg1 arg2))
  (let ((rslt (+ arg1 arg2)))
    (declare (type integer rslt))
    (if (typep rslt 'fixnum)
        rslt
      (error "Fixnum-overflow in (safe-+ ~D ~D)." arg1 arg2))))

#|
  (safe-+ 268435456 268435455)
  (safe-+ 268435456 268435456)
|#


(DEFUN SAFE-- (arg)
  #| This is in 'New-Smith' |#
  (declare (type fixnum arg))
  (let ((rslt (- arg)))
    (declare (type integer rslt))
    (if (typep rslt 'fixnum)
        rslt
      (error "Fixnum-overflow in (SAFE-- ~D)." arg))))


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFUN NEWSMITH-EXTRACT-COLUMN (mtrx icol)
  #| This is in 'New-Smith' |#
  (declare (type matrice mtrx) (type fixnum icol))
  (the list
    (let ((pc0 (basecol mtrx icol)))
      (declare (type t-mat pc0))
      (do ((pc (up pc0) (up pc))
           (rslt +empty-list+ (cons (list (ilig pc) (val pc)) rslt)))
          ((eq pc pc0) rslt)
        (declare (type t-mat pc) (type list rslt))))))


(DEFUN NEWSMITH-EXTRACT-LINE (mtrx ilig)
  #| This is in 'New-Smith' |#
  (declare (type matrice mtrx) (type fixnum ilig))
  (the list
    (let ((pl0 (baselig mtrx ilig)))
      (declare (type t-mat pl0))
      (do ((pl (left pl0) (left pl))
           (rslt +empty-list+ (cons (list (icol pl) (val pl)) rslt)))
          ((eq pl pl0) rslt)
        (declare (type t-mat pl) (type list rslt))))))


#|------------------------------------------------------------------------|#
#|------------------------------------------------------------------------|#


(DEFUN NEWSMITH-LINE-OP-5 (mtrx-list lambda line1 line2)
  #| This is in 'New-Smith' |#
  (let ((slambda (gensym)))
    (let ((slambda lambda))
       (newsmith-column-op (first mtrx-list) (- slambda) line2 line1)
       (newsmith-line-op (second mtrx-list) slambda line1 line2)
       (newsmith-line-op (third mtrx-list) slambda line1 line2))))


(DEFUN NEWSMITH-COLUMN-OP-5 (mtrx-list lambda column1 column2)
  #| This is in 'New-Smith' |#
  (let ((slambda (gensym)))
    (let ((slambda lambda))
       (newsmith-column-op (third mtrx-list) slambda column1 column2)
       (newsmith-column-op (fourth mtrx-list) slambda column1 column2)
       (newsmith-line-op (fifth mtrx-list) (- slambda) column2 column1))))


(DEFUN NEWSMITH-LINE-SWAP-5 (mtrx-list line1 line2)
  #| This is in 'New-Smith' |#
  (progn
     (newsmith-column-swap (first mtrx-list) line1 line2)
     (newsmith-line-swap (second mtrx-list) line1 line2)
     (newsmith-line-swap (third mtrx-list) line1 line2)))


(DEFUN NEWSMITH-COLUMN-SWAP-5 (mtrx-list column1 column2)
  #| This is in 'New-Smith' |#
  (progn
     (newsmith-column-swap (third mtrx-list) column1 column2)
     (newsmith-column-swap (fourth mtrx-list) column1 column2)
     (newsmith-line-swap (fifth mtrx-list) column1 column2)))


(DEFUN NEWSMITH-LINE-MINUS-5 (mtrx-list line)
  #| This is in 'New-Smith' |#
  (progn
     (newsmith-column-minus (first mtrx-list) line)
     (newsmith-line-minus (second mtrx-list) line)
     (newsmith-line-minus (third mtrx-list) line)))


(DEFUN NEWSMITH-COLUMN-MINUS-5 (mtrx-list column)
  #| This is in 'New-Smith' |#
  (progn
     (newsmith-column-minus (third mtrx-list) column)
     (newsmith-column-minus (fourth mtrx-list) column)
     (newsmith-line-minus (fifth mtrx-list) column)))

