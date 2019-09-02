;;;============================================================================
;;;
;;; Subject: Direct Sum
;;;         
;;;
;;;----------------------------------------------------------------------------
;;;
;;; Author: Jónathan Heras Vicente  
;;; 13/10/09
;;;
;;;============================================================================

;;; Adapted by Julián Cuevas, Jose Divasón, Miguel Marco Buzunariz and Ana Romero.
;;; August 2019

(in-package #:cat)

(provide "direct-sum")

(DEFSTRUCT (direct-sum-gsm (:print-function direct-sum-gsm-print) (:conc-name nil))
  (d-sum-ind #+allegro :type #+allegro (member :undefined 0 1))
  (d-sum-old #+allegro :type #+allegro gmsm))

(defmacro direct-sum-gsm (d-sum-ind d-sum-old)
  `(make-direct-sum-gsm :d-sum-ind ,d-sum-ind :d-sum-old ,d-sum-old)    
  )

(DEFMACRO WITH-direct-sum-gsm-2 ((d-sum-ind d-sum-old) direct-sum-gsm . body)
  `(let ((,d-sum-ind (d-sum-ind ,direct-sum-gsm))
         (,d-sum-old (d-sum-old ,direct-sum-gsm)))
     (declare (fixnum ,d-sum-ind) (type gmsm ,d-sum-old))
     ,@body))

(DEFMACRO WITH-direct-sum-gsm (&rest rest)
  (ecase (length (first rest))
    (2 `(WITH-direct-sum-gsm-2 ,@rest))
    ))



(defun direct-sum-gsm-print (direct-sum-gsm stream depth)
  (declare
   (type direct-sum-gsm direct-sum-gsm)
   (stream stream)
   (ignore depth))
  (with-direct-sum-gsm (d-sum-ind d-sum-old) direct-sum-gsm
    (format stream "<DIRECT-SUM-GSM ~A ~A>"
      d-sum-ind
      d-sum-old)
    direct-sum-gsm))



(defun direct-sum-cmpr (cmpr1 cmpr2)
  (declare (type cmprf cmpr1 cmpr2))
  (flet ((rslt (direct-sum-gsm1 direct-sum-gsm2)
               (ecase (d-sum-ind direct-sum-gsm1)
                 (0 (if (equal 0 (d-sum-ind direct-sum-gsm2))
                        (funcall cmpr1 (d-sum-old direct-sum-gsm1) (d-sum-old direct-sum-gsm2))
                      :less
                      ))
                 (1 (if (equal 0 (d-sum-ind direct-sum-gsm2))
                        :greater
                      (funcall cmpr2 (d-sum-old direct-sum-gsm1) (d-sum-old direct-sum-gsm2))))
                 )))
    (the cmprf #'rslt)))



(defun direct-sum-basis (basis1 basis2)
  (declare (type basis basis1 basis2))
  (when (or (eq basis1 :locally-effective)
            (eq basis2 :locally-effective)
            )
    (return-from direct-sum-basis :locally-effective))
  (flet ((rslt (dmns)
               (declare (fixnum dmns))
               (append  (mapcar #'(lambda (old)
                                     (declare (type gmsm old))
                                     (direct-sum-gsm 0 old))
                           (funcall basis1 dmns))
                         (mapcar #'(lambda (old)
                                     (declare (type gmsm old))
                                     (direct-sum-gsm 1 old))
                           (funcall basis2 dmns)))))
    (the basis #'rslt)))


(DEFUN TERM-d-sum-0 (term)
  (declare (type term term))
  (the term
    (term (cffc term)
          (direct-sum-gsm 0 (gnrt term)))))

(DEFUN TERM-d-sum-1 (term)
  (declare (type term term))
  (the term
    (term (cffc term)
          (direct-sum-gsm 1 (gnrt term)))))


(DEFUN direct-sum-cmbn-split (cmbn)
  (declare (type cmbn cmbn))
  (the (values cmbn cmbn)
    (with-cmbn (degr list) cmbn
      (let ((list0 +empty-list+)
            (list1 +empty-list+))
        (declare (type list list0 list1))
        (do ((mark list (cdr mark)))
            ((endp mark)
             (setf list0 (nreverse list0)
                   list1 (nreverse list1)))
          (declare (type list mark))
          (if (equal 0 (d-sum-ind (cdar mark)))
              (push (term (caar mark) (d-sum-old (cdar mark))) list0)
            (push (term (caar mark) (d-sum-old (cdar mark))) list1)))
        (values (make-cmbn :degr degr :list list0)
                (make-cmbn :degr degr :list list1))))))



(DEFUN direct-sum-2cmbn-append (cmbn0 cmbn1)
  (declare (type cmbn cmbn0 cmbn1))
  (the cmbn
    (make-cmbn
     :degr (cmbn-degr cmbn0)
     :list (append (mapcar #'term-d-sum-0 (cmbn-list cmbn0))
                   (mapcar #'term-d-sum-1 (cmbn-list cmbn1))))))



(defun direct-sum-dffr (dffr1 dffr2)
  (declare (type morphism dffr1 dffr2))
  (flet ((rslt (cmbn)
               (declare (type cmbn cmbn))
               (multiple-value-bind (cmbn0 cmbn1)
                       (direct-sum-cmbn-split cmbn)
                     (declare (type cmbn cmbn0 cmbn1))
                     (direct-sum-2cmbn-append (? dffr1 cmbn0)
                                              (? dffr2 cmbn1)))
               ))
    (the intr-mrph #'rslt)))




(defun direct-sum (chcm1 chcm2)
  (declare (type chain-complex chcm1 chcm2))
  (the chain-complex
    (build-chcm
     :cmpr (direct-sum-cmpr (cmpr chcm1) (cmpr chcm2))
     :basis (direct-sum-basis (basis chcm1) (basis chcm2))
     :bsgn (direct-sum-gsm 0 (bsgn chcm1))
     :intr-dffr (direct-sum-dffr (dffr chcm1) (dffr chcm2))
     :strt :cmbn
     :orgn `(direct-sum ,chcm1 ,chcm2))))


(defun direct-sum-mrph (sorc trgt mrph1 mrph2)
  (declare (type morphism mrph1 mrph2))
  (the morphism
    (progn
      (unless (= (degr mrph1) (degr mrph2))
        (error "Non-coherent degrees in direct-sum-mrph"))
      (build-mrph
       :sorc sorc
       :trgt trgt
       :degr (degr mrph1)
       :intr (direct-sum-dffr mrph1 mrph2)
       :strt :cmbn
       :orgn `(direct-sum-mrph ,sorc ,trgt ,mrph1 ,mrph2)))))

(DEFUN direct-sum-efhm (chcm1 chcm2)
  (declare (type chain-complex chcm1 chcm2))
  (the homotopy-equivalence
    (let* ((efhm1 (efhm chcm1)) (efhm2 (efhm chcm2))
           (lf1 (lf efhm1)) (lg1 (lg efhm1)) (lh1 (lh efhm1))
           (rf1 (rf efhm1)) (rg1 (rg efhm1)) (rh1 (rh efhm1))
           (lf2 (lf efhm2)) (lg2 (lg efhm2)) (lh2 (lh efhm2))
           (rf2 (rf efhm2)) (rg2 (rg efhm2)) (rh2 (rh efhm2))
           (tcc1 (tcc efhm1))
           (tcc2 (tcc efhm2))
           (rbcc1 (rbcc efhm1))
           (rbcc2 (rbcc efhm2))
           (lbcc (direct-sum chcm1 chcm2))
           (tcc (direct-sum tcc1 tcc2))
           (rbcc (direct-sum rbcc1 rbcc2))
           (LF (direct-sum-mrph tcc lbcc lf1 lf2))
           (LG (direct-sum-mrph lbcc tcc lg1 lg2))
           (LH (direct-sum-mrph tcc tcc lh1 lh2))
           (RF (direct-sum-mrph tcc rbcc rf1 rf2))
           (RG (direct-sum-mrph rbcc tcc rg1 rg2))
           (RH (direct-sum-mrph tcc tcc rh1 rh2)))
      (build-hmeq
       :lrdct (build-rdct :f LF :g LG :h LH
                          :orgn `(direct-sum-efhm ,chcm1 ,chcm2 lrdct))
       :rrdct (build-rdct :f RF :g RG :h RH
                          :orgn `(direct-sum-efhm ,chcm1 ,chcm2 rrdct))
       :orgn `(direct-sum-efhm ,chcm1 ,chcm2)))))




(DEFMETHOD SEARCH-EFHM (direct-sum (orgn (eql 'direct-sum)))
  (declare (type chain-complex direct-sum))
  (the homotopy-equivalence
    (direct-sum-efhm (second (orgn direct-sum))
                     (third (orgn direct-sum)))))


