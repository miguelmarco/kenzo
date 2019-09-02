;;;============================================================================
;;;
;;; Subject: Cones 2nd definition
;;;         
;;;
;;;----------------------------------------------------------------------------
;;;
;;; Author: Jónathan Heras Vicente
;;;             
;;; 13/10/09
;;;
;;;============================================================================

;;; Adapted by Julián Cuevas, Jose Divasón, Miguel Marco Buzunariz and Ana Romero.
;;; August 2019

(in-package #:cat)

(provide "cones2")


(DEFUN CONE-BASIS2 (basis0 basis1)
  (declare (type basis basis0 basis1))
  (the basis
    (progn
      (when (or (eq basis0 :locally-effective)
                (eq basis1 :locally-effective))
        (return-from cone-basis2 :locally-effective))
      (flet ((rslt (degr)
               (declare (type fixnum degr))
               (the list
                 (append
                  (mapcar #'(lambda (item) (con0 item))
                    (funcall basis0 (1+ degr)))
                  (mapcar #'(lambda (item) (con1 item))
                    (funcall basis1 degr))))))
        (declare (ftype (function (fixnum) list) rslt))
        #'rslt))))


(DEFUN CONE-CMBN-SPLIT2 (cmbn)
  (declare (type cmbn cmbn))
  (the (values cmbn cmbn)
    (with-cmbn (degr list) cmbn
      (let ((list0 +empty-list+)
            (list1 +empty-list+))
        (declare (type list list0 list1))
        (do ((mark list (cdr mark)))
            ((endp mark)
             (setf list0 (nreverse list0)
                   list1 +empty-list+))
          (declare (type list mark))
          (when (eql 1 (conx (cdar mark)))
            (setf list0 (nreverse list0)
              list1 (mapcar #'term-uncon mark))
            (return))
          (push (term-uncon (car mark)) list0))
        (values (make-cmbn :degr (1+ degr) :list list0)
                (make-cmbn :degr degr :list list1))))))



(DEFUN CMBN-CON02 (cmbn)
  (declare (type cmbn cmbn))
  (the cmbn
    (make-cmbn
     :degr (1- (cmbn-degr cmbn))
     :list (mapcar #'term-con0 (cmbn-list cmbn)))))


(DEFUN CONE-2CMBN-APPEND2 (cmbn0 cmbn1)
  (declare (type cmbn cmbn0 cmbn1))
  (the cmbn
    (make-cmbn
     :degr (cmbn-degr cmbn1)
     :list (append (mapcar #'term-con0 (cmbn-list cmbn0))
                   (mapcar #'term-con1 (cmbn-list cmbn1))))))
 
(DEFUN CONE-2MRPH-DIAG-IMPL2 (mrph0 mrph1)
  (declare (type morphism mrph0 mrph1))
  (the intr-mrph
    (flet ((rslt (cmbn)
                 (declare (type cmbn cmbn))
                 (the cmbn
                   (multiple-value-bind (cmbn0 cmbn1)
                       (cone-cmbn-split2 cmbn)
                     (declare (type cmbn cmbn0 cmbn1))
                     (cone-2cmbn-append2 (? mrph0 cmbn0)
                                         (? mrph1 cmbn1))))))
      #'rslt)))



(DEFUN cone-3mrph-triangle-impl2 (cmpr0 mrph0 mrph1 phi)
  (declare (ignore cmpr0) (type morphism mrph0 mrph1 phi))
  (the intr-mrph
    (flet ((rslt (degr gnrt)
             (declare (type fixnum degr) (type gnrt gnrt))
                 (the cmbn
                   (ecase (conx gnrt)
                   (0 (cmbn-con02 (? mrph0 (1+ degr) (icon gnrt))))
                   (1
                    (let ((gnrt (icon gnrt)))
                      (declare (type gnrt gnrt))
                      (cone-2cmbn-append2
                       (? phi degr gnrt)
                       (? mrph1 degr gnrt)))))
                 )))
      #'rslt)))


(DEFUN CONE2 (mrph)
  (declare (type morphism mrph))
  (the chain-complex
    (let ((chcm0 (trgt mrph))
          (chcm1 (sorc mrph)))
      (declare (type chain-complex chcm0 chcm1))
      (build-chcm
       :cmpr (cone-cmpr (cmpr chcm0) (cmpr chcm1))
       :basis (cone-basis2 (basis chcm0) (basis chcm1))
       :bsgn (con0 (bsgn chcm0))
       :intr-dffr (cone-3mrph-triangle-impl2 (cmpr chcm0)
                                            (dffr chcm0)
                                            (n-mrph -1 (dffr chcm1))
                                            mrph)
       :strt :gnrt
       :orgn `(cone2 ,mrph)))))




(DEFUN CONE-2MRPH-DIAG2 (sorc-cone trgt-cone mrph0 mrph1)
  (declare
   (type chain-complex sorc-cone trgt-cone)
   (type morphism mrph0 mrph1))
  (the morphism
    (progn
      (unless (= (degr mrph0) (degr mrph1))
        (error "Non-coherent degrees in CONE-2MRPH-DIAG2."))
      (build-mrph
       :sorc sorc-cone
       :trgt trgt-cone
       :degr (degr mrph0)
       :intr (cone-2mrph-diag-impl2 mrph0 mrph1)
       :strt :cmbn
       :orgn `(cone-2mrph-diag2 ,sorc-cone ,trgt-cone ,mrph0 ,mrph1)))))



(DEFUN CONE-3MRPH-TRIANGLE2 (sorc-cone trgt-cone mrph0 mrph1 phi)
  (declare
   (type chain-complex sorc-cone trgt-cone)
   (type morphism mrph0 mrph1 phi))
    (the morphism
    (progn
      (unless (= (degr mrph0) (degr mrph1))
        (error "Non-coherent degrees in CONE-3MRPH-TRIANGLE."))
      (unless (= (1+ (degr mrph0)) (degr phi))
       (error "Non-coherent-degrees in CONE-3MRPH-TRIANGLE."))
      (build-mrph
       :sorc sorc-cone
       :trgt trgt-cone
       :degr (degr mrph0)
       :intr (cone-3mrph-triangle-impl2
              (cmpr (trgt phi))
              mrph0 mrph1 phi)
       :strt :gnrt
       :orgn `(cone-3mrph-triangle2 ,sorc-cone ,trgt-cone
                                    ,mrph0 ,mrph1 ,phi)))))


(DEFUN CONE2-EFHM (cone)
  (declare (type chain-complex cone))
  (the homotopy-equivalence
    (let* ((phi (second (orgn cone)))
           (chcm0 (trgt phi)) (chcm1 (sorc phi))
           (efhm0 (efhm chcm0)) (efhm1 (efhm chcm1))
           (lf0 (lf efhm0)) (lg0 (lg efhm0)) (lh0 (lh efhm0))
           (rf0 (rf efhm0)) (rg0 (rg efhm0)) (rh0 (rh efhm0))
           (lf1 (lf efhm1)) (lg1 (lg efhm1)) (lh1 (lh efhm1))
           (rf1 (rf efhm1)) (rg1 (rg efhm1)) (rh1 (rh efhm1))
           (hphi (cmps lg0 (cmps phi lf1)))
           (ephi (cmps rf0 (cmps hphi rg1)))
           (hcone (cone2 hphi))
           (econe (cone2 ephi))
           (LF (cone-2mrph-diag2 hcone cone lf0 lf1))
           (LG (cone-2mrph-diag2 cone hcone lg0 lg1))
           (LH (cone-2mrph-diag2 hcone hcone lh0 (n-mrph -1 lh1)))
           (RF (cone-3mrph-triangle2 hcone econe rf0 rf1
                 (cmps rf0 (cmps hphi rh1))))
           (RG (cone-3mrph-triangle2 econe hcone rg0 rg1
                 (n-mrph -1 (cmps rh0 (cmps hphi rg1)))))
           (RH (cone-3mrph-triangle2 hcone hcone rh0 (n-mrph -1 rh1)
                 (cmps rh0 (cmps hphi rh1)))))
      (declare
       (type chain-complex chcm0 chcm1 hcone econe)
       (type morphism phi lf0 lg0 lh0 rf0 rg0 rh0 lf1 lg1 lh1 rf1 rg1 rh1
                      hphi ephi LF LG LH RF RG RH))
      (build-hmeq
       :lrdct (build-rdct :f LF :g LG :h LH
                          :orgn `(cone2-efhm ,cone lrdct))
       :rrdct (build-rdct :f RF :g RG :h RH
                          :orgn `(cone2-efhm ,cone rrdct))
       :orgn `(cone2-efhm ,cone)))))

(DEFMETHOD SEARCH-EFHM (chcm (orgn (eql 'cone2)))
  (declare (type chain-complex chcm))
  (the homotopy-equivalence
    (cone2-efhm chcm)))

