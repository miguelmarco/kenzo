;;;============================================================================
;;;
;;; Subject: Remove Covers
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

(provide "remove-covers")

(defun remove-covers-basis (basis1)
  (declare (type basis basis1))
  (when (eq basis1 :locally-effective)
    (return-from remove-covers-basis :locally-effective))
  (flet ((rslt (dmns)
               (remove-if #'(lambda (item) 
                              (member item '(1 2))) 
                          (funcall basis1 dmns) :key #'gmsm2)
               ))
    (the basis #'rslt)))



(defun remove-covers-dffr (dffr1)
  (declare (type morphism dffr1))
  (flet ((rslt (cmbn)
               (declare (type cmbn cmbn))
               (with-cmbn (degr list) cmbn
                 (declare (ignore list))
                 (make-cmbn :degr (1- degr)
                            :list (remove-if #'(lambda (item) 
                              (member (gnrt2 (gnrt item)) '(1 2))) 
                          (cmbn-list (? dffr1 cmbn))))
                 )
               ))
    (the intr-mrph #'rslt)))



(defun remove-covers (chcm1)
  (declare (type simplicial-set chcm1))
  (let ((XxI (crts-prdc chcm1 (delta 1))))
  (the chain-complex
    (build-chcm
     :cmpr (cmpr XxI)
     :basis (remove-covers-basis (basis XxI))
     :bsgn nil
     :intr-dffr (remove-covers-dffr (dffr XxI))
     :strt :cmbn
     :orgn `(remove-covers ,chcm1)))))






(DEFUN remove-covers-efhm (chcm)
  (let* ((A (direct-sum (second (orgn chcm)) (second (orgn chcm))))
         (B (crts-prdc (second (orgn chcm)) (delta 1)))
         (C chcm)
         (i (build-mrph
             :sorc A :trgt B :degr 0
             :intr
             #'(lambda (cmbn)
                 (with-cmbn (degr list) cmbn
                   (make-cmbn
                      :degr degr
                      :list (mapcar
                                #'(lambda (term)
                                    (with-term (cffc gnrt) term
                                      (if (= (d-sum-ind gnrt) 0)
                                          (term cffc
                                                (crpr 0 (d-sum-old gnrt) (mask degr) 1))
                                        (term cffc
                                                (crpr 0 (d-sum-old gnrt) (mask degr) 2)))))
                              list))))
             :strt :cmbn
             :orgn `(,A -> ,B)))
         (rho (build-mrph
             :sorc B :trgt A :degr 0
             :intr
             #'(lambda (cmbn)
                 (with-cmbn (degr list) cmbn
                   (make-cmbn
                      :degr degr
                    :list (mapcar #'(lambda (mnm)
                                      (if (equal (gmsm2 (gnrt mnm)) 1)
                                          (term (cffc mnm)
                                                (direct-sum-gsm 0 (gmsm1 (gnrt mnm)))
                                                )
                                        (term (cffc mnm)
                                              (direct-sum-gsm 1 (gmsm1 (gnrt mnm))))))
                    (remove-if #'(lambda (item) 
                                               (member (gmsm2 (gnrt item)) '(3))) 
                                           list)))))
                 :strt :cmbn
               :orgn `(,B -> ,A)))
         (sigma (build-mrph
             :sorc C :trgt B :degr 0
             :intr
             #'(lambda (cmbn)
                 (with-cmbn (degr list) cmbn
                   (make-cmbn
                      :degr degr
                      :list list)))
             :strt :cmbn
                 :orgn `(,C -> ,B)))
         (j (build-mrph
             :sorc B :trgt C :degr 0
             :intr
             #'(lambda (cmbn)
                 (with-cmbn (degr list) cmbn
                   (make-cmbn
                      :degr degr
                      :list (remove-if #'(lambda (item) 
                                               (member (gmsm2 (gnrt item)) '(1 2))) 
                                           list))))
             :strt :cmbn
                 :orgn `(,B -> ,C)))
         (cone-i (cone i))
         (efhm-cone (efhm cone-i))
         (rrdct-efhm-cone (rrdct efhm-cone))
         (lrdct-efhm-cone (lrdct efhm-cone))
         (rdct-C-cone-i (aibjc-rdct A i rho B j sigma C))
         (final-lrdct (cmps rdct-C-cone-i lrdct-efhm-cone))
         )
  
  (the homotopy-equivalence (progn                                                   
                              (build-hmeq
                               :lrdct final-lrdct
                               :rrdct rrdct-efhm-cone
                               :orgn `(remove-covers-efhm ,chcm)))))
  )


(DEFMETHOD SEARCH-EFHM (remove-covers (orgn (eql 'remove-covers)))
  (declare (type chain-complex remove-covers))
  (the homotopy-equivalence
    (remove-covers-efhm remove-covers
                     )))

