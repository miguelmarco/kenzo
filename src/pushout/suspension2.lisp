;;;============================================================================
;;;
;;; Subject: Suspensions2
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

(provide "suspension2")

(DEFUN suspension2-basis (basis)
  (declare (type basis basis))
  (when (eq basis :locally-effective)
    (return-from suspension2-basis :locally-effective))
  (flet ((rslt (degr)
	   (declare (fixnum degr))
	   (funcall basis (1- degr))))
    (the basis #'rslt)))



(defun suspension2-intr-dffr (dffr)
  (declare (type morphism dffr))
  (flet ((rslt (cmbn)
	   (declare (type cmbn cmbn))
               (with-cmbn (degr list) cmbn
                 (if (> degr 0)
                     (let ((cmbn-bis (make-cmbn :degr (1- degr) :list list)))
                       (make-cmbn :degr (1- degr)
                                  :list (mapcar #'(lambda (old)
                                                    (term (* -1 (cffc old))
                                                          (gnrt old))
                                                    )
                                          (cmbn-list (? dffr cmbn-bis)) )
                                  ))
                   (zero-cmbn (1- degr))))))
    (the intr-mrph #'rslt))  
  )



(DEFMETHOD SUSPENSION2 ((chcm chain-complex))
  (the chain-complex
    (with-slots (cmpr basis dffr) chcm
      (build-chcm
       :cmpr cmpr
       :basis (suspension2-basis basis)
       :bsgn nil
       :intr-dffr (suspension2-intr-dffr dffr)
       :strt :cmbn
       :orgn `(suspension2 ,chcm)))))





(DEFUN SUSPENSION2-INTR (mrph)
  (declare (type morphism mrph))
  (flet ((rslt (cmbn)
               (declare (type cmbn cmbn))
               (with-cmbn (degr list) cmbn
                 (make-cmbn :degr (+ degr (degr mrph))
                            :list (mapcar #'(lambda (old)
                                     (term (* -1 (cffc old))
                                           (gnrt old))
                                     )
                           (cmbn-list
				 (? mrph (make-cmbn :degr (1- degr)
							 :list list))) )
                            ))))
    (the intr-mrph #'rslt)))




(DEFMETHOD SUSPENSION2 ((mrph morphism))
  (the morphism
    (let ((orgn (orgn mrph)))
      (declare (list orgn))
      (when (eq (first orgn) 'zero-mrph)
	(return-from suspension2
	  (zero-mrph (suspension2 (sorc mrph))
		     (suspension2 (trgt mrph))
		     (degr mrph))))
      (when (eq (first orgn) 'idnt-mrph)
	(return-from suspension2
	  (idnt-mrph (suspension2 (sorc mrph)))))
      (build-mrph
       :sorc (suspension2 (sorc mrph))
       :trgt (suspension2 (trgt mrph))
       :degr (degr mrph)
       :intr (suspension2-intr mrph)
       :strt :cmbn
       :orgn `(suspension2 ,mrph)))))




(DEFMETHOD SUSPENSION2 ((rdct reduction))
  (the reduction (progn
    (when (eq (first (orgn rdct)) 'trivial-rdct)
      (return-from suspension2
	(trivial-rdct (suspension2 (second (orgn rdct))))))
    (build-rdct
      :f (suspension2 (f rdct))
      :g (suspension2 (g rdct))
      :h (suspension2 (h rdct))
     :orgn `(suspension2 ,rdct)))))


(DEFMETHOD SUSPENSION2 ((hmeq homotopy-equivalence))
  (the homotopy-equivalence (progn
    (when (eq (first (orgn hmeq)) 'trivial-hmeq)
      (return-from suspension2
	(trivial-hmeq (suspension2 (second (orgn hmeq))))))
    (build-hmeq
      :lrdct (suspension2 (lrdct hmeq))
      :rrdct (suspension2 (rrdct hmeq))
     :orgn `(suspension2 ,hmeq)))))


(DEFMETHOD SEARCH-EFHM (suspension2 (orgn (eql 'suspension2)))
  (declare (type chain-complex suspension2))
  (suspension2 (efhm (second (orgn suspension2)))))





