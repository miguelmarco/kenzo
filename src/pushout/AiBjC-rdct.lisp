;;;============================================================================
;;;
;;; Subject: AiBjC reduction
;;;         
;;;
;;;----------------------------------------------------------------------------
;;;
;;; Author: Jónathan Heras Vicente
;;; Extracted from Cones of GiftQ program of Francis Sergeraert
;;; 13/10/09
;;;
;;;============================================================================

;;; Adapted by Julián Cuevas, Jose Divasón, Miguel Marco Buzunariz and Ana Romero.
;;; August 2019

(in-package #:cat)
(provide "AiBjC-rdct")

(DEFUN AiBjC-RDCT-F-IMPL (j)
  (declare (type morphism j))
  #'(lambda (cmbn)
        (declare (type cmbn cmbn))
        (the cmbn
          (multiple-value-bind (cmbnB cmbnA)
              (cone-cmbn-split cmbn)
            (declare (type cmbn cmbnB) (ignore cmbnA))
            (? j cmbnB)))))


(DEFUN AiBjC-RDCT-F (A i rho B j sigma C)
  (declare (type chain-complex A B C) (type morphism i rho j sigma))
  (the morphism
    (build-mrph
     :sorc (cone i)
     :trgt C
     :degr 0
     :intr (AiBjC-RDCT-F-IMPL j)
     :strt :cmbn
     :orgn `(AiBjC-rdct-f ,A ,i ,rho ,B ,j ,sigma ,C))))

  
(DEFUN AiBjC-RDCT-G-IMPL (rho B sigma)
  (declare (type morphism sigma) (type chain-complex B))
  
    #'(lambda (cmbn)
        (declare (type cmbn cmbn))
        (the cmbn
          (let ((cmbn-sigma (? sigma cmbn)))
            (declare (type cmbn cmbn-sigma))
            (cone-2cmbn-append
             cmbn-sigma
             (n-cmbn -1 (? rho (? B cmbn-sigma))))))))


(DEFUN AiBjC-RDCT-G (A i rho B j sigma C)
  (declare (type chain-complex A B C) (type morphism i rho j sigma))
  (the morphism
    (build-mrph
     :sorc C :trgt (cone i) :degr 0
     :intr (Aibjc-rdct-g-impl rho B sigma)
     :strt :cmbn
     :orgn `(Aibjc-rdct-g ,A ,i ,rho ,B ,j ,sigma ,C))))


(DEFUN AiBjC-RDCT-H-IMPL (rho)
  (declare (type morphism))
     #'(lambda (cmbn)
        (declare (type cmbn cmbn))
        (the cmbn
          (multiple-value-bind (cmbn0 cmbn1)
              (cone-cmbn-split cmbn)
            (declare (type cmbn cmbn0) (ignore cmbn1))
            (cone-2cmbn-append
             (zero-cmbn (1+ (cmbn-degr cmbn)))
             (? rho cmbn0))))))


(DEFUN AiBjC-RDCT-H (A i rho B j sigma C)
  (declare (type chain-complex A B C) (type morphism i rho j sigma))
  (the morphism
    (build-mrph
     :sorc (cone i) :trgt (cone i) :degr +1
     :intr (AiBjC-rdct-h-impl rho) :strt :cmbn
     :orgn `(AiBjC-RDCT-H ,A ,i ,rho ,B ,j ,sigma ,C))))


(DEFUN AiBjC-RDCT (A i rho B j sigma C)
  (declare (type chain-complex A B C) (type morphism i rho j sigma))
  (the reduction
    (build-rdct
     :f (AiBjC-rdct-f A i rho B j sigma C)
     :g (AiBjC-rdct-g A i rho B j sigma C)
     :h (AiBjC-rdct-h A i rho B j sigma C)
     :orgn `(AiBjC-rdct ,A ,i ,rho ,B ,j ,sigma ,C))))
