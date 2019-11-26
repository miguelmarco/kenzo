;;;============================================================================
;;;
;;; Subject: Pushout
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

(provide "pushout")


;;;
;;; Struct
;;;


(DEFSTRUCT (pushout-gsm (:print-function pushout-gsm-print) (:conc-name nil))
  (p-ind #+allegro :type #+allegro (member :undefined 0 1 2))
  ;;;JJJ Why the possible value ":undefined" ???
  (p-old #+allegro :type #+allegro gmsm))

;;;
;;; MACROS
;;;


(defmacro pushout-gsm (p-ind p-old)
  `(make-pushout-gsm :p-ind ,p-ind :p-old ,p-old)    
  )


(DEFMACRO WITH-pushout-gsm (&rest rest)
  (ecase (length (first rest))
    (2 `(with-pushout-gsm-2 ,@rest))
    ))

;;;JJJ  Why these two macros which in fact are the same ?
;;;JJJ  In fact only the with-pushout-gsm-2 can be used.


(DEFMACRO with-pushout-gsm-2 ((p-ind p-old) pushout-gsm . body)
  `(let ((,p-ind (p-ind ,pushout-gsm))
         (,p-old (p-old ,pushout-gsm)))
     (declare (fixnum ,p-ind) (type gmsm ,p-old))
     ,@body))

(defun pushout-gsm-print (pushout-gsm stream depth)
  (declare
   (type pushout-gsm pushout-gsm)
   (stream stream)
   (ignore depth))
  (with-pushout-gsm (p-ind p-old) pushout-gsm
    (format stream "<PUSHOUT-GSM ~A ~A>"
      p-ind
      p-old)
    pushout-gsm))


;;; Pushout Pushout Pushout Pushout Pushout Pushout Pushout Pushout Pushout Pushout
;;; Pushout Pushout Pushout Pushout Pushout Pushout Pushout Pushout Pushout Pushout
;;; Pushout Pushout Pushout Pushout Pushout Pushout Pushout Pushout Pushout Pushout

#|
(defun pushout-cmpr (cmpr1 cmpr2 cmpr3)
  (declare (type cmprf cmpr1 cmpr2 cmpr3))
  (flet ((rslt (pushout-gsm1 pushout-gsm2)
               (case (p-ind pushout-gsm1)
                 (0 (lexico (f-cmpr (p-ind pushout-gsm1) (p-ind pushout-gsm2))
                            (funcall cmpr1 (p-old pushout-gsm1) (p-old pushout-gsm2))))
                 (1 (lexico (f-cmpr (p-ind pushout-gsm1) (p-ind pushout-gsm2))
                            (funcall cmpr2 (p-old pushout-gsm1) (p-old pushout-gsm2))))
                 (2 (lexico (f-cmpr (p-ind pushout-gsm1) (p-ind pushout-gsm2))
                            (funcall cmpr3 (p-old pushout-gsm1) (p-old pushout-gsm2))))
                 )
               ))
    (the cmprf #'rslt)))
|#

;;;JJJ For the fun...

(defun pushout-cmpr (cmpr1 cmpr2 cmpr3)
  (declare (type cmprf cmpr1 cmpr2 cmpr3))
  (flet ((rslt (pushout-gsm1 pushout-gsm2)
               (ecase (+ (* 10 (p-ind pushout-gsm1)) (p-ind pushout-gsm2))
                 ((01 02 12) :less)
                 ((10 20 21) :greater)
                 (00 (funcall cmpr1 (p-old pushout-gsm1) (p-old pushout-gsm2)))
                 (11 (funcall cmpr2 (p-old pushout-gsm1) (p-old pushout-gsm2)))
                 (22 (funcall cmpr3 (p-old pushout-gsm1) (p-old pushout-gsm2))))))
    (the cmprf #'rslt)))



(defun pushout-basis (basis1 basis2 basis3)
  (declare (type basis basis1 basis2 basis3))
  (when (or (eq basis1 :locally-effective)
            (eq basis2 :locally-effective)
            (eq basis3 :locally-effective))
    (return-from pushout-basis :locally-effective))
  (flet ((rslt (dmns)
               (declare (fixnum dmns))
               (let ((basis-from-1 (remove-if #'(lambda (item) 
                                                  (member item '(1 2))) 
                                              (funcall basis1 dmns) :key #'gmsm2)))
                 (append (mapcar #'(lambda (old)
                                     (declare (type gmsm old))
                                     (pushout-gsm 0 old))
                           basis-from-1)
                         (mapcar #'(lambda (old)
                                     (declare (type gmsm old))
                                     (pushout-gsm 1 old))
                           (funcall basis2 dmns))
                         (mapcar #'(lambda (old)
                                     (declare (type gmsm old))
                                     (pushout-gsm 2 old))
                           (funcall basis3 dmns))))))
    (the basis #'rslt)))



(defun pushout-face (face1 face2 face3 f g)
  (declare (type face face1 face2 face3))
  (flet ((rslt (indx dmns gmsm)
               (declare (fixnum indx dmns) (type pushout-gsm gmsm))
               (with-pushout-gsm (ind old) gmsm
                 (ecase ind
                   (0 (with-absm (dgop crpr)
                        (funcall face1 indx dmns old)
                        (with-crpr (dgop1 gmsm1 dgop2 gmsm2) crpr
                          (declare (ignore dgop1 dgop2))
                          (ecase gmsm2
                            (1 (with-absm (dgop gmsm) (? f (1- dmns) (absm dgop gmsm1))
                                 (absm dgop (pushout-gsm 1 gmsm))))
                            (2 (with-absm (dgop gmsm) (? g (1- dmns) (absm dgop gmsm1))
                                 (absm dgop (pushout-gsm 2 gmsm))))
                            (3 (absm dgop (pushout-gsm 0 crpr)))))))
                   (1 (with-absm (dgop gmsm) (funcall face2 indx dmns old)
                        (absm dgop (pushout-gsm 1 gmsm))))
                   (2 (with-absm (dgop gmsm) (funcall face3 indx dmns old)
                        (absm dgop (pushout-gsm 2 gmsm))))))))
    (the face #'rslt)))



(defmethod pushout ((f morphism) (g morphism))
  (if (equal (sorc f) (sorc g))
      (let* ((X (sorc f))
             (XxI (crts-prdc X (delta 1)))
             (Y (trgt f))
             (Z (trgt g)))
        (the simplicial-set
          (let ((rslt (build-smst
                       :cmpr (pushout-cmpr (cmpr XxI) (cmpr Y) (cmpr Z))
                       :basis (pushout-basis (basis XxI) (basis Y) (basis Z))
                       ;;;JJJ               :bspn nil
                       :bspn (pushout-gsm 1 (bspn Y))
                       :face (pushout-face (face XxI) (face Y) (face Z) f g)
                       :orgn `(pushout ,f ,g))))
            (declare (type simplicial-set rslt))
            rslt)))
    (error "The source of f and g must be the same")))



(defun pushout-efhm (f g)
  (declare (morphism f g))
  (if (equal (sorc f) (sorc g))
      (let* ((sorc (sorc f))
             (trgt-f (trgt f))
             (trgt-g (trgt g))
             (rc (remove-covers sorc))
             (ds (direct-sum trgt-f trgt-g))
             (sds (suspension2 ds))
             (p (pushout f g))
             (sigma 
              (build-mrph
               :sorc rc :trgt p :degr 0
               :intr
               #'(lambda (cmbn)
                   (with-cmbn (degr list) cmbn
                     (make-cmbn
                      :degr degr
                      :list (mapcar
                                #'(lambda (term)
                                    (declare (type term term))
                                    (the term
                                      (with-term (cffc gnrt) term
                                        (term
                                         cffc
                                         (pushout-gsm 0 gnrt)))))
                              list)))
                   )
               :strt :cmbn
               :orgn `(sigma ,rc -> ,p)))
             (rho (build-mrph
                   :sorc p :trgt ds :degr 0
                   :intr
                   #'(lambda (cmbn)
                       (with-cmbn (degr list) cmbn
                         (make-cmbn
                          :degr degr
                          :list (mapcar
                                    #'(lambda (term)
                                        (declare (type term term))
                                        (the term
                                          (with-term (cffc gnrt) term
                                            (ecase (p-ind gnrt)
                                              (1 (term
                                                  cffc
                                                  (direct-sum-gsm 0 (p-old gnrt)))
                                               )
                                              (2 (term
                                                  cffc
                                                  (direct-sum-gsm 1 (p-old gnrt))))
                                              )
                                            )))
                                  (remove-if #'(lambda (item) 
                                                  (equal (p-ind (gnrt item)) 0)) 
                                              list))))
                       )
                   :strt :cmbn
                   :orgn `(rho ,p -> ,ds)))
             (shift (build-mrph
                   :sorc ds :trgt sds :degr 1
                   :intr
                   #'(lambda (cmbn)
                       (with-cmbn (degr list) cmbn
                         (make-cmbn
                          :degr (1+ degr)
                          :list list))
                       )
                   :strt :cmbn
                   :orgn `(shift ,ds -> ,sds)))
             (chi (build-mrph
                   :sorc rc :trgt sds :degr 0
                   :intr
                   #'(lambda (cmbn)
                       (with-cmbn (degr list) (? (cmps shift (cmps rho (cmps p sigma))) cmbn)
                         (make-cmbn
                          :degr degr
                          :list list))
                       )
                   :strt :cmbn
                   :orgn `(chi ,rc -> ,sds)))
             (cone2 (cone2 chi))
             (f-isom (build-mrph :sorc cone2 :trgt p :degr 0 :intr
              #'(lambda (cmbn)
                  (with-cmbn (degr list)
                    cmbn
                    (make-cmbn :degr degr :list
                               (mapcar #'(lambda (term)
                                           (declare (type term term))
                                           (ecase (conx (gnrt term))
                                             (0
                                              (term (cffc term)
                                                    (if(=
                                                         (d-sum-ind (icon (gnrt term)))
                                                         0)
                                                        (pushout-gsm 1 (d-sum-old (icon (gnrt term))))
                                                      (pushout-gsm
                                                       2
                                                       (d-sum-old (icon (gnrt term)))))))
                                              (1
                                               (term (cffc term)
                                                     (pushout-gsm
                                                      0
                                                      (icon (gnrt term)))))))
                                 list))))
                                 :strt :cmbn :orgn `(f-isom ,cone2 -> ,p)))
             (g-isom
              (build-mrph :sorc p :trgt cone2 :degr 0 :intr
                          #'(lambda (cmbn)
                              (with-cmbn (degr list)
                                cmbn
                                (make-cmbn :degr degr :list
                               (mapcar #'(lambda (term)
                                           (declare (type term term))
                                           (ecase (p-ind (gnrt term))
                                             (0
                                              (term (cffc term)
                                                    (con1 (p-old (gnrt term)))))
                                              (1
                                               (term (cffc term)
                                                     (con0
                                                      (direct-sum-gsm
                                                       0
                                                       (p-old (gnrt term))))))
                                              (2
                                               (term (cffc term)
                                                     (con0
                                                      (direct-sum-gsm
                                                       1
                                                       (p-old (gnrt term))))))))
                                 list))))
              :strt :cmbn :orgn `(g-isom ,p -> ,cone2)))
             (cone2-efhm (efhm cone2))
             )
        (the homotopy-equivalence 
          (build-hmeq :lf (cmps f-isom (lf cone2-efhm))
                      :lg (cmps (lg cone2-efhm) g-isom)
                      :lh (lh cone2-efhm)
                      :rf (rf cone2-efhm) :rg (rg cone2-efhm)
                      :rh (rh cone2-efhm))
          )
        )       
    (error "The source of f and g must be the same")))


(DEFMETHOD SEARCH-EFHM (pushout (orgn (eql 'pushout)))
  (declare (type chain-complex pushout))
  (the homotopy-equivalence
    (pushout-efhm (second (orgn pushout))
                  (third (orgn pushout))
                  )))



;;; Wedge


(defun wedge (smst1 smst2)
  (declare (simplicial-set smst1 smst2))
  (let* ((unipunctual (build-finite-ss '(x)))
        (f (build-smmr :sorc unipunctual :trgt smst1 :degr 0
                       :sintr #'(lambda (dmns gmsm)
                                  (declare (ignore gmsm))
                               (if (= dmns 0)
                                   (absm 0 (bspn smst1))
                                 nil
                                 ))
                       :orgn `(triv ,unipunctual ,smst1)))
        (g (build-smmr :sorc unipunctual :trgt smst2 :degr 0
                       :sintr #'(lambda (dmns gmsm)
                                  (declare (ignore gmsm))
                               (if (= dmns 0)
                                   (absm 0 (bspn smst2))
                                 nil
                                 ))
                       :orgn `(triv ,unipunctual ,smst2))))
    (pushout f g)
  ))



;;; Join


(defun join (smst1 smst2)
  (declare (simplicial-set smst1 smst2))
  (let* ((smst1xsmst2 (crts-prdc smst1 smst2))
        (f (build-smmr :sorc smst1xsmst2 :trgt smst1 :degr 0
                       :sintr #'(lambda (dmns gmsm)
                                  (declare (ignore dmns))
                               (absm (dgop1 gmsm) (gmsm1 gmsm)))
                       :orgn `(projection ,smst1xsmst2 ,smst1)))
        (g (build-smmr :sorc smst1xsmst2 :trgt smst2 :degr 0
                       :sintr #'(lambda (dmns gmsm)
                                  (declare (ignore dmns))
                               (absm (dgop2 gmsm) (gmsm2 gmsm)))
                       :orgn `(projection ,smst1xsmst2 ,smst2))))
    (pushout f g)
  ))


;;; Smash product

(defun smash-product (smst1 smst2)
  (declare (simplicial-set smst1 smst2))
  (let* ((smst1wsmst2 (wedge smst1 smst2))
         (smst1xsmst2 (crts-prdc smst1 smst2))
         (unipunctual (build-finite-ss '(x)))
         (f (build-smmr :sorc smst1wsmst2 :trgt smst1xsmst2 :degr 0
                        :sintr #'(lambda (dmns gmsm)
                                   (declare (fixnum dmns))
                                   (with-pushout-gsm (ind old) gmsm
                                     (if (= 1 ind) (absm 0 (crpr 0 old (1- (expt 2 dmns)) (bspn  smst2)))
                                       (if (= 2 ind) (absm 0 (crpr (1- (expt 2 dmns)) (bspn  smst1) 0 old ))
                                         (absm (1- (expt 2 dmns)) (bspn smst1xsmst2))))))                                
                        :orgn `(inclusion ,smst1wsmst2 ,smst1xsmst2)))
         (g (build-smmr :sorc smst1wsmst2 :trgt unipunctual  :degr 0
                        :sintr #'(lambda (dmns gmsm)
                                   (declare (ignore gmsm))
                                   (absm (1- (expt 2 dmns)) (bspn unipunctual)))
                        :orgn `(triv ,smst1wsmst2 ,unipunctual))))
    (declare (type simplicial-set smst1wsmst2 smst1xsmst2 unipunctual)
             (type simplicial-mrph f g))
    (pushout f g)))




