;;;  CYLINDERS  CYLINDERS  CYLINDERS  CYLINDERS  CYLINDERS  CYLINDERS  CYLINDERS
;;;  CYLINDERS  CYLINDERS  CYLINDERS  CYLINDERS  CYLINDERS  CYLINDERS  CYLINDERS
;;;  CYLINDERS  CYLINDERS  CYLINDERS  CYLINDERS  CYLINDERS  CYLINDERS  CYLINDERS


(IN-PACKAGE #:cat)

(provide "cyclinders")


;; MACROS

(DEFMACRO CYLNX (cyln)
  `(car ,CYLN))

(DEFMACRO ICYLN (cyln)
  `(cdr ,cyln))

(DEFMACRO CYLNA1 (gnrt)
  `(cons :CYLNA1 ,gnrt))

(DEFMACRO CYLNB (gnrt)
  `(cons :CYLNB ,gnrt))

(DEFMACRO CYLNA2 (gnrt)
  `(cons :CYLNA2 ,gnrt))

(DEFMACRO WITH-CYLN ((cylnx icyln) cyln . body)
  `(let ((,cylnx (cylnx ,cyln))
         (,icyln (icyln ,cyln)))
     (declare
      (type (member :CYLNA1 :CYLNB :CYLNA2) ,cylnx)
      (type gnrt ,icyln))
     ,@body))


(DEFUN DISPATCH-CYLN-CMBN (cmbn)
  (declare (type cmbn cmbn))
  (the (values cmbn cmbn cmbn)
    (with-cmbn (degr list) cmbn
      (let ((lista1 +empty-list+)
            (listb +empty-list+)
            (lista2 +empty-list+))
        (declare (list lista1 listb lista2))
        (dolist (term list)
          (declare (type term term))
          (with-term (cffc cyln) term
            (with-cyln (cylnx icyln) cyln
              (ecase cylnx
                (:cylna1 (push (term cffc icyln) lista1))
                (:cylnb (push (term cffc icyln) listb))
                (:cylna2 (push (term cffc icyln) lista2))))))
        (values
         (make-cmbn :degr (1- degr) :list (nreverse lista1))
         (make-cmbn :degr degr :list (nreverse listb))
         (make-cmbn :degr degr :list (nreverse lista2)))))))


(DEFUN CYLN-CMBN-CMBNA1 (cmbn)
  (declare (type cmbn cmbn))
  (the cmbn
    (with-cmbn (degr list) cmbn
      (make-cmbn :degr (1- degr)
                 :list (do ((rslt +empty-list+)
                            (mark list (cdr mark)))
                           ((endp mark) (nreverse rslt))
                         (declare (list rslt mark))
                         (with--term (cffc cyln) mark
                           (with-cyln (cylnx gnrta1) cyln
                             (if (eq cylnx :cylna1)
                                 (push (term cffc gnrta1) rslt)
                               (return (nreverse rslt))))))))))


(DEFUN CYLN-CMBN-CMBNB (cmbn)
  (declare (type cmbn cmbn))
  (the cmbn
    (with-cmbn (degr list) cmbn
      (make-cmbn :degr degr
                 :list (do ((rslt +empty-list+)
                            (mark (member :cylnb list :key #'cadr) (cdr mark)))
                           ((endp mark) (nreverse rslt))
                         (declare (list rslt mark))
                         (with--term (cffc cyln) mark
                           (with-cyln (cylnx gnrtb) cyln
                             (if (eq cylnx :cylnb)
                                 (push (term cffc gnrtb) rslt)
                               (return (nreverse rslt))))))))))


(DEFUN CYLN-CMBN-CMBNA2 (cmbn)
  (declare (type cmbn cmbn))
  (the cmbn
    (with-cmbn (degr list) cmbn
      (make-cmbn :degr degr
                 :list (mapcar
                           #'(lambda (term)
                               (declare (type term term))
                               (with-term (cffc cyln) term
                                 (term cffc (icyln cyln))))
                         (member :cylna2 list :key #'cadr))))))


(DEFUN MAKE-CYLN-CMBN (cmbna1 cmbnb cmbna2)
  (declare (type cmbn cmbna1 cmbnb cmbna2))
  (the cmbn
    (with-cmbn (degra1 lista1) cmbna1
      (with-cmbn (degrb listb) cmbnb
        (with-cmbn (degra2 lista2) cmbna2
          (unless (= (1+ degra1) degrb degra2)
            (error "In Make-cyln-cmbn, the degrees are not coherent."))
          (make-cmbn :degr degrb
                     :list (nconc
                            (mapcar #'(lambda (terma1)
                                        (declare (type term terma1))
                                        (with-term (cffc gnrta1) terma1
                                          (term cffc (cylna1 gnrta1))))
                              lista1)
                            (mapcar #'(lambda (termb)
                                        (declare (type term termb))
                                        (with-term (cffc gnrtb) termb
                                          (term cffc (cylnb gnrtb))))
                              listb)
                            (mapcar #'(lambda (terma2)
                                        (declare (type term terma2))
                                        (with-term (cffc gnrta2) terma2
                                          (term cffc (cylna2 gnrta2))))
                              lista2))))))))


(DEFMETHOD PRINT-KEYCONS ((car (eql :CYLNA1)) cdr stream)
  (declare
   (type gnrt cdr)
   (stream stream))
  (the (eql t)
    (progn
      (format stream "<CylnA1 ~A>" cdr)
      t)))


(DEFMETHOD PRINT-KEYCONS ((car (eql :CYLNB)) cdr stream)
  (declare
   (type gnrt cdr)
   (stream stream))
  (the (eql t)
    (progn
      (format stream "<CylnB ~A>" cdr)
      t)))


(DEFMETHOD PRINT-KEYCONS ((car (eql :CYLNA2)) cdr stream)
  (declare
   (type gnrt cdr)
   (stream stream))
  (the (eql t)
    (progn
      (format stream "<CylnA2 ~A>" cdr)
      t)))


(DEFUN CYLNX-CMPR (s1 s2)
  (declare (symbol s1 s2))
  (the cmpr
    (if (string= s1 s2) :equal
      (if (string= s1 :cylnA1) :less
        (if (string= s2 :cylnA1) :greater
          (if (string= s1 :cylnB) :less
            :greater))))))


(DEFUN CYLINDER-CMPR (cmpra cmprb)
  (declare (type cmprf cmpra cmprb))
  (flet ((rslt (cyln1 cyln2)
               (let ((cylnx1 (cylnx cyln1)))
                 (declare (type (member :cylna1 :cylnb :cylna2) cylnx1))
                 (lexico
                  (cylnx-cmpr cylnx1 (cylnx cyln2))
                  (case cylnx1
                    (:cylna1 (funcall cmpra (icyln cyln1) (icyln cyln2)))
                    (:cylnb (funcall cmprb (icyln cyln1) (icyln cyln2)))
                    (:cylna2 (funcall cmpra (icyln cyln1) (icyln cyln2))))))))
    (the cmprf #'rslt)))


(DEFUN CYLINDER-BASIS (basisa basisb)
  (declare (type basis basisa basisb))
  (the basis
    (if (or (eq :locally-effective basisa)
            (eq :locally-effective basisb)
            )
        :locally-effective
      (flet ((rslt (degr)
                   (declare (fixnum degr))
                   (append
                    (mapcar #'(lambda (item)
                                (declare (type gnrt item))
                                (cylna1 item))
                      (funcall basisa (1- degr)))
                    (mapcar #'(lambda (item)
                                (declare (type gnrt item))
                                (cylnb item))
                      (funcall basisb degr))
                    (mapcar #'(lambda (item)
                                (declare (type gnrt item))
                                (cylna2 item))
                      (funcall basisa degr)))))
        (the basis #'rslt)))))


(DEFUN CYLINDER-INTR-DFFR (cmpra cmprb dffra dffrb f)
  (declare
   (type cmprf cmpra cmprb)
   (type morphism dffra dffrb f))
  (flet ((rslt (cmbn)
               (declare (type cmbn cmbn))
               (the cmbn
                 (multiple-value-bind (cmbna1 cmbnb cmbna2) (dispatch-cyln-cmbn cmbn)
                   (declare (type cmbn cmbna1 cmbnb cmbna2))
                   (let ((dffra1-cmbn (cmbn-? dffra cmbna1))
                         (dffrb-cmbn (cmbn-? dffrb cmbnb))
                         (dffra2-cmbn (cmbn-? dffra cmbna2))
                         (f-cmbn (cmbn-? f cmbna1))
                         )
                     (declare (type cmbn dffra1-cmbn dffrb-cmbn dffra2-cmbn f-cmbn))
                     (let ((cmbnb (2cmbn-add cmprb f-cmbn dffrb-cmbn))
                           (cmbna2 (2cmbn-sbtr cmpra dffra2-cmbn cmbna1)))
                       (declare (type cmbn cmbnb cmbna2))
                       (make-cyln-cmbn (cmbn-opps dffra1-cmbn) cmbnb cmbna2)))))))
    (the intr-mrph #'rslt)))


(DEFUN CYLINDER (f)
  (declare (type morphism f))
  (the chain-complex
    (with-slots (sorc trgt) f
      (declare
       (type chain-complex sorc trgt ))
      (let ((rslt (build-chcm
                   :cmpr (cylinder-cmpr (cmpr sorc) (cmpr trgt))
                   :basis (cylinder-basis (basis sorc) (basis trgt))
                   :bsgn :undefined
                   :intr-dffr (cylinder-intr-dffr (cmpr sorc) (cmpr trgt)
                                                  (dffr sorc) (dffr trgt) f)
                   :strt :cmbn
                   :orgn `(cylinder ,f))))
        (declare (type chain-complex rslt))
        (slot-makunbound rslt 'bsgn)
        rslt))))


(DEFUN CYLINDER-RRDCT (f)
  (declare (type morphism f))
  (the reduction
    (with-slots (sorc (bcc trgt)) f
      (with-slots ((cmprb cmpr)) bcc
        (let ((cylinder (cylinder f)))
          (let ((rf (build-mrph
                     :sorc cylinder :trgt bcc :degr 0
                     :intr #'(lambda (cmbn)
                               (declare (type cmbn cmbn))
                               (2cmbn-add cmprb (cyln-cmbn-cmbnb cmbn) (cmbn-? f (cyln-cmbn-cmbna2 cmbn))))
                     :strt :cmbn
                     :orgn `(cylinder ,f rf)))
                (rg (build-mrph
                     :sorc bcc :trgt cylinder :degr 0
                     :intr #'(lambda (cmbn)
                               (declare (type cmbn cmbn))
                               (let* ((cmbna1 (zero-cmbn (1- (cmbn-degr cmbn))))
                                      (cmbna2 (zero-cmbn (cmbn-degr cmbn))))
                                 (declare (type cmbn cmbna1 cmbna2))
                                 (make-cyln-cmbn cmbna1 cmbn cmbna2)))
                     :strt :cmbn
                     :orgn `(cylinder ,f rg)))
                (rh (build-mrph
                     :sorc cylinder :trgt cylinder :degr +1
                     :intr #'(lambda (cmbn)
                               (declare (type cmbn cmbn))
                               (let* ((cmbna1 (cmbn-opps (cyln-cmbn-cmbna2 cmbn)))
                                      (cmbnb (zero-cmbn (1+ (cmbn-degr cmbn))))
                                      (cmbna2 (zero-cmbn (1+ (cmbn-degr cmbn)))))
                                 (declare (type cmbn cmbna1 cmbnb cmbna2))
                                 (make-cyln-cmbn cmbna1 cmbnb cmbna2)))
                     :strt :cmbn
                     :orgn `(cylinder ,f rh))))
            (declare
             (type chain-complex cylinder)
             (type morphism rf rg rh))
            (build-rdct :f rf :g rg :h rh :orgn `(Cylinder right reduction ,f))))))))


(DEFUN CYLINDER-LRDCT (f g h k)
  (declare (type morphism f g h k))
  (the reduction
    (with-slots ((acc sorc) (bcc trgt)) f
      (with-slots ((cmpra cmpr)) acc
        (with-slots ((cmprb cmpr)) bcc
          (let ((cylinder (cylinder f)))
            (let ((lf (build-mrph
                       :sorc cylinder :trgt acc :degr 0
                       :intr #'(lambda (cmbn)
                                 (declare (type cmbn cmbn))
                                 (2cmbn-sbtr cmpra
                                             (2cmbn-add cmpra (cyln-cmbn-cmbna2 cmbn) (cmbn-? g (cyln-cmbn-cmbnb cmbn)))
                                             (cmbn-? h (cyln-cmbn-cmbna1 cmbn))))
                       :strt :cmbn
                       :orgn `(cylinder hmtpeq ,f lf)))
                  (lg (build-mrph
                       :sorc acc :trgt cylinder :degr 0
                       :intr #'(lambda (cmbn)
                                 (declare (type cmbn cmbn))
                                 (let* ((cmbna1 (zero-cmbn (1- (cmbn-degr cmbn))))
                                        (cmbnb (zero-cmbn (cmbn-degr cmbn)))
                                        )
                                   (declare (type cmbn cmbna1 cmbnb ))
                                   (make-cyln-cmbn cmbna1 cmbnb cmbn)))
                       :strt :cmbn
                       :orgn `(cylinder hmtpeq ,f lg)))
                  (lh (build-mrph
                       :sorc cylinder :trgt cylinder :degr +1
                       :intr #'(lambda (cmbn)
                                 (declare (type cmbn cmbn))
                                 (multiple-value-bind (cmbna1 cmbnb cmbna2)
                                     (dispatch-cyln-cmbn cmbn)
                                   (let* ((cmbna11 (cmbn-opps (cmbn-? h cmbna1)))
                                          (cmbna12 (cmbn-opps (cmbn-? g (cmbn-? k (cmbn-? f cmbna1)))))
                                          (cmbna13 (cmbn-? g (cmbn-? f (cmbn-? h cmbna1))))
                                          (cmbna14 (cmbn-? g cmbnb))
                                          (cmbnb1 (cmbn-opps (cmbn-? k (cmbn-? k (cmbn-? f cmbna1)))))
                                          (cmbnb2 (cmbn-? k (cmbn-? f (cmbn-? h cmbna1))))
                                          (cmbnb3 (cmbn-? k cmbnb))
                                          (cmbna21 (cmbn-opps (cmbn-? h (cmbn-? h cmbna1))))
                                          (cmbna22 (cmbn-opps (cmbn-? h (cmbn-? g (cmbn-? k (cmbn-? f cmbna1))))))
                                          (cmbna23 (cmbn-? h (cmbn-? g (cmbn-? f (cmbn-? h cmbna1)))))
                                          (cmbna24 (cmbn-? g (cmbn-? k (cmbn-? k (cmbn-? f cmbna1)))))
                                          (cmbna25 (cmbn-opps (cmbn-? g (cmbn-? k (cmbn-? f (cmbn-? h cmbna1))))))
                                          (cmbna26 (cmbn-? h (cmbn-? g cmbnb)))
                                          (cmbna27 (cmbn-opps (cmbn-? g (cmbn-? k cmbnb)))))
                                     (declare (type cmbn cmbna11 cmbna12 cmbna13 cmbna14
                                                    cmbnb1 cmbnb2 cmbnb3
                                                    cmbna21 cmbna22 cmbna23 cmbna24 cmbna25 cmbna26 cmbna27)
                                              (ignore cmbna2))
                                     (make-cyln-cmbn (ncmbn-add cmpra cmbna11 cmbna12 cmbna13 cmbna14)
                                                     (ncmbn-add cmprb cmbnb1 cmbnb2 cmbnb3)
                                                     (ncmbn-add cmpra cmbna21 cmbna22 cmbna23 cmbna24 cmbna25 cmbna26 cmbna27)))))

                       :strt :cmbn
                       :orgn `(cylinder hmtpeq ,f lh))))
              (declare
               (type chain-complex cylinder)
               (type morphism lf lg lh))
              (build-rdct :f lf :g lg :h lh :orgn `(Cylinder left reduction ,f)))))))))


;; Redefined in file "effective-homotopy"
(DEFMETHOD SEARCH-EFHM (smst (orgn (eql 'k-g-1)))
  (declare (type simplicial-group smst))
  (let* ((group (second (orgn smst)))
         (rsltn1 (bar-rsltn group))
         (rsltn2 (resolution group))
         (f (zgmrph-twi (2rsltn-zgmrph rsltn1 rsltn2)))
         (g (zgmrph-twi (2rsltn-zgmrph rsltn2 rsltn1)))
         (h (zgmrph-twi (2rsltn-hmtpop rsltn1 rsltn2)))
         (k (zgmrph-twi (2rsltn-hmtpop rsltn2 rsltn1)))
         (rrdct (cylinder-rrdct f))
         (lrdct (cylinder-lrdct f g h k)))
    (the homotopy-equivalence
      (build-hmeq :lrdct lrdct :rrdct rrdct :orgn `(efhm of ,smst)))))

