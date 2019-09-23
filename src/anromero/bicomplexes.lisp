;; BICOMPLEXES   BICOMPLEXES   BICOMPLEXES   BICOMPLEXES   BICOMPLEXES   BICOMPLEXES  
;; BICOMPLEXES   BICOMPLEXES   BICOMPLEXES   BICOMPLEXES   BICOMPLEXES   BICOMPLEXES
;; BICOMPLEXES   BICOMPLEXES   BICOMPLEXES   BICOMPLEXES   BICOMPLEXES   BICOMPLEXES

(IN-PACKAGE #:cat)

(PROVIDE "bicomplexes")

(DEFMETHOD PRINT-KEYCONS ((car (eql :BcGnrt)) cdr stream)
  (declare
   (cons cdr)
   (stream stream))
  (the (eql t)
    (let* 
        ((cdr (cons car cdr))
         (degr1 (BcGnrt-Degr1 cdr))
         (degr2 (BcGnrt-Degr2 cdr))
         (gnrt (BcGnrt-gnrt cdr)))
      (progn
        (format stream "<BcGnrt ")
        (format stream "[~D ~D] ~A" degr1 degr2 gnrt)
        (format stream ">"))
      t)))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; TYPES FOR BICOMPLEXES ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; TYPE BcGnrt (a generator in a Bicomplex)
(DEFUN BcGnrt-p (object)
  (declare (type any object))
  (the boolean
    (and (consp object)
         (eql :BcGnrt (car object))
         (consp (cdr object))
         (consp (cadr object))
         (typep (caadr object) 'fixnum)
         (typep (cdadr object) 'fixnum))))
(DEFTYPE BcGnrt () '(satisfies BcGnrt-p))

#|
(typep  `(:BcGnrt ,(cons 5 6) u) 'BcGnrt)
(typep `(:BcGnrt ,(cons 3 2) '(f g)) 'BcGnrt)
(typep `(:BcGnrt ,(cons 4 5.4) u) 'BcGnrt)
|#

 
;; TYPE BcBasis (basis of a Bicomplex, function of the variables p and q or
;; ":locally-effective")
(DEFTYPE BcBASIS () '(or function (eql :locally-effective))) 
;; (function (degr degr) (list gnrt)) 


;; TYPE Bicomplex
(DEFUN OrgnBiCmpl-p (object)
  (declare (type chain-complex object))
  (let* ((orgn (orgn object)))
    (the boolean
      (and (consp orgn)
           (eql (car orgn) 'BiCmpl)))))
(DEFTYPE Bicomplex () '(satisfies OrgnBiCmpl-p))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; FUCTIONS FOR THE REPRESENTATION OF (GENERAL) BICOMPLEXES ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; Function that builds a generator of a Bicomplex from two degrees and 
;; a simple generator.
(DEFUN build-BcGnrt (degr1 degr2 gnrt)
  (declare 
   (type fixnum degr1 degr2)
   (type any gnrt))
  (list :BcGnrt (cons degr1 degr2) gnrt))

#|
(build-BcGnrt 2 3 'c)
(typep  (build-BcGnrt 2 2 'd) 'BcGnrt)
|#
 

;; Function that returns the first degree from a bicomplex' generator.
(DEFUN BcGnrt-Degr1 (BcGnrt)
  (declare (type BcGnrt BcGnrt))
  (caadr BcGnrt))

;; Function that returns the second degree from a bicomplex' generator.
(DEFUN BcGnrt-Degr2 (BcGnrt)
  (declare (type BcGnrt BcGnrt))
  (cdadr BcGnrt))

;; Function that returns the original generator from a bicomplex' generator.
(DEFUN BcGnrt-Gnrt (BcGnrt)
  (declare (type BcGnrt BcGnrt))
  (caddr BcGnrt))

#|
(build-BcGnrt 5 6 8)
(BcGnrt-Degr1 (build-BcGnrt 5 6 8))
(BcGnrt-Degr2 (build-BcGnrt 6 4 '(4 6)))
(BcGnrt-Gnrt (build-BcGnrt 6 4 '(4 6)))
(BcGnrt-Gnrt (build-BcGnrt 2 8 'u))
|# 

 
;; Comparison function of a bicomplex. 
(defun BC-CMPR (cmpr)
  (declare (type cmprf cmpr))
  (flet ((comp (bc1 bc2)
               (declare (type BcGnrt bc1 bc2))
               (let ((degr11 (BcGnrt-Degr1 bc1))
                     (degr21 (BcGnrt-Degr2 bc1))
                     (degr12 (BcGnrt-Degr1 bc2))
                     (degr22 (BcGnrt-Degr2 bc2))
                     (gnrt1 (BcGnrt-Gnrt bc1))
                     (gnrt2 (BcGnrt-Gnrt bc2)))
                 (declare 
                  (type fixnum degr11 degr21 degr12 degr22)
                  (type any gnrt1 gnrt2))
                 (lexico
                  (f-cmpr degr11 degr12)
                  (f-cmpr degr21 degr22)
                  (funcall cmpr gnrt1 gnrt2)))))
    (the cmprf #'comp)))

#|
(setf cmpr (BC-CMPR #'s-cmpr)) 
(funcall cmpr (build-BcGnrt 5 6 'a) (build-BcGnrt 7 8 'b))
(funcall cmpr (build-BcGnrt 8 6 'a) (build-BcGnrt 7 8 'b))
(funcall cmpr (build-BcGnrt 8 5 'c) (build-BcGnrt 8 6 'b))
(funcall cmpr (build-BcGnrt 8 5 'c) (build-BcGnrt 8 5 'b))
|# 
      
;; Basis of a bicomplex (from a function of type BcBasis).
(DEFUN BC-BASIS (bcbasis)
  (declare (type BcBasis bcbasis))
  (when (eq bcbasis :locally-effective)
    (return-from bc-basis :locally-effective))
  (flet ((bas (degr)
              (declare (type fixnum degr))
              (the list
                (progn
                  (when (minusp degr)
                    (return-from bas +empty-list+))
                  (mapcan
                      #'(lambda (degr1)
                          (declare (type fixnum degr1))
                          (let* ((degr2 (- degr degr1))
                                 (basis (funcall bcbasis degr1 degr2)))
                            (declare
                             (fixnum degr2)
                             (list basis))
                            (the list
                              (mapcar
                                  #'(lambda (gnrt)
                                      (declare (type gnrt gnrt))
                                      (the BcGnrt
                                        (build-BcGnrt degr1 degr2 gnrt)))
                                basis)))
                          )
                    (<a-b> 0 degr))))))
    (the basis #'bas)))

#|
(defun bas (degr1 degr2)
  (if (and (= degr1 0) (= degr2 1)) (return-from bas '(a)))
  (if (and (= degr1 1) (= degr2 0)) (return-from bas '(b)))
  (if (and (= degr1 1) (= degr2 1)) (return-from bas '(c)))
  (if (and (= degr1 2) (= degr2 0))  (return-from bas '(d)))
  (return-from bas nil))
(setf basis (BC-BASIS #'bas))
(dotimes (i 5)
  (print (funcall basis i)))
|#

;; Differential function of a Bicomplex (from the differential functions,
;; dffr1=horizontal dffr, dffr2=vertical dffr). 
(DEFUN BC-INTR-DFFR (dffr1 dffr2)
  (declare (type function dffr1 dffr2))
  (flet ((dif (degr bc)
              (declare
               (fixnum degr)
               (type BcGnrt bc))
              (the cmbn
                (let* ((degr1 (BcGnrt-Degr1 bc))
                       (degr2 (BcGnrt-Degr2 bc))
                       (gnrt (BcGnrt-Gnrt bc))
                       (degr1-1 (1- degr1))
                       (degr2-1 (1- degr2))
                       (list1 (funcall dffr1 degr1 degr2 gnrt))
                       (list2 (funcall dffr2 degr1 degr2 gnrt)))
                  (declare (type fixnum degr1 degr2 degr1-1)
                           (list list1 list2))
                  (make-cmbn
                   :degr (1- degr)
                   :list (nconc
                          (mapcar
                              #'(lambda (term1)
                                  (declare (type term term1))
                                  (let* ((cffc1 (cffc term1))
                                         (gnrt1 (gnrt term1)))
                                    (declare (type fixnum cffc1)
                                             (type gnrt gnrt1))
                                    (term cffc1 (build-BcGnrt degr1-1 degr2 gnrt1))))
                            list1)
                          (mapcar 
                              #'(lambda (term2)
                                  (declare (type term term2))
                                  (let* ((cffc2 (cffc term2))
                                         (gnrt2 (gnrt term2)))
                                    (declare (type fixnum cffc2)
                                             (type gnrt gnrt2))
                                    (term cffc2 (build-BcGnrt degr1 degr2-1 gnrt2))))
                            list2)))))))
    (the intr-mrph #'dif)))

#|
(defun dif1 (degr1 degr2 gnrt)
  (if (and (= degr1 1) (= degr2 1) (eql gnrt 'c)) (return-from dif1 (list (cons 2 'a))))
  (if (and (= degr1 2) (= degr2 0) (eql gnrt 'd)) (return-from dif1 (list (cons 2 'b))))
  (return-from dif1 nil))
(defun dif2 (degr1 degr2 gnrt)
  (if (and (= degr1 1) (= degr2 1) (eql gnrt 'c)) (return-from dif2 (list (cons 1 'b))))
  (return-from dif2 nil))
(setf dif (BC-INTR-DFFR #'dif1 #'dif2))
(funcall dif 2 (build-BcGnrt 1 1 'c))
|#

;; Function that builds a bicomplex from a basis function, the two differential
;; functions and the original comparison. 
(DEFUN BUILD-BICM (&key bcbasis dffr1 dffr2 cmpr orgn)
  (declare 
   (type BCbasis bcbasis)
   (type function dffr1 dffr2)
   (type cmprf cmpr)
   (type list orgn))
  (the chain-complex
    (let ((chcm
           (build-chcm
            :cmpr (BC-CMPR cmpr)
            :basis (BC-BASIS bcbasis)
            :intr-dffr (BC-INTR-DFFR dffr1 dffr2)
            :strt :gnrt
            :orgn orgn)))
      (declare (type chain-complex chcm))
      (slot-makunbound chcm 'bsgn)
      chcm)))

#|
(setf bc (Build-Bicm :bcbasis #'bas :dffr1 #'dif1 :dffr2 #'dif2 :cmpr #'s-cmpr 
                     :orgn '(BC-test)))
(dotimes (i 5)
  (print (basis bc i)))
(cmpr bc (build-BcGnrt 0 1 'a) (build-BcGnrt 1 1 'c))

(? bc (? bc 2 (build-BcGnrt 1 1 'c))) 
(? bc 2 (build-BcGnrt 2 0 'd))
(? bc 2 '(:BcGnrt (1 . 1) c))
(? bc (cmbn 2 3 '(:BcGnrt (1 . 1) c)))
(? bc 1 (build-BcGnrt 0  1 'a))
|#


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;; BICOMPLEX FROM A LIST OF MORPHISMS OF CHAIN COMPLEXES ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; EACH CHAIN COMPLEX IS A COLUMN OF THE BICOMPLEX
;; VERTICAL DIFFERENTIAL MAPS ARE THE DIFFERENTIAL MAPS OF THE CHAIN COMPLEXES
;; HORIZONTAL DIFFERENTIAL MAPS ARE THE MAPS INDUCED BY THE MORPHISMS


;; Function that builds the bicomplex basis from a list of morphisms of chain complexes. 
(DEFUN LIST-OF-MRPH-BCBASIS (l)
  (declare (type list l))
  (when (eq (basis (sorc (first l))) :locally-effective)
    (return-from list-of-mrph-bcbasis :locally-effective))
  (dotimes (i (length l))
    (declare (type fixnum i))
    (when (eq (basis (trgt (nth i l))) :locally-effective)
      (return-from list-of-mrph-bcbasis :locally-effective)))
  (flet ((bcbasis (degr1 degr2)
                  (declare (type fixnum degr1 degr2))
                  (if (< degr1 (length l)) 
                      (basis (trgt (nth degr1 l)) degr2)
                    (if (= degr1 (length l))
                        (basis (sorc (nth (1- (length l)) l)) degr2)
                      nil
                      ))))
    (the bcbasis #'bcbasis)))

;; Function that builds the bicomplex horizontal internal differential
;; from a list of morphisms of chain complexes. 
(DEFUN LIST-OF-MRPH-DFFR1 (l)
  (declare (type list l))
  (flet ((dffr1 (degr1 degr2 gnrt)
                (declare (type fixnum degr1 degr2)
                         (type gnrt gnrt))
                (if (and (> degr1 0) (<= degr1 (length l)))
                    (let ((mrph (nth (1- degr1) l)))
                      (declare (type morphism mrph))
                      (cmbn-list (? mrph degr2 gnrt)))
                  nil)))
    (the function #'dffr1)))

;; Function that builds the bicomplex vertical internal differential
;; from a list of morphisms of chain complexes.
(DEFUN LIST-OF-MRPH-DFFR2 (l)
  (declare (type list l))
  (flet ((dffr2 (degr1 degr2 gnrt)
                (declare (type fixnum degr1 degr2)
                         (type gnrt gnrt))
                (if (and (> degr1 0) (<= degr1 (length l)))
                    (let ((mrph (nth (1- degr1) l)))
                      (cmbn-list (? (sorc mrph) degr2 gnrt)))
                  (if (= degr1 0)
                      (cmbn-list (? (trgt (first l)) degr2 gnrt))
                    nil))))
    (the function #'dffr2)))

;; Function that constructs the comparison function of a bicomplex 
;; from a list of morphisms of chain complexes.
(DEFUN LIST-OF-MRPH-BC-CMPR (l)
  (declare (type list l))
  (flet ((comp (bc1 bc2)
               (declare (type BcGnrt bc1 bc2))
               (let ((degr11 (BcGnrt-Degr1 bc1))
                     (degr21 (BcGnrt-Degr2 bc1))
                     (degr12 (BcGnrt-Degr1 bc2))
                     (degr22 (BcGnrt-Degr2 bc2))
                     (gnrt1 (BcGnrt-Gnrt bc1))
                     (gnrt2 (BcGnrt-Gnrt bc2)))
                 (declare 
                  (type fixnum degr11 degr21 degr12 degr22)
                  (type any gnrt1 gnrt2))
                 (lexico
                  (f-cmpr degr11 degr12)
                  (f-cmpr degr21 degr22)
                  (if (and (> degr11 0) (<= degr11 (length l)))
                      (let ((chcm (sorc (nth (1- degr11) l))))
                        (funcall (cmpr chcm) gnrt1 gnrt2))
                    (if (= degr11 0)
                        (let ((chcm (trgt (nth 0 l))))
                          (funcall (cmpr chcm) gnrt1 gnrt2))
                      :equal))
                  ))))
    (the cmprf #'comp)))

;; Function that builds a bicomplex from a list of morphisms of chain complexes.
(DEFUN LIST-OF-MRPH-BUILD-BICM (l)
  (declare (type list l))
  (dotimes (i (1- (length l)))
    (declare (fixnum i))
    (let* ((mrphi (nth i l))
           (mrphi+1 (nth (1+ i) l))           
           (sorci (sorc mrphi))
           (trgti+1 (trgt mrphi+1)))
      (unless (eq (grmd sorci) (grmd trgti+1))
        (error "In LIST-OF-MRPH-BUILD-BICM, the morphisms ~A and ~A may not be composed (cf source and target)."
          mrphi mrphi+1))))
  (the chain-complex
    (let ((chcm
           (build-chcm
            :cmpr (list-of-mrph-BC-CMPR l)
            :basis (BC-BASIS (list-of-mrph-bcbasis l))
            :intr-dffr (BC-INTR-DFFR (list-of-mrph-dffr1 l) (list-of-mrph-dffr2 l))
            :strt :gnrt
            :orgn `(BiCmpl from list  ,l))))
      (declare (type chain-complex chcm))
      (slot-makunbound chcm 'bsgn)
      chcm)))

#|
;; Example of bicomplex from a list of morphisms of two chain complexes

;; First we build the three chain complexes

(cat-init)
(defun chcm1-basis (degr)
  (if (= degr 1) 
      (return-from chcm1-basis (list 'a1 'a2))
    (return-from chcm1-basis nil)))

(defun chcm2-basis (degr)
  (if (= degr 0) (return-from chcm2-basis (list 'b))
    (if (= degr 1) (return-from chcm2-basis (list 'c1 'c2))
      (return-from chcm2-basis nil))))

(defun chcm3-basis (degr)
  (if (= degr 0) 
      (return-from chcm3-basis (list 'd1 'd2))
    (return-from chcm3-basis nil)))

(defun chcm1-intr-dffr (degr gnrt)
  (return-from chcm1-intr-dffr (cmbn (1- degr))))

(defun chcm2-intr-dffr (degr gnrt)
  (if (and (= degr 1) (eql gnrt 'c1))
      (return-from chcm2-intr-dffr (cmbn 0 1 'b))
    (return-from chcm2-intr-dffr (cmbn (1- degr)))))

(defun chcm3-intr-dffr (degr gnrt)
  (return-from chcm3-intr-dffr (cmbn (1- degr))))

(setf chcm1 (build-chcm :cmpr #'s-cmpr :basis #'chcm1-basis :intr-dffr #'chcm1-intr-dffr :strt :gnrt :orgn `(chcm1)))
(setf chcm2 (build-chcm :cmpr #'s-cmpr :basis #'chcm2-basis :intr-dffr #'chcm2-intr-dffr :strt :gnrt :orgn `(chcm2)))
(setf chcm3 (build-chcm :cmpr #'s-cmpr :basis #'chcm3-basis :intr-dffr #'chcm3-intr-dffr :strt :gnrt :orgn `(chcm3)))          


;; Then we build the two morphisms

(defun mrph1-intr (degr gnrt)
  (if (and (= degr 1) (eql gnrt 'c1))
      (return-from mrph1-intr (cmbn 1 2 'a1))
    (if (and (= degr 1) (eql gnrt 'c2))
        (return-from mrph1-intr (cmbn 1 2 'a2))
      (return-from mrph1-intr (cmbn degr)))))

(defun mrph2-intr (degr gnrt)
  (if (and (= degr 0) (eql gnrt 'd1))
      (return-from mrph2-intr (cmbn 0 2 'b))
    (return-from mrph2-intr (cmbn degr))))

(setf mrph1 (build-mrph :sorc chcm2 :trgt chcm1 :intr #'mrph1-intr :strt :gnrt :orgn `(mrph1)))
(setf mrph2 (build-mrph :sorc chcm3 :trgt chcm2 :intr #'mrph2-intr :strt :gnrt :orgn `(mrph2)))

;; Finally we build the list and the bicomplex and we change it into a filtered chain complex

(setf l (list mrph1 mrph2))
(setf bic (list-of-mrph-build-bicm l))
(change-bicm-to-flcc bic)


;; We compute the groups of the associated spectral sequence
(dotimes (r 5)
  (dotimes (n 3)
    (dotimes (p (1+ n))
      (let ((q (- n p)))
        (spsq-group bic r p q)
        ))))


|#
