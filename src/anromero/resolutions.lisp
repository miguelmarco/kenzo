;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*

;; RESOLUTIONS   RESOLUTIONS   RESOLUTIONS   RESOLUTIONS   RESOLUTIONS   RESOLUTIONS
;; RESOLUTIONS   RESOLUTIONS   RESOLUTIONS   RESOLUTIONS   RESOLUTIONS   RESOLUTIONS
;; RESOLUTIONS   RESOLUTIONS   RESOLUTIONS   RESOLUTIONS   RESOLUTIONS   RESOLUTIONS


(IN-PACKAGE #:cat)

(provide "resolutions")


;; CLASSES AND TYPES


(DEFUN ZGGNRT-P (object)
  (declare (type any object))
  (the boolean
    (and (consp object)
         (eql :ZGgnrt (car object))
         (consp (cdr object))
         (typep (cadr object) 'gnrt)
         (typep (cddr object) 'gnrt)
         )))


(DEFTYPE ZGGNRT () '(satisfies zggnrt-p))


(DEFUN ZGTERM-P (object)
  (declare (type any object))
  (the boolean
    (and (consp object)
         (typep (car object) 'fixnum)
         (typep (cdr object) 'zggnrt))))


(DEFTYPE ZGTERM () '(satisfies zgterm-p))


;; IZGCMBN = Internal-ZGCoMBiNation
(DEFUN IZGCMBN-P (object)
  ;; IZGCMBN = Internal ZGCoMBiNation,
  ;; without the keyword and the degree
  (declare (type any object))
  (the boolean
    (and (listp object)
         (every #'zgterm-p object))))


(DEFTYPE IZGCMBN () '(satisfies izgcmbn-p))


;; ZGCMBN = ZGCoMBiNation
(DEFUN ZGCMBN-P (object)
  (declare (type any object))
  (the boolean
    (and (consp object)
         (eql (car object) :cmbn)
         (consp (cdr object))
         (typep (second object) 'fixnum)
         (typep (cddr object) 'izgcmbn))))


(DEFTYPE ZGCMBN () '(satisfies zgcmbn-p))


(DEFTYPE GROUP-BASIS () '(or list (eql :locally-effective)))


;;(DEFCLASS GROUP-MRPH () ()) ;; to be redefined later


(DEFCLASS GROUP ()
  ((elements :type group-basis :initarg :elements :reader elements)
   (cmpr :type cmprf :initarg :cmpr :reader cmpr1)
   (mult :type function :initarg :mult :reader mult1)
   (inv :type function :initarg :inv :reader inv1)
   (nullel :type gnrt :initarg :nullel :reader nullel)
   (idnm :type fixnum :initform (incf *idnm-counter*) :reader idnm)
   (orgn :type list :initarg :orgn :reader orgn)
   (resolution :type reduction :initarg :resolution :reader resolution))
  )


(DEFMETHOD INITIALIZE-INSTANCE :after ((group group) &rest rest)
  (set (intern (format nil "K~D" (idnm group))) group))


(DEFCLASS AB-GROUP (GROUP)
  ())


(DEFCLASS GROUP-MRPH ()
  ;; SOuRCe
  ((sorc :type group :initarg :sorc :reader sorc)
   ;; TaRGeT
   (trgt :type group :initarg :trgt :reader trgt)
   ;; INTeRnal
   (intr :type function :initarg :intr :reader intr)
   ;; IDentification NuMber
   (idnm :type fixnum :initform (incf *idnm-counter*) :reader idnm)
   ;; ORiGiN
   (orgn :type list :initarg :orgn :reader orgn)))


(DEFMETHOD INITIALIZE-INSTANCE :after ((group-mrph group-mrph) &rest rest)
  (set (intern (format nil "K~D" (idnm group-mrph))) group-mrph))


(DEFVAR *GROUP-LIST*
    "The variable *GROUP-LIST* is bound to a list of user created groups.")
(SETF *GROUP-LIST* +EMPTY-LIST+)
(PUSHNEW '*GROUP-LIST* *LIST-LIST*)

(DEFVAR *AB-GROUP-LIST*
    "The variable *AB-GROUP-LIST* is bound to a list of user created abelian groups.")
(SETF *AB-GROUP-LIST* +EMPTY-LIST+)
(PUSHNEW '*AB-GROUP-LIST* *LIST-LIST*)

(DEFVAR *GROUP-MRPH-LIST*
    "The variable *GROUP-MRPH-LIST* is bound to a list of user created group morphisms.")
(SETF *GROUP-MRPH-LIST* +EMPTY-LIST+)
(PUSHNEW '*GROUP-MRPH-LIST* *LIST-LIST*)

#+clisp(eval-when (:compile-toplevel :load-toplevel :execute)
         (setf (ext:package-lock :clos) nil))
(DEFMETHOD PRINT-OBJECT ((group group) stream)
  (the group
    (progn
      (format stream "[K~D Group]" (idnm group))
      group)))

(DEFMETHOD PRINT-OBJECT ((group ab-group) stream)
  (the ab-group
    (progn
      (format stream "[K~D Abelian-Group]" (idnm group))
      group)))

(DEFMETHOD PRINT-OBJECT ((group-mrph group-mrph) stream)
  (the group-mrph
    (progn
      (format stream "[K~D Group-morphism]" (idnm group-mrph))
      group-mrph)))
#+clisp(eval-when (:compile-toplevel :load-toplevel :execute)
         (setf (ext:package-lock :clos) t))


(DEFUN GROUP (idnm)
  (declare (type fixnum idnm))
  (the (or group null)
    (find idnm *group-list* :key #'idnm)))

(DEFUN AB-GROUP (idnm)
  (declare (type fixnum idnm))
  (the (or ab-group null)
    (find idnm *ab-group-list* :key #'idnm)))

(DEFUN GROUP-MRPH (idnm)
  (declare (type fixnum idnm))
  (the (or group-mrph null)
    (find idnm *group-mrph-list* :key #'idnm)))



;;(DEFCLASS ZG-MORPHISM () ())  ;; will be redefined later.


(DEFCLASS ZG-CHAIN-COMPLEX (chain-complex)
  ((group :type group :initarg :group :reader group1)
   (zgcmpr :type cmprf :initarg :zgcmpr :reader zgcmpr1)
   (zgbasis :type basis :initarg :zgbasis :reader zgbasis1)
   (zgbsgn :type gnrt :initarg :zgbsgn :reader zgbsgn)
   (zgbndr :type ZG-morphism :initarg :zgbndr :reader zgbndr1)
   ))


(DEFVAR *ZGCC-LIST*
    "The variable *ZGCC-LIST* is bound to a list of user created zg-chain complexes.")
(SETF *ZGCC-LIST* +empty-list+)
(PUSHNEW '*ZGCC-LIST* *LIST-LIST*)

#+clisp(eval-when (:compile-toplevel :load-toplevel :execute)
         (setf (ext:package-lock :clos) nil))
(DEFMETHOD PRINT-OBJECT ((zgcm ZG-CHAIN-COMPLEX) stream)
  (the ZG-CHAIN-COMPLEX
    (progn
      (format stream "[K~D ZG-Chain-Complex]" (idnm zgcm))
      zgcm)))
#+clisp(eval-when (:compile-toplevel :load-toplevel :execute)
         (setf (ext:package-lock :clos) t))


(DEFUN ZGCC (idnm)
  (declare (type fixnum idnm))
  (the (or zg-chain-complex null)
    (find idnm *zgcc-list* :key #'idnm)))


;; ZGSTRT = ZGSTRaTegy
(DEFTYPE ZGSTRT () '(member :gnrt :zggnrt :zgcmbn))

;; INTR-ZGMRPH = INTeRnal-ZGMoRPHism
(DEFTYPE INTR-ZGMRPH () 'function)
;; (or (function (degr gnrt) zgcmbn)    ;; if :gnrt strategy
;;     (function (degr zggnrt) zgcmbn)  ;; if :zggnrt strategy
;;     (function (zgcmbn) zgcmbn)     ;; if :zgcmbn strategy

(DEFCLASS ZG-MORPHISM (morphism)
  ((zgintr :type intr-zgmrph :initarg :zg-intr :reader intr-zgmrph)
   (zgstrt :type zgstrt :initarg :zgstrt :reader zgstrt)
   ))



(DEFVAR *zgmrph-list*
    "The variable *ZGMRPH-LIST* is bound to a list of user created morphisms of zg-chain complexes.")
(SETF *zgmrph-list* +empty-list+)
(PUSHNEW '*zgmrph-list* *list-list*)

#+clisp(eval-when (:compile-toplevel :load-toplevel :execute)
         (setf (ext:package-lock :clos) nil))
(DEFMETHOD PRINT-OBJECT ((mrph ZG-MORPHISM) stream)
  (the ZG-MORPHISM
    (progn
      (format stream "[K~D ZG-Morphism (degree ~D): K~D -> K~D]"
        (idnm mrph) (degr mrph)
        (idnm (sorc mrph)) (idnm (trgt mrph)))
      mrph)))
#+clisp(eval-when (:compile-toplevel :load-toplevel :execute)
         (setf (ext:package-lock :clos) t))

(DEFUN ZGMRPH (idnm)
  (declare (type fixnum idnm))
  (the (or zg-morphism null)
    (find idnm *zgmrph-list* :key #'idnm)))


;; MACROS


(DEFMACRO ZGGNRT (gnrt1 gnrt2)
  `(cons  :zggnrt
         (cons ,gnrt1 ,gnrt2)))


(DEFMACRO ZGGNRT1 (ZGgnrt)
  `(car (cdr ,zggnrt)))

(DEFMACRO ZGGNRT2 (ZGgnrt)
  `(cddr ,zggnrt))


(DEFMACRO WITH-ZGGNRT ((zggnrt1 zggnrt2) zggnrt . body)
  `(let (,@(if zggnrt1 `((,zggnrt1 (zggnrt1 ,zggnrt))) nil)
            ,@(if zggnrt2 `((,zggnrt2 (zggnrt2 ,zggnrt))) nil))
     (declare
      (type gnrt ,@(if zggnrt1 `(,zggnrt1) nil) ,@(if zggnrt2 `(,zggnrt2) nil)))
     ,@body))


(DEFMACRO ZGTERM (cffc zggnrt)
  `(cons ,cffc ,zggnrt))


(DEFMACRO WITH-ZGTERM ((cffc zggnrt) zgterm . body)
  `(let (,@(if cffc `((,cffc (cffc ,zgterm))) nil)
            ,@(if zggnrt `((,zggnrt (gnrt ,zgterm))) nil))
     (declare
      (fixnum ,@(if cffc `(,cffc) nil))
      (type zggnrt ,@(if zggnrt `(,zggnrt) nil)))
     ,@body))


(DEFUN BUILD-GROUP
    (&key elements cmpr mult inv nullel orgn)
  (declare
   (type group-basis elements)
   (type cmprf cmpr)
   (type function mult)
   (type function inv)
   (type gnrt nullel)
   (list orgn))
  (the group
    (progn
      (let ((already (find orgn *group-list* :test #'equal :key #'orgn)))
        (declare (type (or group null) already))
        (when already
          (return-from build-group already)))
      (unless elements
        (setf elements :locally-effective))
      (let ((group (make-instance 'group
                     :cmpr cmpr
                     :elements elements
                     :mult mult
                     :inv inv
                     :nullel nullel
                     :orgn orgn)))
        (declare (type group group))
        (push group *group-list*)
        group))))

(DEFUN BUILD-GRMR (&key sorc trgt intr orgn)
  (declare
   (type group sorc trgt)
   (type (or intr-mrph null) intr)
   (list orgn))
  (the group-mrph
    (progn
      (let ((already (find orgn *group-mrph-list* :test #'equal :key #'orgn)))
        (declare (type (or group-mrph null) already))
        (when already
          (return-from build-grmr already)))
      (let ((mrph (make-instance 'group-mrph
                    :sorc sorc :trgt trgt
                    :intr intr
                    :orgn orgn)))
        (declare (type group-mrph mrph))
        (push mrph *group-mrph-list*)
        mrph))))


(DEFMACRO BUILD-AB-GROUP (&rest rest)
  `(change-class (build-group ,@rest) 'ab-group))


(DEFMETHOD ?2 ((grmrph group-mrph) grelt)
  (declare (type gnrt grelt))
  (the gnrt
    (with-slots (intr) grmrph
      (funcall intr grelt))))


(DEFMETHOD ?3 ((group group) g1 g2)
  (declare
   (type gnrt g1 g2))
  (the gnrt
    (with-slots (mult) group
      (funcall mult g1 g2))))


(DEFUN ZGGNRT-CMPR (grcmpr zgcmpr)
  (declare (type cmprf grcmpr zgcmpr))
  (flet ((rslt (gnrt1 gnrt2)
               (declare (type zggnrt gnrt1 gnrt2))
               (the cmpr
                 (let ((gnrt12 (zggnrt2 gnrt1))
                       (gnrt22 (zggnrt2 gnrt2)))
                   (declare (type gnrt gnrt12 gnrt22))
                   (lexico
                    (funcall zgcmpr gnrt12 gnrt22)
                    (let ((gnrt11 (zggnrt1 gnrt1))
                          (gnrt21 (zggnrt1 gnrt2)))
                      (declare (type gnrt gnrt11 gnrt21))
                      (funcall grcmpr gnrt11 gnrt21)))))))
    (the cmprf #'rslt)))


(DEFUN ZGCHCM-BASIS (grbasis zgbasis)
  (declare (type group-basis grbasis)
           (type basis zgbasis))
  (when (or (eq grbasis :locally-effective)
            (eq zgbasis :locally-effective))
    (return-from ZGCHCM-BASIS :locally-effective))
  (flet ((rslt (degr)
               (declare (fixnum degr))
               (the list
                 (progn
                   (when (minusp degr)
                     (return-from rslt +empty-list+))
                   (the list
                     (mapcan
                         #'(lambda (gnrt2)
                             (declare (type gnrt gnrt2))
                             (the list
                               (mapcar
                                   #'(lambda (gnrt1)
                                       (declare (type gnrt gnrt1))
                                       (the zggnrt
                                         (zggnrt gnrt1 gnrt2)))
                                 grbasis)))
                       (funcall zgbasis degr)))))))
    (the basis #'rslt)))


(DEFUN ZGMRPH-GNRT-GNRT (intr degr gnrt)
  (declare
   (type intr-zgmrph intr)
   (fixnum degr)
   (type gnrt gnrt))
  (funcall intr degr gnrt))

(DEFUN ZGMRPH-GNRT-ZGGNRT (grcmpr zgcmpr grmult intr degr zggnrt)
  (declare
   (type cmprf grcmpr zgcmpr)
   (type function grmult)
   (type intr-zgmrph intr)
   (fixnum degr)
   (type zggnrt zggnrt))
  (let* ((zggnrt1 (zggnrt1 zggnrt))
         (zggnrt2 (zggnrt2 zggnrt))
         (zgcmbn (ZGMRPH-GNRT-GNRT intr degr zggnrt2))
         (degr (cmbn-degr zgcmbn))
         (list (cmbn-list zgcmbn))
         (n-cmbn-list +empty-list+))
    (declare (list n-cmbn-list))
    (if list
        (progn
          (mapcar #'(lambda (zgterm)
                      (declare (type zgterm zgterm))
                      (with-zgterm (cffc gnrt) zgterm
                        (declare (type fixnum cffc)
                                 (type zggnrt gnrt))
                        (with-zggnrt (gnrt1 gnrt2) gnrt
                          (declare (type gnrt gnrt1 gnrt2))
                          (push (cons 1
                                      (term-cmbn degr cffc (zggnrt (funcall grmult zggnrt1 gnrt1) gnrt2)))
                                n-cmbn-list))))
            list)
          (cmbn-cmbn (zggnrt-cmpr grcmpr zgcmpr) n-cmbn-list))
      (zero-cmbn degr))))


(DEFUN ZGMRPH-GNRT-ZGCMBN (grcmpr zgcmpr grmult intr zgcmbn)
  (declare
   (type cmprf grcmpr zgcmpr)
   (type function grmult)
   (type intr-zgmrph intr)
   ;;(type zgcmbn zgcmbn)
)
  (the cmbn
    (with-cmbn (degr list) zgcmbn
      (if list
          (let ((n-cmbn-list +empty-list+))
            (declare (list n-cmbn-list))
            (do ((mark list (cdr mark)))
                ((endp mark))
              (declare (list mark))
              (with-term (cffc gnrt) (car mark)
                (push
                 (cons cffc
                       (zgmrph-gnrt-zggnrt grcmpr zgcmpr grmult intr degr gnrt))
                 n-cmbn-list)))
            (cmbn-cmbn (zggnrt-cmpr grcmpr zgcmpr) n-cmbn-list))
        (cmbn degr)))))


(DEFUN ZGMRPH-ZGGNRT-ZGGNRT (intr degr zggnrt)
  (declare
   (type intr-zgmrph intr)
   (fixnum degr)
   (type zggnrt zggnrt))
  (funcall intr degr zggnrt))


(DEFUN ZGMRPH-ZGGNRT-ZGCMBN (grcmpr zgcmpr intr zgcmbn)
  (declare
   (type intr-zgmrph intr)
   (type cmprf grcmpr zgcmpr)
   ;;(type zgcmbn zgcmbn)
)

  (the cmbn
    (with-cmbn (degr list) zgcmbn
      (let ((n-cmbn-list +empty-list+))
        (declare (list n-cmbn-list))
        (do ((mark list (cdr mark)))
            ((endp mark))
          (declare (list mark))
          (with-term (cffc gnrt) (car mark)
            (push
             (cons cffc
                   (zgmrph-zggnrt-zggnrt intr degr gnrt))
             n-cmbn-list)))
        (cmbn-cmbn (zggnrt-cmpr grcmpr zgcmpr) n-cmbn-list)))))


(DEFUN ZGMRPH-GNRT-? (ZGmrph degr gnrt)
  (declare
   (type zg-morphism ZGmrph)
   (fixnum degr)
   (type gnrt gnrt))
  (the cmbn
    (with-slots (zgintr zgstrt trgt) zgmrph
      (declare
       (type zg-chain-complex trgt)
       (type intr-zgmrph zgintr)
       (type zgstrt zgstrt))
      (with-slots (group) trgt
        (declare
         (type group group))

        (with-slots (nullel) group
          (declare
           (type gnrt nullel))
          (prog1
              (ecase zgstrt
                (:gnrt
                 (zgmrph-gnrt-gnrt zgintr degr gnrt))
                (:zggnrt
                 (funcall zgintr degr (zggnrt nullel gnrt)))
                (:cmbn
                 (funcall zgintr (term-cmbn degr 1 (zggnrt nullel gnrt)))))
            ))))))


(DEFUN ZGMRPH-ZGGNRT-? (ZGmrph degr zggnrt)
  (declare
   (type ZG-morphism ZGmrph)
   (fixnum degr)
   (type zggnrt zggnrt))
  (the cmbn
    (with-slots (trgt zgintr zgstrt ) zgmrph
      (declare
       (type zg-chain-complex trgt)
       (type intr-zgmrph zgintr)
       (type zgstrt zgstrt))
      (with-slots (group) trgt
        (declare
         (type group group))
        (with-slots ((grcmpr cmpr) mult nullel) group
          (declare (type cmprf grcmpr))
          (with-slots (zgcmpr)  trgt
            (declare (type cmprf zgcmpr))
            (prog1
                (ecase zgstrt
                  (:gnrt
                   (zgmrph-gnrt-zggnrt grcmpr zgcmpr mult zgintr degr zggnrt))
                  (:zggnrt
                   (funcall zgintr degr zggnrt))
                  (:cmbn
                   (funcall zgintr (term-cmbn degr 1 zggnrt)))))))))))


(DEFUN ZGMRPH-ZGCMBN-? (ZGmrph zgcmbn)
  (declare
   (type zg-morphism zgmrph)
   ;;(type zgcmbn zgcmbn)
)
  (the cmbn
    (with-slots (trgt zgintr zgstrt ) zgmrph
      (declare
       (type zg-chain-complex trgt)
       (type intr-zgmrph zgintr)
       (type zgstrt zgstrt))
      (with-slots (group zgcmpr) trgt
        (declare
         (type group group)
         (type cmprf zgcmpr))
        (with-slots ((grcmpr cmpr) mult nullel) group
          (declare
           (type cmprf grcmpr)
           (type function mult)
           (type gnrt nullel))
          (prog1
              (ecase zgstrt
                (:gnrt
                 (zgmrph-gnrt-zgcmbn grcmpr zgcmpr mult zgintr zgcmbn))
                (:zggnrt
                 (zgmrph-zggnrt-zgcmbn grcmpr zgcmpr zgintr zgcmbn))
                (:cmbn
                 (funcall zgintr zgcmbn)))
            ))))))

(DEFUN IZGMRPH-IMRPH (grcmpr zgcmpr mult intr strt)
  (declare (type intr-zgmrph intr)
           (type zgstrt strt))
  (flet ((rslt (degr zggnrt)
               (declare
                (fixnum degr)
                (type zggnrt zggnrt))
               (ecase strt
                 (:gnrt
                  (zgmrph-gnrt-zggnrt grcmpr zgcmpr mult intr degr zggnrt))
                 (:zggnrt
                  (funcall intr degr zggnrt))
                 (:cmbn
                  (funcall intr (term-cmbn degr 1 zggnrt))))))
    (the intr-mrph #'rslt)))


(DEFUN INTR-ZGCHCM-MRPH (zgmrph)
  (declare (type zg-morphism zgmrph))
  (flet ((rslt (degr zggnrt)
               (declare
                (fixnum degr)
                (type zggnrt zggnrt))
               (zgmrph-zggnrt-? zgmrph degr zggnrt)))
    (the intr-mrph #'rslt)))


(DEFUN BUILD-ZGCHCM
    (&key group zgcmpr zgbasis zgbsgn intr-zgbndr zgbndr-strt orgn)
  (declare
   (type group group)
   (type cmprf zgcmpr)
   (type basis zgbasis)
   (type gnrt zgbsgn)
   (type intr-zgmrph intr-zgbndr)
   (type zgstrt zgbndr-strt)
   (list orgn))
  (let ((already (find orgn *zgcc-list* :test #'equal :key #'orgn)))
    (declare (type (or zg-chain-complex null) already))
    (when already
      (return-from build-zgchcm already)))
  (if intr-zgbndr
      (unless zgbndr-strt
        (error "In BUILD-ZGCHCM, an intr-zgbndr is given but not its strategy."))
    )
  (the zg-chain-complex
    (with-slots ((grcmpr cmpr) (grbasis elements) mult nullel) group
      (let* ((cmpr (zggnrt-cmpr grcmpr zgcmpr))
             (basis (ZGCHCM-BASIS grbasis zgbasis))
             (bsgn (zggnrt nullel zgbsgn))
             (intr-dffr (izgmrph-imrph grcmpr zgcmpr mult intr-zgbndr zgbndr-strt))
             (rslt (build-chcm :cmpr cmpr :basis basis :bsgn bsgn
                               :intr-dffr intr-dffr :strt :gnrt :orgn orgn)))
        (change-class rslt 'zg-chain-complex)
        (setf (slot-value rslt 'group) group)
        (setf (slot-value rslt 'zgcmpr) zgcmpr)
        (setf (slot-value rslt 'zgbasis) zgbasis)
        (setf (slot-value rslt 'zgbsgn) zgbsgn)
        (setf (slot-value rslt 'zgbndr) (build-zgmrph :sorc rslt :trgt rslt :degr -1
                                                      :zgintr intr-zgbndr :zgstrt zgbndr-strt :orgn `(zgbndr ,rslt)))
        (push rslt *zgcc-list*)
        rslt))))


(DEFUN BUILD-ZGMRPH (&key sorc trgt degr zgintr zgstrt orgn)
  (declare
   (type zg-chain-complex sorc trgt)
   (fixnum degr)
   (type intr-zgmrph zgintr)
   (type zgstrt zgstrt)
   (list orgn))
  (the zg-morphism
    (with-slots (group zgcmpr) trgt
      (with-slots ((grcmpr cmpr) (grbasis elements) mult nullel ) group
        (let* ((intr (izgmrph-imrph grcmpr zgcmpr mult zgintr zgstrt))
               (rslt (build-mrph :sorc sorc :trgt trgt :degr degr
                                 :intr intr :strt :gnrt :orgn orgn)))
          (change-class rslt 'zg-morphism)
          (setf (slot-value rslt 'zgintr) zgintr)
          (setf (slot-value rslt 'zgstrt) zgstrt)
          (push rslt *zgcc-list*)
          rslt)))))


;; TENSOR WITH INTEGERS



(DEFUN ZGCMBN-TWI (zgcmpr zgcmbn)
  ;;(declare (type zgcmbn zgcmbn))
  (with-cmbn (degr list) zgcmbn
    (if list
        (let* ((cmbn-list +empty-list+))
          (mapcar #'(lambda (zgterm)
                      (with-zgterm (cffc zggnrt) zgterm
                        (let ((gnrt2 (zggnrt2 zggnrt)))
                          (push (cons 1
                                      (term-cmbn degr cffc gnrt2))
                                cmbn-list))))
            list)
          (cmbn-cmbn zgcmpr cmbn-list))
      (zero-cmbn degr))))


(DEFUN ZGMRPH-TWI-INTR (mrph)
  (declare (type zg-morphism mrph))
  (flet ((rslt (degr gnrt)
               (declare
                (fixnum degr)
                (type gnrt gnrt))
               (with-slots (trgt) mrph
                 (with-slots (zgcmpr) trgt
                   (zgcmbn-twi zgcmpr (zgmrph-gnrt-? mrph degr gnrt))))))
    (the intr-mrph #'rslt)))

(DEFUN ZGMRPH-TWI (mrph)
  (declare (type zg-morphism mrph))
  (the morphism
    (with-slots (sorc trgt degr zgintr zgstrt orgn) mrph
      (build-mrph :sorc (chcm-twi sorc) :trgt (chcm-twi trgt) :degr degr
                  :intr (zgmrph-twi-intr mrph) :strt :gnrt :orgn `(tensor-with-integers ,orgn))
      )))

(DEFUN CHCM-TWI (chcm)
  (declare (type zg-chain-complex chcm))
  (the chain-complex
    (if (string= (car (orgn chcm)) 'bar-zgchcm)
        (with-slots (group) chcm
          (k-g-1 group))
      (with-slots (zgcmpr zgbasis zgbsgn zgbndr orgn) chcm
        (build-chcm :cmpr zgcmpr :basis zgbasis :bsgn zgbsgn :intr-dffr (zgmrph-twi-intr zgbndr)
                    :strt :gnrt :orgn `(tensor-with-integers ,orgn))))))

(DEFUN CYCLICGROUP (n)
  (declare (fixnum n))
  (build-ab-group :elements (<a-b> 0 (1- n)) :cmpr #'f-cmpr :mult #'(lambda (g1 g2)
                                                                      (mod (+ g1 g2) n))
                  :inv #'(lambda (g)
                           (mod (- n g) n)) :nullel 0 :orgn `(Cyclic-group of order ,n)))


(DEFUN Z-GROUP ()
  (build-ab-group :elements :locally-effective :cmpr #'f-cmpr :mult #'(lambda (g1 g2)
                                                                        (+ g1 g2) )
                  :inv #'(lambda (g)
                           (- g)) :nullel 0 :orgn `(z-group)))


(DEFMETHOD SLOT-UNBOUND (class (group group) (slot-name (eql 'resolution)))
  (declare (ignore class))
  (the reduction
    (let ((rsltn (search-resolution group (first (orgn group)))))
      (setf (slot-value group 'resolution)
        (or rsltn
            (error "I don't know how to determine ~
                         a resolution for
                         ~A (Origin: ~A)." group (orgn group)))))))

(DEFGENERIC SEARCH-RESOLUTION (group orgn))

(DEFMETHOD SEARCH-RESOLUTION (group orgn)
  (declare (ignore group orgn))
  nil)


(DEFUN CYCLICGROUP-ZGCHCM (n)
  (let* ((orgn `(Cyclic-group of order ,n))
         (already (find orgn *group-list* :test #'equal :key #'orgn))
         (group (or already (cyclicgroup n)))
         (zgbasis #'(lambda (degr) (list degr)))
         (zgcmpr #'f-cmpr)
         (intr #'(lambda(degr gnrt)
                   (let ((degr- (1- degr)))
                     (if (and gnrt (zerop degr)) (zero-cmbn -1)
                       (if (oddp degr)
                           (cmbn degr- -1 (zggnrt 0 degr-)  1 (zggnrt 1 degr-))
                         (let ((list +empty-list+))
                           (mapcar #'(lambda (i)
                                       (push (cons 1 (zggnrt (- n i) degr-)) list))
                             (<a-b> 1  n))
                           (make-cmbn :degr degr- :list list))))))))
    (BUILD-ZGCHCM
     :group group :zgcmpr zgcmpr :zgbasis zgbasis :zgbsgn 0
     :intr-zgbndr intr :zgbndr-strt :gnrt :orgn  `(ZG-Chain complex for ,group))))


(DEFUN  CYCLICGROUP-RSLTN (n)
  (let* ((orgn `(Cyclic-group of order ,n))
         (already (find orgn *group-list* :test #'equal :key #'orgn))
         (group (or already (cyclicgroup n)))
         (zgchcm (cyclicgroup-zgchcm n))
         (f-intr #'(lambda (degr zggnrt)
                     (with-zggnrt (gnrt1 gnrt2) zggnrt
                       (declare (ignore gnrt1 gnrt2))
                       (if (zerop degr) (cmbn degr 1 :zgnrt) (zero-cmbn degr)))))
         (g-intr #'(lambda (degr gnrt)
                     (declare (ignore gnrt))
                     (if (zerop degr)  (cmbn degr 1 (zggnrt 0 0)) (zero-cmbn degr))))
         (h-intr #'(lambda (degr zggnrt)
                     (with-zggnrt (gnrt1 gnrt2) zggnrt
                       (declare (ignore gnrt2))
                       (if (oddp degr)
                           (if (= gnrt1 (1- n)) (cmbn (1+ degr) 1 (zggnrt 0 (1+ degr)))
                             (zero-cmbn (1+ degr)))
                         (let ((list +empty-list+))
                           (dotimes (i gnrt1)
                             (push (cons 1 (zggnrt i (1+ degr))) list ))
                           (make-cmbn :degr (1+ degr) :list (nreverse list)))))))
         (f (build-mrph :sorc zgchcm :trgt (z-chcm) :degr 0 :intr f-intr :strt :gnrt :orgn `(f-mrph for ,group)))
         (g (build-mrph :sorc (z-chcm) :trgt zgchcm :degr 0 :intr g-intr :strt :gnrt :orgn  `(g-mrph for ,group)))
         (h (build-mrph :sorc zgchcm :trgt zgchcm :degr 1 :intr h-intr :strt :gnrt :orgn `(h-mrph for ,group))))
    (build-rdct :f f :g g :h h :orgn `(Resolution of ,group))))


(DEFMETHOD SEARCH-RESOLUTION (group (orgn (eql 'Cyclic-group)))
  (declare
   (type group group))
  (the reduction
    (with-slots (elements) group
      (let ((n (length elements)))
        (cyclicgroup-rsltn n)))))


(DEFMETHOD GROUP-HOMOLOGY ((group group) degr1 &optional (degr2 (1+ degr1)))
  (declare (fixnum degr1 degr2))
  (with-slots (resolution) group
    (let ((chcm (chcm-twi (tcc resolution))))
      (homology chcm degr1 degr2))))


;;;; BAR RESOLUTION


(DEFUN BAR-ZGCMPR (grcmpr)
  (declare (type cmprf grcmpr))
  (labels ((rslt (gnrt1 gnrt2)
                 (declare (list gnrt1 gnrt2))
                 (unless gnrt1
                   (return-from rslt (if gnrt2 :less :equal)))
                 (unless gnrt2
                   (return-from rslt :greater))
                 (lexico
                  (funcall grcmpr (first gnrt1) (first gnrt2))
                  (rslt (rest gnrt1) (rest gnrt2)))
                 ))
    (the cmprf #'rslt)))


(DEFUN BAR-ZGBASIS (elements nullel)
  (declare (type group-basis elements)
           (type gnrt nullel))
  (when (eq elements :locally-effective)
    (return-from bar-zgbasis :locally-effective))
  (labels ((rslt (degr)
                 (declare (fixnum degr))
                 ;;(the list
                 (progn
                   (when (minusp degr)
                     (return-from rslt +empty-list+))
                   (when (zerop degr)
                     (return-from rslt (list (list))))
                   (let ((basis-1 (rslt (1- degr)))
                         (elements-0 (remove nullel elements)))
                     (mapcan #'(lambda (gbar)
                                 (mapcar #'(lambda (g)
                                             (append gbar (list g)))
                                   elements-0))
                       basis-1)))))
    ;;)
    (the basis #'rslt)))


(DEFUN BAR-BNDR-I (grmult nullel indx degr gnrt)
  (declare
   (fixnum indx degr)
   (list gnrt))
  (the cmbn
    (cond ((zerop indx) (cmbn (1- degr) 1 (zggnrt (car gnrt) (rest gnrt))))
          ((= indx degr) (cmbn (1- degr) 1 (zggnrt nullel (butlast gnrt))))
          (t
           (do ((rslt +empty-list+ (cons (car mark) rslt))
                (mark gnrt (cdr mark))
                (i 1 (1+ i)))
               ((= i indx) (let ((new-k (funcall grmult (car mark) (cadr mark))))
                             (declare (type gnrt new-k))
                             (if (equal nullel new-k)
                                 (zero-cmbn (1- degr))
                               (cmbn (1- degr) 1 (zggnrt nullel (nreconc rslt (cons new-k (cddr mark))))))))
             (declare
              (list rslt mark)
              (fixnum i)))))))


(DEFUN BAR-BNDR-INTR (grcmpr grmult nullel)
  #'(lambda (degr gnrt)
      (if (zerop degr)
          (zero-cmbn -1)
        (let ((cmbn-list
               (mapcar #'(lambda (indx)
                           (let ((cmbn (bar-bndr-i grmult nullel indx degr gnrt)))
                             (if (oddp indx)
                                 (cons -1 cmbn)
                               (cons 1 cmbn))))
                 (<a-b> 0 degr))))
          (cmbn-cmbn (zggnrt-cmpr grcmpr (bar-zgcmpr grcmpr)) cmbn-list)))))


(DEFUN BAR-ZGCHCM (group)
  (declare (type group group))
  (the zg-chain-complex
    (with-slots (elements cmpr mult nullel) group
      (build-zgchcm
       :group group
       :zgcmpr (bar-zgcmpr cmpr)
       :zgbasis (bar-zgbasis elements nullel)
       :zgbsgn +empty-list+
       :intr-zgbndr (bar-bndr-intr cmpr mult nullel)
       :zgbndr-strt :gnrt
       :orgn  `(bar-zgchcm ,group)))))


(DEFUN BAR-CNTRHMTP-INTR (grcmpr grmult nullel)
  #'(lambda (degr zggnrt)
      (declare (ignore grcmpr grmult))
      (with-zggnrt (gnrt1 gnrt2) zggnrt
        (if (equal nullel gnrt1)
            (zero-cmbn (1+ degr))
          (if (zerop degr)
              (cmbn 1 1 (zggnrt nullel (list gnrt1)))
            (cmbn (1+ degr) 1 (zggnrt nullel (cons gnrt1 gnrt2))))))))


(DEFUN BAR-F-INTR ()
  #'(lambda (degr zggnrt)
      (declare (ignore zggnrt))
      (if (zerop degr) (cmbn degr 1 :zgnrt) (zero-cmbn degr))))


(DEFUN BAR-G-INTR (nullel)
  #'(lambda (degr gnrt)
      (declare (ignore gnrt))
      (if (zerop degr) (cmbn degr 1 (zggnrt nullel +empty-list+)) (zero-cmbn degr))))

(DEFUN BAR-RSLTN (group)
  (declare (type group group))
  (the reduction
    (with-slots (cmpr mult nullel) group
      (let* ((zgchcm (bar-zgchcm group))
             (f (build-mrph :sorc zgchcm :trgt (z-chcm) :degr 0 :intr (bar-f-intr) :strt :gnrt :orgn `(f-mrph in bar-resolution for ,group)))
             (g (build-mrph :sorc (z-chcm) :trgt zgchcm :degr 0 :intr (bar-g-intr nullel) :strt :gnrt :orgn  `(g-mrph in bar-resolution for ,group)))
             (h (build-mrph :sorc zgchcm :trgt zgchcm :degr 1 :intr (bar-cntrhmtp-intr cmpr mult nullel)
                            :strt :gnrt :orgn `(h-mrph in bar-resolution for ,group))))
        (build-rdct :f f :g g :h h :orgn `(Bar-resolution of ,group))))))


(DEFUN K-G-1-CMPR (grcmpr)
  (declare (type cmprf grcmpr))
  (labels ((rslt (gnrt1 gnrt2)
                 (declare (list gnrt1 gnrt2))
                 (unless gnrt1
                   (return-from rslt (if gnrt2 :less :equal)))
                 (unless gnrt2
                   (return-from rslt :greater))
                 (lexico
                  (funcall grcmpr (first gnrt1) (first gnrt2))
                  (rslt (rest gnrt1) (rest gnrt2)))))
    (the cmprf #'rslt)))


(DEFUN K-G-1-BASIS (elements nullel)
  (declare (type group-basis elements)
           (type gnrt nullel))
  (when (eq elements :locally-effective)
    (return-from K-G-1-basis :locally-effective))
  (labels ((rslt (degr)
                 (declare (fixnum degr))
                 ;;(the list
                 (progn
                   (when (minusp degr)
                     (return-from rslt +empty-list+))
                   (when (zerop degr)
                     (return-from rslt (list (list))))
                   (let ((basis-1 (rslt (1- degr)))
                         (elements-0

                          (do ((indx 0 (1+ indx))
                               (aux +empty-list+))
                              ((>= indx (length elements)) (nreverse aux))
                            (if (not (equal nullel (nth indx elements)))
                                (push (nth indx elements) aux)))))
                     (mapcan #'(lambda (gbar)
                                 (mapcar #'(lambda (g)
                                             (append gbar (list g)))
                                   elements-0))
                       basis-1)))))
    ;;)
    (the basis #'rslt)))



(DEFUN K-G-1-FACE (grmult nullel)
  (flet ((face (indx dmns gmsm)
               (declare
                (fixnum indx dmns)
                (list gmsm))
               (the absm
                 (cond ((zerop indx) (absm 0 (rest gmsm)))
                       ((= indx dmns) (absm 0 (butlast gmsm)))
                       (t
                        (do ((rslt +empty-list+ (cons (car mark) rslt))
                             (mark gmsm (cdr mark))
                             (i 1 (1+ i)))
                            ((= i indx) (let ((new-k (funcall grmult (car mark) (cadr mark))))
                                          (declare (type gnrt new-k))
                                          (if (equal nullel new-k)
                                              (absm (2-exp (1- indx)) (nreconc rslt (cddr mark)))
                                            (absm 0 (nreconc rslt (cons new-k (cddr mark)))))))
                          (declare
                           (list rslt mark)
                           (fixnum i))))))))
    (the face #'face)))


(DEFUN G-ABSM-BAR (nullel absm)
  (declare (type absm absm))
  (the list
    (with-absm (dgop bar1) absm
      (do ((dgop dgop (ash dgop -1))
           (rslt +empty-list+))
          ((and (zerop dgop) (endp bar1)) (nreverse rslt))
        (declare
         (fixnum dgop)
         (list rslt))
        (push (if (oddp dgop)
                  nullel
                (pop bar1))
              rslt)))))


(DEFUN G-BAR-ABSM (nullel bar)
  (declare (list bar))
  (the absm
    (do ((rslt-dgop 0)
         (rslt-bar +empty-list+)
         (mark bar (cdr mark))
         (bark 1 (+ bark bark)))
        ((endp mark) (absm rslt-dgop (nreverse rslt-bar)))
      (declare
       (fixnum rslt-dgop bark)
       (list rslt-bar mark))
      (let ((k-i (car mark)))
        (declare (fixnum k-i))
        (if (equal nullel k-i)
            (incf rslt-dgop bark)
          (push k-i rslt-bar))))))


(DEFUN K-G-1-GRML (grmult nullel)
  (flet ((grml (dmns crpr)
               (declare
                (fixnum dmns)
                (type crpr crpr))
               (the absm
                 (with-crpr (dgop1 bar1 dgop2 bar2) crpr
                   (do ((indx 1 (1+ indx))
                        (dgop1 dgop1 (ash dgop1 -1))
                        (dgop2 dgop2 (ash dgop2 -1))
                        (bark 1 (ash bark +1))
                        (rslt-dgop 0)
                        (rslt-bar +empty-list+))
                       ((> indx dmns) (absm rslt-dgop (nreverse rslt-bar)))
                     (declare
                      (fixnum indx dgop1 dgop2 bark rslt-dgop)
                      (list rslt-bar))
                     (let ((item (if (evenp dgop1)
                                     (if (evenp dgop2)
                                         (funcall grmult (pop bar1) (pop bar2))
                                       (pop bar1))
                                   (if (evenp dgop2)
                                       (pop bar2)
                                     nullel))))
                       (declare (type gnrt item))
                       (if (equal nullel item)
                           (incf rslt-dgop bark)
                         (push item rslt-bar))))))))
    (the intr-mrph #'grml)))


(DEFUN K-G-1-GRIN (grin)
  (flet ((grin (dmns bar)
               (declare
                (ignore dmns)
                (list bar))
               (the absm
                 (absm 0 (mapcar #'(lambda (g) (funcall grin g)) bar)

                       ))))
    (the intr-mrph #'grin)))


(DEFMETHOD K-G-1 ((group GROUP))
  (declare (type group group))
  (the simplicial-group
    (if (string= (car (orgn group)) 'z-group)
        (k-z-1)

      (with-slots (elements cmpr mult nullel inv) group
        (build-smst
         :cmpr (k-g-1-cmpr cmpr)
         :basis (k-g-1-basis elements nullel)
         :bspn +empty-list+
         :face (k-g-1-face mult nullel)
         :orgn  `(k-g-1 ,group)))
      )))

(DEFMETHOD K-G-1 ((group AB-GROUP))
  (declare (type AB-group group))
  (the AB-simplicial-group
    (if (string= (car (orgn group)) 'z-group)
        (k-z-1)
      (with-slots (elements cmpr mult nullel inv) group
        (build-AB-smgr
         :cmpr (k-g-1-cmpr cmpr)
         :basis (k-g-1-basis elements nullel)
         :bspn +empty-list+
         :face (k-g-1-face mult nullel)
         :sintr-grml (k-g-1-grml mult nullel)
         :sintr-grin (k-g-1-grin inv)
         :orgn  `(k-g-1 ,group)))
      )))


(defun K-G-0-BASIS (elements)
  (declare (type group-basis elements))
  (when (eq elements :locally-effective)
    (return-from K-G-0-basis :locally-effective))
  (labels ((rslt (degr)
                 (declare (fixnum degr))
                 ;;(the list
                 (progn
                   (when (not (equal degr 0))
                     (return-from rslt +empty-list+))
                   (return-from rslt elements)
                   )))
    ;;)
    (the basis #'rslt)))


(DEFUN K-G-0-FACE ()
  (flet ((face (indx dmns gmsm)
               (declare
                (ignore indx dmns)
                (gnrt gmsm))
               (the absm
                 (absm 0 gmsm))))
    (the face #'face)))

(DEFUN K-G-0-GRML (grmult)
  (flet ((grml (dmns crpr)
               (declare
                (fixnum dmns)
                (type crpr crpr))
               (the absm
                 (with-crpr (dgop1 g1 dgop2 g2) crpr
                   (declare (ignore dgop1 dgop2))
                   (absm (dgop-ext-int (nreverse (<a-b> 0 (1- dmns)))) (funcall grmult g1 g2)))
                 )))
    (the intr-mrph #'grml)))


(DEFUN K-G-0-GRIN (grin)
  (flet ((grin (dmns g)
               (declare
                (ignore dmns)
                (gnrt g))
               (the absm
                 (absm 0 (funcall grin g)))))
    (the intr-mrph #'grin)))


(DEFMETHOD K-G-0 ((group GROUP))
  (declare (type group group))
  (the simplicial-group
    (with-slots (elements cmpr mult nullel inv) group
      (build-smgr
       :cmpr cmpr
       :basis (k-g-0-basis elements)
       :bspn nullel
       :face (k-g-0-face)
       :sintr-grml (k-g-0-grml mult)
       :sintr-grin (k-g-0-grin inv)
       :orgn  `(k-g-0 ,group)))
    ))

(DEFMETHOD K-G-0 ((group AB-GROUP))
  (declare (type AB-group group))
  (the AB-simplicial-group
    (with-slots (elements cmpr mult nullel inv) group
      (build-AB-smgr
       :cmpr cmpr
       :basis (k-g-0-basis elements)
       :bspn nullel
       :face (k-g-0-face)
       :sintr-grml (k-g-0-grml mult)
       :sintr-grin (k-g-0-grin inv)
       :orgn  `(k-g-0 ,group)))))


(DEFMETHOD K-G ((group ab-group) n)
  (declare (fixnum n)
           (type AB-group group))
  (the ab-simplicial-group
    (if (= n 0)
        (k-g-0 group)
      (if (= n 1)
          (k-g-1 group)
        (classifying-space (k-g group (1- n)))))))


(DEFUN K-ZP-1 (p)
  (k-g-1 (cyclicgroup p)))

(DEFUN K-ZP (p n)
  (declare (fixnum p n)
           )
  (the ab-simplicial-group
    (if (= n 0)
        (k-g-0 (cyclicgroup p))
      (if (= n 1)
          (k-Zp-1 p)
        (classifying-space (k-Zp p (1- n)))))))


(DEFUN UNIV-FBRT-INTR (G n)
  (declare
   (ignore G)
   (fixnum n))
  (flet ((rslt (degr gmsm)
               (declare
                (fixnum degr)
                (type gnrt gmsm))
               (if (equal n 1)
                   (absm (dgop-ext-int (nreverse (<a-b> 0 (- degr 2)))) (first (last gmsm)))
                 (with-gbar (d list) gmsm
                   (declare (ignore d))
                   (first list)))))
    (the sintr #'rslt)))


(DEFUN UNIV-FBRT-TW (G n)
  (the fibration
    (build-smmr
     :sorc (k-g g n) :trgt (k-g g (1- n)) :degr -1
     :sintr (univ-fbrt-intr G n)
     :orgn `(universal fibration ,G ,n))))


(DEFUN RSLTN-ZGCHCM (rsltn)
  (declare (type reduction rsltn))
  (tcc rsltn))


(DEFUN RSLTN-CNTRHMTP (rsltn)
  (declare (type reduction rsltn))
  (h rsltn))


(DEFUN 2RSLTN-ZGMRPH-INTR (rsltn1 rsltn2)
  (let ((zgchcm2 (rsltn-zgchcm rsltn2))
        (h2 (rsltn-cntrhmtp rsltn2))
        (f1 (f rsltn1))
        (g2 (g rsltn2))
        (d1 (zgbndr1 (rsltn-zgchcm rsltn1))))
    (with-slots (group (zgcmpr2 zgcmpr)) zgchcm2
      (with-slots ((grcmpr cmpr) mult nullel) group
        (labels ((rslt (degr gnrt)
                       (if (zerop degr)
                           (? g2 (? f1 degr (zggnrt nullel gnrt)))
                         (? h2

                            (with-cmbn (degr- list) (zgmrph-gnrt-? d1 degr gnrt)
                              (let ((n-cmbn-list +empty-list+))
                                (declare (list n-cmbn-list))
                                (do ((mark list (cdr mark))
                                     )
                                    ((endp mark))
                                  (declare (list mark))
                                  (with-term (cffc ZGgnrt) (car mark)
                                    (push
                                     (cons cffc

                                           (let* ((gnrt1 (zggnrt1 zggnrt))
                                                  (gnrt2 (zggnrt2 zggnrt))
                                                  (zgcmbn (rslt degr- gnrt2))
                                                  (list (cmbn-list zgcmbn))
                                                  (n-cmbn-list2 +empty-list+))
                                             (declare (list n-cmbn-list2))
                                             (if list
                                                 (progn
                                                   (mapcar #'(lambda (zgterm)
                                                               (declare (type zgterm zgterm))
                                                               (with-zgterm (cffc2 zggnrt2) zgterm
                                                                 (declare (type fixnum cffc2)
                                                                          (type zggnrt zggnrt2))
                                                                 (with-zggnrt (gnrt21 gnrt22) zggnrt2
                                                                   (declare (type gnrt gnrt21 gnrt22))
                                                                   (push (cons 1
                                                                               (term-cmbn degr- cffc2 (zggnrt (funcall mult gnrt1 gnrt21) gnrt22)))
                                                                         n-cmbn-list2))))
                                                     list)
                                                   (cmbn-cmbn (zggnrt-cmpr grcmpr zgcmpr2) n-cmbn-list2))
                                               (zero-cmbn degr-))))
                                     n-cmbn-list)))
                                (cmbn-cmbn (zggnrt-cmpr grcmpr zgcmpr2) n-cmbn-list)))))))
          #'rslt)))))


(DEFUN 2RSLTN-ZGMRPH (rsltn1 rsltn2)
  (let ((zgchcm1 (rsltn-zgchcm rsltn1))
        (zgchcm2 (rsltn-zgchcm rsltn2)))
    (build-zgmrph :sorc zgchcm1 :trgt zgchcm2 :degr 0 :zgintr (2rsltn-zgmrph-intr rsltn1 rsltn2)
                  :zgstrt :gnrt :orgn `(ZG-morphism from ,rsltn1 to ,rsltn2))))


(DEFUN 2RSLTN-HMTPOP-INTR (rsltn1 rsltn2)
  (let ((zgchcm1 (rsltn-zgchcm rsltn1))
        (h1 (rsltn-cntrhmtp rsltn1))
        (m (2rsltn-zgmrph rsltn1 rsltn2))
        (n (2rsltn-zgmrph rsltn2 rsltn1))
        (d1 (zgbndr1 (rsltn-zgchcm rsltn1))))
    (with-slots (group (zgcmpr1 zgcmpr)) zgchcm1
      (with-slots ((grcmpr cmpr) mult nullel) group
        (let ((cmpr (zggnrt-cmpr grcmpr zgcmpr1)))

          (labels ((rslt (degr gnrt)
                         (if (zerop degr)
                             (2cmbn-sbtr cmpr (? h1 degr (zggnrt nullel gnrt))
                                         (? h1 (zgmrph-zgcmbn-? n (zgmrph-gnrt-? m degr gnrt))))
                           (2cmbn-sbtr cmpr (2cmbn-sbtr cmpr (? h1 degr (zggnrt nullel gnrt))
                                                        (? h1 (zgmrph-zgcmbn-? n (zgmrph-gnrt-? m degr gnrt))))
                                       (? h1
                                          (with-cmbn (degr- list) (zgmrph-gnrt-? d1 degr gnrt)
                                            (let ((n-cmbn-list +empty-list+))
                                              (declare (list n-cmbn-list))
                                              (do ((mark list (cdr mark))
                                                   )
                                                  ((endp mark))
                                                (declare (list mark))
                                                (with-term (cffc ZGgnrt) (car mark)
                                                  (push
                                                   (cons cffc

                                                         (let* ((gnrt1 (zggnrt1 zggnrt))
                                                                (gnrt2 (zggnrt2 zggnrt))
                                                                (zgcmbn (rslt degr- gnrt2))
                                                                (list (cmbn-list zgcmbn))
                                                                (n-cmbn-list2 +empty-list+))
                                                           (declare (list n-cmbn-list2))
                                                           (if list
                                                               (progn
                                                                 (mapcar #'(lambda (zgterm)
                                                                             (declare (type zgterm zgterm))
                                                                             (with-zgterm (cffc2 zggnrt2) zgterm
                                                                               (declare (type fixnum cffc2)
                                                                                        (type zggnrt zggnrt2))
                                                                               (with-zggnrt (gnrt21 gnrt22) zggnrt2
                                                                                 (declare (type gnrt gnrt21 gnrt22))
                                                                                 (push (cons 1
                                                                                             (term-cmbn degr cffc2 (zggnrt (funcall mult gnrt1 gnrt21) gnrt22)))
                                                                                       n-cmbn-list2))))
                                                                   list)
                                                                 (cmbn-cmbn cmpr n-cmbn-list2))
                                                             (zero-cmbn degr))))
                                                   n-cmbn-list)))
                                              (cmbn-cmbn cmpr n-cmbn-list))))))))
            #'rslt))))))

(DEFUN 2RSLTN-HMTPOP (rsltn1 rsltn2)
  (let ((zgchcm1 (rsltn-zgchcm rsltn1)))
    (build-zgmrph :sorc zgchcm1 :trgt zgchcm1 :degr 1 :zgintr (2rsltn-hmtpop-intr rsltn1 rsltn2)
                  :zgstrt :gnrt :orgn `(ZG-homotopy-operator from ,rsltn1 to ,rsltn2))))


(DEFUN K-CM-N (m n)
  (K-Zp m n))

