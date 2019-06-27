;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*

;;;  FILES  FILES  FILES  FILES  FILES  FILES  FILES  FILES  FILES  FILES
;;;  FILES  FILES  FILES  FILES  FILES  FILES  FILES  FILES  FILES  FILES
;;;  FILES  FILES  FILES  FILES  FILES  FILES  FILES  FILES  FILES  FILES

(DEFPACKAGE :CAT
  (:USE CL))

(IN-PACKAGE "COMMON-LISP-USER")

(DEFUN KENZO-VERSION ()
  (format t "~%*** Kenzo-Version 1.1.7 ***~2%")
  (values))

(PROCLAIM '(OPTIMIZE (speed 3) (safety 1) (space 0) (debug 0)))

(DEFCONSTANT +FILE-LIST+
  '("kenzo"
    "macros"
    "various"
    "classes"
    "combinations"
    "chain-complexes"
    "chcm-elementary-op"
    "effective-homology"
    "homology-groups"
    "searching-homology"
    "cones"
    "bicones"
    "tensor-products"
    "coalgebras"
    "cobar"
    "algebras"
    "bar"
    "simplicial-sets"
    "simplicial-mrphs"
    "delta"
    "special-smsts"
    "suspensions"
    "disk-pasting"
    "cartesian-products"
    "eilenberg-zilber"
    "kan"
    "simplicial-groups"
    "fibrations"
    "loop-spaces"
    "ls-twisted-products"
    "lp-space-efhm"
    "classifying-spaces"
    "k-pi-n"
    "serre"
    "cs-twisted-products"
    "cl-space-efhm"
    "whitehead"
    "smith"
    "anromero\\resolutions"
    "anromero\\cylinders"
    "anromero\\fundamental-classes"
    "anromero\\homotopy"
    "anromero\\bicomplexes"
    "anromero\\filtered-complexes"
    "anromero\\spectral-sequences"
    "anromero\\persistent-homology"
    ))

(DO ((i 1 (1+ i))
     (mark +file-list+ (cdr mark)))
    ((endp mark) (terpri))
  (format t "~% FILE ~2D: ~A" i (car mark)))

(DEFCONSTANT +SOURCE-EXTENSION+
  #+(or lispworks) "cl"
  #+(or allegro clisp ccl ecl sbcl) "lisp"
  #-(or allegro ccl clisp ecl lispworks sbcl)
  (error "Not an Allegro or CCL or CLisp or LispWorks or SBCL environment."))

(DEFCONSTANT +COMPILED-EXTENSION+
  #+(or allegro sbcl) "fasl"
  #+ccl "lafsl"
  #+ (or clisp ecl) "fas"
  #+lispworks "ofasl")

(DEFUN LOAD-SFILES ()
  (mapc #'(lambda (file-name)
            (load (concatenate 'string file-name "." +source-extension+)))
        +file-list+))

(DEFVAR *CMBN-CONTROL*)
(SETF *CMBN-CONTROL* T)

#+SBCL (DECLAIM (SB-EXT:MUFFLE-CONDITIONS style-warning compiler-note))

(DEFUN COMPILE-FILES ()
  (format t "~%*CMBN-CONTROL*  = ~A~2%" *cmbn-control*)
  (mapc #'(lambda (file-name)
            (compile-file (concatenate 'string file-name "." +source-extension+))
            (load (concatenate 'string file-name "." +compiled-extension+)))
        +file-list+))

(DEFUN LOAD-CFILES ()
  (mapc #'(lambda (file-name)
            (load (concatenate 'string file-name "." +compiled-extension+)))
        +file-list+))
