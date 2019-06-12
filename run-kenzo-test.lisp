;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp; Base: 10 -*

(DECLAIM (OPTIMIZE (speed 3) (space 0) (debug 3)))

(require :asdf)
(require :kenzo)
(require :ecl-quicklisp)
(ql:quickload "fiveam")
(push #P"./" asdf:*central-registry*)
(require :kenzo-test)
(setf fiveam:*debug-on-error* t
      fiveam:*debug-on-failure* t)
(setf *debugger-hook*
      (lambda (c h)
              (declare (ignore c h))
              (uiop:quit -1)))
(fiveam:run! :kenzo)

#+(or ecl ccl) (quit)
#+sbcl(sb-ext:exit)
