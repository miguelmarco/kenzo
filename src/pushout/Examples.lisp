;;;============================================================================
;;;
;;; Subject: Examples to test the pushout
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
;; August 2019
;;; Most of the examples have been collected from the rest of the development 
;;; and unified in this file.

;; ====================================
;; Examples cones2.lisp
;; ====================================

#|
()
(setf basis (cone-basis2
             #'(lambda (n) (list n))
             #'(lambda (n) (list n))))
(funcall basis 4)
|#

#|
()
(cone-cmbn-split2
 (cmbn 3 4 (con0 'a) 5 (con1 'b)))
(cone-cmbn-split2 (cmbn 3 4 (con0 'a)))
(cone-cmbn-split2 (cmbn 3 4 (con1 'a)))
|#

#|
()
(cmbn-con02 (cmbn 3 1 'a 2 'b))
|#

#|
()
(setf mrph (dffr (delta 4)))
(setf rslt (cone-2mrph-diag-impl2 mrph mrph))
(funcall rslt (cmbn 3 4 (con0 15) 5 (con1 7)))
|#


#|
()
(setf rslt (cone-3mrph-triangle-impl2 #'f-cmpr
                                     (dffr (delta 4))
                                     (dffr (delta 4))
                                     (idnt-mrph (delta 4))))
(funcall rslt 4 (con0 '31))
(funcall rslt 5 (con1 '31))
(funcall rslt 4 (con1 '15))
; (funcall rslt (cmbn 4 1 (con0 '31) 1 (con1 '15)))
|#


#|
()
(cat-init)
(setf k (k-z 1))
(setf u (idnt-mrph k))
(setf c (cone2 u))
(setf *tc* (cmbn 3 1 (con0 '(1 2 3 4))))
(? c 3 (con0 '(1 2 3 4)))
(? c *tc*)
|#


#|
()
(setf idnt (idnt-mrph (delta 4)))
(setf cone (cone2 idnt))
(setf ff (cone-2mrph-diag2 cone cone idnt idnt))
(? ff (cmbn 4 1 (con0 '31) 10 (con1 '15)))
|#

#|
()
(cat-init)
(setf d (delta 4))
(setf df (dffr (delta 4)))
(setf n (idnt-mrph d))
(setf c (cone n))
(setf z (sbtr (dffr c)
              (cone-3mrph-triangle c c df (n-mrph -1 df) n)))
(? z 0 (con0 '1))
(? z 1 (con1 '1))
(? z 1 (con0 '3))
(? z 2 (con1 '3))
(? z 4 (con0 '31))
(? z 5 (con1 '31))
|#


;; ====================================
;; Examples direct-sum.lisp
;; ====================================


#|
()
(direct-sum-gsm 0 'a)
(direct-sum-gsm 1 'b)
|#

#|
()
(setf cmpr (direct-sum-cmpr #'s-cmpr #'f-cmpr))
(funcall cmpr (direct-sum-gsm 0 'a) (direct-sum-gsm 0 'b))
(funcall cmpr (direct-sum-gsm 0 'b) (direct-sum-gsm 0 'b))
(funcall cmpr (direct-sum-gsm 0 'b) (direct-sum-gsm 0 'a))
(funcall cmpr (direct-sum-gsm 0 'b) (direct-sum-gsm 1 '2))
(funcall cmpr (direct-sum-gsm 1 '2) (direct-sum-gsm 0 'b))
(funcall cmpr (direct-sum-gsm 1 '1) (direct-sum-gsm 1 '2))
(funcall cmpr (direct-sum-gsm 1 '2) (direct-sum-gsm 1 '2))
(funcall cmpr (direct-sum-gsm 1 '2) (direct-sum-gsm 1 '1))
|#


#|
(setf basis (direct-sum-basis
             #'(lambda (n) (list n))
             #'(lambda (n) (list n))))
(funcall basis 4)
|#


#|
()
(TERM-d-sum-0 (term 2 'a))
(TERM-d-sum-1 (term 2 'a))
|#

#|
()
(direct-sum-cmbn-split
 (cmbn 3 4 (direct-sum-gsm 0 'a) 5 (direct-sum-gsm 1 'b)))
(cone-cmbn-split (cmbn 3 4 (direct-sum-gsm 0 'a)))
(cone-cmbn-split (cmbn 3 4 (direct-sum-gsm 1 'a)))
|#

#|
()
(direct-sum-2cmbn-append (cmbn 3 4 'a) (cmbn 2 5 'b))
|#


#|
()
(setf mrph (dffr (delta 4)))
(setf rslt (direct-sum-dffr mrph mrph))
(funcall rslt (cmbn 3 4 (direct-sum-gsm 0 15) 5 (direct-sum-gsm 1 7)))
|#

#|
()
(setf s3 (sphere 3))
(setf s3+s3 (direct-sum s3 s3))

|#


#|

(cat-init)
(setf s3xs3 (crts-prdc (sphere 3) (sphere 3)))
(setf s3xs3+s3xs3-efhm (direct-sum-efhm s3xs3 s3xs3))
(setf lrdct (lrdct s3xs3+s3xs3-efhm))

(pre-check-rdct lrdct)
(setf *tc* (cmbn 3 1 (first (basis (k 40) 3)) -2 (second (basis (k 40) 3)) 3 (third (basis (k 40) 3))))
(setf *bc* (cmbn 4 -1 (first (basis (k 38) 4)) -5 (second (basis (k 38) 4)) 7 (third (basis (k 38) 4)) -4 (fourth (basis (k 38) 4))))
(check-rdct)

(setf rrdct (rrdct s3xs3+s3xs3-efhm))

(pre-check-rdct rrdct)
(setf *tc* (cmbn 2 1 (first (basis (k 40) 2)) -2 (second (basis (k 40) 2)) 3 (third (basis (k 40) 2))))
(setf *bc* (cmbn 3 -1 (first (basis (k 42) 3)) -5 (second (basis (k 42) 3)) 7 (third (basis (k 42) 3)) -4 (fourth (basis (k 42) 3))))
(check-rdct)

|#

#|
()
(setf ls3 (loop-space (sphere 3)))
(setf s3 (sphere 3))
(setf ls3+s3 (direct-sum ls3 s3))

(homology ls3+s3 0 6)

|#



;; ====================================
;; Examples remove-covers.lisp
;; ====================================

#|
(setf s3xI (crts-prdc (sphere 3) (delta 1)))

(setf basis (remove-covers-basis
             (basis s3xI)))
(funcall basis 4)
|#

#|


(setf rcs3 (remove-covers (sphere 3)))

(homology rcs3 0 6)

|#


#|

(cat-init)
(setf s3xs3 (crts-prdc (sphere 3) (sphere 3)))
(setf rcs3xs3-efhm (remove-covers-efhm (remove-covers s3xs3)))
(setf lrdct (lrdct rcs3xs3-efhm))

(pre-check-rdct lrdct)
(setf *tc* (cmbn 3 1 (first (basis (k 117) 3)) -2 (second (basis (k 117) 3)) 4 (third (basis (k 117) 3))))
(setf *bc* (cmbn 4 -1 (first (basis (k 21) 4)) -5 (second (basis (k 21) 4)) 7 (third (basis (k 21) 4))))
(check-rdct)

(setf rrdct (rrdct rcs3xs3-efhm))

(pre-check-rdct rrdct)
(setf *tc* (cmbn 3 1 (first (basis (k 117) 3)) -2 (second (basis (k 117) 3)) 4 (third (basis (k 117) 3))))
(setf *bc* (cmbn 4 -1 (first (basis (k 120) 4)) -5 (second (basis (k 120) 4)) 3 (car (last (basis (k 120) 4)))))
(check-rdct)

|#


#|

(setf rcls3 (remove-covers (loop-space (sphere 3))))

(homology rcls3 0 6)

|#



;; ====================================
;; Examples suspension2.lisp
;; ====================================



#|
  (setf cc (sphere 2))
  (setf basis (suspension2-basis (basis cc)))
  (dotimes (i 6)
    (print (funcall basis i)))
  (suspension-basis :locally-effective)
|#

#|
  (setf cc (deltab))
  (setf intr (suspension2-intr-dffr (dffr cc)))
  (funcall intr (cmbn 3 66 7))
|#


#|
  (cat-init)
  (setf cc (deltab))
  (setf scc (suspension2 cc))
  (cmpr scc 7 11)
  (basis scc 3)  ;; error
  (bsgn scc)
  (? scc 3 11)
|#

#|
  (setf f (idnt-mrph (deltab)))
  (setf sf (suspension2-intr f))
  (funcall sf (cmbn 2 4 3))
  (setf d (dffr (deltab)))
  (setf sd (suspension2-intr d))
  (funcall sd (cmbn 2 4 3))
  (funcall sd (cmbn 3 4 7))
|#

#|
  (setf f (idnt-mrph (deltab)))
  (setf sf (suspension2 f))
  (? sf (cmbn 2 4 3))
  (setf d (dffr (deltab)))
  (setf sd (suspension2 d))
  (? sd (cmbn 2 4 3))
  (? sd (cmbn 3 4 7))
|#


#|

(cat-init)
(setf s3xs3 (crts-prdc (sphere 3) (sphere 3)))
(setf s3xs3-efhm (efhm s3xs3))
(setf s_2-s3xs3-efhm (suspension2 s3xs3-efhm))
(setf lrdct (lrdct s_2-s3xs3-efhm))

(pre-check-rdct lrdct)
(setf *tc* (cmbn 3 1 (first (basis (k 38) 3)) -2 (second (basis (k 38) 3))))
(setf *bc* (cmbn 4 -1 (first (basis (k 40) 4)) -5 (second (basis (k 40) 4)) 7 (third (basis (k 40) 4))))
(check-rdct)

(setf rrdct (rrdct s_2-s3xs3-efhm))

(pre-check-rdct rrdct)
(setf *tc* (cmbn 3 1 (first (basis (k 38) 3)) -2 (second (basis (k 38) 3))))
(setf *bc* (cmbn 4 -1 (first (basis (k 46) 4)) -5 (second (basis (k 46) 4))))
(check-rdct)

|#



#|
  (setf sk (suspension2 (k-z 2)))
  (homology sk 0 10)
|#



;; ====================================
;; Examples Pushout.lisp
;; ====================================




#|
(setf s3 (sphere 3))
(setf s3xs3 (crts-prdc (sphere 3) (sphere 3)))

(setf f (build-smmr :sorc s3xs3 :trgt s3 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               (absm (dgop1 gmsm) (gmsm1 gmsm)))
                    :orgn `(first ,s3xs3)))

(setf g (build-smmr :sorc s3xs3 :trgt s3 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               (absm (dgop2 gmsm) (gmsm2 gmsm)))
                    :orgn `(second ,s3xs3)))


(setf p (pushout-cmpr (cmpr (crts-prdc s3xs3 (delta 1)))
                      (cmpr s3)
                      (cmpr s3)))

;;;JJJ Coherent examples.

(funcall p (pushout-gsm 0 (crpr 0 (crpr 0 's3 7 '*) 7 1))
         (pushout-gsm 0 (crpr 0 (crpr 0 's3 7 '*) 7 1)))
(funcall p (pushout-gsm 0 (crpr 0 (crpr 0 's3 7 '*) 7 1))
         (pushout-gsm 0 (crpr 0 (crpr 7 '* 0 's3) 7 1)))
(funcall p (pushout-gsm 0 (crpr 0 (crpr 7 '* 0 's3) 7 1))
         (pushout-gsm 0 (crpr 0 (crpr 0 's3 7 '*) 7 1)))
(funcall p (pushout-gsm 0 (crpr 0 (crpr 0 's3 7 '*) 7 1))
         (pushout-gsm 1 's3))
(funcall p (pushout-gsm 1 's3) (pushout-gsm 1 's3))
(funcall p (pushout-gsm 1 's3) (pushout-gsm 2 's3))
(funcall p (pushout-gsm 2 's3) (pushout-gsm 1 's3))
|#


#|
(setf s3 (sphere 3))
(setf s3xs3 (crts-prdc (sphere 3) (sphere 3)))

(setf f (build-smmr :sorc s3xs3 :trgt s3 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               (absm (dgop1 gmsm) (gmsm1 gmsm)))
                    :orgn `(first ,s3xs3)))

(setf g (build-smmr :sorc s3xs3 :trgt s3 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               (absm (dgop2 gmsm) (gmsm2 gmsm)))
                    :orgn `(second ,s3xs3)))


(setf p-basis (pushout-basis (basis (crts-prdc s3xs3 (delta 1)))
                             (basis s3)
                             (basis s3)))

(dolist (item (funcall p-basis 3))
  (print item))


|#


#|
(setf s3 (sphere 3))
(setf s3xs3 (crts-prdc (sphere 3) (sphere 3)))

(setf f (build-smmr :sorc s3xs3 :trgt s3 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               (absm (dgop1 gmsm) (gmsm1 gmsm)))
                    :orgn `(first ,s3xs3)))

(setf g (build-smmr :sorc s3xs3 :trgt s3 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               (absm (dgop2 gmsm) (gmsm2 gmsm)))
                    :orgn `(second ,s3xs3)))


(setf p-face (pushout-face (face (crts-prdc s3xs3 (delta 1)))
                           (face s3)
                           (face s3) f g))

(setf p-basis (pushout-basis (basis (crts-prdc s3xs3 (delta 1)))
                             (basis s3) (basis s3)))

(map nil #'print (funcall p-basis 3))

(dolist (item (funcall p-basis 3))
  (terpri)
  (print item)
  (dotimes (i 4)
    (print (funcall p-face i 3 item))))

|#


#|
(cat-init)
(setf s3 (sphere 3))
(setf s3xs3 (crts-prdc (sphere 3) (sphere 3)))

(setf f (build-smmr :sorc s3xs3 :trgt s3 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               (absm (dgop1 gmsm) (gmsm1 gmsm)))
                    :orgn `(first ,s3xs3)))

(setf g (build-smmr :sorc s3xs3 :trgt s3 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               (absm (dgop2 gmsm) (gmsm2 gmsm)))
                    :orgn `(second ,s3xs3)))
(setf p (pushout f g))


;;;JJJ Very useful: several bugs so discovered.

(check-smst p 0 8)

|#


#|
()
(cat-init)
(setf unipunctual (build-finite-ss '(x)))
(setf k-z (k-z 1))


(setf f (build-smmr :sorc unipunctual :trgt k-z :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               (if (= dmns 0)
                                   (absm 0 (bspn k-z))
                                 nil
                                 ))
                    :orgn `(triv ,unipunctual ,k-z)))


(setf p (pushout f f))

(homology p 0 3)

|#


#|

  (setf wkk (wedge (k-z 2) (k-z 2)))
  (homology wkk 0 10)
  (setf wedge-loop (wedge (loop-space (sphere 3)) (loop-space (sphere 3))))
  (homology wedge-loop 0 10)

|#



#|
  (cat-init)
  (setf s3*s4 (join (sphere 3) (sphere 4)))
  (homology s3*s4 0 9)
  (setf k-z2*k-z2 (join (k-z2 1) (k-z2 1)))  
  (homology k-z2*k-z2 0 10)

|#



#|
  (cat-init)
  (setf s3*s4 (join (sphere 3) (sphere 4)))
  (setf hom-equ (efhm s3*s4))
  (setf lrdct (lrdct hom-equ))
  (pre-check-rdct lrdct)
  (setf *bc* (CMBN 3 2 (CAR (basis (k 28) 3)) -8 (CADR (basis (k 28) 3)) 7 (CADDR (basis (k 28) 3))))
  (SETF *tc* (CMBN 2 -1 (CAR (basis (k 195) 2)) 4 (CADR (basis (k 195) 2)) -5 (CADDR (basis (k 195) 2))))
  (check-rdct)

  (setf rrdct (rrdct hom-equ))
  (pre-check-rdct rrdct)
  (setf *bc* (CMBN 3 2 (CAR (basis (k 198) 3)) -8 (CADR (basis (k 198) 3)) 7 (CADDR (basis (k 198) 3))))
  (SETF *tc* (CMBN 2 -1 (CAR (basis (k 195) 2)) 4 (CADR (basis (k 195) 2)) -5 (CADDR (basis (k 195) 2))))
  (check-rdct)
|#


;;;----------------------------------------------------------------------------
;;;
;;; Subject: Disjoint Union (Direct Sum)
;;; 
;;; X = \emptyset        
;;;
;;;----------------------------------------------------------------------------


(IN-PACKAGE #:cat)
(setf s3 (sphere 3))
(setf s4 (sphere 4))
(setf nil-ss (build-finite-ss nil))

(setf f (build-smmr :sorc nil-ss :trgt s3 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               nil)
                    :orgn `(triv ,nil-ss ,s3)))

(setf g (build-smmr :sorc nil-ss :trgt s4 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               nil)
                    :orgn `(triv ,nil-ss ,s4)))

(setf p (pushout f g))

(homology p 0 6)

(cat-init)
(setf l2s3 (loop-space (sphere 3) 2))
(setf l3s4 (loop-space (sphere 4) 3))
(setf nil-ss (build-finite-ss nil))

(setf f (build-smmr :sorc nil-ss :trgt l2s3 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               nil)
                    :orgn `(triv ,nil-ss ,l2s3)))

(setf g (build-smmr :sorc nil-ss :trgt l3s4 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               nil)
                    :orgn `(triv ,nil-ss ,l3s4)))

(setf p (pushout f g))
(setf ds (direct-sum l2s3 l3s4))

(homology p 0 6)
(homology ds 0 6)



;;;----------------------------------------------------------------------------
;;;
;;; Subject: Suspension
;;; 
;;; Arbitrary X, Y=Z={*}
;;;
;;;----------------------------------------------------------------------------

(cat-init)
(setf ls3 (loop-space (sphere 3)))
(setf unipunctual (build-finite-ss '(x)))

(setf f (build-smmr :sorc ls3 :trgt unipunctual :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               (if (and (equal dmns 0)
                                        (equal gmsm (bsgn ls3)))
                                   'x
                                 nil))
                    :orgn `(triv ,ls3 ,unipunctual)))

(setf p (pushout f f))

(homology p 0 6)

(homology (suspension2 (loop-space (sphere 3))) 0 6)


;;;----------------------------------------------------------------------------
;;;
;;; Subject: P^2(C)
;;; 
;;; X = X3, Y = S^2, Z = {x}
;;;
;;; Homology groups of P^2(C) -> (Z,0,Z,0,Z)
;;;
;;;----------------------------------------------------------------------------

(cat-init)
(setf s2 (sphere 2))
(setf ch2 ( chml-clss s2 2))
(setf f2 (z-whitehead s2 ch2))
(setf x3 (fibration-total f2))
(setf unipunctual (build-finite-ss '(x)))

(setf f (build-smmr :sorc x3 :trgt s2 :degr 0
                    :sintr #'(lambda (dmns gmsm)
                                  (declare (ignore dmns))
                               (absm (dgop1 gmsm) (gmsm1 gmsm)))
                    :orgn `(proj ,x3 ,s2)))

(setf g (build-smmr :sorc x3 :trgt unipunctual :degr 0
                    :sintr #'(lambda (dmns gmsm)
                               (if (and (equal dmns 0)
                                        (equal gmsm (bsgn x3)))
                                   'x
                                 nil))
                    :orgn `(proj ,x3 ,unipunctual)))

(setf p (pushout f g))

(homology p 0 10)

