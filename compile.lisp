;;; This script compiles kenzo into a standalone .fasl file
;;; It relies on asdf
;;;
;;; The output is a file called kenzoo--all-systems.fasl


(require :asdf)
(push #P"./" asdf:*central-registry*)
(require :kenzo)
(asdf:make-build :kenzo :type :fasl :monolithic t :move-here #P".")
(quit)
