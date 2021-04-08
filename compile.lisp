;;; This script compiles kenzo into a standalone .fasl file
;;; It relies on asdf
;;;
;;; The output is a file called kenzoo--all-systems.fasl


(require :asdf)
(asdf:disable-output-translations)
(push #P"./" asdf:*central-registry*)
(require :kenzo)
(asdf:operate 'asdf:monolithic-compile-bundle-op :kenzo)
(quit)
