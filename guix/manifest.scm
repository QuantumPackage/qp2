;; This "manifest" file can be passed to 'guix package -m' to reproduce
;; the content of your profile.  This is "symbolic": it only specifies
;; package names.  To reproduce the exact same profile, you also need to
;; capture the channels being used, as returned by "guix describe".
;; See the "Replicating Guix" section in the manual.
;;
;; guix shell  --load-path="guix/"   \
;;             --manifest="manifest.scm"

(specifications->manifest
  (list "nss-certs"
        "ocaml-cryptokit"
        "ocaml-zmq"
        "git"
        "ocaml-sexplib"
        "ocaml-ppx-sexp-conv"
        "ocaml-ppx-deriving"
        "ocaml-getopt"
        "ocamlbuild"
        "python-resultsfile"
        "f77_zmq"
        "bubblewrap"
        "ocaml-findlib"
        "ocaml"
        "bash"
        "glibc-locales"
        "python-docopt"
        "bats"
        "python"
        "gmp"
        "zeromq"
        "zlib"
        "ninja"
        "openblas"
        "gfortran-toolchain"
        "gcc-toolchain"))
