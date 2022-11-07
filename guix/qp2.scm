;; Custom Guix module definitions for dependencies of QP

(define-module (qp2)
  #:use-module (guix download)
  #:use-module (guix build-system dune)
  #:use-module (guix build-system gnu)
  #:use-module (guix build-system ocaml)
  #:use-module (guix build-system python)
  #:use-module (guix build-system copy)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix packages)
  #:use-module (gnu packages base)
  #:use-module (gnu packages bash)
  #:use-module (gnu packages certs)
  #:use-module (gnu packages commencement)
  #:use-module (gnu packages compression)
  #:use-module (gnu packages gcc)
  #:use-module (gnu packages maths)
  #:use-module (gnu packages multiprecision)
  #:use-module (gnu packages networking)
  #:use-module (gnu packages ninja)
  #:use-module (gnu packages ocaml)
  #:use-module (gnu packages python)
  #:use-module (gnu packages python-xyz)
  #:use-module (gnu packages version-control)
  #:use-module (gnu packages virtualization))

(define-public ocaml-cryptokit
  (package
    (name "ocaml-cryptokit")
    (version "1.17")
    (source (origin
              (method url-fetch)
              (uri
               "https://github.com/xavierleroy/cryptokit/archive/release117.tar.gz")
              (sha256
               (base32
                "0gsgb4bwy389ndkribgk5fmzjp3fgml0xi8sd9amx3krwd1k1qak"))))
    (build-system dune-build-system)
    (propagated-inputs (list zlib ocaml-zarith gmp))
    (home-page "https://github.com/xavierleroy/cryptokit")
    (synopsis "A library of cryptographic primitives")
    (description
     "Cryptokit includes authenticated encryption (AES-GCM, Chacha20-Poly1305), block
ciphers (AES, DES, 3DES), stream ciphers (Chacha20, ARCfour), public-key
cryptography (RSA, DH), hashes (SHA-256, SHA-512, SHA-3, Blake2), MACs,
compression, random number generation -- all presented with a compositional,
extensible interface.")
    (license license:lgpl2.0+)))

(define-public ocaml-getopt
  (package
    (name "ocaml-getopt")
    (version "20120615")
    (source (origin
              (method url-fetch)
              (uri
               "https://download.ocamlcore.org/ocaml-getopt/ocaml-getopt/20120615/ocaml-getopt-20120615.tar.gz")
              (sha256
               (base32
                "1rz2mi3gddwpd249bvz6h897swiajk4d6cczrsscibwpkmdvrfwa"))))
    (build-system ocaml-build-system)
    (native-inputs (list ocamlbuild))
    (home-page #f)
    (synopsis
     "Parsing of command line arguments (similar to GNU GetOpt) for OCaml")
    (description
     "General command line syntax of GNU getopt and getopt_long, but is close to the
spirit of the Arg module.")
    (license license:expat)))


(define-public ocaml-stdint
  (package
    (name "ocaml-stdint")
    (version "0.7.2")
    (source (origin
              (method url-fetch)
              (uri
               "https://github.com/andrenth/ocaml-stdint/releases/download/0.7.2/stdint-0.7.2.tbz")
              (sha256
               (base32
                "1rkpx48k5d7ypgwzvhzc2116jri6x2mhvi2jm4zaziy9if6ijq0m"))))
    (build-system dune-build-system)
    (arguments `(#:tests? #f))
    (propagated-inputs (list ocaml-odoc ocaml-base))
    (native-inputs (list ocaml-qcheck))
    (home-page "https://github.com/andrenth/ocaml-stdint")
    (synopsis "Signed and unsigned integer types having specified widths")
    (description
     "The stdint library provides signed and unsigned integer types of various fixed
widths: 8, 16, 24, 32, 40, 48, 56, 64 and 128 bit.  This interface is similar to
Int32 and Int64 from the base library but provides more functions and constants
like arithmetic and bit-wise operations, constants like maximum and minimum
values, infix operators conversion to and from every other integer type
(including int, float and nativeint), parsing from and conversion to readable
strings (binary, octal, decimal, hexademical), conversion to and from buffers in
both big endian and little endian byte order.")
    (license license:expat)))


(define-public ocaml-zmq
  (package
    (name "ocaml-zmq")
    (version "5.1.5")
    (source (origin
              (method url-fetch)
              (uri
               "https://github.com/issuu/ocaml-zmq/releases/download/5.1.5/zmq-lwt-5.1.5.tbz")
              (sha256
               (base32
                "0mkvlrr9ykmaidwch39gb7jpqnb90n50iw30rn95fg2bmcyx2iwr"))))
    (build-system dune-build-system)
    (arguments
     `(#:package "zmq"
       #:test-target "."))
    (propagated-inputs (list zeromq ocaml-stdint))
    (native-inputs (list ocaml-ounit2))
    (home-page "https://github.com/issuu/ocaml-zmq")
    (synopsis "OCaml bindings for ZeroMQ 4.x")
    (description
     "This library contains basic bindings for ZMQ.")
    (license license:expat)))


(define-public f77_zmq
  (package
   (name "f77_zmq")
   (version "4.3.3")
   (source (origin
            (method url-fetch)
            (uri (string-append
                  "https://github.com/zeromq/f77_zmq/releases/download/v"
                  version
                  "/f77-zmq-"
                  version
                  ".tar.gz"))
            (sha256
             (base32
              "0k70riil3fczymp17184rzfvzfvy0k8a6s9yfwyrrh2qyrz797hf"))))
   (build-system gnu-build-system)
   (arguments '(#:configure-flags '("--enable-silent-rules")))
   (inputs `(
             ("gcc", gcc)
             ("gfortran", gfortran)
             ("python", python)
             ("zeromq", zeromq)
             ))
   (synopsis "ZeroMQ Fortran77 binding")
   (description "Fortran77 binding for the ZeroMQ lightweight messaging library.")
   (home-page "https://github.com/zeromq/f77_zmq")
   (license license:lgpl2.1+)))


(define-public python-resultsfile
  (package
    (name "python-resultsfile")
    (version "2.4")
    (source (origin
              (method url-fetch)
              (uri (pypi-uri "resultsFile" version))
              (sha256
              (base32
                "199agrkihv3vh0y4h0gpsa5w2d24q991k4sgxmscs75nl3z8sx90"))))
    (build-system python-build-system)
    (home-page "https://gitlab.com/scemama/resultsFile")
    (synopsis "Module for reading output files of quantum chemistry codes.")
    (description "Module for reading output files of quantum chemistry codes.")
    (license license:gpl2+)))



;(define-public qp2-environment
;  (package
;    (name "qp2-environment")
;    (version "2.0")
;    (source #f)
;    (build-system copy-build-system)
;    (arguments
;     '(#:install-plan #f))
;    (inputs
;      (list nss-certs
;        ocaml-cryptokit
;        ocaml-zmq
;        git
;        ocaml-sexplib
;        ocaml-ppx-sexp-conv
;        ocaml-ppx-deriving
;        ocaml-getopt
;        ocamlbuild
;        python-resultsfile
;        f77_zmq
;        bubblewrap
;        ocaml-findlib
;        ocaml
;        bash
;        glibc-locales
;        python-docopt
;        bats
;        python
;        gmp
;        zeromq
;        zlib
;        ninja
;        openblas
;        gfortran-toolchain
;        gcc-toolchain))
;    (home-page "https://quantumpackage.github.io/qp2/")
;    (synopsis "A programming environment for wave function methods")
;    (description "Quantum Package is an open-source programming environment for
;    quantum chemistry specially designed for wave function methods. Its main
;    goal is the development of determinant-driven selected configuration
;    interaction (sCI) methods and multi-reference second-order perturbation
;    theory (PT2).")
;    (license license:agpl3+)))

