.NOPARALLEL:

# Check if QP_ROOT is defined
ifndef QP_ROOT
$(info  -------------------- Error --------------------)
$(info  QP_ROOT undefined. Source the quantum_package.rc script)
$(info  -----------------------------------------------)
$(error )
endif


LIBS=
PKGS=
OCAMLCFLAGS="-g"
OCAMLOPTFLAGS="opt -O3 -remove-unused-arguments -rounds 16 -inline 100 -inline-max-unroll 100"
OCAMLBUILD=ocamlbuild -j 0 -cflags $(OCAMLCFLAGS) -lflags $(OCAMLCFLAGS) -ocamlopt $(OCAMLOPTFLAGS)
MLLFILES=$(wildcard *.mll) 
MLFILES=$(wildcard *.ml) ezfio.ml Qptypes.ml Input_auto_generated.ml qp_edit.ml 
MLIFILES=$(wildcard *.mli) git
ALL_TESTS=$(patsubst %.ml,%.byte,$(wildcard test_*.ml))
ALL_EXE=$(patsubst %.ml,%.native,$(wildcard qp_*.ml)) qp_edit.native element_create_db.byte

.PHONY: executables default remake_executables


default: $(ALL_EXE) .gitignore
tests: $(ALL_TESTS) 

.gitignore: $(MLFILES) $(MLIFILES)
	@for i in .gitignore ezfio.ml element_create_db Qptypes.ml Git.ml *.byte *.native _build $(ALL_EXE) $(ALL_TESTS) \
	$(patsubst %.ml,%,$(wildcard test_*.ml)) $(patsubst %.ml,%,$(wildcard qp_*.ml)) \
	Input_*.ml \
	qp_edit.ml qp_edit qp_edit.native Input_auto_generated.ml;\
	do \
	   echo $$i ; \
	done |sort | uniq > .gitignore

executables: $(QP_ROOT)/data/executables


$(QP_ROOT)/data/executables: remake_executables element_create_db.byte Qptypes.ml
	$(QP_ROOT)/scripts/module/create_executables_list.sh
	$(QP_ROOT)/ocaml/element_create_db.byte

external_libs:
	opam install cryptokit sexplib

qpackage.odocl: $(MLIFILES)
	ls $(MLIFILES) | sed "s/\.mli//" > qpackage.odocl

doc: qpackage.odocl
	$(OCAMLBUILD) qpackage.docdir/index.html -use-ocamlfind  $(PKGS) 

%.inferred.mli: $(MLFILES)
	$(OCAMLBUILD) $*.inferred.mli -use-ocamlfind  $(PKGS) 
	mv _build/$*.inferred.mli .

%.byte: $(MLFILES) $(MLIFILES)
	rm -f -- $*
	$(OCAMLBUILD) $*.byte  -use-ocamlfind  $(PKGS)
	ln -s $*.byte $*

qp_run.native: $(MLFILES) $(MLIFILES) executables

%.native: $(MLFILES) $(MLIFILES)
	rm -f -- $*
	$(OCAMLBUILD) $*.native -use-ocamlfind $(PKGS)
	ln -s $*.native $*

ezfio.ml: ${QP_ROOT}/external/ezfio/Ocaml/ezfio.ml
	cp ${QP_ROOT}/external/ezfio/Ocaml/ezfio.ml .

qptypes_generator.byte: qptypes_generator.ml
	$(OCAMLBUILD) qptypes_generator.byte  -use-ocamlfind 

Qptypes.ml: qptypes_generator.byte
	./qptypes_generator.byte > Qptypes.ml

git: 
	./create_git_sha1.sh

${QP_EZFIO}/Ocaml/ezfio.ml: 
	ninja -C ${QP_ROOT}/config ${QP_ROOT}/lib/libezfio_irp.a

qp_edit.ml: ../scripts/ezfio_interface/qp_edit_template

Input_auto_generated.ml qp_edit.ml: $(filter-out Input_auto_generated.ml, $(wildcard Input_*.ml))
	ei_handler.py ocaml_global

clean: 
	rm -rf _build Qptypes.ml Git.ml Input_auto_generated.ml $(ALL_EXE) $(ALL_TESTS)

