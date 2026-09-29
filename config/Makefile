#
#  _______                     _____
#  __  __ \___  _______ _________  /____  ________ ___
#  _  / / /  / / /  __ `/_  __ \  __/  / / /_  __ `__ \
#  / /_/ // /_/ // /_/ /_  / / / /_ / /_/ /_  / / / / /
#  \___\_\\__,_/ \__,_/ /_/ /_/\__/ \__,_/ /_/ /_/ /_/
#
#               ________             ______
#               ___  __ \_____ _________  /_______ _______ _____
#               __  /_/ /  __ `/  ___/_  //_/  __ `/_  __ `/  _ \
#               _  ____// /_/ // /__ _  ,<  / /_/ /_  /_/ //  __/
#               /_/     \__,_/ \___/ /_/|_| \__,_/ _\__, / \___/
#                                                  /____/
#
# https://github.com/QuantumPackage/qp2,
#
# Generated automatically by $(QP_ROOT)/scripts/compilation/qp_create_makefile
#
#

QP_ROOT = /home/runner/work/qp2/qp2

# Build configuration

FC = gfortran -g -ffree-line-length-none -I . -fPIC -march=native -std=legacy -fopenmp
FCFLAGS = -Ofast
IRPF90 = irpf90
IRPF90_FLAGS = --make --align=32 --assert -DSET_NESTED --openmp
LIB = -lblas -llapack $(QP_ROOT)/lib/libezfio_irp.a  -lz  -ltrexio -lf77zmq -lzmq
CONFIG_FILE = config/gfortran.cfg

# Main targets
.PHONY: all clean modules install

all: $(EZFIO_LIB) modules

modules: module_gpu_arch module_ao_basis module_ao_one_e_ints module_ao_two_e_ints module_aux_quantities module_basis module_becke_numerical_grid module_bitmask module_cas_based_on_top module_casscf_cipsi module_ccsd module_cipsi module_cipsi_utils module_cis module_cisd module_csf module_dav_general_mat module_davidson module_davidson_dressed module_davidson_keywords module_davidson_undressed module_density_for_dft module_determinants module_dft_keywords module_dft_one_e module_dft_utils_func module_dft_utils_in_r module_dressing module_dummy module_electrons module_ezfio_files module_fci module_functionals module_generators_cas module_generators_full module_generators_full_tc module_gpu module_hamiltonian module_hartree_fock module_iterations module_json module_kohn_sham module_kohn_sham_rs module_mo_basis module_mo_guess module_mo_one_e_ints module_mo_optimization module_mo_optimization_utils module_mo_two_e_ints module_mol_properties module_mp2 module_mpi module_mu_of_r module_nuclei module_perturbation module_pseudo module_psiref_cas module_psiref_utils module_scf_utils module_selectors_cassd module_selectors_full module_selectors_utils module_single_ref_method module_tools module_trexio module_two_body_rdm module_two_rdm_routines module_utils module_utils_cc module_utils_periodic module_utils_trust_region module_zmq

# Rule to build EZFIO interface files
%.ezfio_interface_config %ezfio_interface.irp.f: %EZFIO.cfg
	@echo 'Building EZFIO interface for $<'
	ei_handler.py --path_module $(dir $<)

# Rule to run IRPF90 preprocessing
%/IRPF90_temp/Makefile:
	@echo 'Running IRPF90 for $(dir $@)'
	@cd $(dir $@)/.. && \
	export FC='$(FC)' && \
	export FCFLAGS='$(FCFLAGS)' && \
	export LIB='$(LIB)' && \
	$(IRPF90) $(IRPF90_FLAGS)

# Rule to build module binaries
%.bin: %/IRPF90_temp/Makefile $(EZFIO_LIB)
	@echo 'Building binaries for $*'
	@cd $*/IRPF90_temp && $(MAKE) && \
	for i in $$(find . -type f -executable); do [ -x $$i ] && touch $$i; done

# EZFIO library
$(EZFIO_LIB):
	@echo 'Building EZFIO library'
	@cd $(QP_ROOT)/external/ezfio && $(MAKE)

module_gpu_arch: $(QP_ROOT)/plugins/local/gpu_x86.bin

module_ao_basis: $(QP_ROOT)/src/ao_basis.bin

module_ao_one_e_ints: $(QP_ROOT)/src/ao_one_e_ints.bin

module_ao_two_e_ints: $(QP_ROOT)/src/ao_two_e_ints.bin

module_aux_quantities: $(QP_ROOT)/src/aux_quantities.bin

module_basis: $(QP_ROOT)/src/basis.bin

module_becke_numerical_grid: $(QP_ROOT)/src/becke_numerical_grid.bin

module_bitmask: $(QP_ROOT)/src/bitmask.bin

module_cas_based_on_top: $(QP_ROOT)/src/cas_based_on_top.bin

module_casscf_cipsi: $(QP_ROOT)/src/casscf_cipsi.bin

module_ccsd: $(QP_ROOT)/src/ccsd.bin

module_cipsi: $(QP_ROOT)/src/cipsi.bin

module_cipsi_utils: $(QP_ROOT)/src/cipsi_utils.bin

module_cis: $(QP_ROOT)/src/cis.bin

module_cisd: $(QP_ROOT)/src/cisd.bin

module_csf: $(QP_ROOT)/src/csf.bin

module_dav_general_mat: $(QP_ROOT)/src/dav_general_mat.bin

module_davidson: $(QP_ROOT)/src/davidson.bin

module_davidson_dressed: $(QP_ROOT)/src/davidson_dressed.bin

module_davidson_keywords: $(QP_ROOT)/src/davidson_keywords.bin

module_davidson_undressed: $(QP_ROOT)/src/davidson_undressed.bin

module_density_for_dft: $(QP_ROOT)/src/density_for_dft.bin

module_determinants: $(QP_ROOT)/src/determinants.bin

module_dft_keywords: $(QP_ROOT)/src/dft_keywords.bin

module_dft_one_e: $(QP_ROOT)/src/dft_one_e.bin

module_dft_utils_func: $(QP_ROOT)/src/dft_utils_func.bin

module_dft_utils_in_r: $(QP_ROOT)/src/dft_utils_in_r.bin

module_dressing: $(QP_ROOT)/src/dressing.bin

module_electrons: $(QP_ROOT)/src/electrons.bin

module_ezfio_files: $(QP_ROOT)/src/ezfio_files.bin

module_fci: $(QP_ROOT)/src/fci.bin

module_functionals: $(QP_ROOT)/src/functionals.bin

module_generators_cas: $(QP_ROOT)/src/generators_cas.bin

module_generators_full: $(QP_ROOT)/src/generators_full.bin

module_generators_full_tc: $(QP_ROOT)/src/generators_full_tc.bin

module_gpu: $(QP_ROOT)/src/gpu.bin

module_hamiltonian: $(QP_ROOT)/src/hamiltonian.bin

module_hartree_fock: $(QP_ROOT)/src/hartree_fock.bin

module_iterations: $(QP_ROOT)/src/iterations.bin

module_json: $(QP_ROOT)/src/json.bin

module_kohn_sham: $(QP_ROOT)/src/kohn_sham.bin

module_kohn_sham_rs: $(QP_ROOT)/src/kohn_sham_rs.bin

module_mo_basis: $(QP_ROOT)/src/mo_basis.bin

module_mo_guess: $(QP_ROOT)/src/mo_guess.bin

module_mo_one_e_ints: $(QP_ROOT)/src/mo_one_e_ints.bin

module_mo_optimization: $(QP_ROOT)/src/mo_optimization.bin

module_mo_optimization_utils: $(QP_ROOT)/src/mo_optimization_utils.bin

module_mo_two_e_ints: $(QP_ROOT)/src/mo_two_e_ints.bin

module_mol_properties: $(QP_ROOT)/src/mol_properties.bin

module_mp2: $(QP_ROOT)/src/mp2.bin

module_mpi: $(QP_ROOT)/src/mpi.bin

module_mu_of_r: $(QP_ROOT)/src/mu_of_r.bin

module_nuclei: $(QP_ROOT)/src/nuclei.bin

module_perturbation: $(QP_ROOT)/src/perturbation.bin

module_pseudo: $(QP_ROOT)/src/pseudo.bin

module_psiref_cas: $(QP_ROOT)/src/psiref_cas.bin

module_psiref_utils: $(QP_ROOT)/src/psiref_utils.bin

module_scf_utils: $(QP_ROOT)/src/scf_utils.bin

module_selectors_cassd: $(QP_ROOT)/src/selectors_cassd.bin

module_selectors_full: $(QP_ROOT)/src/selectors_full.bin

module_selectors_utils: $(QP_ROOT)/src/selectors_utils.bin

module_single_ref_method: $(QP_ROOT)/src/single_ref_method.bin

module_tools: $(QP_ROOT)/src/tools.bin

module_trexio: $(QP_ROOT)/src/trexio.bin

module_two_body_rdm: $(QP_ROOT)/src/two_body_rdm.bin

module_two_rdm_routines: $(QP_ROOT)/src/two_rdm_routines.bin

module_utils: $(QP_ROOT)/src/utils.bin

module_utils_cc: $(QP_ROOT)/src/utils_cc.bin

module_utils_periodic: $(QP_ROOT)/src/utils_periodic.bin

module_utils_trust_region: $(QP_ROOT)/src/utils_trust_region.bin

module_zmq: $(QP_ROOT)/src/zmq.bin

clean:
	@echo 'Cleaning build files'
	@find $(QP_ROOT)/src -name 'IRPF90_temp' -type d -exec rm -rf {} + 2>/dev/null || true
	@find $(QP_ROOT)/src -name '*.bin' -delete 2>/dev/null || true
