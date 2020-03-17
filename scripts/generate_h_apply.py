#!/usr/bin/env python3

import os

keywords = """
check_double_excitation
copy_buffer
filter_only_connected_to_hf_single
filter_only_connected_to_hf_double
declarations
decls_main
deinit_thread
init_main
filter_integrals
filter2p
filter2h2p_double
filter2h2p_single
filter1h
filter1p
only_2p_single
only_2p_double
only_2h_single
only_2h_double
only_1h_single
only_1h_double
only_1p_single
only_1p_double
only_2h1p_single
only_2h1p_double
filter_only_1h1p_single
filter_only_1h1p_double
filter_only_1h2p_single
filter_only_1h2p_double
filter_only_2h2p_single
filter_only_2h2p_double
filterhole
filter_only_1h1p_double
filter_only_1h1p_single
filterparticle
filter_vvvv_excitation
finalization
generate_psi_guess
initialization
init_main
init_thread
keys_work
omp_barrier
omp_do
omp_enddo
omp_end_master
omp_end_parallel
omp_master
omp_parallel
only_2p_double
only_2p_single
parameters
params_main
printout_always
printout_now
subroutine
""".split()

class H_apply(object):

  def read_template(self):
    file = open(os.environ["QP_ROOT"]+'/src/determinants/h_apply.template.f','r')
    self.template = file.read()
    file.close()
    file = open(os.environ["QP_ROOT"]+'/src/determinants/h_apply_nozmq.template.f','r')
    self.template += file.read()
    file.close()

  def __init__(self,sub,SingleRef=False,do_mono_exc=True, do_double_exc=True):
    self.read_template()
    s = {}
    for k in keywords:
      s[k] = ""
    s["subroutine"] = "H_apply_%s"%(sub)
    s["params_post"] = ""

    self.selection_pt2 = None
    self.energy = "CI_electronic_energy"
    self.perturbation = None
    self.do_double_exc = do_double_exc
#    s["omp_parallel"]     = """ PROVIDE elec_num_tab
#        !$OMP PARALLEL DEFAULT(SHARED)        &
#        !$OMP PRIVATE(i,j,k,l,keys_out,hole,particle,                &
#        !$OMP  occ_particle,occ_hole,j_a,k_a,other_spin,             &
#        !$OMP  hole_save,ispin,jj,l_a,ib_jb_pairs,array_pairs,       &
#        !$OMP  accu,i_a,hole_tmp,particle_tmp,occ_particle_tmp,      &
#        !$OMP  occ_hole_tmp,key_idx,i_b,j_b,key,N_elec_in_key_part_1,&
#        !$OMP  N_elec_in_key_hole_1,N_elec_in_key_part_2,            &
#        !$OMP  N_elec_in_key_hole_2,ia_ja_pairs,key_union_hole_part) &
#        !$OMP SHARED(key_in,N_int,elec_num_tab,mo_num,           &
#        !$OMP  hole_1, particl_1, hole_2, particl_2,                 &
#        !$OMP  elec_alpha_num,i_generator) FIRSTPRIVATE(iproc)"""
#    s["omp_end_parallel"] = "!$OMP END PARALLEL"
#    s["omp_master"]       = "!$OMP MASTER"
#    s["omp_end_master"]   = "!$OMP END MASTER"
#    s["omp_barrier"]      = "!$OMP BARRIER"
#    s["omp_do"]           = "!$OMP DO SCHEDULE (static,1)"
#    s["omp_enddo"]        = "!$OMP ENDDO"

    d = { True : '.True.', False : '.False.'}
    s["do_mono_excitations"] = d[do_mono_exc]
    s["do_double_excitations"] = d[do_double_exc]
    s["keys_work"]  += "call fill_H_apply_buffer_no_selection(key_idx,keys_out,N_int,iproc)"
    s["filter_integrals"] = "array_pairs = .True."
    if SingleRef:
      s["filter_integrals"] = """
      call get_mo_bielec_integrals_existing_ik(i_a,j_a,mo_num,array_pairs,mo_integrals_map)
      """

    s["generate_psi_guess"]  = """
  ! Sort H_jj to find the N_states lowest states
  integer                        :: i
  integer, allocatable           :: iorder(:)
  double precision, allocatable  :: H_jj(:)
  double precision, external     :: diag_h_mat_elem
  allocate(H_jj(N_det),iorder(N_det))
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP SHARED(psi_det,N_int,H_jj,iorder,N_det)                  &
      !$OMP PRIVATE(i)
  !$OMP DO
  do i = 1, N_det
    H_jj(i) = diag_h_mat_elem(psi_det(1,1,i),N_int)
    iorder(i) = i
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dsort(H_jj,iorder,N_det)
  do k=1,N_states
    psi_coef(iorder(k),k) = 1.d0
  enddo
  deallocate(H_jj,iorder)
    """

    s["size_max"] = "8192"
    s["copy_buffer"] = """call copy_H_apply_buffer_to_wf
  if (s2_eig) then
    call make_s2_eigenfunction
  endif
  SOFT_TOUCH psi_det psi_coef N_det
"""
    s["printout_now"]   = """write(6,*)  &
       100.*float(i_generator)/float(N_det_generators), '% in ', wall_1-wall_0, 's'"""
    self.data = s

  def __setitem__(self,key,value):
    self.data[key] = value

  def __getitem__(self,key):
    return self.data[key]

  def __repr__(self):
    buffer = self.template
    for key,value in list(self.data.items()):
      buffer = buffer.replace('$'+key, value)
    return buffer

  def unset_double_excitations(self):
    self["do_double_excitations"] = ".False."
    self["check_double_excitation"] = """
     check_double_excitation = .False.
    """

  def filter_vvvv_excitation(self):
    self["filter_vvvv_excitation"] = """
     key_union_hole_part = 0_bit_kind
     call set_bit_to_integer(i_a,key_union_hole_part,N_int)
     call set_bit_to_integer(j_a,key_union_hole_part,N_int)
     call set_bit_to_integer(i_b,key_union_hole_part,N_int)
     call set_bit_to_integer(j_b,key_union_hole_part,N_int)
     do jtest_vvvv = 1, N_int
      if(iand(key_union_hole_part(jtest_vvvv),virt_bitmask(jtest_vvvv,1).ne.key_union_hole_part(jtest_vvvv)))then
       b_cycle = .False.
      endif
     enddo
     if(b_cycle) cycle
    """
  def set_filter_holes(self):
    self["filterhole"] = """
     if(iand(ibset(0_bit_kind,j),hole(k,other_spin)).eq.0_bit_kind )cycle
    """
  def set_filter_particl(self):
    self["filterparticle"] = """
     if(iand(ibset(0_bit_kind,j_a),hole(k_a,other_spin)).eq.0_bit_kind )cycle
    """
  def filter_1h(self):
    self["filter1h"] = """
!    ! DIR$ FORCEINLINE
     if (is_a_1h(hole)) cycle
    """
  def filter_2p(self):
    self["filter2p"] = """
!    ! DIR$ FORCEINLINE
     if (is_a_2p(hole)) cycle
    """
  def filter_1p(self):
    self["filter1p"] = """
!    ! DIR$ FORCEINLINE
     if (is_a_1p(hole)) cycle
    """

  def filter_only_2h(self):
    self["only_2h_single"] = """
!    ! DIR$ FORCEINLINE
     if (.not.is_a_2h(hole)) cycle
    """
    self["only_2h_double"] = """
!    ! DIR$ FORCEINLINE
     if (.not.is_a_2h(key))cycle
    """

  def filter_only_1h(self):
    self["only_1h_single"] = """
!    ! DIR$ FORCEINLINE
     if (.not.is_a_1h(hole)) cycle
    """
    self["only_1h_double"] = """
!    ! DIR$ FORCEINLINE
     if (.not.is_a_1h(key)  ) cycle
    """

  def filter_only_1p(self):
    self["only_1p_single"] = """
!    ! DIR$ FORCEINLINE
     if (.not. is_a_1p(hole) ) cycle
    """
    self["only_1p_double"] = """
!    ! DIR$ FORCEINLINE
     if (.not. is_a_1p(key) ) cycle
    """

  def filter_only_2h1p(self):
    self["only_2h1p_single"] = """
!    ! DIR$ FORCEINLINE
     if (.not. is_a_2h1p(hole) ) cycle
    """
    self["only_2h1p_double"] = """
!    ! DIR$ FORCEINLINE
     if (.not.is_a_2h1p(key) ) cycle
    """


  def filter_only_2p(self):
    self["only_2p_single"] = """
!    ! DIR$ FORCEINLINE
     if (.not.is_a_2p(hole)) cycle
    """
    self["only_2p_double"] = """
!    ! DIR$ FORCEINLINE
     if (.not.is_a_2p(key)) cycle
    """


  def filter_only_1h1p(self):
    self["filter_only_1h1p_single"] = """
     if (.not.is_a_1h1p(hole)) cycle
    """
    self["filter_only_1h1p_double"] = """
     if (.not.is_a_1h1p(key)) cycle
    """



  def filter_only_2h2p(self):
    self["filter_only_2h2p_single"] = """
!    ! DIR$ FORCEINLINE
     if (.not.is_a_two_holes_two_particles(hole)) cycle
    """
    self["filter_only_2h2p_double"] = """
!    ! DIR$ FORCEINLINE
     if (.not.is_a_two_holes_two_particles(key)) cycle
    """


  def filter_only_1h2p(self):
    self["filter_only_1h2p_single"] = """
!    ! DIR$ FORCEINLINE
     if (.not.is_a_1h2p(hole)) cycle
    """
    self["filter_only_1h2p_double"] = """
!    ! DIR$ FORCEINLINE
     if (.not.is_a_1h2p(key)) cycle
    """


  def set_filter_2h_2p(self):
    self["filter2h2p_double"] = """
     if (is_a_two_holes_two_particles(key)) cycle
    """
    self["filter2h2p_single"] = """
     if (is_a_two_holes_two_particles(hole)) cycle
    """

  def filter_only_connected_to_hf(self):
    self["filter_only_connected_to_hf_single"] = """
     call connected_to_hf(hole,yes_no)
     if (.not.yes_no) cycle
    """
    self["filter_only_connected_to_hf_double"] = """
     call connected_to_hf(key,yes_no)
     if (.not.yes_no) cycle
    """


  def set_perturbation(self,pert):
    if self.perturbation is not None:
        raise
    self.perturbation = pert
    if pert is not None:
      self.data["parameters"] = ",sum_e_2_pert_in,sum_norm_pert_in,sum_H_pert_diag_in,N_st,Nint"
      self.data["declarations"] = """
      integer, intent(in)             :: N_st,Nint
      double precision, intent(inout) :: sum_e_2_pert_in(N_st)
      double precision, intent(inout) :: sum_norm_pert_in(N_st)
      double precision, intent(inout) :: sum_H_pert_diag_in(N_st)
      double precision                :: sum_e_2_pert(N_st)
      double precision                :: sum_norm_pert(N_st)
      double precision                :: sum_H_pert_diag(N_st)
      double precision, allocatable   :: e_2_pert_buffer(:,:)
      double precision, allocatable   :: coef_pert_buffer(:,:)
      ASSERT (Nint == N_int)
      """
      self.data["init_thread"] = """
      allocate (e_2_pert_buffer(N_st,size_max), coef_pert_buffer(N_st,size_max))
      do k=1,N_st
        sum_e_2_pert(k) = 0.d0
        sum_norm_pert(k) = 0.d0
        sum_H_pert_diag(k) = 0.d0
      enddo
      """

      self.data["deinit_thread"] = """
      ! OMP CRITICAL
      do k=1,N_st
        sum_e_2_pert_in(k) = sum_e_2_pert_in(k) + sum_e_2_pert(k)
        sum_norm_pert_in(k) = sum_norm_pert_in(k) + sum_norm_pert(k)
        sum_H_pert_diag_in(k) = sum_H_pert_diag_in(k) + sum_H_pert_diag(k)
      enddo
      ! OMP END CRITICAL
      deallocate (e_2_pert_buffer, coef_pert_buffer)
      """
      self.data["size_max"] = "8192"
      self.data["initialization"] = """
      PROVIDE psi_selectors_coef psi_selectors E_corr_per_selectors psi_det_sorted_bit
      """
      if self.do_double_exc == True:
          self.data["keys_work"] = """
!          if(check_double_excitation)then
            call perturb_buffer_%s(i_generator,keys_out,key_idx,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert, &
             sum_norm_pert,sum_H_pert_diag,N_st,N_int,key_mask,fock_diag_tmp,%s)
          """%(pert,self.energy)
      else:
          self.data["keys_work"] = """
            call perturb_buffer_by_mono_%s(i_generator,keys_out,key_idx,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert, &
             sum_norm_pert,sum_H_pert_diag,N_st,N_int,key_mask,fock_diag_tmp,%s)
          """%(pert,self.energy)


      self.data["finalization"] = """
      """
      self.data["copy_buffer"] = ""
      self.data["generate_psi_guess"] = ""

      self.data["params_main"] = "pt2, norm_pert, H_pert_diag, N_st"
      self.data["params_post"] = ","+self.data["params_main"] +", N_int"
      self.data["decls_main"] = """  integer, intent(in)            :: N_st
  double precision, intent(inout):: pt2(N_st)
  double precision, intent(inout):: norm_pert(N_st)
  double precision, intent(inout):: H_pert_diag(N_st)
  double precision               :: delta_pt2(N_st), norm_psi(N_st), pt2_old(N_st)
  PROVIDE N_det_generators
  do k=1,N_st
    pt2(k) = 0.d0
    norm_pert(k) = 0.d0
    H_pert_diag(k) = 0.d0
    norm_psi(k) = 0.d0
    delta_pt2(k) = 0.d0
    pt2_old(k) = 0.d0
  enddo
        write(6,'(A12, 1X, A8, 3(2X, A9), 2X, A8, 2X, A8, 2X, A8)') &
                 'N_generators', 'Norm', 'Delta PT2', 'PT2', 'Est. PT2', 'secs'
        write(6,'(A12, 1X, A8, 3(2X, A9), 2X, A8, 2X, A8, 2X, A8)') &
                 '============', '========', '=========', '=========', '=========', &
                 '========='
      """

      self.data["printout_always"] = """
      do k=1,N_st
          norm_psi(k) = norm_psi(k) + psi_coef_generators(i_generator,k)*psi_coef_generators(i_generator,k)
        delta_pt2(k) = pt2(k) - pt2_old(k)
      enddo
      """
      self.data["printout_now"] = """
      do k=1,N_st
        write(6,'(I10, 4(2X, F9.6), 2X, F8.1)') &
                 i_generator, norm_psi(k), delta_pt2(k), pt2(k), &
                 pt2(k)/(norm_psi(k)*norm_psi(k)), &
                 wall_1-wall_0
         pt2_old(k) = pt2(k)
      enddo
      """
#      self.data["omp_parallel"]    += """&
# !$OMP SHARED(N_st) PRIVATE(e_2_pert_buffer,coef_pert_buffer) &
# !$OMP PRIVATE(sum_e_2_pert, sum_norm_pert, sum_H_pert_diag)"""

  def set_selection_pt2(self,pert):
    if self.selection_pt2 is not None:
        raise
    self.set_perturbation(pert)
    self.selection_pt2 = pert
    if pert is not None:
      self.data["parameters"] += ",select_max_out"
      self.data["declarations"] += """
      double precision, intent(inout) :: select_max_out"""

      self.data["params_post"] += ", select_max(min(i_generator,size(select_max,1)))"
      self.data["size_max"] = "8192"
      self.data["copy_buffer"] = """
      call copy_H_apply_buffer_to_wf
      if (s2_eig) then
        call make_s2_eigenfunction
      endif
      SOFT_TOUCH psi_det psi_coef N_det
      selection_criterion_min = min(selection_criterion_min, maxval(select_max))*0.1d0
      selection_criterion = selection_criterion_min
      call write_double(6,selection_criterion,'Selection criterion')
      """
      self.data["keys_work"] = """
      e_2_pert_buffer = 0.d0
      coef_pert_buffer = 0.d0
      """ + self.data["keys_work"]
      self.data["keys_work"] += """
      call fill_H_apply_buffer_selection(key_idx,keys_out,e_2_pert_buffer, &
        coef_pert_buffer,N_st,N_int,iproc,select_max_out)
      """
#      self.data["omp_parallel"]    += """&
# !$OMP REDUCTION (max:select_max_out)"""


  def unset_openmp(self):
    for k in keywords:
      if k.startswith("omp_"):
        self[k] = ""


class H_apply_zmq(H_apply):

  def read_template(self):
    file = open(os.environ["QP_ROOT"]+'/src/determinants/h_apply.template.f','r')
    self.template = file.read()
    file.close()
    file = open(os.environ["QP_ROOT"]+'/src/determinants/h_apply_zmq.template.f','r')
    self.template += file.read()
    file.close()

  def set_perturbation(self,pert):
     H_apply.set_perturbation(self,pert)
     self.data["printout_now"] = ""
     self.data["printout_always"] = ""
     self.data["decls_main"] = """  integer, intent(in)            :: N_st
  double precision, intent(inout):: pt2(N_st)
  double precision, intent(inout):: norm_pert(N_st)
  double precision, intent(inout):: H_pert_diag(N_st)
  double precision               :: delta_pt2(N_st), norm_psi(N_st), pt2_old(N_st)
  PROVIDE N_det_generators
  do k=1,N_st
    pt2(k) = 0.d0
    norm_pert(k) = 0.d0
    H_pert_diag(k) = 0.d0
    norm_psi(k) = 0.d0
    energy(k) = %s(k)
  enddo
     """ % (self.energy)
     self.data["copy_buffer"] = """
    do i=1,N_det_generators
      do k=1,N_st
        pt2(k) = pt2(k) + pt2_generators(k,i)
        norm_pert(k) = norm_pert(k) + norm_pert_generators(k,i)
        H_pert_diag(k) = H_pert_diag(k) + H_pert_diag_generators(k,i)
      enddo
    enddo
     """

  def set_selection_pt2(self,pert):
     H_apply.set_selection_pt2(self,pert)

