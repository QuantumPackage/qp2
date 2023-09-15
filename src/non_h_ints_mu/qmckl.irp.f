BEGIN_PROVIDER [ integer*8, qmckl_ctx_jastrow ]
  use qmckl
  implicit none
  BEGIN_DOC
  ! Context for the QMCKL library
  END_DOC
  integer(qmckl_exit_code) :: rc

  qmckl_ctx_jastrow = qmckl_context_create()

  rc = qmckl_set_nucleus_num(qmckl_ctx_jastrow, nucl_num*1_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc = qmckl_set_nucleus_charge(qmckl_ctx_jastrow, nucl_charge, nucl_num*1_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc = qmckl_set_nucleus_coord(qmckl_ctx_jastrow, 'T', nucl_coord, nucl_num*3_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc = qmckl_set_electron_num(qmckl_ctx_jastrow, 1_8, 1_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1


  ! Jastrow parameters
  rc =  qmckl_set_jastrow_champ_type_nucl_num     (qmckl_ctx_jastrow, 2_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_type_nucl_vector  (qmckl_ctx_jastrow, (/0_8,1_8,1_8/), 1_8*nucl_num)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_rescale_factor_ee (qmckl_ctx_jastrow, 0.6d0)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_rescale_factor_en (qmckl_ctx_jastrow, (/0.6d0, 0.6d0 /), 2_8 )
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_aord_num          (qmckl_ctx_jastrow, 5_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_bord_num          (qmckl_ctx_jastrow, 5_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_cord_num          (qmckl_ctx_jastrow, 0_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

!  double precision :: a_vector(12) = dble(&
!           (/ 0.00000000,  0.00000000, -0.71168405, -0.44415699, -0.13865109,  0.07002267 , &
!              0.00000000,  0.00000000, -0.11379992,  0.04542846,  0.01696997, -0.01809299 /) )

!  double precision :: b_vector(6) = dble(&
!           (/  0.00000000,  0.65603311,  0.14581988,  0.03138163,  0.00153156, -0.00447302 /) )

!  double precision :: c_vector(46) = &
!           (/ 1.06384279d0, -1.44303973d0, -0.92409833d0,  0.11845356d0, -0.02980776d0, &
!              1.07048863d0,  0.06009623d0, -0.01854872d0, -0.00915398d0,  0.01324198d0, &
!             -0.00504959d0, -0.01202497d0, -0.00531644d0,  0.15101629d0, -0.00723831d0, &
!             -0.00384182d0, -0.00295036d0, -0.00114583d0,  0.00158107d0, -0.00078107d0, &
!             -0.00080000d0, -0.14140576d0, -0.00237271d0, -0.03006706d0,  0.01537009d0, &
!             -0.02327226d0,  0.16502789d0, -0.01458259d0, -0.09946065d0,  0.00850029d0, &
!             -0.02969361d0, -0.01159547d0,  0.00516313d0,  0.00405247d0, -0.02200886d0, &
!              0.03376709d0,  0.01277767d0, -0.01523013d0, -0.00739224d0, -0.00463953d0, &
!              0.00003174d0, -0.01421128d0,  0.00808140d0,  0.00612988d0, -0.00610632d0, &
!              0.01926215d0 /)

!     a_vector = 0.d0
!     b_vector = 0.d0
!     c_vector = 0.d0

   double precision :: a_vector(12) = dble(&
            (/ 0.00000000 , 0.00000000, -0.45105821, -0.23519218, -0.03825391,  0.10072866, &
               0.00000000 , 0.00000000, -0.06930592, -0.02909224, -0.00134650,  0.01477242 /) )

   double precision :: b_vector(6) = dble(&
            (/  0.00000000,  0.00000000,  0.29217862, -0.00450671, -0.02925982, -0.01381532 /) )

   double precision :: c_vector(46)
   c_vector = 0.d0

  rc =  qmckl_set_jastrow_champ_a_vector(qmckl_ctx_jastrow, a_vector, 12_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_b_vector(qmckl_ctx_jastrow, b_vector, 6_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

! rc =  qmckl_set_jastrow_champ_c_vector(qmckl_ctx_jastrow, c_vector, 46_8)
! rc = qmckl_check(qmckl_ctx_jastrow, rc)
! if (rc /= QMCKL_SUCCESS) stop -1

END_PROVIDER
