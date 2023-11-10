BEGIN_PROVIDER [ integer*8, qmckl_ctx_jastrow ]
  use qmckl
  use iso_c_binding
  implicit none
  BEGIN_DOC
  ! Context for the QMCKL library
  END_DOC
  integer(qmckl_exit_code) :: rc

  qmckl_ctx_jastrow = qmckl_context_create()

  rc =  qmckl_set_jastrow_champ_spin_independent(qmckl_ctx_jastrow, 1)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

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
  rc =  qmckl_set_jastrow_champ_type_nucl_num(qmckl_ctx_jastrow, 1_8*jast_qmckl_type_nucl_num)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_type_nucl_vector(qmckl_ctx_jastrow, 1_8*jast_qmckl_type_nucl_vector-1_8, 1_8*nucl_num)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_rescale_factor_ee(qmckl_ctx_jastrow, jast_qmckl_rescale_ee)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_rescale_factor_en(qmckl_ctx_jastrow, jast_qmckl_rescale_en, 1_8*jast_qmckl_type_nucl_num)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_aord_num(qmckl_ctx_jastrow, jast_qmckl_aord_num*1_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_a_vector(qmckl_ctx_jastrow, jast_qmckl_a_vector, 1_8*size(jast_qmckl_a_vector))
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_bord_num(qmckl_ctx_jastrow, jast_qmckl_bord_num*1_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_b_vector(qmckl_ctx_jastrow, jast_qmckl_b_vector, 1_8*size(jast_qmckl_b_vector))
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1


  rc =  qmckl_set_jastrow_champ_cord_num(qmckl_ctx_jastrow, jast_qmckl_cord_num*1_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  if (jast_qmckl_cord_num > 0) then
    rc =  qmckl_set_jastrow_champ_c_vector(qmckl_ctx_jastrow, jast_qmckl_c_vector, 1_8*jast_qmckl_c_vector_size)
    rc = qmckl_check(qmckl_ctx_jastrow, rc)
    if (rc /= QMCKL_SUCCESS) stop -1
  endif

END_PROVIDER
