subroutine qp_bug(from, code, message)
  implicit none
  BEGIN_DOC
! This routine prints a bug report
  END_DOC
  character*(*) :: from
  integer :: code
  character*(*) :: message

  print *, ''
  print *, '======================='
  print *, 'Bug in Quantum Package!'
  print *, '======================='
  print *, ''
  print *, '   from: ', trim(from)
  print *, '   code: ', code
  print *, '   info: ', trim(message)
  print *, ''
  print *, 'Please report this bug at https://github.com/QuantumPackage/qp2/issues'
  print *, 'with your output file attached.'
  print *, ''
  stop -1
end subroutine qp_bug
