BEGIN_PROVIDER [ double precision, inertia_tensor, (3,3) ]
  implicit none
  BEGIN_DOC
  ! Inertia tensor
  END_DOC
  integer                        :: i,j,k
  inertia_tensor = 0.d0
  do k=1,nucl_num
    inertia_tensor(1,1) += element_mass(int(nucl_charge(k))) * ((nucl_coord(k,2)-center_of_mass(2))**2 + (nucl_coord(k,3)-center_of_mass(3))**2)
    inertia_tensor(2,2) += element_mass(int(nucl_charge(k))) * ((nucl_coord(k,1)-center_of_mass(1))**2 + (nucl_coord(k,3)-center_of_mass(3))**2)
    inertia_tensor(3,3) += element_mass(int(nucl_charge(k))) * ((nucl_coord(k,1)-center_of_mass(1))**2 + (nucl_coord(k,2)-center_of_mass(2))**2)
    inertia_tensor(1,2) -= element_mass(int(nucl_charge(k))) * ((nucl_coord(k,1)-center_of_mass(1))    * (nucl_coord(k,2)-center_of_mass(2))   )
    inertia_tensor(1,3) -= element_mass(int(nucl_charge(k))) * ((nucl_coord(k,1)-center_of_mass(1))    * (nucl_coord(k,3)-center_of_mass(3))   )
    inertia_tensor(2,3) -= element_mass(int(nucl_charge(k))) * ((nucl_coord(k,2)-center_of_mass(2))    * (nucl_coord(k,3)-center_of_mass(3))   )
  enddo
  inertia_tensor(2,1) = inertia_tensor(1,2)
  inertia_tensor(3,1) = inertia_tensor(1,3)
  inertia_tensor(3,2) = inertia_tensor(2,3)
END_PROVIDER

 BEGIN_PROVIDER [ double precision, inertia_tensor_eigenvectors, (3,3) ]
&BEGIN_PROVIDER [ double precision, inertia_tensor_eigenvalues , (3) ]
 implicit none
 BEGIN_DOC
 ! Eigenvectors/eigenvalues of the inertia_tensor. Used to find normal orientation.
 END_DOC
 integer :: k
 call lapack_diagd(inertia_tensor_eigenvalues,inertia_tensor_eigenvectors,-inertia_tensor,3,3)
 inertia_tensor_eigenvalues = -inertia_tensor_eigenvalues
 print *, 'Rotational constants (GHZ):'
 print *, (1805.65468542d0/(inertia_tensor_eigenvalues(k)+1.d-32), k=1,3)
END_PROVIDER

