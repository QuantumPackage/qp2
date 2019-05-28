BEGIN_PROVIDER [ integer, grid_atomic_number, (nucl_num) ]
 implicit none
 BEGIN_DOC
 ! Atomic number used to adjust the grid
 END_DOC
 grid_atomic_number(:) = max(1,int(nucl_charge(:)))

END_PROVIDER

