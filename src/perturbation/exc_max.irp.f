BEGIN_PROVIDER [integer, max_exc_pert]
 implicit none
 BEGIN_DOC
 ! Constant integer provider, currently fixed to 0.
 ! No code in the tree currently reads it.
 ! TODO(doc): is this meant to limit the excitation degree of the
 ! perturber determinants used in the perturbation?
 END_DOC
 max_exc_pert = 0
END_PROVIDER
