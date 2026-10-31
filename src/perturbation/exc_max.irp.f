BEGIN_PROVIDER [integer, max_exc_pert]
 implicit none
 BEGIN_DOC
 ! Constant integer provider, currently fixed to 0.
 ! No code in the tree currently reads it.
 ! No in-tree caller; probably unused (author, 2026-10-31); kept
 ! because external plugins may use it.
 END_DOC
 max_exc_pert = 0
END_PROVIDER
