BEGIN_PROVIDER [ character*(32), DFT_TYPE]
 implicit none
 BEGIN_DOC
! defines the type of DFT applied: LDA, GGA etc ...
 END_DOC
 logical :: is_lda
 if(correlation_functional.eq."None")then
  is_lda             = (index(exchange_functional,"LDA")  .ne. 0)
 else if(exchange_functional.eq."None")then
  is_lda             = (index(correlation_functional,"LDA")  .ne. 0)
 else
  is_lda             = (index(correlation_functional,"LDA")  .ne. 0) .and. (index(exchange_functional,"LDA")  .ne. 0)
 endif
 if(is_lda)then
  DFT_TYPE = "LDA"
 else
  DFT_TYPE = "GGA"
 endif
END_PROVIDER
