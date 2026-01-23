#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

/* Generic functions */

int gpu_ndevices() {
  return 0;
}

void gpu_set_device(int32_t i) {
  return;
}

void gpu_get_memory(size_t* free, size_t* total) {
  *free = 0;
  *total = 0;
}


/* Allocation functions */

void gpu_allocate(void** ptr, const int64_t n) {
  *ptr = malloc((size_t) n);
  if (*ptr == NULL) {
    perror("Allocation failed");
  }
}

void gpu_deallocate(void** ptr) {
  free(*ptr);
  *ptr = NULL;
}


/* Memory transfer functions */

void gpu_upload(const void* cpu_ptr, void* gpu_ptr, const int64_t n) {
  memcpy(gpu_ptr, cpu_ptr, n);
}

void gpu_download(const void* gpu_ptr, void* cpu_ptr, const int64_t n) {
  memcpy(cpu_ptr, gpu_ptr, n);
}

void gpu_copy(const void* gpu_ptr_src, void* gpu_ptr_dest, const int64_t n) {
  memcpy(gpu_ptr_dest, gpu_ptr_src, n);
}


/* Streams */

void gpu_stream_create(void** ptr) {
  *ptr = (void*) malloc(sizeof(char));
}

void gpu_stream_destroy(void** ptr) {
  free(*ptr);
  *ptr = NULL;
}

void gpu_set_stream(void* handle, void* stream) {
  return;
}

void gpu_synchronize() {
  return;
}

void gpu_stream_synchronize(void* stream) {
  return;
}


/* BLAS functions */

void gpu_blas_create(void** handle) {
  *handle = (void*) malloc(sizeof(char));
}


void gpu_blas_destroy(void** handle) {
  free(*handle);
  *handle = NULL;
}


double ddot_(const int32_t* n, const double* x, const int32_t* incx, const double* y, const int32_t* incy);

void gpu_ddot(void* handle, const int64_t n, const double* x, const int64_t incx, const double* y, const int64_t incy, double* result) {
  assert (handle != NULL);

  /* Convert to int32_t */
  int32_t n_, incx_, incy_;

  n_    = (int32_t) n;
  incx_ = (int32_t) incx;
  incy_ = (int32_t) incy;

  /* Check for integer overflows */
  assert ( (int64_t)    n_ == n   );
  assert ( (int64_t) incx_ == incx);
  assert ( (int64_t) incy_ == incy);

  *result = ddot_(&n_, x, &incx_, y, &incy_);
}


float sdot_(const int32_t* n, const float* x, const int32_t* incx, const float* y, const int32_t* incy);

void gpu_sdot(void* handle, const int64_t n, const float* x, const int64_t incx, const float* y, const int64_t incy, float* result) {
  assert (handle != NULL);

  /* Convert to int32_t */
  int32_t n_, incx_, incy_;

  n_    = (int32_t) n;
  incx_ = (int32_t) incx;
  incy_ = (int32_t) incy;

  /* Check for integer overflows */
  assert ( (int64_t)    n_ == n   );
  assert ( (int64_t) incx_ == incx);
  assert ( (int64_t) incy_ == incy);

  *result = sdot_(&n_, x, &incx_, y, &incy_);
}


void dgemv_(const char* transa, const int32_t* m, const int32_t* n, const double* alpha,
            const double* a, const int32_t* lda, const double* x, const int32_t* incx, const double* beta, double* y, const int32_t* incy);

void gpu_dgemv(void* handle, const char* transa, const int64_t m, const int64_t n, const double* alpha,
               const double* a, const int64_t lda, const double* x, const int64_t incx, const double* beta, double* y, const int64_t incy) {

  assert (handle != NULL);

  /* Convert to int32_t */
  int32_t m_, n_, lda_, incx_, incy_;

  m_    = (int32_t) m;
  n_    = (int32_t) n;
  lda_  = (int32_t) lda;
  incx_ = (int32_t) incx;
  incy_ = (int32_t) incy;

  /* Check for integer overflows */
  assert ( (int64_t)    m_ == m   );
  assert ( (int64_t)    n_ == n   );
  assert ( (int64_t)  lda_ == lda );
  assert ( (int64_t) incx_ == incx);
  assert ( (int64_t) incy_ == incy);

  dgemv_(transa, &m_, &n_, alpha, a, &lda_, x, &incx_, beta, y, &incy_);
}


void sgemv_(const char* transa, const int32_t* m, const int32_t* n, const float* alpha,
               const float* a, const int32_t* lda, const float* x, const int32_t* incx, const float* beta, float* y, const int32_t* incy);

void gpu_sgemv(void* handle, const char* transa, const int64_t m, const int64_t n, const float* alpha,
               const float* a, const int64_t lda, const float* x, const int64_t incx, const float* beta, float* y, const int64_t incy) {

  assert (handle != NULL);

  /* Convert to int32_t */
  int32_t m_, n_, lda_, incx_, incy_;

  m_    = (int32_t) m;
  n_    = (int32_t) n;
  lda_  = (int32_t) lda;
  incx_ = (int32_t) incx;
  incy_ = (int32_t) incy;

  /* Check for integer overflows */
  assert ( (int64_t)    m_ == m   );
  assert ( (int64_t)    n_ == n   );
  assert ( (int64_t)  lda_ == lda );
  assert ( (int64_t) incx_ == incx);
  assert ( (int64_t) incy_ == incy);

  sgemv_(transa, &m_, &n_, alpha, a, &lda_, x, &incx_, beta, y, &incy_);
}


void dgemm_(const char* transa, const char* transb, const int32_t* m, const int32_t* n, const int32_t* k, const double* alpha,
            const double* a, const int32_t* lda, const double* b, const int32_t* ldb, const double* beta, double* c, const int32_t* ldc);

void gpu_dgemm(void* handle, const char* transa, const char* transb, const int64_t m, const int64_t n, const int64_t k, const double* alpha,
               const double* a, const int64_t lda, const double* b, const int64_t ldb, const double* beta, double* c, const int64_t ldc) {

  assert (handle != NULL);

  /* Convert to int32_t */
  int32_t m_, n_, k_, lda_, ldb_, ldc_;

  m_   = (int32_t) m;
  n_   = (int32_t) n;
  k_   = (int32_t) k;
  lda_ = (int32_t) lda;
  ldb_ = (int32_t) ldb;
  ldc_ = (int32_t) ldc;

  /* Check for integer overflows */
  assert ( (int64_t)   m_ == m  );
  assert ( (int64_t)   n_ == n  );
  assert ( (int64_t)   k_ == k  );
  assert ( (int64_t) lda_ == lda);
  assert ( (int64_t) ldb_ == ldb);
  assert ( (int64_t) ldc_ == ldc);

  dgemm_(transa, transb, &m_, &n_, &k_, alpha, a, &lda_, b, &ldb_, beta, c, &ldc_);
}



void sgemm_(const char* transa, const char* transb, const int32_t* m, const int32_t* n, const int32_t* k, const float* alpha,
            const float* a, const int32_t* lda, const float* b, const int32_t* ldb, const float* beta, float* c, const int32_t* ldc);

void gpu_sgemm(void* handle, const char* transa, const char* transb, const int64_t m, const int64_t n, const int64_t k, const float* alpha,
               const float* a, const int64_t lda, const float* b, const int64_t ldb, const float* beta, float* c, const int64_t ldc) {

  assert (handle != NULL);

  /* Convert to int32_t */
  int32_t m_, n_, k_, lda_, ldb_, ldc_;

  m_   = (int32_t) m;
  n_   = (int32_t) n;
  k_   = (int32_t) k;
  lda_ = (int32_t) lda;
  ldb_ = (int32_t) ldb;
  ldc_ = (int32_t) ldc;

  /* Check for integer overflows */
  assert ( (int64_t)   m_ == m  );
  assert ( (int64_t)   n_ == n  );
  assert ( (int64_t)   k_ == k  );
  assert ( (int64_t) lda_ == lda);
  assert ( (int64_t) ldb_ == ldb);
  assert ( (int64_t) ldc_ == ldc);

  sgemm_(transa, transb, &m_, &n_, &k_, alpha, a, &lda_, b, &ldb_, beta, c, &ldc_);
}


void gpu_dgeam(void* handle, const char* transa, const char* transb, const int64_t m, const int64_t n, const double* alpha,
               const double* a, const int64_t lda, const double* beta, const double* b, const int64_t ldb, double* c, const int64_t ldc) {
  assert (handle != NULL);

  if ( (*transa == 'N' && *transb == 'N') ||
       (*transa == 'n' && *transb == 'N') ||
       (*transa == 'N' && *transb == 'n') ||
       (*transa == 'n' && *transb == 'n') ) {

     if (*alpha == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *beta * b[j*ldb+i];
         }
       }

     } else if (*beta == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[j*lda+i];
         }
       }

     } else {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[j*lda+i] + *beta * b[j*ldb+i];
         }
       }

     }

  } else if ( (*transa == 'N' && *transb == 'T') ||
              (*transa == 'n' && *transb == 'T') ||
              (*transa == 'N' && *transb == 't') ||
              (*transa == 'n' && *transb == 't') ) {

     if (*alpha == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *beta * b[i*ldb+j];
         }
       }

     } else if (*beta == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[j*lda+i];
         }
       }

     } else {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[j*lda+i] + *beta * b[i*ldb+j];
         }
       }

     }

  } else if ( (*transa == 'T' && *transb == 'N') ||
              (*transa == 't' && *transb == 'N') ||
              (*transa == 'T' && *transb == 'n') ||
              (*transa == 't' && *transb == 'n') ) {

     if (*alpha == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *beta * b[j*ldb+i];
         }
       }

     } else if (*beta == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[i*lda+j];
         }
       }

     } else {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[i*lda+j] + *beta * b[j*ldb+i];
         }
       }

     }

  } else if ( (*transa == 'T' && *transb == 'T') ||
              (*transa == 't' && *transb == 'T') ||
              (*transa == 'T' && *transb == 't') ||
              (*transa == 't' && *transb == 't') ) {

     if (*alpha == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *beta * b[i*ldb+j];
         }
       }

     } else if (*beta == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[i*lda+j];
         }
       }

     } else {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[i*lda+j] + *beta * b[i*ldb+j];
         }
       }

     }

  }
}


void gpu_sgeam(void* handle, const char* transa, const char* transb, const int64_t m, const int64_t n, const float* alpha,
               const float* a, const int64_t lda, const float* beta, const float* b, const int64_t ldb, float* c, const int64_t ldc) {
  assert (handle != NULL);

  if ( (*transa == 'N' && *transb == 'N') ||
       (*transa == 'n' && *transb == 'N') ||
       (*transa == 'N' && *transb == 'n') ||
       (*transa == 'n' && *transb == 'n') ) {

     if (*alpha == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *beta * b[j*ldb+i];
         }
       }

     } else if (*beta == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[j*lda+i];
         }
       }

     } else {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[j*lda+i] + *beta * b[j*ldb+i];
         }
       }

     }

  } else if ( (*transa == 'N' && *transb == 'T') ||
              (*transa == 'n' && *transb == 'T') ||
              (*transa == 'N' && *transb == 't') ||
              (*transa == 'n' && *transb == 't') ) {

     if (*alpha == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *beta * b[i*ldb+j];
         }
       }

     } else if (*beta == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[j*lda+i];
         }
       }

     } else {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[j*lda+i] + *beta * b[i*ldb+j];
         }
       }

     }

  } else if ( (*transa == 'T' && *transb == 'N') ||
              (*transa == 't' && *transb == 'N') ||
              (*transa == 'T' && *transb == 'n') ||
              (*transa == 't' && *transb == 'n') ) {

     if (*alpha == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *beta * b[j*ldb+i];
         }
       }

     } else if (*beta == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[i*lda+j];
         }
       }

     } else {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[i*lda+j] + *beta * b[j*ldb+i];
         }
       }

     }

  } else if ( (*transa == 'T' && *transb == 'T') ||
              (*transa == 't' && *transb == 'T') ||
              (*transa == 'T' && *transb == 't') ||
              (*transa == 't' && *transb == 't') ) {

     if (*alpha == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *beta * b[i*ldb+j];
         }
       }

     } else if (*beta == 0.) {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[i*lda+j];
         }
       }

     } else {

       for (int64_t j=0 ; j<n ; ++j) {
         for (int64_t i=0 ; i<m ; ++i) {
           c[j*ldc+i] = *alpha * a[i*lda+j] + *beta * b[i*ldb+j];
         }
       }

     }

  }
}
