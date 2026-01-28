#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <cublas_v2.h>
#include <cuda_runtime.h>


/* Generic functions */

int gpu_ndevices() {
  int ngpus;
  cudaGetDeviceCount(&ngpus);
  return ngpus;
}

void gpu_set_device(int32_t igpu) {
 cudaError_t rc = cudaSetDevice((int) igpu);
 if (rc != cudaSuccess) {
    fprintf(stderr,"cudaSetDevice(%d) failed: %s\n", igpu, cudaGetErrorString(rc));
    assert (rc == cudaSuccess);
 }
}

void gpu_get_memory(size_t* free, size_t* total) {
    cudaError_t rc = cudaMemGetInfo( free, total );
    if (rc != cudaSuccess) {
      *free = 0;
      *total = 0;
    }
}

/* Allocation functions */

void gpu_allocate(void** ptr, const int64_t size) {
    size_t free, total;
    cudaError_t rc = cudaMemGetInfo( &free, &total );
    if (rc != cudaSuccess) {
      free = INT64_MAX;
    }

//    rc = cudaMallocManaged(ptr, size, cudaMemAttachGlobal);
      rc= cudaMalloc(ptr, size);

//    /* Use managed memory if it does not fit on the GPU */
//    if (size < free && size < total/2) {
//      rc= cudaMalloc(ptr, size);
//    } else {
//      rc = cudaMallocManaged(ptr, size, cudaMemAttachGlobal);
//    }
    if (rc != cudaSuccess) {
      fprintf(stderr,"cudaMallocManaged failed: %s\n", cudaGetErrorString(rc));
      assert (rc == cudaSuccess);
    }
}

void gpu_deallocate(void** ptr) {
  assert (*ptr != NULL);
  cudaFree(*ptr);
  *ptr = NULL;
}


/* Memory transfer functions */

void gpu_upload(const void* cpu_ptr, void* gpu_ptr, const int64_t n) {
 cudaError_t rc = cudaMemcpy (gpu_ptr, cpu_ptr, n, cudaMemcpyHostToDevice);
 if (rc != cudaSuccess) {
    fprintf(stderr,"cudaMemcpy (upload) failed: %s\n", cudaGetErrorString(rc));
    assert (rc == cudaSuccess);
 }
}

void gpu_download(const void* gpu_ptr, void* cpu_ptr, const int64_t n) {
 cudaError_t rc = cudaMemcpy (cpu_ptr, gpu_ptr, n, cudaMemcpyDeviceToHost);
 if (rc != cudaSuccess) {
    fprintf(stderr,"cudaMemcpy (download) failed: %s\n", cudaGetErrorString(rc));
    assert (rc == cudaSuccess);
 }
}

void gpu_copy(const void* gpu_ptr_src, void* gpu_ptr_dest, const int64_t n) {
 cudaError_t rc = cudaMemcpy (gpu_ptr_dest, gpu_ptr_src, n, cudaMemcpyDeviceToDevice);
 if (rc != cudaSuccess) {
   fprintf(stderr,"cudaMemcpy (copy) failed: %s\n", cudaGetErrorString(rc));
   assert (rc == cudaSuccess);
 }
}


/* Streams */

void gpu_stream_create(cudaStream_t* ptr) {
  cudaError_t rc = cudaStreamCreate(ptr);
  if (rc != cudaSuccess) {
    fprintf(stderr,"cudaStreamCreate failed: %s\n", cudaGetErrorString(rc));
    assert (rc == cudaSuccess);
  }
}

void gpu_stream_destroy(cudaStream_t* ptr) {
  assert (ptr != NULL);
  cudaError_t rc = cudaStreamDestroy(*ptr);
  if (rc != cudaSuccess) {
    fprintf(stderr,"cudaStreamDestroy failed: %s\n", cudaGetErrorString(rc));
    assert (rc == cudaSuccess);
  }
  *ptr = NULL;
}

void gpu_set_stream(cublasHandle_t handle, cudaStream_t stream) {
  cudaError_t rc = cublasSetStream(handle, stream);
  if (rc != cudaSuccess) {
    fprintf(stderr,"cudaSetStream failed: %s\n", cudaGetErrorString(rc));
    assert (rc == cudaSuccess);
  }
}

void gpu_synchronize() {
  cudaError_t rc = cudaDeviceSynchronize();
  if (rc != cudaSuccess) {
    fprintf(stderr,"cudaDeviceSynchronize failed: %s\n", cudaGetErrorString(rc));
    assert (rc == cudaSuccess);
  }
}

void gpu_stream_synchronize(void* stream) {
  cudaError_t rc = cudaStreamSynchronize(stream);
  if (rc != cudaSuccess) {
    fprintf(stderr,"cudaStreamSynchronize failed: %s\n", cudaGetErrorString(rc));
    assert (rc == cudaSuccess);
  }
}


/* BLAS functions */

void gpu_blas_create(cublasHandle_t* ptr) {
  cublasStatus_t rc = cublasCreate(ptr);
  assert (rc == CUBLAS_STATUS_SUCCESS);
}


void gpu_blas_destroy(cublasHandle_t* ptr) {
  assert (ptr != NULL);
  cublasStatus_t rc = cublasDestroy(*ptr);
  assert (rc == CUBLAS_STATUS_SUCCESS);
  ptr = NULL;
}


void gpu_ddot(cublasHandle_t handle, const int64_t n, const double* x, const int64_t incx, const double* y, const int64_t incy, double* result) {
  assert (handle != NULL);
  /* Convert to int */
  int n_, incx_, incy_;

  n_    = (int) n;
  incx_ = (int) incx;
  incy_ = (int) incy;

  assert ( (int64_t)    n_ == n   );
  assert ( (int64_t) incx_ == incx);
  assert ( (int64_t) incy_ == incy);

  cublasStatus_t rc = cublasDdot(handle, n_, x, incx_, y, incy_, result);
  assert (rc == CUBLAS_STATUS_SUCCESS);
}



void gpu_sdot(cublasHandle_t handle, const int64_t n, const float* x, const int64_t incx, const float* y, const int64_t incy, float* result) {
  assert (handle != NULL);

  /* Convert to int */
  int n_, incx_, incy_;

  n_    = (int) n;
  incx_ = (int) incx;
  incy_ = (int) incy;

  /* Check for integer overflows */
  assert ( (int64_t)    n_ == n   );
  assert ( (int64_t) incx_ == incx);
  assert ( (int64_t) incy_ == incy);

  float result_ = 0.;
  cublasStatus_t rc = cublasSdot(handle, n_, x, incx_, y, incy_, &result_);
  assert (rc == CUBLAS_STATUS_SUCCESS);
  *result = result_;
}



void gpu_dgemv(cublasHandle_t handle, const char* transa, const int64_t m, const int64_t n, const double* alpha,
               const double* a, const int64_t lda, const double* x, const int64_t incx, const double* beta, double* y, const int64_t incy) {

  assert (handle != NULL);

  /* Convert to int */
  int m_, n_, lda_, incx_, incy_;

  m_    = (int) m;
  n_    = (int) n;
  lda_  = (int) lda;
  incx_ = (int) incx;
  incy_ = (int) incy;

  /* Check for integer overflows */
  assert ( (int64_t)    m_ == m   );
  assert ( (int64_t)    n_ == n   );
  assert ( (int64_t)  lda_ == lda );
  assert ( (int64_t) incx_ == incx);
  assert ( (int64_t) incy_ == incy);

  cublasOperation_t transa_ = CUBLAS_OP_N;
  if (*transa == 'T' || *transa == 't') transa_ = CUBLAS_OP_T;

  cublasDgemv(handle, transa_, m_, n_, alpha, a, lda_, x, incx_, beta, y, incy_);
}



void gpu_sgemv(cublasHandle_t handle, const char* transa, const int64_t m, const int64_t n, const float* alpha,
               const float* a, const int64_t lda, const float* x, const int64_t incx, const float* beta, float* y, const int64_t incy) {

  assert (handle != NULL);

  /* Convert to int */
  int m_, n_, lda_, incx_, incy_;

  m_    = (int) m;
  n_    = (int) n;
  lda_  = (int) lda;
  incx_ = (int) incx;
  incy_ = (int) incy;

  /* Check for integer overflows */
  assert ( (int64_t)    m_ == m   );
  assert ( (int64_t)    n_ == n   );
  assert ( (int64_t)  lda_ == lda );
  assert ( (int64_t) incx_ == incx);
  assert ( (int64_t) incy_ == incy);

  cublasOperation_t transa_ = CUBLAS_OP_N;
  if (*transa == 'T' || *transa == 't') transa_ = CUBLAS_OP_T;

  cublasSgemv(handle, transa_, m_, n_, alpha, a, lda_, x, incx_, beta, y, incy_);
}


void gpu_dgemm(cublasHandle_t handle, const char* transa, const char* transb, const int64_t m, const int64_t n, const int64_t k, const double* alpha,
               const double* a, const int64_t lda, const double* b, const int64_t ldb, const double* beta, double* c, const int64_t ldc) {

  assert (handle != NULL);

  /* Convert to int */
  int m_, n_, k_, lda_, ldb_, ldc_;

  m_   = (int) m;
  n_   = (int) n;
  k_   = (int) k;
  lda_ = (int) lda;
  ldb_ = (int) ldb;
  ldc_ = (int) ldc;

  /* Check for integer overflows */
  assert ( (int64_t)   m_ == m  );
  assert ( (int64_t)   n_ == n  );
  assert ( (int64_t)   k_ == k  );
  assert ( (int64_t) lda_ == lda);
  assert ( (int64_t) ldb_ == ldb);
  assert ( (int64_t) ldc_ == ldc);

  cublasOperation_t transa_ = CUBLAS_OP_N;
  cublasOperation_t transb_ = CUBLAS_OP_N;
  if (*transa == 'T' || *transa == 't') transa_ = CUBLAS_OP_T;
  if (*transb == 'T' || *transb == 't') transb_ = CUBLAS_OP_T;

  cublasDgemm(handle, transa_, transb_, m_, n_, k_, alpha, a, lda_, b, ldb_, beta, c, ldc_);
}



void gpu_sgemm(cublasHandle_t handle, const char* transa, const char* transb, const int64_t m, const int64_t n, const int64_t k, const float* alpha,
               const float* a, const int64_t lda, const float* b, const int64_t ldb, const float* beta, float* c, const int64_t ldc) {

  assert (handle != NULL);

  /* Convert to int */
  int m_, n_, k_, lda_, ldb_, ldc_;

  m_   = (int) m;
  n_   = (int) n;
  k_   = (int) k;
  lda_ = (int) lda;
  ldb_ = (int) ldb;
  ldc_ = (int) ldc;

  /* Check for integer overflows */
  assert ( (int64_t)   m_ == m  );
  assert ( (int64_t)   n_ == n  );
  assert ( (int64_t)   k_ == k  );
  assert ( (int64_t) lda_ == lda);
  assert ( (int64_t) ldb_ == ldb);
  assert ( (int64_t) ldc_ == ldc);

  cublasOperation_t transa_ = CUBLAS_OP_N;
  cublasOperation_t transb_ = CUBLAS_OP_N;
  if (*transa == 'T' || *transa == 't') transa_ = CUBLAS_OP_T;
  if (*transb == 'T' || *transb == 't') transb_ = CUBLAS_OP_T;

  cublasSgemm(handle, transa_, transb_, m_, n_, k_, alpha, a, lda_, b, ldb_, beta, c, ldc_);
}


void gpu_dgeam(cublasHandle_t handle, const char* transa, const char* transb, const int64_t m, const int64_t n, const double* alpha,
               const double* a, const int64_t lda, const double* beta, const double* b, const int64_t ldb, double* c, const int64_t ldc) {
  assert (handle != NULL);

  /* Convert to int */
  int m_, n_, lda_, ldb_, ldc_;

  m_   = (int) m;
  n_   = (int) n;
  lda_ = (int) lda;
  ldb_ = (int) ldb;
  ldc_ = (int) ldc;

  /* Check for integer overflows */
  assert ( (int64_t)   m_ == m  );
  assert ( (int64_t)   n_ == n  );
  assert ( (int64_t) lda_ == lda);
  assert ( (int64_t) ldb_ == ldb);
  assert ( (int64_t) ldc_ == ldc);

  cublasOperation_t transa_ = CUBLAS_OP_N;
  cublasOperation_t transb_ = CUBLAS_OP_N;
  if (*transa == 'T' || *transa == 't') transa_ = CUBLAS_OP_T;
  if (*transb == 'T' || *transb == 't') transb_ = CUBLAS_OP_T;

  cublasDgeam(handle, transa_, transb_, m_, n_, alpha, a, lda_, beta, b, ldb_, c, ldc_);

}


void gpu_sgeam(cublasHandle_t handle, const char* transa, const char* transb, const int64_t m, const int64_t n, const float* alpha,
               const float* a, const int64_t lda, const float* beta, const float* b, const int64_t ldb, float* c, const int64_t ldc) {
  assert (handle != NULL);

  /* Convert to int */
  int m_, n_, lda_, ldb_, ldc_;

  m_   = (int) m;
  n_   = (int) n;
  lda_ = (int) lda;
  ldb_ = (int) ldb;
  ldc_ = (int) ldc;

  /* Check for integer overflows */
  assert ( (int64_t)   m_ == m  );
  assert ( (int64_t)   n_ == n  );
  assert ( (int64_t) lda_ == lda);
  assert ( (int64_t) ldb_ == ldb);
  assert ( (int64_t) ldc_ == ldc);

  cublasOperation_t transa_ = CUBLAS_OP_N;
  cublasOperation_t transb_ = CUBLAS_OP_N;
  if (*transa == 'T' || *transa == 't') transa_ = CUBLAS_OP_T;
  if (*transb == 'T' || *transb == 't') transb_ = CUBLAS_OP_T;

  cublasSgeam(handle, transa_, transb_, m_, n_, alpha, a, lda_, beta, b, ldb_, c, ldc_);

}
