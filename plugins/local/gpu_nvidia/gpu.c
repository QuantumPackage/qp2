#include <stdint.h>
#include <stdio.h>
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
  cudaSetDevice(igpu);
}


/* Allocation functions */

void gpu_allocate(void** ptr, const int64_t size) {
    size_t free, total;
    cudaError_t rc = cudaMemGetInfo( &free, &total );
    if (rc != cudaSuccess) {
      free = INT64_MAX;
    }

    /* Use managed memory if it does not fit on the GPU */
    if (size < free && size < total/2) {
//      rc= cudaMalloc(ptr, size);
      rc = cudaMallocManaged(ptr, size, cudaMemAttachGlobal);
    } else {
      rc = cudaMallocManaged(ptr, size, cudaMemAttachGlobal);
    }
    assert (rc == cudaSuccess);
}

void gpu_deallocate(void** ptr) {
  assert (*ptr != NULL);
  cudaFree(*ptr);
  *ptr = NULL;
}


/* Memory transfer functions */

void gpu_upload(const void* cpu_ptr, void* gpu_ptr, const int64_t n) {
  cudaMemcpy (gpu_ptr, cpu_ptr, n, cudaMemcpyHostToDevice);
}

void gpu_download(const void* gpu_ptr, void* cpu_ptr, const int64_t n) {
  cudaMemcpy (cpu_ptr, gpu_ptr, n, cudaMemcpyDeviceToHost);
}

void gpu_copy(const void* gpu_ptr_src, void* gpu_ptr_dest, const int64_t n) {
  cudaMemcpy (gpu_ptr_dest, gpu_ptr_src, n, cudaMemcpyDeviceToDevice);
}


/* Streams */

void gpu_stream_create(void** ptr) {
  cudaStream_t stream;
  cudaError_t rc = cudaStreamCreate(&stream);
  assert (rc == cudaSuccess);
  *ptr = (void*) stream;
}

void gpu_stream_destroy(void** ptr) {
  assert (*ptr != NULL);
  cudaError_t rc = cudaStreamDestroy( (cudaStream_t) *ptr);
  assert (rc == cudaSuccess);
  *ptr = NULL;
}

void gpu_set_stream(void** handle, void** stream) {
  cublasSetStream( (cublasHandle_t) *handle, (cudaStream_t) *stream);
}

void gpu_synchronize() {
  cudaDeviceSynchronize();
}


/* BLAS functions */

void gpu_blas_create(void** handle) {
  cublasHandle_t cublas_handle;
  cublasStatus_t rc = cublasCreate(&cublas_handle);
  assert (rc == CUBLAS_STATUS_SUCCESS);
  *handle = (void*) cublas_handle;
}


void gpu_blas_destroy(void** handle) {
  assert (*handle != NULL);
  cublasStatus_t rc = cublasDestroy( (cublasHandle_t) *handle);
  assert (rc == CUBLAS_STATUS_SUCCESS);
  *handle = NULL;
}


void gpu_ddot(void** handle, const int64_t n, const double* x, const int64_t incx, const double* y, const int64_t incy, double* result) {
  assert (*handle != NULL);

  /* Convert to int32_t */
  int32_t n_, incx_, incy_;

  n_    = (int32_t) n;
  incx_ = (int32_t) incx;
  incy_ = (int32_t) incy;

  /* Check for integer overflows */
  assert ( (int64_t)    n_ == n   );
  assert ( (int64_t) incx_ == incx);
  assert ( (int64_t) incy_ == incy);

  cublasDdot((cublasHandle_t) *handle, n_, x, incx_, y, incy_, result);
}



void gpu_sdot(void** handle, const int64_t n, const float* x, const int64_t incx, const float* y, const int64_t incy, float* result) {
  assert (*handle != NULL);

  /* Convert to int32_t */
  int32_t n_, incx_, incy_;

  n_    = (int32_t) n;
  incx_ = (int32_t) incx;
  incy_ = (int32_t) incy;

  /* Check for integer overflows */
  assert ( (int64_t)    n_ == n   );
  assert ( (int64_t) incx_ == incx);
  assert ( (int64_t) incy_ == incy);

  cublasSdot((cublasHandle_t) *handle, n_, x, incx_, y, incy_, result);
}



void gpu_dgemv(void** handle, const char transa, const int64_t m, const int64_t n, const double alpha,
               const double* a, const int64_t lda, const double* x, const int64_t incx, const double beta, double* y, const int64_t incy) {

  assert (*handle != NULL);

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

  cublasOperation_t transa_ = CUBLAS_OP_N;
  if (transa == 'T' || transa == 't') transa_ = CUBLAS_OP_T;

  cublasDgemv((cublasHandle_t) *handle, transa_, m_, n_, &alpha, a, lda_, x, incx_, &beta, y, incy_);
}



void gpu_sgemv(void** handle, const char transa, const int64_t m, const int64_t n, const float alpha,
               const float* a, const int64_t lda, const float* x, const int64_t incx, const float beta, float* y, const int64_t incy) {

  assert (*handle != NULL);

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

  cublasOperation_t transa_ = CUBLAS_OP_N;
  if (transa == 'T' || transa == 't') transa_ = CUBLAS_OP_T;

  cublasSgemv((cublasHandle_t) *handle, transa_, m_, n_, &alpha, a, lda_, x, incx_, &beta, y, incy_);
}


void gpu_dgemm(void** handle, const char transa, const char transb, const int64_t m, const int64_t n, const int64_t k, const double alpha,
               const double* a, const int64_t lda, const double* b, const int64_t ldb, const double beta, double* c, const int64_t ldc) {

  assert (*handle != NULL);

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

  cublasOperation_t transa_ = CUBLAS_OP_N;
  cublasOperation_t transb_ = CUBLAS_OP_N;
  if (transa == 'T' || transa == 't') transa_ = CUBLAS_OP_T;
  if (transb == 'T' || transb == 't') transb_ = CUBLAS_OP_T;

  cublasDgemm((cublasHandle_t) *handle, transa_, transb_, m_, n_, k_, &alpha, a, lda_, b, ldb_, &beta, c, ldc_);
}



void gpu_sgemm(void** handle, const char transa, const char transb, const int64_t m, const int64_t n, const int64_t k, const float alpha,
               const float* a, const int64_t lda, const float* b, const int64_t ldb, const float beta, float* c, const int64_t ldc) {

  assert (*handle != NULL);

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

  cublasOperation_t transa_ = CUBLAS_OP_N;
  cublasOperation_t transb_ = CUBLAS_OP_N;
  if (transa == 'T' || transa == 't') transa_ = CUBLAS_OP_T;
  if (transb == 'T' || transb == 't') transb_ = CUBLAS_OP_T;

  cublasSgemm((cublasHandle_t) *handle, transa_, transb_, m_, n_, k_, &alpha, a, lda_, b, ldb_, &beta, c, ldc_);
}


void gpu_dgeam(void** handle, const char transa, const char transb, const int64_t m, const int64_t n, const double alpha,
               const double* a, const int64_t lda, const double beta, const double* b, const int64_t ldb, double* c, const int64_t ldc) {
  assert (*handle != NULL);

  /* Convert to int32_t */
  int32_t m_, n_, lda_, ldb_, ldc_;

  m_   = (int32_t) m;
  n_   = (int32_t) n;
  lda_ = (int32_t) lda;
  ldb_ = (int32_t) ldb;
  ldc_ = (int32_t) ldc;

  /* Check for integer overflows */
  assert ( (int64_t)   m_ == m  );
  assert ( (int64_t)   n_ == n  );
  assert ( (int64_t) lda_ == lda);
  assert ( (int64_t) ldb_ == ldb);
  assert ( (int64_t) ldc_ == ldc);

  cublasOperation_t transa_ = CUBLAS_OP_N;
  cublasOperation_t transb_ = CUBLAS_OP_N;
  if (transa == 'T' || transa == 't') transa_ = CUBLAS_OP_T;
  if (transb == 'T' || transb == 't') transb_ = CUBLAS_OP_T;

  cublasDgeam((cublasHandle_t) *handle, transa_, transb_, m_, n_, &alpha, a, lda_, &beta, b, ldb_, c, ldc_);

}


void gpu_sgeam(void** handle, const char transa, const char transb, const int64_t m, const int64_t n, const float alpha,
               const float* a, const int64_t lda, const float beta, const float* b, const int64_t ldb, float* c, const int64_t ldc) {
  assert (*handle != NULL);

  /* Convert to int32_t */
  int32_t m_, n_, lda_, ldb_, ldc_;

  m_   = (int32_t) m;
  n_   = (int32_t) n;
  lda_ = (int32_t) lda;
  ldb_ = (int32_t) ldb;
  ldc_ = (int32_t) ldc;

  /* Check for integer overflows */
  assert ( (int64_t)   m_ == m  );
  assert ( (int64_t)   n_ == n  );
  assert ( (int64_t) lda_ == lda);
  assert ( (int64_t) ldb_ == ldb);
  assert ( (int64_t) ldc_ == ldc);

  cublasOperation_t transa_ = CUBLAS_OP_N;
  cublasOperation_t transb_ = CUBLAS_OP_N;
  if (transa == 'T' || transa == 't') transa_ = CUBLAS_OP_T;
  if (transb == 'T' || transb == 't') transb_ = CUBLAS_OP_T;

  cublasSgeam((cublasHandle_t) *handle, transa_, transb_, m_, n_, &alpha, a, lda_, &beta, b, ldb_, c, ldc_);

}
