#include <stdint.h>

int  gpu_ndevices();
void gpu_set_device(int32_t i);

void gpu_allocate(void** ptr, const int64_t n);
void gpu_free(void** ptr);

void gpu_upload(const void* cpu_ptr, void* gpu_ptr, const int64_t n);
void gpu_download(const void* gpu_ptr, void* cpu_ptr, const int64_t n);
void gpu_copy(const void* gpu_ptr_src, void* gpu_ptr_dest, const int64_t n);

void gpu_stream_create(void** ptr);
void gpu_stream_destroy(void** ptr);
void gpu_set_stream(void* handle, void* stream);
void gpu_synchronize();

void gpu_blas_create(void** handle);
void gpu_blas_destroy(void** handle);

void gpu_ddot(const void* handle, const int64_t n, const double* x, const int64_t incx, const double* y, const int64_t incy, double* result);

void gpu_sdot(const void* handle, const int64_t n, const float* x, const int64_t incx, const float* y, const int64_t incy, float* result);

void gpu_dgemv(const void* handle, const char transa, const int64_t m, const int64_t n, const double* alpha,
               const double* a, const int64_t lda, const double* x, const int64_t incx, const double* beta, double* y, const int64_t incy);

void gpu_sgemv(const void* handle, const char transa, const int64_t m, const int64_t n, const float* alpha,
               const float* a, const int64_t lda, const float* x, const int64_t incx, const float* beta, float* y, const int64_t incy);

void gpu_dgemm(const void* handle, const char transa, const char transb, const int64_t m, const int64_t n, const int64_t k, const double* alpha,
               const double* a, const int64_t lda, const double* b, const int64_t ldb, const double* beta, double* c, const int64_t ldc);

void gpu_sgemm(const void* handle, const char transa, const char transb, const int64_t m, const int64_t n, const int64_t k, const float* alpha,
               const float* a, const int64_t lda, const float* b, const int64_t ldb, const float* beta, float* c, const int64_t ldc);

void gpu_dgeam(const void* handle, const char transa, const char transb, const int64_t m, const int64_t n, const double* alpha,
               const double* a, const int64_t lda, const double* beta, const double* b, const int64_t ldb, double* c, const int64_t ldc);

void gpu_sgeam(const void* handle, const char transa, const char transb, const int64_t m, const int64_t n, const float* alpha,
               const float* a, const int64_t lda, const float* beta, const float* b, const int64_t ldb, float* c, const int64_t ldc);
