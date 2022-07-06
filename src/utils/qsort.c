/* [[file:~/qp2/src/utils/qsort.org::*Generated%20C%20file][Generated C file:1]] */
#include <stdlib.h>
#include <stdint.h>

struct int16_t_comp {
  int16_t    x;
  int32_t i;
};

int compare_int16_t( const void * l, const void * r )
{
  const int16_t * restrict _l= l;
  const int16_t * restrict _r= r;
  if( *_l > *_r ) return 1;
  if( *_l < *_r ) return -1;
  return 0;
}

void qsort_int16_t(int16_t* restrict A_in, int32_t* restrict iorder, int32_t isize) {
  struct int16_t_comp* A = malloc(isize * sizeof(struct int16_t_comp));
  if (A == NULL) return;

  for (int i=0 ; i<isize ; ++i) {
    A[i].x = A_in[i];
    A[i].i = iorder[i];
  }

  qsort( (void*) A, (size_t) isize, sizeof(struct int16_t_comp), compare_int16_t);

  for (int i=0 ; i<isize ; ++i) {
    A_in[i] = A[i].x;
    iorder[i] = A[i].i;
  }
  free(A);
}

void qsort_int16_t_noidx(int16_t* A, int32_t isize) {
  qsort( (void*) A, (size_t) isize, sizeof(int16_t), compare_int16_t);
}


struct int16_t_comp_big {
  int16_t    x;
  int64_t i;
};

int compare_int16_t_big( const void * l, const void * r )
{
  const int16_t * restrict _l= l;
  const int16_t * restrict _r= r;
  if( *_l > *_r ) return 1;
  if( *_l < *_r ) return -1;
  return 0;
}

void qsort_int16_t_big(int16_t* restrict A_in, int64_t* restrict iorder, int64_t isize) {
  struct int16_t_comp_big* A = malloc(isize * sizeof(struct int16_t_comp_big));
  if (A == NULL) return;

  for (int i=0 ; i<isize ; ++i) {
    A[i].x = A_in[i];
    A[i].i = iorder[i];
  }

  qsort( (void*) A, (size_t) isize, sizeof(struct int16_t_comp_big), compare_int16_t_big);

  for (int i=0 ; i<isize ; ++i) {
    A_in[i] = A[i].x;
    iorder[i] = A[i].i;
  }
  free(A);
}

void qsort_int16_t_noidx_big(int16_t* A, int64_t isize) {
  qsort( (void*) A, (size_t) isize, sizeof(int16_t), compare_int16_t_big);
}


struct int32_t_comp {
  int32_t    x;
  int32_t i;
};

int compare_int32_t( const void * l, const void * r )
{
  const int32_t * restrict _l= l;
  const int32_t * restrict _r= r;
  if( *_l > *_r ) return 1;
  if( *_l < *_r ) return -1;
  return 0;
}

void qsort_int32_t(int32_t* restrict A_in, int32_t* restrict iorder, int32_t isize) {
  struct int32_t_comp* A = malloc(isize * sizeof(struct int32_t_comp));
  if (A == NULL) return;

  for (int i=0 ; i<isize ; ++i) {
    A[i].x = A_in[i];
    A[i].i = iorder[i];
  }

  qsort( (void*) A, (size_t) isize, sizeof(struct int32_t_comp), compare_int32_t);

  for (int i=0 ; i<isize ; ++i) {
    A_in[i] = A[i].x;
    iorder[i] = A[i].i;
  }
  free(A);
}

void qsort_int32_t_noidx(int32_t* A, int32_t isize) {
  qsort( (void*) A, (size_t) isize, sizeof(int32_t), compare_int32_t);
}


struct int32_t_comp_big {
  int32_t    x;
  int64_t i;
};

int compare_int32_t_big( const void * l, const void * r )
{
  const int32_t * restrict _l= l;
  const int32_t * restrict _r= r;
  if( *_l > *_r ) return 1;
  if( *_l < *_r ) return -1;
  return 0;
}

void qsort_int32_t_big(int32_t* restrict A_in, int64_t* restrict iorder, int64_t isize) {
  struct int32_t_comp_big* A = malloc(isize * sizeof(struct int32_t_comp_big));
  if (A == NULL) return;

  for (int i=0 ; i<isize ; ++i) {
    A[i].x = A_in[i];
    A[i].i = iorder[i];
  }

  qsort( (void*) A, (size_t) isize, sizeof(struct int32_t_comp_big), compare_int32_t_big);

  for (int i=0 ; i<isize ; ++i) {
    A_in[i] = A[i].x;
    iorder[i] = A[i].i;
  }
  free(A);
}

void qsort_int32_t_noidx_big(int32_t* A, int64_t isize) {
  qsort( (void*) A, (size_t) isize, sizeof(int32_t), compare_int32_t_big);
}


struct int64_t_comp {
  int64_t    x;
  int32_t i;
};

int compare_int64_t( const void * l, const void * r )
{
  const int64_t * restrict _l= l;
  const int64_t * restrict _r= r;
  if( *_l > *_r ) return 1;
  if( *_l < *_r ) return -1;
  return 0;
}

void qsort_int64_t(int64_t* restrict A_in, int32_t* restrict iorder, int32_t isize) {
  struct int64_t_comp* A = malloc(isize * sizeof(struct int64_t_comp));
  if (A == NULL) return;

  for (int i=0 ; i<isize ; ++i) {
    A[i].x = A_in[i];
    A[i].i = iorder[i];
  }

  qsort( (void*) A, (size_t) isize, sizeof(struct int64_t_comp), compare_int64_t);

  for (int i=0 ; i<isize ; ++i) {
    A_in[i] = A[i].x;
    iorder[i] = A[i].i;
  }
  free(A);
}

void qsort_int64_t_noidx(int64_t* A, int32_t isize) {
  qsort( (void*) A, (size_t) isize, sizeof(int64_t), compare_int64_t);
}


struct int64_t_comp_big {
  int64_t    x;
  int64_t i;
};

int compare_int64_t_big( const void * l, const void * r )
{
  const int64_t * restrict _l= l;
  const int64_t * restrict _r= r;
  if( *_l > *_r ) return 1;
  if( *_l < *_r ) return -1;
  return 0;
}

void qsort_int64_t_big(int64_t* restrict A_in, int64_t* restrict iorder, int64_t isize) {
  struct int64_t_comp_big* A = malloc(isize * sizeof(struct int64_t_comp_big));
  if (A == NULL) return;

  for (int i=0 ; i<isize ; ++i) {
    A[i].x = A_in[i];
    A[i].i = iorder[i];
  }

  qsort( (void*) A, (size_t) isize, sizeof(struct int64_t_comp_big), compare_int64_t_big);

  for (int i=0 ; i<isize ; ++i) {
    A_in[i] = A[i].x;
    iorder[i] = A[i].i;
  }
  free(A);
}

void qsort_int64_t_noidx_big(int64_t* A, int64_t isize) {
  qsort( (void*) A, (size_t) isize, sizeof(int64_t), compare_int64_t_big);
}


struct double_comp {
  double    x;
  int32_t i;
};

int compare_double( const void * l, const void * r )
{
  const double * restrict _l= l;
  const double * restrict _r= r;
  if( *_l > *_r ) return 1;
  if( *_l < *_r ) return -1;
  return 0;
}

void qsort_double(double* restrict A_in, int32_t* restrict iorder, int32_t isize) {
  struct double_comp* A = malloc(isize * sizeof(struct double_comp));
  if (A == NULL) return;

  for (int i=0 ; i<isize ; ++i) {
    A[i].x = A_in[i];
    A[i].i = iorder[i];
  }

  qsort( (void*) A, (size_t) isize, sizeof(struct double_comp), compare_double);

  for (int i=0 ; i<isize ; ++i) {
    A_in[i] = A[i].x;
    iorder[i] = A[i].i;
  }
  free(A);
}

void qsort_double_noidx(double* A, int32_t isize) {
  qsort( (void*) A, (size_t) isize, sizeof(double), compare_double);
}


struct double_comp_big {
  double    x;
  int64_t i;
};

int compare_double_big( const void * l, const void * r )
{
  const double * restrict _l= l;
  const double * restrict _r= r;
  if( *_l > *_r ) return 1;
  if( *_l < *_r ) return -1;
  return 0;
}

void qsort_double_big(double* restrict A_in, int64_t* restrict iorder, int64_t isize) {
  struct double_comp_big* A = malloc(isize * sizeof(struct double_comp_big));
  if (A == NULL) return;

  for (int i=0 ; i<isize ; ++i) {
    A[i].x = A_in[i];
    A[i].i = iorder[i];
  }

  qsort( (void*) A, (size_t) isize, sizeof(struct double_comp_big), compare_double_big);

  for (int i=0 ; i<isize ; ++i) {
    A_in[i] = A[i].x;
    iorder[i] = A[i].i;
  }
  free(A);
}

void qsort_double_noidx_big(double* A, int64_t isize) {
  qsort( (void*) A, (size_t) isize, sizeof(double), compare_double_big);
}


struct float_comp {
  float    x;
  int32_t i;
};

int compare_float( const void * l, const void * r )
{
  const float * restrict _l= l;
  const float * restrict _r= r;
  if( *_l > *_r ) return 1;
  if( *_l < *_r ) return -1;
  return 0;
}

void qsort_float(float* restrict A_in, int32_t* restrict iorder, int32_t isize) {
  struct float_comp* A = malloc(isize * sizeof(struct float_comp));
  if (A == NULL) return;

  for (int i=0 ; i<isize ; ++i) {
    A[i].x = A_in[i];
    A[i].i = iorder[i];
  }

  qsort( (void*) A, (size_t) isize, sizeof(struct float_comp), compare_float);

  for (int i=0 ; i<isize ; ++i) {
    A_in[i] = A[i].x;
    iorder[i] = A[i].i;
  }
  free(A);
}

void qsort_float_noidx(float* A, int32_t isize) {
  qsort( (void*) A, (size_t) isize, sizeof(float), compare_float);
}


struct float_comp_big {
  float    x;
  int64_t i;
};

int compare_float_big( const void * l, const void * r )
{
  const float * restrict _l= l;
  const float * restrict _r= r;
  if( *_l > *_r ) return 1;
  if( *_l < *_r ) return -1;
  return 0;
}

void qsort_float_big(float* restrict A_in, int64_t* restrict iorder, int64_t isize) {
  struct float_comp_big* A = malloc(isize * sizeof(struct float_comp_big));
  if (A == NULL) return;

  for (int i=0 ; i<isize ; ++i) {
    A[i].x = A_in[i];
    A[i].i = iorder[i];
  }

  qsort( (void*) A, (size_t) isize, sizeof(struct float_comp_big), compare_float_big);

  for (int i=0 ; i<isize ; ++i) {
    A_in[i] = A[i].x;
    iorder[i] = A[i].i;
  }
  free(A);
}

void qsort_float_noidx_big(float* A, int64_t isize) {
  qsort( (void*) A, (size_t) isize, sizeof(float), compare_float_big);
}
/* Generated C file:1 ends here */
