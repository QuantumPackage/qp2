#pragma once


#if defined(__ARM_FEATURE_SVE)
#  include <arm_sve.h>
#  define USE_SVE
#  warning "SVE detected"
#elif defined(__ARM_NEON)
#  include <arm_neon.h>
#  define USE_NEON
#  warning "Neon detected"
#elif defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__)
#  define USE_AVX512
#  warning "AVX512 detected"
#elif defined(__AVX2__)
#  define USE_AVX2
#  warning "AVX2 detected"
#endif

#ifdef USE_NEON
#define HAVE_NEON 1
#else
#define HAVE_NEON 0
#endif

#ifdef USE_SVE
#define HAVE_SVE  1
#else
#define HAVE_SVE  0
#endif

#ifdef USE_AVX2
#define HAVE_AVX2 1
#else
#define HAVE_AVX2 0
#endif

#ifdef USE_AVX512
#define HAVE_AVX512 1
#else
#define HAVE_AVX512 0
#endif

