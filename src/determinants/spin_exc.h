#pragma once
#include <stdint.h>
#include "local_cpu.h"

/* 
 * Memory layout
 * The Fortran buffer is column-major: word k of determinant i lives at
 * buffer[k + N_int*i]. The C code indexes it as buffer[N_int*i + k] —
 * identical in practice, and the compiler sees it as a stride-N_int pattern,
 * which it vectorizes well.
 *
 * __builtin_popcountll for the scalar path. GCC/Clang/ICC all lower this to a
 * single POPCNT instruction on x86-64 with -march=native. No need for -mpopcnt
 * separately.
 *
 * AVX2 path — the Muła nibble-LUT. Since AVX2 has no native 64-bit VPOPCNTQ, the
 * best general approach is the byte-level nibble lookup (VPSHUFB) followed by
 * VPSADBW to accumulate byte counts into per-lane 64-bit sums, then a horizontal
 * add. This is typically 2–4× faster than four scalar POPCNT calls for N≥4
 * because it hides instruction latency and eliminates the scalar-vector
 * round-trips.
 *
 * AVX-512 path — _mm256_popcnt_epi64 / _mm512_popcnt_epi64. On CPUs with
 * AVX512VPOPCNTDQ (Intel Icelake+, AMD Zen4+) this is a single instruction per
 * 256/512-bit register. The N=4 case fits in one 256-bit register and needs only
 * a 4-way horizontal reduce. The general-N path processes 8 words at a time with
 * 512-bit registers.
 *
 */

/* Singles only */
void get_all_spin_singles_1(const uint64_t *buffer, const int *idx,
    uint64_t spindet0, int size_buffer,
    int *singles, int *n_singles);

void get_all_spin_singles_2(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int size_buffer,
    int *singles, int *n_singles);

void get_all_spin_singles_3(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int size_buffer,
    int *singles, int *n_singles);

void get_all_spin_singles_4(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int size_buffer,
    int *singles, int *n_singles);

void get_all_spin_singles_N_int(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int N_int, int size_buffer,
    int *singles, int *n_singles);

void get_all_spin_singles(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int N_int, int size_buffer,
    int *singles, int *n_singles);


/* Doubles only */
void get_all_spin_doubles_1(const uint64_t *buffer, const int *idx,
    uint64_t spindet0, int size_buffer,
    int *doubles, int *n_doubles);

void get_all_spin_doubles_2(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int size_buffer,
    int *doubles, int *n_doubles);

void get_all_spin_doubles_3(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int size_buffer,
    int *doubles, int *n_doubles);

void get_all_spin_doubles_4(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int size_buffer,
    int *doubles, int *n_doubles);

void get_all_spin_doubles_N_int(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int N_int, int size_buffer,
    int *doubles, int *n_doubles);

void get_all_spin_doubles(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int N_int, int size_buffer,
    int *doubles, int *n_doubles);


/* Singles + doubles */
void get_all_spin_singles_and_doubles_1(const uint64_t *buffer, const int *idx,
    uint64_t spindet0, int size_buffer,
    int *singles, int *doubles, int *n_singles, int *n_doubles);

void get_all_spin_singles_and_doubles_2(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int size_buffer,
    int *singles, int *doubles, int *n_singles, int *n_doubles);

void get_all_spin_singles_and_doubles_3(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int size_buffer,
    int *singles, int *doubles, int *n_singles, int *n_doubles);

void get_all_spin_singles_and_doubles_4(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int size_buffer,
    int *singles, int *doubles, int *n_singles, int *n_doubles);

void get_all_spin_singles_and_doubles_N_int(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int N_int, int size_buffer,
    int *singles, int *doubles, int *n_singles, int *n_doubles);

void get_all_spin_singles_and_doubles(const uint64_t *buffer, const int *idx,
    const uint64_t *spindet, int N_int, int size_buffer,
    int *singles, int *doubles, int *n_singles, int *n_doubles);





/* ================================================================== */
/* Scalar popcount helpers                                              */
/* ================================================================== */

static inline int pc1(uint64_t a)
    { return __builtin_popcountll(a); }

static inline int pc2(uint64_t a, uint64_t b)
    { return __builtin_popcountll(a) + __builtin_popcountll(b); }

static inline int pc3(uint64_t a, uint64_t b, uint64_t c)
    { return __builtin_popcountll(a) + __builtin_popcountll(b)
           + __builtin_popcountll(c); }

static inline int pc4(uint64_t a, uint64_t b, uint64_t c, uint64_t d)
    { return __builtin_popcountll(a) + __builtin_popcountll(b)
           + __builtin_popcountll(c) + __builtin_popcountll(d); }


#if HAVE_AVX2 || HAVE_AVX512
#  include <immintrin.h>


/* ================================================================== */
/* AVX2 horizontal popcount via Muła lookup-table method               */
/* ================================================================== */

static inline int popcnt_avx2_256(__m256i v)
{
    /*
     * Classic Muła/Wilkes-Wilson nibble-LUT algorithm.
     * 1. Split each byte into low / high nibble.
     * 2. Lookup bit-count per nibble in a 16-entry table.
     * 3. Sum bytes within each 64-bit lane via PSADBW.
     * 4. Horizontal-add the 4 lane sums.
     */
    const __m256i LO_MASK = _mm256_set1_epi8(0x0f);
    const __m256i LUT     = _mm256_setr_epi8(
        /* popcount of nibbles 0..15, duplicated for both 128-bit lanes */
        0,1,1,2, 1,2,2,3, 1,2,2,3, 2,3,3,4,
        0,1,1,2, 1,2,2,3, 1,2,2,3, 2,3,3,4);

    __m256i lo  = _mm256_and_si256(v, LO_MASK);
    __m256i hi  = _mm256_and_si256(_mm256_srli_epi16(v, 4), LO_MASK);
    __m256i cnt = _mm256_add_epi8(_mm256_shuffle_epi8(LUT, lo),
                                   _mm256_shuffle_epi8(LUT, hi));
    /* sum 8-bit counts into four 64-bit lane sums */
    __m256i sad = _mm256_sad_epu8(cnt, _mm256_setzero_si256());
    /* add the two 128-bit halves, then the two 64-bit lanes */
    __m128i lo128 = _mm256_castsi256_si128(sad);
    __m128i hi128 = _mm256_extracti128_si256(sad, 1);
    __m128i sum   = _mm_add_epi64(lo128, hi128);
    sum = _mm_add_epi64(sum, _mm_shuffle_epi32(sum, _MM_SHUFFLE(1,0,3,2)));
    return (int)_mm_cvtsi128_si64(sum);
}
#endif /* HAVE_AVX2 || HAVE_AVX512 */


/* ================================================================== */
/* NEON helper: popcount of a 128-bit vector (2 x uint64)              */
/*                                                                      */
/* vcntq_u8  counts bits in each byte (16 results).                    */
/* vaddlvq_u8 horizontally adds all 16 bytes into one uint64.          */
/* This is equivalent to popcount(lo) + popcount(hi) in one sequence.  */
/* ================================================================== */

#if HAVE_NEON
static inline int popcnt_neon_128(uint64x2_t v)
{
    /* Reinterpret as bytes, count bits per byte, sum all bytes */
    uint8x16_t bytes = vreinterpretq_u8_u64(v);
    uint8x16_t cnt   = vcntq_u8(bytes);
    return (int)vaddlvq_u8(cnt);   /* vaddlv: unsigned add-long across vector */
}
#endif /* HAVE_NEON */


/* ================================================================== */
/* SVE helper: popcount of an SVE uint64 vector                        */
/*                                                                      */
/* svcnt_u64_z  : per-lane popcount (VCNT on SVE = bit-count per lane) */
/* svaddv_u64   : horizontal reduction to scalar u64                   */
/*                                                                      */
/* On Neoverse V1, svcntd() == 4 (256-bit vectors), so one register    */
/* holds exactly the N=4 case. The general path processes svcntd()      */
/* words per iteration, making it correct on wider SVE implementations. */
/* ================================================================== */

#if HAVE_SVE
static inline int popcnt_sve(svuint64_t v, svbool_t pg)
{
    svuint64_t cnt = svcnt_u64_z(pg, v);
    return (int)svaddv_u64(pg, cnt);
}
#endif /* HAVE_SVE */
