/*
 * For each determinant d in buffer[N_int * size_buffer] (column-major / Fortran order):
 *   degree = popcount( spindet XOR d )   (summed over N_int 64-bit words)
 *   degree == 2  =>  single excitation
 *   degree == 4  =>  double excitation
 *
 * Build flags:
 *   -DUSE_AVX512   VPOPCNTDQ  (Icelake+, Zen4+)   - best on supported CPUs
 *   -DUSE_AVX2     LUT method  (Haswell+)         - good everywhere
 *   (default)      scalar __builtin_popcountll    - universal fallback
 *
 * gcc -O3 -march=native [-DUSE_AVX512|-DUSE_AVX2] spin_singles.c -c
 */

#include "spin_exc.h"


/* ================================================================== */
/* SINGLES ONLY                                                         */
/* ================================================================== */

/* ---- N=1 --------------------------------------------------------- */
void get_all_spin_singles_1(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        uint64_t                  spindet0,
        int                       size_buffer,
        int            * restrict singles,
        int            * restrict n_singles)
{
    int ns = 0;
    for (int i = 0; i < size_buffer; ++i)
        if (pc1(spindet0 ^ buffer[i]) == 2)
            singles[ns++] = idx[i];
    *n_singles = ns;
}

/* ---- N=2 --------------------------------------------------------- */
void get_all_spin_singles_2(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       size_buffer,
        int            * restrict singles,
        int            * restrict n_singles)
{
    const uint64_t s0=spindet[0], s1=spindet[1];
    int ns = 0;
    for (int i = 0; i < size_buffer; ++i)
        if (pc2(s0^buffer[2*i], s1^buffer[2*i+1]) == 2)
            singles[ns++] = idx[i];
    *n_singles = ns;
}

/* ---- N=3 --------------------------------------------------------- */
void get_all_spin_singles_3(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       size_buffer,
        int            * restrict singles,
        int            * restrict n_singles)
{
    const uint64_t s0=spindet[0], s1=spindet[1], s2=spindet[2];
    int ns = 0;
    for (int i = 0; i < size_buffer; ++i)
        if (pc3(s0^buffer[3*i], s1^buffer[3*i+1], s2^buffer[3*i+2]) == 2)
            singles[ns++] = idx[i];
    *n_singles = ns;
}

/* ---- N=4 --------------------------------------------------------- */
void get_all_spin_singles_4(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       size_buffer,
        int            * restrict singles,
        int            * restrict n_singles)
{
    int ns = 0;

#if defined(USE_AVX512)
    /* VPOPCNTDQ gives per-lane 64-bit popcnt directly */
    const __m256i sv = _mm256_loadu_si256((const __m256i *)spindet);
    for (int i = 0; i < size_buffer; ++i) {
        __m256i xv = _mm256_xor_si256(sv,
            _mm256_loadu_si256((const __m256i *)(buffer + 4*i)));
        __m256i pc = _mm256_popcnt_epi64(xv);   /* AVX-512VPOPCNTDQ */
        /* horizontal sum of 4 x 64-bit lanes */
        __m128i lo = _mm256_castsi256_si128(pc);
        __m128i hi = _mm256_extracti128_si256(pc, 1);
        __m128i s  = _mm_add_epi64(lo, hi);
        s = _mm_add_epi64(s, _mm_srli_si128(s, 8));
        if ((int)_mm_cvtsi128_si64(s) == 2) singles[ns++] = idx[i];
    }

#elif defined(USE_AVX2)
    const __m256i sv = _mm256_loadu_si256((const __m256i *)spindet);
    for (int i = 0; i < size_buffer; ++i) {
        __m256i xv = _mm256_xor_si256(sv,
            _mm256_loadu_si256((const __m256i *)(buffer + 4*i)));
        if (popcnt_avx2_256(xv) == 2) singles[ns++] = idx[i];
    }

#else
    const uint64_t s0=spindet[0],s1=spindet[1],s2=spindet[2],s3=spindet[3];
    for (int i = 0; i < size_buffer; ++i)
        if (pc4(s0^buffer[4*i],s1^buffer[4*i+1],
                s2^buffer[4*i+2],s3^buffer[4*i+3]) == 2)
            singles[ns++] = idx[i];
#endif

    *n_singles = ns;
}

/* ---- N=arbitrary ------------------------------------------------- */
void get_all_spin_singles_N_int(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       N_int,
        int                       size_buffer,
        int            * restrict singles,
        int            * restrict n_singles)
{
    int ns = 0;
    for (int i = 0; i < size_buffer; ++i) {
        const uint64_t *col = buffer + N_int * i;
        int deg = 0;

#if defined(USE_AVX512)
        int k = 0;
        for (; k + 8 <= N_int; k += 8) {
            __m512i xv = _mm512_xor_si512(
                _mm512_loadu_si512((const __m512i *)(spindet + k)),
                _mm512_loadu_si512((const __m512i *)(col + k)));
            deg += (int)_mm512_reduce_add_epi64(_mm512_popcnt_epi64(xv));
        }
        for (; k < N_int; ++k) deg += __builtin_popcountll(spindet[k] ^ col[k]);

#elif defined(USE_AVX2)
        int k = 0;
        for (; k + 4 <= N_int; k += 4) {
            __m256i xv = _mm256_xor_si256(
                _mm256_loadu_si256((const __m256i *)(spindet + k)),
                _mm256_loadu_si256((const __m256i *)(col + k)));
            deg += popcnt_avx2_256(xv);
        }
        for (; k < N_int; ++k) deg += __builtin_popcountll(spindet[k] ^ col[k]);

#else
        for (int k = 0; k < N_int; ++k) deg += __builtin_popcountll(spindet[k] ^ col[k]);
#endif

        if (deg == 2) singles[ns++] = idx[i];
    }
    *n_singles = ns;
}

/* ---- Dispatcher -------------------------------------------------- */
void get_all_spin_singles(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int N_int, int size_buffer,
        int * restrict singles, int * restrict n_singles)
{
    switch (N_int) {
        case 1: get_all_spin_singles_1(buffer,idx,spindet[0],size_buffer,singles,n_singles); return;
        case 2: get_all_spin_singles_2(buffer,idx,spindet,size_buffer,singles,n_singles); return;
        case 3: get_all_spin_singles_3(buffer,idx,spindet,size_buffer,singles,n_singles); return;
        case 4: get_all_spin_singles_4(buffer,idx,spindet,size_buffer,singles,n_singles); return;
        default: get_all_spin_singles_N_int(buffer,idx,spindet,N_int,size_buffer,singles,n_singles);
    }
}

/*
 *
 * For each determinant d in buffer[N_int * size_buffer] (column-major / Fortran order):
 *   degree = popcount( spindet XOR d )   (summed over N_int 64-bit words)
 *   degree == 2  =>  double excitation
 *   degree == 4  =>  double excitation
 *
 * Build flags:
 *   -DUSE_AVX512   VPOPCNTDQ  (Icelake+, Zen4+)   - best on supported CPUs
 *   -DUSE_AVX2     LUT method  (Haswell+)         - good everywhere
 *   (default)      scalar __builtin_popcountll    - universal fallback
 *
 * gcc -O3 -march=native [-DUSE_AVX512|-DUSE_AVX2] spin_doubles.c -c
 */

#include "spin_exc.h"


/* ================================================================== */
/* DOUBLES ONLY                                                         */
/* ================================================================== */

/* ---- N=1 --------------------------------------------------------- */
void get_all_spin_doubles_1(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        uint64_t                  spindet0,
        int                       size_buffer,
        int            * restrict doubles,
        int            * restrict n_doubles)
{
    int ns = 0;
    for (int i = 0; i < size_buffer; ++i)
        if (pc1(spindet0 ^ buffer[i]) == 4)
            doubles[ns++] = idx[i];
    *n_doubles = ns;
}

/* ---- N=2 --------------------------------------------------------- */
void get_all_spin_doubles_2(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       size_buffer,
        int            * restrict doubles,
        int            * restrict n_doubles)
{
    const uint64_t s0=spindet[0], s1=spindet[1];
    int ns = 0;
    for (int i = 0; i < size_buffer; ++i)
        if (pc2(s0^buffer[2*i], s1^buffer[2*i+1]) == 4)
            doubles[ns++] = idx[i];
    *n_doubles = ns;
}

/* ---- N=3 --------------------------------------------------------- */
void get_all_spin_doubles_3(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       size_buffer,
        int            * restrict doubles,
        int            * restrict n_doubles)
{
    const uint64_t s0=spindet[0], s1=spindet[1], s2=spindet[2];
    int ns = 0;
    for (int i = 0; i < size_buffer; ++i)
        if (pc3(s0^buffer[3*i], s1^buffer[3*i+1], s2^buffer[3*i+2]) == 4)
            doubles[ns++] = idx[i];
    *n_doubles = ns;
}

/* ---- N=4 --------------------------------------------------------- */
void get_all_spin_doubles_4(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       size_buffer,
        int            * restrict doubles,
        int            * restrict n_doubles)
{
    int ns = 0;

#if defined(USE_AVX512)
    /* VPOPCNTDQ gives per-lane 64-bit popcnt directly */
    const __m256i sv = _mm256_loadu_si256((const __m256i *)spindet);
    for (int i = 0; i < size_buffer; ++i) {
        __m256i xv = _mm256_xor_si256(sv,
            _mm256_loadu_si256((const __m256i *)(buffer + 4*i)));
        __m256i pc = _mm256_popcnt_epi64(xv);   /* AVX-512VPOPCNTDQ */
        /* horizontal sum of 4 x 64-bit lanes */
        __m128i lo = _mm256_castsi256_si128(pc);
        __m128i hi = _mm256_extracti128_si256(pc, 1);
        __m128i s  = _mm_add_epi64(lo, hi);
        s = _mm_add_epi64(s, _mm_srli_si128(s, 8));
        if ((int)_mm_cvtsi128_si64(s) == 4) doubles[ns++] = idx[i];
    }

#elif defined(USE_AVX2)
    const __m256i sv = _mm256_loadu_si256((const __m256i *)spindet);
    for (int i = 0; i < size_buffer; ++i) {
        __m256i xv = _mm256_xor_si256(sv,
            _mm256_loadu_si256((const __m256i *)(buffer + 4*i)));
        if (popcnt_avx2_256(xv) == 4) doubles[ns++] = idx[i];
    }

#else
    const uint64_t s0=spindet[0],s1=spindet[1],s2=spindet[2],s3=spindet[3];
    for (int i = 0; i < size_buffer; ++i)
        if (pc4(s0^buffer[4*i],s1^buffer[4*i+1],
                s2^buffer[4*i+2],s3^buffer[4*i+3]) == 4)
            doubles[ns++] = idx[i];
#endif

    *n_doubles = ns;
}

/* ---- N=arbitrary ------------------------------------------------- */
void get_all_spin_doubles_N_int(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       N_int,
        int                       size_buffer,
        int            * restrict doubles,
        int            * restrict n_doubles)
{
    int ns = 0;
    for (int i = 0; i < size_buffer; ++i) {
        const uint64_t *col = buffer + N_int * i;
        int deg = 0;

#if defined(USE_AVX512)
        int k = 0;
        for (; k + 8 <= N_int; k += 8) {
            __m512i xv = _mm512_xor_si512(
                _mm512_loadu_si512((const __m512i *)(spindet + k)),
                _mm512_loadu_si512((const __m512i *)(col + k)));
            deg += (int)_mm512_reduce_add_epi64(_mm512_popcnt_epi64(xv));
        }
        for (; k < N_int; ++k) deg += __builtin_popcountll(spindet[k] ^ col[k]);

#elif defined(USE_AVX2)
        int k = 0;
        for (; k + 4 <= N_int; k += 4) {
            __m256i xv = _mm256_xor_si256(
                _mm256_loadu_si256((const __m256i *)(spindet + k)),
                _mm256_loadu_si256((const __m256i *)(col + k)));
            deg += popcnt_avx2_256(xv);
        }
        for (; k < N_int; ++k) deg += __builtin_popcountll(spindet[k] ^ col[k]);

#else
        for (int k = 0; k < N_int; ++k) deg += __builtin_popcountll(spindet[k] ^ col[k]);
#endif

        if (deg == 4) doubles[ns++] = idx[i];
    }
    *n_doubles = ns;
}

/* ---- Dispatcher -------------------------------------------------- */
void get_all_spin_doubles(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int N_int, int size_buffer,
        int * restrict doubles, int * restrict n_doubles)
{
    switch (N_int) {
        case 1: get_all_spin_doubles_1(buffer,idx,spindet[0],size_buffer,doubles,n_doubles); return;
        case 2: get_all_spin_doubles_2(buffer,idx,spindet,size_buffer,doubles,n_doubles); return;
        case 3: get_all_spin_doubles_3(buffer,idx,spindet,size_buffer,doubles,n_doubles); return;
        case 4: get_all_spin_doubles_4(buffer,idx,spindet,size_buffer,doubles,n_doubles); return;
        default: get_all_spin_doubles_N_int(buffer,idx,spindet,N_int,size_buffer,doubles,n_doubles);
    }
}

/*
 *
 * For each determinant d in buffer[N_int * size_buffer] (column-major / Fortran order):
 *   degree = popcount( spindet XOR d )   (summed over N_int 64-bit words)
 *   degree == 2  =>  single excitation
 *   degree == 4  =>  double excitation
 *
 * Build flags:
 *   -DUSE_AVX512   VPOPCNTDQ  (Icelake+, Zen4+)   - best on supported CPUs
 *   -DUSE_AVX2     LUT method  (Haswell+)         - good everywhere
 *   (default)      scalar __builtin_popcountll    - universal fallback
 *
 * gcc -O3 -march=native [-DUSE_AVX512|-DUSE_AVX2] spin_singles.c -c
 */

#include "spin_exc.h"
#include <stdint.h>


/* ================================================================== */
/* SINGLES + DOUBLES                                                    */
/* ================================================================== */

/* Shared inner-loop macro to avoid copy-paste between N=1..4 variants */
#define CLASSIFY(deg)                                  \
    do {                                               \
        if      ((deg) == 2) singles[ns++] = idx[i];  \
        else if ((deg) == 4) doubles[nd++] = idx[i];  \
    } while (0)

/* ---- N=1 --------------------------------------------------------- */
void get_all_spin_singles_and_doubles_1(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        uint64_t                  spindet0,
        int                       size_buffer,
        int * restrict singles, int * restrict doubles,
        int * restrict n_singles, int * restrict n_doubles)
{
    int ns=0, nd=0;
    for (int i=0; i<size_buffer; ++i) { CLASSIFY(pc1(spindet0^buffer[i])); }
    *n_singles=ns; *n_doubles=nd;
}

/* ---- N=2 --------------------------------------------------------- */
void get_all_spin_singles_and_doubles_2(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       size_buffer,
        int * restrict singles, int * restrict doubles,
        int * restrict n_singles, int * restrict n_doubles)
{
    const uint64_t s0=spindet[0], s1=spindet[1];
    int ns=0, nd=0;
    for (int i=0; i<size_buffer; ++i) { CLASSIFY(pc2(s0^buffer[2*i], s1^buffer[2*i+1])); }
    *n_singles=ns; *n_doubles=nd;
}

/* ---- N=3 --------------------------------------------------------- */
void get_all_spin_singles_and_doubles_3(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       size_buffer,
        int * restrict singles, int * restrict doubles,
        int * restrict n_singles, int * restrict n_doubles)
{
    const uint64_t s0=spindet[0], s1=spindet[1], s2=spindet[2];
    int ns=0, nd=0;
    for (int i=0; i<size_buffer; ++i)
        { CLASSIFY(pc3(s0^buffer[3*i], s1^buffer[3*i+1], s2^buffer[3*i+2])); }
    *n_singles=ns; *n_doubles=nd;
}

/* ---- N=4 --------------------------------------------------------- */
void get_all_spin_singles_and_doubles_4(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       size_buffer,
        int * restrict singles, int * restrict doubles,
        int * restrict n_singles, int * restrict n_doubles)
{
    int ns=0, nd=0;

#if defined(USE_AVX512)
    const __m256i sv = _mm256_loadu_si256((const __m256i *)spindet);
    for (int i=0; i<size_buffer; ++i) {
        __m256i xv = _mm256_xor_si256(sv,
            _mm256_loadu_si256((const __m256i *)(buffer + 4*i)));
        __m256i pc = _mm256_popcnt_epi64(xv);
        __m128i lo = _mm256_castsi256_si128(pc);
        __m128i hi = _mm256_extracti128_si256(pc, 1);
        __m128i s  = _mm_add_epi64(lo, hi);
        s = _mm_add_epi64(s, _mm_srli_si128(s, 8));
        CLASSIFY((int)_mm_cvtsi128_si64(s));
    }

#elif defined(USE_AVX2)
    const __m256i sv = _mm256_loadu_si256((const __m256i *)spindet);
    for (int i=0; i<size_buffer; ++i) {
        __m256i xv = _mm256_xor_si256(sv,
            _mm256_loadu_si256((const __m256i *)(buffer + 4*i)));
        CLASSIFY(popcnt_avx2_256(xv));
    }

#else
    const uint64_t s0=spindet[0],s1=spindet[1],s2=spindet[2],s3=spindet[3];
    for (int i=0; i<size_buffer; ++i)
        CLASSIFY(pc4(s0^buffer[4*i],s1^buffer[4*i+1],
                     s2^buffer[4*i+2],s3^buffer[4*i+3]));
#endif

    *n_singles=ns; *n_doubles=nd;
}

/* ---- N=arbitrary ------------------------------------------------- */
void get_all_spin_singles_and_doubles_N_int(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       N_int,
        int                       size_buffer,
        int * restrict singles, int * restrict doubles,
        int * restrict n_singles, int * restrict n_doubles)
{
    int ns=0, nd=0;
    for (int i=0; i<size_buffer; ++i) {
        const uint64_t *col = buffer + N_int * i;
        int deg = 0;

#if defined(USE_AVX512)
        int k=0;
        for (; k+8<=N_int; k+=8) {
            __m512i xv = _mm512_xor_si512(
                _mm512_loadu_si512((const __m512i *)(spindet+k)),
                _mm512_loadu_si512((const __m512i *)(col+k)));
            deg += (int)_mm512_reduce_add_epi64(_mm512_popcnt_epi64(xv));
        }
        for (; k<N_int; ++k) deg += __builtin_popcountll(spindet[k]^col[k]);

#elif defined(USE_AVX2)
        int k=0;
        for (; k+4<=N_int; k+=4) {
            __m256i xv = _mm256_xor_si256(
                _mm256_loadu_si256((const __m256i *)(spindet+k)),
                _mm256_loadu_si256((const __m256i *)(col+k)));
            deg += popcnt_avx2_256(xv);
        }
        for (; k<N_int; ++k) deg += __builtin_popcountll(spindet[k]^col[k]);

#else
        for (int k=0; k<N_int; ++k) deg += __builtin_popcountll(spindet[k]^col[k]);
#endif

        CLASSIFY(deg);
    }
    *n_singles=ns; *n_doubles=nd;
}

/* ---- Dispatcher -------------------------------------------------- */
void get_all_spin_singles_and_doubles(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int N_int, int size_buffer,
        int * restrict singles, int * restrict doubles,
        int * restrict n_singles, int * restrict n_doubles)
{
    switch (N_int) {
        case 1: get_all_spin_singles_and_doubles_1(buffer,idx,spindet[0],size_buffer,singles,doubles,n_singles,n_doubles); return;
        case 2: get_all_spin_singles_and_doubles_2(buffer,idx,spindet,size_buffer,singles,doubles,n_singles,n_doubles); return;
        case 3: get_all_spin_singles_and_doubles_3(buffer,idx,spindet,size_buffer,singles,doubles,n_singles,n_doubles); return;
        case 4: get_all_spin_singles_and_doubles_4(buffer,idx,spindet,size_buffer,singles,doubles,n_singles,n_doubles); return;
        default: get_all_spin_singles_and_doubles_N_int(buffer,idx,spindet,N_int,size_buffer,singles,doubles,n_singles,n_doubles);
    }
}
