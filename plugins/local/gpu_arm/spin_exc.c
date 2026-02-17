/*
 * spin_singles_arm.c
 *
 * ARM Neoverse V1 implementation of the spin excitation detection functions.
 *
 * Neoverse V1 capabilities exploited:
 *
 *   NEON (always available, 128-bit / 2 x uint64):
 *     - vcntq_u8   : byte-level popcount, then vaddlvq_u8 to reduce
 *     - Used for N=2 fixed path and as the 2-lane building block for N=3
 *
 *   SVE  (-march=neoverse-v1 enables SVE with 256-bit vectors = 4 x uint64):
 *     - svcnt_u64  : per-lane 64-bit popcount in one instruction
 *     - svaddv_u64 : horizontal reduction to scalar
 *     - svcntd()   : runtime lane count (= 4 on V1, but written portably)
 *     - Used for N=4 fixed path and the general N_int path
 *
 * Build:
 *   gcc  -O3 -march=neoverse-v1 spin_singles_arm.c -c -o spin_singles_arm.o
 *   clang -O3 -march=neoverse-v1 spin_singles_arm.c -c -o spin_singles_arm.o
 *
 * Alternatively, choose the ISA level explicitly:
 *   gcc  -O3 -march=armv8.4-a+sve+sve2+profile spin_singles_arm.c -c
 *
 * Linking: replace spin_singles.c with spin_singles_arm.c in your build.
 * The public API (function names and signatures) is identical.
 *
 * Fallback: if neither USE_SVE nor USE_NEON is defined, the file compiles
 * to portable scalar code using __builtin_popcountll (which GCC/Clang lower
 * to the AArch64 CNT+FMOV pair automatically).
 */

#include "spin_exc.h"

#include <stdint.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Feature selection                                                    */
/* ------------------------------------------------------------------ */

/*
 * -march=neoverse-v1 (or armv8.4-a+sve) defines __ARM_FEATURE_SVE.
 * -march=armv8-a or later defines __ARM_NEON.
 * Both imply __ARM_ARCH_ISA_A64 which gives __builtin_popcountll -> CNT.
 */

#if defined(__ARM_FEATURE_SVE)
#  include <arm_sve.h>
#  define HAVE_SVE  1
#else
#  define HAVE_SVE  0
#endif

#if defined(__ARM_NEON)
#  include <arm_neon.h>
#  define HAVE_NEON 1
#else
#  define HAVE_NEON 0
#endif


/* SINGLES */

/* ================================================================== */
/* N=1 : scalar (no SIMD benefit for a single 64-bit word)             */
/* ================================================================== */

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


/* ================================================================== */
/* N=2 : NEON path (128-bit, 2 x uint64)                               */
/* ================================================================== */

void get_all_spin_singles_2(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       size_buffer,
        int            * restrict singles,
        int            * restrict n_singles)
{
    int ns = 0;

#if HAVE_NEON
    /*
     * Load spindet into a 128-bit NEON register once.
     * For each determinant, load 2 x uint64, XOR, popcount.
     */
    const uint64x2_t sv = vld1q_u64(spindet);
    for (int i = 0; i < size_buffer; ++i) {
        uint64x2_t bv = vld1q_u64(buffer + 2*i);
        uint64x2_t xv = veorq_u64(sv, bv);
        if (popcnt_neon_128(xv) == 2)
            singles[ns++] = idx[i];
    }
#else
    const uint64_t s0=spindet[0], s1=spindet[1];
    for (int i = 0; i < size_buffer; ++i)
        if (pc2(s0^buffer[2*i], s1^buffer[2*i+1]) == 2)
            singles[ns++] = idx[i];
#endif

    *n_singles = ns;
}


/* ================================================================== */
/* N=3 : NEON path (128-bit for 2 words + scalar for the third)        */
/*                                                                      */
/* There is no clean 3-word SIMD container. The best approach on NEON  */
/* is to process words 0-1 in a 128-bit register and word 2 scalarly.  */
/* On SVE with svcntd()>=4 we could fit all 3 words, but the overhead  */
/* of predicated loads for an odd count is not worth it for N=3.       */
/* ================================================================== */

void get_all_spin_singles_3(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       size_buffer,
        int            * restrict singles,
        int            * restrict n_singles)
{
    int ns = 0;

#if HAVE_NEON
    const uint64x2_t sv01 = vld1q_u64(spindet);        /* words 0,1 */
    const uint64_t   s2   = spindet[2];                 /* word  2   */
    for (int i = 0; i < size_buffer; ++i) {
        uint64x2_t bv01 = vld1q_u64(buffer + 3*i);
        uint64x2_t xv01 = veorq_u64(sv01, bv01);
        int deg = popcnt_neon_128(xv01)
                + pc1(s2 ^ buffer[3*i + 2]);
        if (deg == 2) singles[ns++] = idx[i];
    }
#else
    const uint64_t s0=spindet[0], s1=spindet[1], s2=spindet[2];
    for (int i = 0; i < size_buffer; ++i)
        if (pc3(s0^buffer[3*i], s1^buffer[3*i+1], s2^buffer[3*i+2]) == 2)
            singles[ns++] = idx[i];
#endif

    *n_singles = ns;
}


/* ================================================================== */
/* N=4 : SVE path (256-bit on V1, exactly 4 x uint64)                  */
/* ================================================================== */

void get_all_spin_singles_4(
        const uint64_t * restrict buffer,
        const int      * restrict idx,
        const uint64_t * restrict spindet,
        int                       size_buffer,
        int            * restrict singles,
        int            * restrict n_singles)
{
    int ns = 0;

#if HAVE_SVE
    /*
     * On Neoverse V1, svcntd() == 4, so svld1_u64 loads all 4 words
     * of spindet into a single SVE register.
     *
     * svptrue_b64() activates all lanes. The result is a single
     * svcnt_u64 + svaddv_u64 sequence per determinant.
     */
    const svbool_t   pg = svptrue_b64();
    const svuint64_t sv = svld1_u64(pg, spindet);
    for (int i = 0; i < size_buffer; ++i) {
        svuint64_t bv  = svld1_u64(pg, buffer + 4*i);
        svuint64_t xv  = sveor_u64_z(pg, sv, bv);
        if (popcnt_sve(xv, pg) == 2)
            singles[ns++] = idx[i];
    }

#elif HAVE_NEON
    /*
     * NEON fallback for N=4: two 128-bit registers, add their
     * popcounts.
     */
    const uint64x2_t sv0 = vld1q_u64(spindet);
    const uint64x2_t sv1 = vld1q_u64(spindet + 2);
    for (int i = 0; i < size_buffer; ++i) {
        uint64x2_t bv0 = vld1q_u64(buffer + 4*i);
        uint64x2_t bv1 = vld1q_u64(buffer + 4*i + 2);
        int deg = popcnt_neon_128(veorq_u64(sv0, bv0))
                + popcnt_neon_128(veorq_u64(sv1, bv1));
        if (deg == 2) singles[ns++] = idx[i];
    }

#else
    const uint64_t s0=spindet[0],s1=spindet[1],s2=spindet[2],s3=spindet[3];
    for (int i = 0; i < size_buffer; ++i)
        if (pc4(s0^buffer[4*i], s1^buffer[4*i+1],
                s2^buffer[4*i+2], s3^buffer[4*i+3]) == 2)
            singles[ns++] = idx[i];
#endif

    *n_singles = ns;
}


/* ================================================================== */
/* N=arbitrary : SVE path processing svcntd() words per iteration      */
/* ================================================================== */

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

#if HAVE_SVE
    /*
     * sve_step is the number of uint64 lanes per SVE register.
     * On Neoverse V1 this is 4; the code is correct for any width.
     *
     * For the tail (N_int % sve_step != 0) we use a predicated load
     * with svwhilelt_b64 which sets exactly the remaining lanes active.
     */
    const int sve_step = (int)svcntd();

    for (int i = 0; i < size_buffer; ++i) {
        const uint64_t *col = buffer + N_int * i;
        int deg = 0;
        int k = 0;
        /* Full SVE-width chunks */
        for (; k + sve_step <= N_int; k += sve_step) {
            svbool_t   pg = svptrue_b64();
            svuint64_t sv = svld1_u64(pg, spindet + k);
            svuint64_t bv = svld1_u64(pg, col + k);
            deg += popcnt_sve(sveor_u64_z(pg, sv, bv), pg);
        }
        /* Tail: fewer than sve_step words remain */
        if (k < N_int) {
            svbool_t   pg = svwhilelt_b64((uint64_t)k, (uint64_t)N_int);
            svuint64_t sv = svld1_u64(pg, spindet + k);
            svuint64_t bv = svld1_u64(pg, col + k);
            deg += popcnt_sve(sveor_u64_z(pg, sv, bv), pg);
        }
        if (deg == 2) singles[ns++] = idx[i];
    }

#elif HAVE_NEON
    const int neon_step = 2; /* 2 x uint64 per NEON register */
    for (int i = 0; i < size_buffer; ++i) {
        const uint64_t *col = buffer + N_int * i;
        int deg = 0;
        int k = 0;
        for (; k + neon_step <= N_int; k += neon_step) {
            uint64x2_t sv = vld1q_u64(spindet + k);
            uint64x2_t bv = vld1q_u64(col + k);
            deg += popcnt_neon_128(veorq_u64(sv, bv));
        }
        if (k < N_int)   /* one word tail */
            deg += pc1(spindet[k] ^ col[k]);
        if (deg == 2) singles[ns++] = idx[i];
    }

#else
    for (int i = 0; i < size_buffer; ++i) {
        const uint64_t *col = buffer + N_int * i;
        int deg = 0;
        for (int k = 0; k < N_int; ++k)
            deg += __builtin_popcountll(spindet[k] ^ col[k]);
        if (deg == 2) singles[ns++] = idx[i];
    }
#endif

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



/* ================================================================== */
/* SINGLES + DOUBLES                                                    */
/*                                                                      */
/* Same SIMD structure as singles-only; just classify against both      */
/* thresholds. The CLASSIFY macro updates both output arrays.           */
/* ================================================================== */

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
    for (int i=0; i<size_buffer; ++i) CLASSIFY(pc1(spindet0 ^ buffer[i]));
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
    int ns=0, nd=0;

#if HAVE_NEON
    const uint64x2_t sv = vld1q_u64(spindet);
    for (int i=0; i<size_buffer; ++i) {
        uint64x2_t bv = vld1q_u64(buffer + 2*i);
        CLASSIFY(popcnt_neon_128(veorq_u64(sv, bv)));
    }
#else
    const uint64_t s0=spindet[0], s1=spindet[1];
    for (int i=0; i<size_buffer; ++i)
        CLASSIFY(pc2(s0^buffer[2*i], s1^buffer[2*i+1]));
#endif

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
    int ns=0, nd=0;

#if HAVE_NEON
    const uint64x2_t sv01 = vld1q_u64(spindet);
    const uint64_t   s2   = spindet[2];
    for (int i=0; i<size_buffer; ++i) {
        uint64x2_t bv01 = vld1q_u64(buffer + 3*i);
        int deg = popcnt_neon_128(veorq_u64(sv01, bv01))
                + pc1(s2 ^ buffer[3*i+2]);
        CLASSIFY(deg);
    }
#else
    const uint64_t s0=spindet[0], s1=spindet[1], s2=spindet[2];
    for (int i=0; i<size_buffer; ++i)
        CLASSIFY(pc3(s0^buffer[3*i], s1^buffer[3*i+1], s2^buffer[3*i+2]));
#endif

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

#if HAVE_SVE
    const svbool_t   pg = svptrue_b64();
    const svuint64_t sv = svld1_u64(pg, spindet);
    for (int i=0; i<size_buffer; ++i) {
        svuint64_t bv = svld1_u64(pg, buffer + 4*i);
        CLASSIFY(popcnt_sve(sveor_u64_z(pg, sv, bv), pg));
    }

#elif HAVE_NEON
    const uint64x2_t sv0 = vld1q_u64(spindet);
    const uint64x2_t sv1 = vld1q_u64(spindet + 2);
    for (int i=0; i<size_buffer; ++i) {
        int deg = popcnt_neon_128(veorq_u64(sv0, vld1q_u64(buffer + 4*i)))
                + popcnt_neon_128(veorq_u64(sv1, vld1q_u64(buffer + 4*i + 2)));
        CLASSIFY(deg);
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

#if HAVE_SVE
    const int sve_step = (int)svcntd();
    for (int i=0; i<size_buffer; ++i) {
        const uint64_t *col = buffer + N_int * i;
        int deg = 0, k = 0;
        for (; k + sve_step <= N_int; k += sve_step) {
            svbool_t   pg = svptrue_b64();
            svuint64_t xv = sveor_u64_z(pg, svld1_u64(pg, spindet+k),
                                             svld1_u64(pg, col+k));
            deg += popcnt_sve(xv, pg);
        }
        if (k < N_int) {
            svbool_t   pg = svwhilelt_b64((uint64_t)k, (uint64_t)N_int);
            svuint64_t xv = sveor_u64_z(pg, svld1_u64(pg, spindet+k),
                                             svld1_u64(pg, col+k));
            deg += popcnt_sve(xv, pg);
        }
        CLASSIFY(deg);
    }

#elif HAVE_NEON
    for (int i=0; i<size_buffer; ++i) {
        const uint64_t *col = buffer + N_int * i;
        int deg = 0, k = 0;
        for (; k + 2 <= N_int; k += 2)
            deg += popcnt_neon_128(veorq_u64(vld1q_u64(spindet+k),
                                             vld1q_u64(col+k)));
        if (k < N_int)
            deg += pc1(spindet[k] ^ col[k]);
        CLASSIFY(deg);
    }

#else
    for (int i=0; i<size_buffer; ++i) {
        const uint64_t *col = buffer + N_int * i;
        int deg = 0;
        for (int k=0; k<N_int; ++k)
            deg += __builtin_popcountll(spindet[k] ^ col[k]);
        CLASSIFY(deg);
    }
#endif

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


/* ================================================================== */
/* DOUBLES ONLY                                                          */
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
    int nd=0;
    for (int i=0; i<size_buffer; ++i)
        if (pc1(spindet0 ^ buffer[i]) == 4)
            doubles[nd++] = idx[i];
    *n_doubles = nd;
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
    int nd=0;

#if HAVE_NEON
    const uint64x2_t sv = vld1q_u64(spindet);
    for (int i=0; i<size_buffer; ++i) {
        uint64x2_t bv = vld1q_u64(buffer + 2*i);
        if (popcnt_neon_128(veorq_u64(sv, bv)) == 4)
            doubles[nd++] = idx[i];
    }
#else
    const uint64_t s0=spindet[0], s1=spindet[1];
    for (int i=0; i<size_buffer; ++i)
        if (pc2(s0^buffer[2*i], s1^buffer[2*i+1]) == 4)
            doubles[nd++] = idx[i];
#endif

    *n_doubles = nd;
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
    int nd=0;

#if HAVE_NEON
    const uint64x2_t sv01 = vld1q_u64(spindet);
    const uint64_t   s2   = spindet[2];
    for (int i=0; i<size_buffer; ++i) {
        uint64x2_t bv01 = vld1q_u64(buffer + 3*i);
        int deg = popcnt_neon_128(veorq_u64(sv01, bv01))
                + pc1(s2 ^ buffer[3*i+2]);
        if (deg == 4) doubles[nd++] = idx[i];
    }
#else
    const uint64_t s0=spindet[0], s1=spindet[1], s2=spindet[2];
    for (int i=0; i<size_buffer; ++i)
        if (pc3(s0^buffer[3*i], s1^buffer[3*i+1], s2^buffer[3*i+2]) == 4)
            doubles[nd++] = idx[i];
#endif

    *n_doubles = nd;
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
    int nd=0;

#if HAVE_SVE
    const svbool_t   pg = svptrue_b64();
    const svuint64_t sv = svld1_u64(pg, spindet);
    for (int i=0; i<size_buffer; ++i) {
        svuint64_t bv = svld1_u64(pg, buffer + 4*i);
        if (popcnt_sve(sveor_u64_z(pg, sv, bv), pg) == 4)
            doubles[nd++] = idx[i];
    }

#elif HAVE_NEON
    const uint64x2_t sv0 = vld1q_u64(spindet);
    const uint64x2_t sv1 = vld1q_u64(spindet + 2);
    for (int i=0; i<size_buffer; ++i) {
        int deg = popcnt_neon_128(veorq_u64(sv0, vld1q_u64(buffer + 4*i)))
                + popcnt_neon_128(veorq_u64(sv1, vld1q_u64(buffer + 4*i + 2)));
        if (deg == 4) doubles[nd++] = idx[i];
    }

#else
    const uint64_t s0=spindet[0],s1=spindet[1],s2=spindet[2],s3=spindet[3];
    for (int i=0; i<size_buffer; ++i)
        if (pc4(s0^buffer[4*i],s1^buffer[4*i+1],
                s2^buffer[4*i+2],s3^buffer[4*i+3]) == 4)
            doubles[nd++] = idx[i];
#endif

    *n_doubles = nd;
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
    int nd=0;

#if HAVE_SVE
    const int sve_step = (int)svcntd();
    for (int i=0; i<size_buffer; ++i) {
        const uint64_t *col = buffer + N_int * i;
        int deg = 0, k = 0;
        for (; k + sve_step <= N_int; k += sve_step) {
            svbool_t   pg = svptrue_b64();
            svuint64_t xv = sveor_u64_z(pg, svld1_u64(pg, spindet+k),
                                             svld1_u64(pg, col+k));
            deg += popcnt_sve(xv, pg);
        }
        if (k < N_int) {
            svbool_t   pg = svwhilelt_b64((uint64_t)k, (uint64_t)N_int);
            svuint64_t xv = sveor_u64_z(pg, svld1_u64(pg, spindet+k),
                                             svld1_u64(pg, col+k));
            deg += popcnt_sve(xv, pg);
        }
        if (deg == 4) doubles[nd++] = idx[i];
    }

#elif HAVE_NEON
    for (int i=0; i<size_buffer; ++i) {
        const uint64_t *col = buffer + N_int * i;
        int deg = 0, k = 0;
        for (; k + 2 <= N_int; k += 2)
            deg += popcnt_neon_128(veorq_u64(vld1q_u64(spindet+k),
                                             vld1q_u64(col+k)));
        if (k < N_int)
            deg += pc1(spindet[k] ^ col[k]);
        if (deg == 4) doubles[nd++] = idx[i];
    }

#else
    for (int i=0; i<size_buffer; ++i) {
        const uint64_t *col = buffer + N_int * i;
        int deg = 0;
        for (int k=0; k<N_int; ++k)
            deg += __builtin_popcountll(spindet[k] ^ col[k]);
        if (deg == 4) doubles[nd++] = idx[i];
    }
#endif

    *n_doubles = nd;
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

