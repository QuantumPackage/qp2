#include <stdint.h>
#include <stddef.h>
#include "local_cpu.h"

#if defined(USE_AVX512) || defined(USE_AVX2)
#include <immintrin.h>
#elif defined(USE_SVE)
#include <arm_sve.h>
#elif defined(USE_NEON)
#include <arm_neon.h>
#endif

void get_excitation_degree_c(const uint64_t *restrict key1,
                             const uint64_t *restrict key2,
                             int32_t *restrict degree_ptr,
                             const int32_t nint) {

  const size_t n = (size_t)nint << 1; // n = Nint * 2
  uint64_t total_pop = 0;
  size_t i = 0;

#if defined(USE_AVX512)
  /* --- AVX-512 Path ---
   * High efficiency for small Nint using masking to avoid remainder loops.
   */
  __m512i vsum = _mm512_setzero_si512();
  for (; i <= (n >= 8 ? n - 8 : -1); i += 8) {
    __m512i x = _mm512_xor_si512(_mm512_loadu_si512(&key1[i]),
                                 _mm512_loadu_si512(&key2[i]));
    vsum = _mm512_add_epi64(vsum, _mm512_popcnt_epi64(x));
  }
  // Handle remainder (n < 8 or leftover) with a single masked operation
  if (i < n) {
    __mmask8 mask = (__mmask8)((1U << (n - i)) - 1);
    __m512i x = _mm512_maskz_xor_epi64(mask,
                                       _mm512_maskz_loadu_epi64(mask, &key1[i]),
                                       _mm512_maskz_loadu_epi64(mask, &key2[i]));
    vsum = _mm512_add_epi64(vsum, _mm512_popcnt_epi64(x));
  }
  total_pop = _mm512_reduce_add_epi64(vsum);

#elif defined(USE_SVE)
  /* --- ARM SVE Path ---
   * Naturally handles small/variable Nint via predication.
   */
  svbool_t pg = svwhilelt_b64(i, n);
  while (svptest_any(svptrue_b64(), pg)) {
    svuint64_t x = sveor_u64_z(pg, svld1_u64(pg, &key1[i]), svld1_u64(pg, &key2[i]));
    total_pop += svaddv_u64(pg, svcnt_u64_z(pg, x));
    i += svcntd();
    pg = svwhilelt_b64(i, n);
  }

#elif defined(USE_AVX2)
  /* --- AVX2 Path ---
   * For Nint < 4, scalar popcnt is usually faster than Muła LUT setup.
   */
  if (n >= 8) {
    const __m256i LO_MASK = _mm256_set1_epi8(0x0f);
    const __m256i LUT     = _mm256_setr_epi8(0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
                                             0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4);
    __m256i vsum_sad = _mm256_setzero_si256();

    for (; i <= n - 4; i += 4) {
      __m256i v   = _mm256_xor_si256(_mm256_loadu_si256((void*)&key1[i]),
                                     _mm256_loadu_si256((void*)&key2[i]));
      __m256i lo  = _mm256_and_si256(v, LO_MASK);
      __m256i hi  = _mm256_and_si256(_mm256_srli_epi16(v, 4), LO_MASK);
      __m256i cnt = _mm256_add_epi8(_mm256_shuffle_epi8(LUT, lo),
                                    _mm256_shuffle_epi8(LUT, hi));
      vsum_sad = _mm256_add_epi64(vsum_sad, _mm256_sad_epu8(cnt, _mm256_setzero_si256()));
    }
    __m128i final_v = _mm_add_epi64(_mm256_castsi256_si128(vsum_sad),
                                    _mm256_extracti128_si256(vsum_sad, 1));
    final_v = _mm_add_epi64(final_v, _mm_unpackhi_epi64(final_v, final_v));
    total_pop = (uint64_t)_mm_cvtsi128_si64(final_v);
  }
  // Fallback for n < 8 or remainder
  for (; i < n; i++) {
    total_pop += __builtin_popcountll(key1[i] ^ key2[i]);
  }

#elif defined(USE_NEON)
  /* --- ARM NEON Path --- */
  for (; i <= n - 2; i += 2) {
    uint8x16_t x = vreinterpretq_u8_u64(veorq_u64(vld1q_u64(&key1[i]), vld1q_u64(&key2[i])));
    total_pop += vaddlvq_u8(vcntq_u8(x));
  }
  if (i < n) total_pop += __builtin_popcountll(key1[i] ^ key2[i]);

#else
  /* --- Pure Scalar Path --- */
  for (int i=0; i < n; i++) {
    total_pop += __builtin_popcountll(key1[i] ^ key2[i]);
  }
#endif

  *degree_ptr = (int32_t)(total_pop >> 1);
}
