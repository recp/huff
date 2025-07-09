/*
 * Copyright (C) 2024 Recep Aslantas
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef huff_read_h
#define huff_read_h
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__ARM_NEON)
#include <arm_neon.h>
#elif defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

HUFF_INLINE
int
huff_read_scalar(const uint8_t ** __restrict buff,
                 bitstream_t    * __restrict bits,
                 const uint8_t  * __restrict end) {
  const uint8_t *p;
  bitstream_t    result;
  size_t         remb;
  int            i, s, maxb, n;

  p       = *buff;
  remb    = end - p;
  maxb    = (int)sizeof(bitstream_t);
  n       = ((int)remb < maxb) ? (int)remb : maxb;
  result  = 0;

  for (i = s = 0; i < n; i++, s += 8) {
    result |= ((bitstream_t)p[i]) << s;
  }

  *buff   += n;
  *bits    = result;

  return n * 8;
}

#if defined(__ARM_NEON)
HUFF_INLINE
int
huff_read_neon(const uint8_t ** __restrict buff,
               bitstream_t    * __restrict bits,
               const uint8_t  * __restrict end) {
  const uint8_t *p;
  bitstream_t    result;
  size_t         remb;
  int            maxb, n;

  p       = *buff;
  remb    = end - p;
  maxb    = (int)sizeof(bitstream_t);
  n       = ((int)remb < maxb) ? (int)remb : maxb;

#ifdef ENABLE_BIG_BITSTREAM
  // For 128-bit mode, load two 64-bit chunks
  uint8x16_t bytes = vld1q_u8(p);
  uint64x2_t chunks = vreinterpretq_u64_u8(bytes);

  // Extract both 64-bit halves
  uint64_t low = vgetq_lane_u64(chunks, 0);
  uint64_t high = vgetq_lane_u64(chunks, 1);

  // Combine into 128-bit result
  result = (((bitstream_t)high) << 64) | low;
#else
  // For 64-bit mode, load single 64-bit chunk
  uint8x8_t bytes = vld1_u8(p);
  uint64x1_t chunk = vreinterpret_u64_u8(bytes);
  result = vget_lane_u64(chunk, 0);
#endif

  *buff   += n;
  *bits    = result;

  return n * 8;
}
#elif defined(__x86_64__) || defined(_M_X64)
HUFF_INLINE
int
huff_read_sse(const uint8_t ** __restrict buff,
              bitstream_t    * __restrict bits,
              const uint8_t  * __restrict end) {
  const uint8_t *p;
  bitstream_t    result;
  size_t         remb;
  int            maxb, n;

  p       = *buff;
  remb    = end - p;
  maxb    = (int)sizeof(bitstream_t);
  n       = ((int)remb < maxb) ? (int)remb : maxb;
  result  = 0;

#ifdef ENABLE_BIG_BITSTREAM
  // For 128-bit mode, load entire 128 bits at once
  __m128i bytes = _mm_loadu_si128((__m128i*)p);

#if defined(__SIZEOF_INT128__)
  // Extract both 64-bit halves and combine
  uint64_t low = _mm_cvtsi128_si64(bytes);
  __m128i high_vec = _mm_srli_si128(bytes, 8);
  uint64_t high = _mm_cvtsi128_si64(high_vec);

  result = (((bitstream_t)high) << 64) | low;
#else
  // If 128-bit integers aren't supported, fall back to scalar
  unsigned char *bytes_ptr = (unsigned char *)&bytes;
  for (int i = 0; i < n; i++) {
    result |= ((bitstream_t)bytes_ptr[i]) << (i * 8);
  }
#endif

#else
  // For 64-bit mode, just extract lower 64 bits
  result = _mm_cvtsi128_si64(_mm_loadu_si128((__m128i*)p));
#endif

  *buff   += n;
  *bits    = result;

  return n * 8;
}

#ifdef __AVX2__
HUFF_INLINE
int
huff_read_avx2(const uint8_t ** __restrict buff,
               bitstream_t    * __restrict bits,
               const uint8_t  * __restrict end) {
  const uint8_t *p;
  bitstream_t    result;
  size_t         remb;
  int            maxb, n;

  p       = *buff;
  remb    = end - p;
  maxb    = (int)sizeof(bitstream_t);
  n       = ((int)remb < maxb) ? (int)remb : maxb;
  result  = 0;

#ifdef ENABLE_BIG_BITSTREAM
  // For 128-bit mode, use AVX2 for loading
  __m256i bytes = _mm256_loadu_si256((__m256i*)p);
  __m128i lower = _mm256_extracti128_si256(bytes, 0);

#if defined(__SIZEOF_INT128__)
  // Extract both 64-bit halves and combine
  uint64_t low = _mm_cvtsi128_si64(lower);
  __m128i high_vec = _mm_srli_si128(lower, 8);
  uint64_t high = _mm_cvtsi128_si64(high_vec);

  result = (((bitstream_t)high) << 64) | low;
#else
  // If 128-bit integers aren't supported, fall back to scalar
  unsigned char *bytes_ptr = (unsigned char *)&lower;
  for (int i = 0; i < n; i++) {
    result |= ((bitstream_t)bytes_ptr[i]) << (i * 8);
  }
#endif

#else
  // For 64-bit mode, extract lower 64 bits
  __m128i lower = _mm256_extracti128_si256(_mm256_loadu_si256((__m256i*)p), 0);
  result = _mm_cvtsi128_si64(lower);
#endif

  *buff   += n;
  *bits    = result;

  return n * 8;
}
#endif
#endif

HUFF_INLINE
int
huff_read(const uint8_t ** __restrict buff,
          bitstream_t    * __restrict bits,
          const uint8_t  * __restrict end) {
#if defined(__AVX2__)
  return huff_read_avx2(buff, bits, end);
#elif defined(__ARM_NEON)
  return huff_read_neon(buff, bits, end);
#elif defined(__x86_64__) || defined(_M_X64)
  return huff_read_sse(buff, bits, end);
#else
  return huff_read_scalar(buff, bits, end);
#endif
}

#ifdef __cplusplus
}
#endif
#endif /* huff_read_h */
