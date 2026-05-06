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
  if (unlikely(p >= end)) {
    *bits = 0;
    return 0;
  }

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
  if (unlikely(p >= end))
    return huff_read_scalar(buff, bits, end);

  remb    = end - p;
  maxb    = (int)sizeof(bitstream_t);
  n       = ((int)remb < maxb) ? (int)remb : maxb;

#if defined(HUFF_ENABLE_BIG_BITSTREAM) && defined(__SIZEOF_INT128__)
  if (likely(n == maxb)) {
    uint8x16_t bytes;
    uint64x2_t chunks;
    uint64_t   low, high;

    bytes  = vld1q_u8(p);
    chunks = vreinterpretq_u64_u8(bytes);
    low    = vgetq_lane_u64(chunks, 0);
    high   = vgetq_lane_u64(chunks, 1);
    result = (((bitstream_t)high) << 64) | low;
  } else {
    return huff_read_scalar(buff, bits, end);
  }
#else
  if (likely(n == maxb)) {
    uint8x8_t  bytes;
    uint64x1_t chunk;

    bytes  = vld1_u8(p);
    chunk  = vreinterpret_u64_u8(bytes);
    result = vget_lane_u64(chunk, 0);
  } else {
    return huff_read_scalar(buff, bits, end);
  }
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
  if (unlikely(p >= end))
    return huff_read_scalar(buff, bits, end);

  remb    = end - p;
  maxb    = (int)sizeof(bitstream_t);
  n       = ((int)remb < maxb) ? (int)remb : maxb;
  result  = 0;

#if defined(HUFF_ENABLE_BIG_BITSTREAM) && defined(__SIZEOF_INT128__)
  if (likely(n == maxb)) {
    __m128i  bytes, high_vec;
    uint64_t low, high;

    bytes    = _mm_loadu_si128((const __m128i*)p);
    low      = (uint64_t)_mm_cvtsi128_si64(bytes);
    high_vec = _mm_srli_si128(bytes, 8);
    high     = (uint64_t)_mm_cvtsi128_si64(high_vec);
    result   = (((bitstream_t)high) << 64) | low;
  } else {
    return huff_read_scalar(buff, bits, end);
  }
#else
  if (likely(n == maxb)) {
    result = (uint64_t)_mm_cvtsi128_si64(_mm_loadl_epi64((const __m128i*)p));
  } else {
    return huff_read_scalar(buff, bits, end);
  }
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
  return huff_read_sse(buff, bits, end);
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
