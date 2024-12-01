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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <stdint.h>
#include <stddef.h>

#if defined(__AVX2__)
#  include <immintrin.h> /* AVX2 intrinsics */
#  define SIMD_WIDTH 32
#elif defined(__AVX__)
#  include <immintrin.h> /* AVX intrinsics */
#  define SIMD_WIDTH 32
#elif defined(__SSE2__)
#  include <emmintrin.h> /* SSE2 intrinsics */
#  define SIMD_WIDTH 16
#elif defined(__ARM_NEON)
#  include <arm_neon.h>  /* NEON intrinsics */
#  define SIMD_WIDTH 16
#else
#  define SIMD_WIDTH 8 /* Fallback to scalar */
#endif

#include "../include/huff.h"

/* 15: DEFLATE, 16: JPEG */
#define MAX_CODE_LENGTH 16  /* Maximum length of Huffman codes */
#define FAST_TABLE_BITS  8  /* Number of bits used for fast lookup */

#define BITSTREAM_T_IS_128 (sizeof(bitstream_t) > 8)

static const uint8_t bit_reverse_table[256] = {
  0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
  0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
  0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
  0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
  0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
  0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
  0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
  0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
  0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
  0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
  0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
  0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
  0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
  0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
  0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
  0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
};

HUFF_EXPORT
void
huff_init(huff_table_t   * __restrict table,
          const uint8_t  * __restrict lengths,
          const uint16_t * __restrict symbols,
          uint16_t                    n) {
  uint_fast16_t sym_idx, code, i, start, end, fast_idx, sym;
  uint_fast8_t  len;

  /* initialize fast table with invalid values */
  memset(table->fast_table, 0xFF, sizeof(table->fast_table));

  sym_idx = code = 0;

  /* process each code length (1 to MAX_CODE_LENGTH) */
  for (len = 1; len <= MAX_CODE_LENGTH; len++) {
    table->sym_offset[len] = sym_idx;

    /* iterate over all symbols to find those with the current length */
    for (i = 0; i < n; i++) {
      if (lengths[i] == len) {
        /* handle sequential symbols */
        sym                    = (symbols != NULL) ? symbols[i] : i;
        table->syms[sym_idx++] = sym;

        /* precompute fast table for short codes (≤ FAST_TABLE_BITS) */
        if (len <= FAST_TABLE_BITS) {
          start = code << (FAST_TABLE_BITS - len);
          end   = (code + 1) << (FAST_TABLE_BITS - len);

          for (fast_idx = start; fast_idx < end; fast_idx++) {
            table->fast_table[fast_idx]  = sym;
            table->fast_length[fast_idx] = len;
          }
        }

        /* increment the canonical code */
        code++;
      }
    }

    /* compute sentinel bits for the current length */
    table->sentinel_bits[len] = code << (MAX_CODE_LENGTH - len);
  }
}

HUFF_EXPORT
big_int_t
huff_rev_bits(big_int_t x) {
#if defined(__AVX2__)
#if sizeof(big_int_t) > 8
  // AVX2 Implementation for 128-bit big_int_t
  __m256i input = _mm256_set_epi64x(x >> 64, x & 0xFFFFFFFFFFFFFFFFULL, 0, 0);
  __m256i lookup = _mm256_loadu_si256((__m256i *)bit_reverse_table);
  __m256i reversed = _mm256_shuffle_epi8(lookup, input);
  uint64_t lower = _mm256_extract_epi64(reversed, 0);
  uint64_t upper = _mm256_extract_epi64(reversed, 1);
  return ((big_int_t)upper << 64) | lower;
#else
  // AVX2 Implementation for 64-bit big_int_t
  __m128i input = _mm_set_epi64x(0, x);
  __m128i lookup = _mm_loadu_si128((__m128i *)bit_reverse_table);
  __m128i reversed = _mm_shuffle_epi8(lookup, input);
  return (big_int_t)_mm_cvtsi128_si64(reversed);
#endif
#elif defined(__SSE2__)
  if (sizeof(big_int_t) == 8) {
    __m128i input = _mm_set_epi64x(0, x);

    // Load the bit_reverse_table into an SSE register
    __m128i lookup = _mm_loadu_si128((__m128i *)bit_reverse_table);

    // Perform byte-wise lookup manually
    __m128i reversed = _mm_shuffle_epi8(lookup, input);

    return (big_int_t)_mm_cvtsi128_si64(reversed);
  }
#elif defined(__ARM_NEON)
  if (sizeof(big_int_t) > 8) {
    uint8x16_t input = vreinterpretq_u8_u64((uint64x2_t){x >> 64, x});
    uint8x16_t lookup = vld1q_u8(bit_reverse_table);
    uint8x16_t reversed = vqtbl1q_u8(lookup, input);
    uint64_t lower = vgetq_lane_u64(vreinterpretq_u64_u8(reversed), 0);
    uint64_t upper = vgetq_lane_u64(vreinterpretq_u64_u8(reversed), 1);
    return ((big_int_t)upper << 64) | lower;
  } else if (sizeof(big_int_t) == 8) {
    uint8x8_t input = vreinterpret_u8_u64(vdup_n_u64(x));
    uint8x8_t lookup = vld1_u8(bit_reverse_table);
    uint8x8_t reversed = vtbl1_u8(lookup, input);
    return (big_int_t)vget_lane_u64(vreinterpret_u64_u8(reversed), 0);
  }
#else
  // Scalar Fallback for all sizes

  size_t    bit_count = sizeof(big_int_t) * 8; // Dynamically determine bit width
  big_int_t result    = 0;

  for (size_t i = 0; i < bit_count; i++) {
    result |= ((x >> i) & 1) << (bit_count - 1 - i);
  }
  return result;

#endif
}

HUFF_EXPORT
bitstream_t
huff_read_lsb(const uint8_t *stream, size_t *bit_offset, size_t stream_size) {
  size_t byte_offset = *bit_offset / 8;
  size_t bit_in_byte = *bit_offset % 8;

  // Ensure we don't read beyond the end of the stream
  size_t remaining_bytes = (stream_size > byte_offset) ? (stream_size - byte_offset) : 0;

  bitstream_t result = 0;

#if defined(__AVX2__) || defined(__AVX__)
  if (remaining_bytes >= 32) {
    // Load 32 bytes using AVX or AVX2
    __m256i data = _mm256_loadu_si256((__m256i *)&stream[byte_offset]);
    uint64_t lower = _mm256_extract_epi64(data, 0);
    uint64_t upper = _mm256_extract_epi64(data, 1);
    result = ((bitstream_t)upper << 64) | lower;
  } else
#endif
#if defined(__SSE2__)
    if (remaining_bytes >= 16) {
      // Load 16 bytes using SSE2
      __m128i data = _mm_loadu_si128((__m128i *)&stream[byte_offset]);
      uint64_t lower = _mm_cvtsi128_si64(data);
      uint64_t upper = _mm_extract_epi64(data, 1);
      result = ((bitstream_t)upper << 64) | lower;
    } else
#endif
#if defined(__ARM_NEON)
      if (remaining_bytes >= 16) {
        // Load 16 bytes using NEON
        uint8x16_t data = vld1q_u8(&stream[byte_offset]);
        uint64_t lower = vgetq_lane_u64(vreinterpretq_u64_u8(data), 0);
        uint64_t upper = vgetq_lane_u64(vreinterpretq_u64_u8(data), 1);
        result = sizeof(bitstream_t) > 8 ? (((bitstream_t)upper << 64) | lower) : (bitstream_t)lower;
      } else
#endif
        {
        // Scalar fallback
        for (size_t i = 0; i < sizeof(bitstream_t) && i < remaining_bytes; i++) {
          result |= (bitstream_t)stream[byte_offset + i] << (i * 8);
        }
        }

  // Align the result with the bit offset
  result >>= bit_in_byte;

  // Reverse bits for LSB-first order
  result = huff_rev_bits(result);

  // Advance the bit offset
  *bit_offset += remaining_bytes * 8 - bit_in_byte;

  return result;
}

HUFF_EXPORT
bitstream_t
huff_read_msb(const uint8_t *stream, size_t *bit_offset, size_t stream_size) {
  size_t byte_offset = *bit_offset / 8;
  size_t bit_in_byte = *bit_offset % 8;

  // Ensure we don't read beyond the end of the stream
  size_t remaining_bytes = (stream_size > byte_offset) ? (stream_size - byte_offset) : 0;

  bitstream_t result = 0;

#if defined(__AVX2__) || defined(__AVX__)
  if (remaining_bytes >= 32) {
    // Load 32 bytes using AVX or AVX2
    __m256i data = _mm256_loadu_si256((__m256i *)&stream[byte_offset]);
    uint64_t lower = _mm256_extract_epi64(data, 0);
    uint64_t upper = _mm256_extract_epi64(data, 1);
    result = ((bitstream_t)upper << 64) | lower;
  } else
#endif
#if defined(__SSE2__)
    if (remaining_bytes >= 16) {
      // Load 16 bytes using SSE2
      __m128i data = _mm_loadu_si128((__m128i *)&stream[byte_offset]);
      uint64_t lower = _mm_cvtsi128_si64(data);
      uint64_t upper = _mm_extract_epi64(data, 1);
      result = ((bitstream_t)upper << 64) | lower;
    } else
#endif
#if defined(__ARM_NEON)
      if (remaining_bytes >= 16) {
        // Load 16 bytes using NEON
        uint8x16_t data = vld1q_u8(&stream[byte_offset]);
        uint64_t lower = vgetq_lane_u64(vreinterpretq_u64_u8(data), 0);
        uint64_t upper = vgetq_lane_u64(vreinterpretq_u64_u8(data), 1);
        result = sizeof(bitstream_t) > 8 ? (((bitstream_t)upper << 64) | lower) : (bitstream_t)lower;
      } else
#endif
        {
        // Scalar fallback
        for (size_t i = 0; i < sizeof(bitstream_t) && i < remaining_bytes; i++) {
          result |= (bitstream_t)stream[byte_offset + i] << (i * 8);
        }
        }

  // Align the result with the bit offset
  result >>= bit_in_byte;

  // Advance the bit offset
  *bit_offset += remaining_bytes * 8 - bit_in_byte;

  return result;
}

HUFF_EXPORT
uint_fast16_t
huff_decode(const huff_table_t *table, bitstream_t bitstream, uint8_t bit_length) {
  bitstream_t   code;
  uint_fast16_t fast_idx, idx;
  uint_fast8_t  l, fast_len;

  /* fast lookup for short codes (≤ FAST_TABLE_BITS) */
  fast_idx = bitstream & ((1 << FAST_TABLE_BITS) - 1);
  fast_len = table->fast_length[fast_idx];

  if (fast_len <= bit_length) {
    return table->fast_table[fast_idx];
  }

  /* fallback for longer codes */
  for (l = FAST_TABLE_BITS + 1; l <= MAX_CODE_LENGTH; l++) {
    if (bitstream < table->sentinel_bits[l]) {
      code = bitstream >> (MAX_CODE_LENGTH - l); /* Extract `l`-bit code */
      idx  = table->sym_offset[l] + (code - table->sentinel_bits[l - 1]);
      return table->syms[idx];
    }
  }

  /* decoding failed */
  return (uint_fast16_t)-1;
}
