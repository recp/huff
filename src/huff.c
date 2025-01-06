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

#include "../include/huff/huff.h"

#define FAST_MASK        ((1U << FAST_TABLE_BITS) - 1)
#define BIT_MASK(l)      (((bitstream_t)1 << (l)) - 1)

#define BITSTREAM_T_IS_128 (sizeof(bitstream_t) > 8)

static const uint8_t bit_reverse_table[256] HUFF_ALIGN(32) = {
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

#if DEBUG
static void
huff_debug_print(const huff_table_t *table) {
  printf("=== Huffman Table Debug Info ===\n");

  // Print sentinel_bits
  printf("Sentinel Bits:\n");
  for (uint_fast8_t len = 1; len <= MAX_CODE_LENGTH; len++) {
    printf("Length %u: Sentinel Bits = %u\n", len, table->sentinel_bits[len]);
  }

  // Print sym_offset
  printf("\nSymbol Offsets:\n");
  for (uint_fast8_t len = 1; len <= MAX_CODE_LENGTH; len++) {
    printf("Length %u: Sym Offset = %u\n", len, table->sym_offset[len]);
  }

  // Print syms array
  printf("\nSymbols Array:\n");
  for (uint_fast16_t i = 0; i < table->num_symbols; i++) {
    printf("Symbol[%u] = %u\n", i, table->syms[i]);
  }

  // Print fast table
  printf("\nFast Table:\n");
  for (uint_fast16_t i = 0; i < FAST_TABLE_SIZE; i++) {
    printf("Fast Table[%u]: Sym = %u, Len = %u\n",
           i,
           table->fast_table[i].sym,
           table->fast_table[i].len);
  }

  printf("================================\n");
}
#endif

HUFF_EXPORT
void
huff_init_lsb(huff_table_t   * __restrict table,
              const uint8_t  * __restrict lengths,
              const uint16_t * __restrict symbols,
              uint16_t                    n) {
  uint_fast16_t sym_idx, code, i, fast_idx, sym, blksz, ext;
  uint_fast8_t  len;

  /* Initialize fast table with invalid entries */
  for (fast_idx = 0; fast_idx < FAST_TABLE_SIZE; fast_idx++) {
    table->fast_table[fast_idx].sym = 0x1FF; /* invalid symbol (max 9 bits) */
    table->fast_table[fast_idx].len = 0;     /* invalid length */
  }

  if (!table->syms) {
    table->syms = calloc(n, sizeof(uint16_t));
  }

  table->num_symbols      = n;
  table->sentinel_bits[0] = 0; /* No codes of length 0 */
  sym_idx                 = 0;
  code                    = 0; /* LSB-based code counter */

  /* Process each code length (1 to MAX_CODE_LENGTH) */
  for (len = 1; len <= MAX_CODE_LENGTH; len++) {
    table->sym_offset[len] = sym_idx;

    /* Iterate over all symbols with the current length */
    for (i = 0; i < n; i++) {
      if (lengths[i] == len) {
        /* Sequential symbols */
        sym = (symbols != NULL) ? symbols[i] : i;
        table->syms[sym_idx++] = sym;

        /* Precompute fast table for short codes (≤ FAST_TABLE_BITS) */
        if (len <= FAST_TABLE_BITS) {
          blksz = 1U << (FAST_TABLE_BITS - len);
          for (ext = 0; ext < blksz; ext++) {
            fast_idx                        = code + (ext << len);
            table->fast_table[fast_idx].sym = sym;
            table->fast_table[fast_idx].len = len;
          }
        }

        /* Increment the canonical LSB-first code */
        code++;
      }
    }

    /* Sentinel bits for current length */
    table->sentinel_bits[len] = code;
  }
}
HUFF_EXPORT
void
huff_init_lsb3(huff_table_t   * __restrict table,
              const uint8_t  * __restrict lengths,
              const uint16_t * __restrict symbols,
              uint16_t                    n) {
  uint_fast16_t sym_idx, code, i, fast_idx, sym, blksz, ext;
  uint_fast8_t  len;

  /* initialize fast table with invalid entries */
  for (fast_idx = 0; fast_idx < FAST_TABLE_SIZE; fast_idx++) {
    table->fast_table[fast_idx].sym = 0x1FF; /* invalid symbol (max 9 bits) */
    table->fast_table[fast_idx].len = 0;     /* invalid length */
  }

  if (!table->syms) {
    table->syms = calloc(n, sizeof(uint16_t));
  }

  table->num_symbols      = n;
  table->sentinel_bits[0] = 0; /* no codes of length 0 */
  sym_idx                 = 0; /* where we store next symbol in table->syms[] */
  code                    = 0; /* LSB-based code counter */

  /*
   * 3) For each code length from 1..MAX_CODE_LENGTH:
   *    Assign the current 'code' to all symbols that have this length.
   *    Then increment 'code' for each such symbol, continuing
   *    in plain ascending integer order (no shifts).
   */
  for (len = 1; len <= MAX_CODE_LENGTH; len++) {
    table->sym_offset[len] = sym_idx;

    /* for all symbols, check if they have 'lengths[i] == len'. */
    for (i = 0; i < n; i++) {
      if (lengths[i] == len) {
        /* handle sequential symbols */
        sym                    = (symbols != NULL) ? symbols[i] : i;
        table->syms[sym_idx++] = sym;

        /* precompute fast table for short codes (≤ FAST_TABLE_BITS) */
        if (len <= FAST_TABLE_BITS) {
          /* We'll replicate 'code' across all possible higher bits:
           *   index = code + (ext << len)
           * for ext in [0 .. (1 << (FAST_TABLE_BITS - len)) - 1].
           */
          blksz = 1U << (FAST_TABLE_BITS - len);
          for (ext = 0; ext < blksz; ext++) {
            fast_idx                        = code + (ext << len);
            table->fast_table[fast_idx].sym = sym;
            table->fast_table[fast_idx].len = len;
          }
        }

        /* increment the canonical code */
        code++;
      }
    }

    /* sentinel_bits[len] = total # of codes up through length 'len' */
    table->sentinel_bits[len] = code;
  }

  // huff_debug_print(table);
}

HUFF_EXPORT
void
huff_init_msb(huff_table_t   * __restrict table,
              const uint8_t  * __restrict lengths,
              const uint16_t * __restrict symbols,
              uint16_t                    n) {
  uint_fast16_t sym_idx, code, i, start, end, fast_idx, sym;
  uint_fast8_t  len, len_shift;

  /* initialize fast table with invalid entries */
  for (fast_idx = 0; fast_idx < FAST_TABLE_SIZE; fast_idx++) {
    table->fast_table[fast_idx].sym = 0x1FF; /* invalid symbol (max 9 bits) */
    table->fast_table[fast_idx].len = 0;     /* invalid length */
  }


  if (!table->syms) {
    table->syms = calloc(n, sizeof(uint16_t));
  }

  table->num_symbols = n;
  sym_idx = code = 0;

  /* process each code length (1 to MAX_CODE_LENGTH) */
  for (len = 1; len <= MAX_CODE_LENGTH; len++) {
    table->sym_offset[len] = sym_idx;
    len_shift              = MAX_CODE_LENGTH - len;

    /* iterate over all symbols to find those with the current length */
    for (i = 0; i < n; i++) {
      if (lengths[i] == len) {
        /* handle sequential symbols */
        sym                    = (symbols != NULL) ? symbols[i] : i;
        table->syms[sym_idx++] = sym;

        /* precompute fast table for short codes (≤ FAST_TABLE_BITS) */
        if (len <= FAST_TABLE_BITS) {
          start = (code << len_shift) >> FAST_SHIFT;
          end   = ((code + 1) << len_shift) >> FAST_SHIFT;

          for (fast_idx = start; fast_idx < end; fast_idx++) {
            table->fast_table[fast_idx].sym = sym;
            table->fast_table[fast_idx].len = len;
          }
        }

        /* increment the canonical code */
        code++;
      }
    }

    /* compute sentinel bits for the current length (MSB-first) */
    table->sentinel_bits[len] = code << len_shift;
  }
}

HUFF_EXPORT
bitstream_t
huff_rev_bits(bitstream_t x) {
#if defined(__AVX2__)
  if (sizeof(bitstream_t) > 8) {
    uint64_t lower    = x & 0xFFFFFFFFFFFFFFFFULL;
    uint64_t upper    = x >> 64;
    __m256i lookup    = _mm256_loadu_si256((__m256i *)bit_reverse_table);
    __m128i lower_vec = _mm_cvtsi64_si128(lower);
    __m128i upper_vec = _mm_cvtsi64_si128(upper);
    lower_vec         = _mm_shuffle_epi8(_mm256_castsi256_si128(lookup), lower_vec);
    upper_vec         = _mm_shuffle_epi8(_mm256_castsi256_si128(lookup), upper_vec);
    lower             = _mm_cvtsi128_si64(lower_vec);
    upper             = _mm_cvtsi128_si64(upper_vec);
    return ((bitstream_t)upper << 64) | lower;
  } else {
    __m128i input    = _mm_cvtsi64_si128(x);
    __m128i lookup   = _mm_loadu_si128((__m128i *)bit_reverse_table);
    __m128i reversed = _mm_shuffle_epi8(lookup, input);
    return (bitstream_t)_mm_cvtsi128_si64(reversed);
  }
#elif defined(__SSSE3__)
  if (sizeof(bitstream_t) == 8) {
    __m128i input    = _mm_cvtsi64_si128(x);
    __m128i lookup   = _mm_loadu_si128((__m128i *)bit_reverse_table);
    __m128i reversed = _mm_shuffle_epi8(lookup, input);
    return (bitstream_t)_mm_cvtsi128_si64(reversed);
  }
#elif defined(__ARM_NEON)
  if (sizeof(bitstream_t) > 8) {
    uint8x16_t input    = vreinterpretq_u8_u64((uint64x2_t){x >> 64, x});
    uint8x16_t lookup   = vld1q_u8(bit_reverse_table);
    uint8x16_t reversed = vqtbl1q_u8(lookup, input);
    uint64_t lower      = vgetq_lane_u64(vreinterpretq_u64_u8(reversed), 0);
    uint64_t upper      = vgetq_lane_u64(vreinterpretq_u64_u8(reversed), 1);
    return ((bitstream_t)upper << 64) | lower;
  } else {
    uint8x8_t input    = vreinterpret_u8_u64(vdup_n_u64(x));
    uint8x8_t lookup   = vld1_u8(bit_reverse_table);
    uint8x8_t reversed = vtbl1_u8(lookup, input);
    return (bitstream_t)vget_lane_u64(vreinterpret_u64_u8(reversed), 0);
  }
#else
  size_t      bit_count = sizeof(bitstream_t) * 8;
  bitstream_t result    = 0;

  for (size_t i = 0; i < bit_count; i++) {
    result |= ((x >> i) & 1) << (bit_count - 1 - i);
  }
  return result;
#endif
}
HUFF_EXPORT
bitstream_t
huff_read(const uint8_t ** __restrict p,
          size_t        * __restrict bit_offset,
          uint8_t       * __restrict nbits,
          const uint8_t * __restrict end)
{
  // Ensure we don't read beyond the end of the buffer
  if (*p >= end) {
    *nbits = 0; // No more bits available
    return 0;
  }

  size_t bit_in_byte = *bit_offset % 8; // Bits already consumed in the current byte
  size_t remaining_bytes = (size_t)(end - *p); // Bytes left to read

  const size_t max_bytes = sizeof(bitstream_t); // Maximum bytes to load
  size_t bytes_to_load = (remaining_bytes < max_bytes) ? remaining_bytes : max_bytes;

  bitstream_t result = 0;

  // Load bytes into result, LSB-first
  for (size_t i = 0; i < bytes_to_load; i++) {
    result |= (bitstream_t)(*(*p + i)) << (i * 8);
  }

  // Adjust for bits already consumed in the first byte
  result >>= bit_in_byte;

  // Update the number of bits read
  *nbits = (uint8_t)(bytes_to_load * 8 - bit_in_byte);

  // Advance the bit offset and pointer
  *bit_offset += *nbits;
  *p += *nbits / 8;       // Advance by full bytes
//  *p += *bit_offset / 8;       // Advance by full bytes
  //*bit_offset %= 8;            // Keep leftover bits in bit_offset

  return result;
}

HUFF_EXPORT
bitstream_t
huff_read2(const uint8_t * __restrict stream,
          size_t        * __restrict bit_offset,
          uint8_t       * __restrict nbits,
          size_t                     stream_size)
{
  /* figure out how many bytes remain from bit_offset up to stream_size */
  size_t byte_offset = *bit_offset / 8;
  size_t bit_in_byte = *bit_offset % 8;

  size_t remaining_bytes = (stream_size > byte_offset)
  ? (stream_size - byte_offset)
  : 0;

  /* We'll read at most `max_bytes = sizeof(bitstream_t)` into `result`. */
  const size_t max_bytes = sizeof(bitstream_t);
  size_t bytes_to_load   = (remaining_bytes < max_bytes)
  ? remaining_bytes
  : max_bytes;

  bitstream_t result = 0;

  /*
   * If we have a 256-bit bitstream_t (sizeof(bitstream_t) >= 32)
   * and we have at least 32 bytes left, we can do a 32-byte AVX2/AVX load.
   */
//#if (defined(__AVX2__) || defined(__AVX__)) && !defined(__clang__)  /* _mm256_extract_epi64 not always in clang */
//  if (sizeof(bitstream_t) >= 32 && remaining_bytes >= 32) {
//    __m256i data   = _mm256_loadu_si256((__m256i *)&stream[byte_offset]);
//    /* If bitstream_t is 256 bits, you can store all 256 bits from `data`.
//     For simplicity, assume it's exactly 256 bits: */
//    result = (bitstream_t)data; /* pseudocode if you have a direct cast or store */
//    bytes_to_load = 32;
//  }
//  else
//#endif
//
//    /*
//     * If we have a 128-bit bitstream_t (sizeof(bitstream_t) >= 16)
//     * and we have at least 16 bytes left, we can do a 16-byte SSE2 / NEON load.
//     */
//#if defined(__SSE2__) || defined(__ARM_NEON)
//    if (sizeof(bitstream_t) >= 16 && remaining_bytes >= 16) {
//#  if defined(__SSE2__)
//      __m128i data   = _mm_loadu_si128((__m128i *)&stream[byte_offset]);
//      /* store 128 bits into result (assuming bitstream_t is exactly 128 bits) */
//      result = (bitstream_t)data;  /* pseudocode; you'll need a union or similar */
//#  elif defined(__ARM_NEON)
//      uint8x16_t data = vld1q_u8(&stream[byte_offset]);
//      /* store 128 bits into result (again, pseudocode if you have 128-bit type) */
//      result = *((bitstream_t *)&data);  /* or some union trick */
//#  endif
//      bytes_to_load = 16;
//    }
//    else
//#endif
      {
      /*
       * Fallback: read up to `bytes_to_load` (which is min(remaining_bytes, max_bytes))
       * in a scalar loop. This covers the case of a 64-bit bitstream or insufficient data
       * for a 16/32-byte load.
       */
      for (size_t i = 0; i < bytes_to_load; i++) {
        result |= (bitstream_t)stream[byte_offset + i] << (i * 8);
      }
      }

  /*
   * Now shift out the bits we already used in the partial byte (`bit_in_byte`).
   * e.g. if bit_in_byte=3, shift off the lowest 3 bits from `result`.
   */
  result >>= bit_in_byte;

  /*
   * Indicate how many bits we loaded total: `bytes_to_load * 8`.
   * Then update `bit_offset` to skip those bits minus the partial shift.
   */
  *nbits       = (uint8_t)(bytes_to_load * 8);
  *bit_offset += *nbits - bit_in_byte;

  return result;
}

HUFF_EXPORT
uint_fast16_t
huff_decode_lsb(const huff_table_t * __restrict table,
                bitstream_t                     bitstream,
                uint8_t                         bit_length,
                uint8_t            * __restrict used_bits) {
  uint_fast16_t code, fast_idx, idx;
  uint_fast8_t  l;

  if (bit_length > MAX_CODE_LENGTH) {
    bit_length = MAX_CODE_LENGTH;
  }

  /* fast lookup for short codes (≤ FAST_TABLE_BITS) */
  fast_idx = bitstream & FAST_MASK;
  if (table->fast_table[fast_idx].len <= bit_length) {
    *used_bits = table->fast_table[fast_idx].len;
    return table->fast_table[fast_idx].sym;
  }

  /* fallback for longer codes (LSB-first) */
  for (l = FAST_TABLE_BITS + 1; l <= bit_length; l++) {
    code = bitstream & ((1ULL << l) - 1);
    if (code < table->sentinel_bits[l]) {
      idx        = table->sym_offset[l] + (code - table->sentinel_bits[l - 1]);
      *used_bits = l;
      return table->syms[idx];
    }
  }

  /* decoding failed */
  *used_bits = 0;
  return (uint_fast16_t)-1;
}

HUFF_EXPORT
uint_fast16_t
huff_decode_msb(const huff_table_t * __restrict table,
                bitstream_t                     bitstream,
                uint8_t                         bit_length,
                uint8_t            * __restrict used_bits) {
  uint_fast16_t code, fast_idx, idx;
  uint_fast8_t  l;

  if (bit_length > MAX_CODE_LENGTH) {
    bit_length = MAX_CODE_LENGTH;
  }

  /* fast lookup for short codes (≤ FAST_TABLE_BITS) */
  fast_idx = (bitstream >> FAST_SHIFT) & FAST_MASK;
  if (table->fast_table[fast_idx].len <= bit_length) {
    *used_bits = table->fast_table[fast_idx].len;
    return table->fast_table[fast_idx].sym;
  }

  /* fallback for longer codes (MSB-first) */
  for (l = FAST_TABLE_BITS + 1; l <= bit_length; l++) {
    code = bitstream >> (MAX_CODE_LENGTH - l);
    if (code < table->sentinel_bits[l]) {
      idx        = table->sym_offset[l] + (code - table->sentinel_bits[l - 1]);
      *used_bits = l;
      return table->syms[idx];
    }
  }

  /* decoding failed */
  *used_bits = 0;
  return (uint_fast16_t)-1;
}
