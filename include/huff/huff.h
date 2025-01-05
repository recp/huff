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

#ifndef huff_h
#define huff_h
#ifdef __cplusplus
extern "C" {
#endif

#ifndef _USE_MATH_DEFINES
#  define _USE_MATH_DEFINES       /* for windows */
#endif

#ifndef _CRT_SECURE_NO_WARNINGS
#  define _CRT_SECURE_NO_WARNINGS /* for windows */
#endif

#ifndef _CRT_NONSTDC_NO_DEPRECATE
#  define _CRT_NONSTDC_NO_DEPRECATE /* for windows */
#endif

/* since C99 or compiler ext */
#include <stdint.h>
#include <stddef.h>
#include <float.h>
#include <stdbool.h>
#include <errno.h>
#include <stdlib.h>

#ifdef DEBUG
#  include <assert.h>
#  include <stdio.h>
#endif

#if defined(_MSC_VER) || defined(__MINGW32__) || defined(__MINGW64__)
#  define HUFF_WINAPI
#  pragma warning (disable : 4068) /* disable unknown pragma warnings */
#endif

#if defined(_MSC_VER) || defined(__MINGW32__) || defined(__MINGW64__)
#  ifdef HUFF_STATIC
#    define HUFF_EXPORT
#  elif defined(HUFF_EXPORTS)
#    define HUFF_EXPORT __declspec(dllexport)
#  else
#    define HUFF_EXPORT __declspec(dllimport)
#  endif
#  define HUFF_HIDE
#else
#  define HUFF_EXPORT   __attribute__((visibility("default")))
#  define HUFF_HIDE     __attribute__((visibility("hidden")))
#endif

#if defined(_MSC_VER)
#  define HUFF_INLINE   __forceinline
#  define HUFF_ALIGN(X) __declspec(align(X))
#else
#  define HUFF_ALIGN(X) __attribute((aligned(X)))
#  define HUFF_INLINE   static inline __attribute((always_inline))
#endif

#ifndef __has_builtin
#  define __has_builtin(x) 0
#endif

#if defined(__SIZEOF_INT128__)
typedef __uint128_t big_int_t;
#else
typedef uintmax_t   big_int_t;
#endif

#ifdef ENABLE_BIG_BITSTREAM
typedef big_int_t   bitstream_t;
#else
typedef uint64_t    bitstream_t;
#endif

/* 15: DEFLATE, 16: JPEG */
#define MAX_CODE_LENGTH  16  /* maximum length (bits) of Huffman codes     */
#define FAST_TABLE_BITS  8   /* number of bits used for fast lookup        */
#define FAST_TABLE_SIZE  (1U << FAST_TABLE_BITS)
#define FAST_SHIFT       (MAX_CODE_LENGTH - FAST_TABLE_BITS)

typedef struct huff_table_t {
  /* fast lookup table for short codes (up to 8 bits) */
  struct {
    uint16_t sym : 9;  /* Symbol (up to 288 symbols for DEFLATE) */
    uint16_t len : 7;  /* Codeword length (0 means invalid) */
  } fast_table[1U << FAST_TABLE_BITS];

  uint32_t sentinel_bits[MAX_CODE_LENGTH + 1];  /* boundaries for identifying code lengths */
  uint16_t sym_offset[MAX_CODE_LENGTH + 1];     /* start offsets for symbols by code length */

  uint16_t *syms;                               /* array of decoded symbols */
  size_t    num_symbols;                        /* number of symbols in the table */

  /* uint8_t max_bitlen; to repalce MAX_CODE_LENGTH */
} huff_table_t;

/*!
 * @brief Initializes a Huffman table for decoding LSB-first bitstreams.
 *
 * This function initializes a Huffman table using the provided symbol lengths
 * and symbols. The table is designed to decode bitstreams in LSB-first order,
 * where the least significant bits are processed first.
 *
 * If the `symbols` array is `NULL`, the symbols will be generated as a sequential
 * range starting from 0, which is common for formats like DEFLATE.
 *
 * @param[in, out]  table    Pointer to the Huffman table to be initialized.
 * @param[in]       lengths  Array of bit lengths for each symbol. The length of
 *                           each symbol must not exceed `MAX_CODE_LENGTH`.
 * @param[in]       symbols  Array of symbols corresponding to the provided lengths.
 *                           Pass `NULL` for sequential symbols.
 * @param[in]       n        Number of symbols in the `lengths` (and `symbols`, if provided) array.
 *
 * @note This function supports LSB-first bitstreams. For MSB-first bitstreams,
 *       use `huff_init_msb()`.
 */
HUFF_EXPORT
void
huff_init_lsb(huff_table_t   * __restrict table,
              const uint8_t  * __restrict lengths,
              const uint16_t * __restrict symbols,
              uint16_t                    n);

/*!
 * @brief Initializes a Huffman table for decoding MSB-first bitstreams.
 *
 * This function initializes a Huffman table using the provided symbol lengths
 * and symbols. The table is designed to decode bitstreams in MSB-first order,
 * where the most significant bits are processed first.
 *
 * If the `symbols` array is `NULL`, the symbols will be generated as a sequential
 * range starting from 0, which is common for formats like DEFLATE.
 *
 * @param[in, out]  table    Pointer to the Huffman table to be initialized.
 * @param[in]       lengths  Array of bit lengths for each symbol. The length of
 *                           each symbol must not exceed `MAX_CODE_LENGTH`.
 * @param[in]       symbols  Array of symbols corresponding to the provided lengths.
 *                           Pass `NULL` for sequential symbols.
 * @param[in]       n        Number of symbols in the `lengths` (and `symbols`, if provided) array.
 *
 * @note This function supports MSB-first bitstreams. For LSB-first bitstreams,
 *       use `huff_init_lsb()`.
 */
HUFF_EXPORT
void
huff_init_msb(huff_table_t   * __restrict table,
              const uint8_t  * __restrict lengths,
              const uint16_t * __restrict symbols,
              uint16_t                    n);

/**
 * @brief Reads a bitstream_t worth of bits from the input stream.
 *
 * @param[in, out] stream     Pointer to the input stream.
 * @param[in, out] bit_offset Current bit offset within the stream.
 * @param[out]     nbits      Number of loaded bits
 * @return                    A `bitstream_t` containing the requested bits.
 */
HUFF_EXPORT
bitstream_t
huff_read_lsb(const uint8_t * __restrict stream,
              size_t        * __restrict bit_offset,
              uint8_t       * __restrict nbits,
              size_t                     stream_size);

/**
 * @brief Reads a bitstream_t worth of bits from the input stream.
 *
 * @param[in, out] stream     Pointer to the input stream.
 * @param[in, out] bit_offset Current bit offset within the stream.
 * @param[out]     nbits      Number of loaded bits
 * @return                    A `bitstream_t` containing the requested bits.
 */
//HUFF_EXPORT
//bitstream_t
//huff_read(const uint8_t * __restrict stream,
//          size_t        * __restrict bit_offset,
//          uint8_t       * __restrict nbits,
//          size_t                     stream_size);
HUFF_EXPORT
bitstream_t
huff_read(const uint8_t ** __restrict p,
          size_t        * __restrict bit_offset,
          uint8_t       * __restrict nbits,
          const uint8_t * __restrict end);
/**
 * @brief Reverses the bit order of the input value.
 *
 * Reverses the bit order of a `big_int_t` value, useful for converting
 * an LSB-first bitstream into MSB-first order.
 *
 * @param[in] x The input value whose bits are to be reversed.
 * @return The input value with its bit order reversed.
 *
 * @example
 * big_int_t original = 0b10110010; // Example 8-bit value
 * big_int_t reversed = huff_rev_bits(original);
 * // reversed is now 0b01001101
 */
HUFF_EXPORT
bitstream_t
huff_rev_bits(bitstream_t x);

/**
 * @brief Decodes a single symbol from a Huffman-encoded bitstream (LSB-first).
 *
 * This function decodes a symbol using a pre-initialized Huffman table. The
 * bitstream is expected to be in LSB-first order, where the least significant
 * bits are processed first. If the bitstream is in MSB-first order, it should
 * be reversed before calling this function. You can use `huff_rev_bits()` to
 * reverse the bitstream if needed.
 *
 * @param[in]     table       Pointer to the initialized Huffman table.
 * @param[in]     bitstream   The bitstream to decode, in LSB-first order.
 * @param[in]     bit_length  The number of valid bits in the bitstream.
 * @param[out]    used_bits   Pointer to store the number of bits used to decode.
 *
 * @return The decoded symbol, or `(uint_fast16_t)-1` if decoding fails.
 *
 * @note The caller is responsible for ensuring the bitstream contains
 *       enough valid bits for decoding a symbol.
 */
HUFF_EXPORT
uint_fast16_t
huff_decode_lsb(const huff_table_t * __restrict table,
                bitstream_t                     bitstream,
                uint8_t                         bit_length,
                uint8_t            * __restrict used_bits);

/**
 * @brief Decodes a single symbol from a Huffman-encoded bitstream (MSB-first).
 *
 * This function decodes a symbol using a pre-initialized Huffman table. The
 * bitstream is expected to be in MSB-first order, where the most significant
 * bits are processed first. If the bitstream is in LSB-first order, it should
 * be reversed before calling this function. You can use `huff_rev_bits()` to
 * reverse the bitstream if needed.
 *
 * @param[in]     table       Pointer to the initialized Huffman table.
 * @param[in]     bitstream   The bitstream to decode, in MSB-first order.
 * @param[in]     bit_length  The number of valid bits in the bitstream.
 * @param[out]    used_bits   Pointer to store the number of bits used to decode.
 *
 * @return The decoded symbol, or `(uint_fast16_t)-1` if decoding fails.
 *
 * @note The caller is responsible for ensuring the bitstream contains
 *       enough valid bits for decoding a symbol.
 */
HUFF_EXPORT
uint_fast16_t
huff_decode_msb(const huff_table_t * __restrict table,
                bitstream_t                     bitstream,
                uint8_t                         bit_length,
                uint8_t            * __restrict used_bits);

#ifdef __cplusplus
}
#endif
#endif /* huff_h */
