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

typedef struct huff_table_t {
  uint16_t fast_table[256];   /* Fast lookup table for short codes (up to 8 bits) */
  uint8_t  fast_length[256];  /* Code lengths for fast decoding table */

  uint32_t sentinel_bits[16]; /* Boundaries for identifying code lengths (up to 15 bits) */
  uint16_t sym_offset[16];    /* Start offsets for symbols by code length */

  uint16_t *syms;             /* Array of decoded symbols */
  size_t    num_symbols;      /* Number of symbols in the table */
} huff_table_t;

/*!
 * @brief initializes huffman table
 *
 *  NOTE: you can pass NULL to symbols for sequential symbols e.g. DEFLATE
 *
 * @param[in, out]  table    huffman table to be initialized
 * @param[in]       lengths  array of bit lengths for each symbol
 * @param[in]       symbols  array of symbols corresponding to the provided lengths
 * @param[in]       n        number of symbols in the lengths and symbols arrays
 */
HUFF_EXPORT
void
huff_init(huff_table_t   * __restrict table,
          const uint8_t  * __restrict lengths,
          const uint16_t * __restrict symbols,
          uint16_t                    n);

#ifdef __cplusplus
}
#endif
#endif /* huff_h */
