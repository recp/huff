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

#ifdef __GNUC__
#  define unlikely(expr) __builtin_expect(!!(expr), 0)
#  define likely(expr)   __builtin_expect(!!(expr), 1)
#else
#  define unlikely(expr) (expr)
#  define likely(expr)   (expr)
#endif

#if defined(__SIZEOF_INT128__)
typedef __uint128_t big_int_t;
#else
typedef uintmax_t   big_int_t;
#endif

#ifdef HUFF_ENABLE_BIG_BITSTREAM
typedef big_int_t     bitstream_t;
#else
typedef uint_fast64_t bitstream_t;
#endif

/* 15: DEFLATE, 16: JPEG */
#define HUFF_MAX_CODE_LENGTH  16  /* maximum length (bits) of Huffman codes     */
#define HUFF_FAST_TABLE_BITS  8   /* number of bits used for fast lookup        */
#define HUFF_FAST_TABLE_SIZE  (1U << HUFF_FAST_TABLE_BITS)
#define HUFF_FAST_SHIFT       (HUFF_MAX_CODE_LENGTH - HUFF_FAST_TABLE_BITS)
#define HUFF_MAX_CODES        288

typedef struct {uint64_t base:16,bits:8,mask:24;} huff_ext_t;

typedef struct huff_fast_entry_t {
  uint8_t  len;
  uint8_t  rev;
  uint16_t sym;
} huff_fast_entry_t;

/* extended fast table entry with extra bits info */
typedef struct huff_fast_entry_ext_t {
  uint8_t  len;
  uint8_t  rev;
  uint16_t sym;
  uint32_t value;
  uint32_t mask;
  uint8_t  total;
} huff_fast_entry_ext_t;

typedef struct huff_table_t {
  HUFF_ALIGN(32) huff_fast_entry_t fast[HUFF_FAST_TABLE_SIZE];

  union {
    uint16_t sentinels[HUFF_MAX_CODE_LENGTH + 1];
    uint16_t maxcode[HUFF_MAX_CODE_LENGTH   + 1];
  } HUFF_ALIGN(32);

  union {
    uint16_t offsets[HUFF_MAX_CODE_LENGTH + 1];
    uint16_t mincode[HUFF_MAX_CODE_LENGTH + 1];
  } HUFF_ALIGN(32);

  /* uint16_t *sym; to allow dynmaic syms if no overhead? */
  HUFF_ALIGN(32) uint16_t syms[HUFF_MAX_CODES];
} huff_table_t;

/* extended table for extra bits (e.g length/distance in deflate) */
typedef struct huff_table_ext_t {
  HUFF_ALIGN(32) huff_fast_entry_ext_t fast[HUFF_FAST_TABLE_SIZE];

  union {
    uint16_t sentinels[HUFF_MAX_CODE_LENGTH + 1];
    uint16_t maxcode[HUFF_MAX_CODE_LENGTH   + 1];
  } HUFF_ALIGN(32);

  union {
    uint16_t offsets[HUFF_MAX_CODE_LENGTH + 1];
    uint16_t mincode[HUFF_MAX_CODE_LENGTH + 1];
  } HUFF_ALIGN(32);

  /* uint16_t *sym; to allow dynmaic syms if no overhead? */
  HUFF_ALIGN(32) uint16_t syms[HUFF_MAX_CODES];
  const huff_ext_t       *extras; /* extra bits info                        */
  int                     offset; /* 257 for lit/len, 0 for dist in deflate */
} huff_table_ext_t;

#include "read.h"
#include "rev.h"
#include "lsb.h"
#include "msb.h"

#ifdef __cplusplus
}
#endif
#endif /* huff_h */
