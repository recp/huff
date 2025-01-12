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

#ifndef huff_msb_h
#define huff_msb_h
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__ARM_NEON)
#include <arm_neon.h>
#elif defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

static inline uint8_t huff_rev8(uint8_t b, int len) {
#if defined(__x86_64__) || defined(__i386__)
  return (uint8_t)__builtin_bitreverse8(b) >> (8 - len);
#elif defined(__ARM_NEON) &&  defined(__aarch64__)
  return vget_lane_u8(vrbit_u8(vdup_n_u8(b)), 0) >> (8 - len);
#else
  b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;
  b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
  b = (b & 0xAA) >> 1 | (b & 0x55) << 1;
  return b >> (8 - len);
#endif
}

static inline uint8_t huff_rev8full(uint8_t b) {
#if defined(__x86_64__) || defined(__i386__)
  return (uint8_t)__builtin_bitreverse8(b);
#elif defined(__ARM_NEON) &&  defined(__aarch64__)
  return vget_lane_u8(vrbit_u8(vdup_n_u8(b)), 0);
#else
  b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;
  b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
  return (b & 0xAA) >> 1 | (b & 0x55) << 1;
#endif
}

#ifdef __cplusplus
}
#endif
#endif /* huff_msb_h */
