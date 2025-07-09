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

#ifndef huff_lsb_h
#define huff_lsb_h
#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>

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
 * @param[out]    used        Pointer to store the number of bits used to decode.
 *
 * @return The decoded symbol, or `(uint_fast16_t)-1` if decoding fails.
 *
 * @note The caller is responsible for ensuring the bitstream contains
 *       enough valid bits for decoding a symbol.
 */

HUFF_INLINE
uint_fast16_t
huff_decode_lsb(const huff_table_t * __restrict table,
                bitstream_t                     bitstream,
                uint8_t                         bit_length,
                uint8_t            * __restrict used) {
  huff_fast_entry_t fe;
  uint16_t          code, bits;
  uint8_t           bits8, l;

  (void)bit_length;

  /* align bits so LSB is always in the first position */
  bits8 = (uint8_t)bitstream;
  fe    = table->fast[bits8];

  if (likely(fe.len)) {
    *used = fe.len;
    return fe.sym;
  }

  bits = (uint16_t)(bitstream >> 8);
  code = fe.rev; /* huff_rev8full(bits8); */

  /* check length and add next bit from LSB to MSB position of our code */
#define CHECK_LENGTH(l)                                                       \
  code = (code << 1) | (bits & 1);                                            \
  if (code < table->sentinels[l]) {                                           \
    *used=l;return table->syms[(uint16_t)(table->offsets[l]+code)];           \
  }                                                                           \
  bits>>=1;

  /*
#define CHECK_LENGTH(l)                                                       \
  if (code < table->sentinels[l]) {                                           \
    *used=l;return table->syms[(uint16_t)(table->offsets[l]+code)];           \
  } code=(code<<1)|((bits>>=1)&1);
  */

  CHECK_LENGTH(9);
  CHECK_LENGTH(10);
  CHECK_LENGTH(11);
  CHECK_LENGTH(12);

  for (l = 13; l <= HUFF_MAX_CODE_LENGTH; l++) { CHECK_LENGTH(l); }

#undef CHECK_LENGTH
  *used = 0;
  return -1;
}

#define HUFF_DECODE_LSB(table, bitstream, bit_length, used, result) do {     \
  huff_fast_entry_t fe_;                                                      \
  uint16_t l_, code_, bits_;                                                  \
  uint8_t bits8_;                                                            \
                                                                             \
  /* align bits so LSB is always in the first position */                    \
  bits8_ = (bitstream);                                                      \
  fe_    = (table)->fast[bits8_];                                      \
                                                                             \
  if (likely(fe_.len)) {                                                     \
    *(used) = fe_.len;                                                       \
    result = fe_.sym;                                                        \
  } else {                                                                   \
    bits_ = (bitstream) >> 8;                                               \
    code_ = fe_.rev;                                                         \
    /* Handle remaining lengths 13+ */                                      \
    for (l_ = 9; l_ <= HUFF_MAX_CODE_LENGTH; l_++) {                           \
      code_ = (code_ << 1) | (bits_ & 1);                                   \
      if (code_ < (table)->sentinels[l_]) {                                \
        *(used) = l_;                                                       \
        result = (table)->syms[(uint16_t)((table)->offsets[l_]+code_)];    \
        break;                                                              \
      }                                                                     \
      bits_ >>= 1;                                                          \
    }                                                                       \
                                                                            \
    if (l_ > HUFF_MAX_CODE_LENGTH) {                                            \
      *(used) = 0;                                                          \
      result = -1;                                                          \
    }                                                                       \
  }                                                                         \
} while(0)

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
 *                           each symbol must not exceed `HUFF_MAX_CODE_LENGTH`.
 * @param[in]       symbols  Array of symbols corresponding to the provided lengths.
 *                           Pass `NULL` for sequential symbols.
 * @param[in]       n        Number of symbols in the `lengths` (and `symbols`, if provided) array.
 *
 * @note This function supports LSB-first bitstreams. For MSB-first bitstreams,
 *       use `huff_init_msb()`.
 */
HUFF_INLINE
bool
huff_init_lsb(huff_table_t   * restrict table,
              const uint8_t  * restrict lengths,
              const uint16_t * restrict symbols,
              uint16_t                  n) {
  uint_fast16_t l, i, prev_code = 0, prev_sym_idx = 0;
  uint_fast16_t count[HUFF_MAX_CODE_LENGTH   + 1] = {0};
  uint_fast16_t code[HUFF_MAX_CODE_LENGTH    + 1];
  uint_fast16_t sym_idx[HUFF_MAX_CODE_LENGTH + 1];

  (void)symbols; /* symbols auto-generated when NULL */

  /**
   * currently table->syms is fixed size array
   *
   * if (!table->syms) {
   *   table->syms = calloc(n, sizeof(uint16_t));
   *   if (!table->syms) return false;
   * }
   */

  /* mark fast table as invalid */
  for (i = 0; i < (1U << HUFF_FAST_TABLE_BITS); i++) {
    table->fast[i].len = 0;
  }

  for (i = 0; i < n; i++) {
    count[lengths[i]]++;
  }

  count[0] = code[0] = sym_idx[0] = 0;

  for (l = 1; l <= HUFF_MAX_CODE_LENGTH; l++) {
    code[l]             = (prev_code   + count[l - 1]) << 1;
    sym_idx[l]          = prev_sym_idx + count[l - 1];
    table->sentinels[l] = code[l]      + count[l];
    table->offsets[l]   = sym_idx[l]   - code[l];

    prev_code           = code[l];
    prev_sym_idx        = sym_idx[l];
  }

  /* fast table and symbol array */
  for (i = 0; i < n; i++) {
    if ((l = lengths[i])) {
      table->syms[sym_idx[l]++] = i;

      /* fill the fast table for short codes */
      if (l <= HUFF_FAST_TABLE_BITS) {
        uint8_t  padlen = HUFF_FAST_TABLE_BITS - l;
        uint16_t pad;
        uint8_t  code8, index;

        code8 = huff_rev8((uint8_t)code[l]++, l);

        for (pad = 0; pad < (1U << padlen); pad++) {
          index = (uint8_t)(code8 | (pad << l));
          table->fast[index].sym = i;
          table->fast[index].len = (uint8_t)l;
        }
      }
    }
  }

  /* reverse bits only for slow path entries */
  for (i = 0; i < (1U << HUFF_FAST_TABLE_BITS); i++) {
    if (table->fast[i].len == 0) {
      table->fast[i].rev = huff_rev8full((uint8_t)i);
    }
  }
  return true;
}

HUFF_INLINE
bool
huff_init_fast_lsb(huff_fast_entry_t         fast[HUFF_FAST_TABLE_SIZE],
                   const uint8_t  * restrict lengths,
                   const uint16_t * restrict symbols,
                   uint16_t                  n) {
  uint_fast16_t l, i, prev_code = 0;
  uint_fast16_t count[HUFF_FAST_TABLE_BITS + 1] = {0};
  uint_fast16_t code[HUFF_FAST_TABLE_BITS  + 1];

  (void)symbols; /* symbols auto-generated when NULL */

  /* mark fast table as invalid */
  for (i = 0; i < (1U << HUFF_FAST_TABLE_BITS); i++) {
    fast[i].len = 0;
  }

  /* only count lengths <= HUFF_FAST_TABLE_BITS */
  for (i = 0; i < n; i++) {
    l = lengths[i];
    if (l && l <= HUFF_FAST_TABLE_BITS) {
      count[l]++;
    }
  }

  /* compute first code for each length */
  code[0] = count[0] = 0;
  for (l = 1; l <= HUFF_FAST_TABLE_BITS; l++) {
    code[l]   = (prev_code + count[l - 1]) << 1;
    prev_code = code[l];
  }

  /* fill fast table */
  for (i = 0; i < n; i++) {
    if ((l = lengths[i]) && l <= HUFF_FAST_TABLE_BITS) {
      uint8_t  code8 = huff_rev8((uint8_t)code[l]++, l);
      uint8_t  padlen = HUFF_FAST_TABLE_BITS - l;
      uint16_t pad;

      for (pad = 0; pad < (1U << padlen); pad++) {
        uint8_t index = code8 | (pad << l);
        fast[index].sym = i;
        fast[index].len = l;
      }
    }
  }

  return true;
}

HUFF_INLINE
bool
huff_init_lsb_ext(huff_table_ext_t   * __restrict table,
                  const uint8_t      * __restrict lengths,
                  const uint16_t     * __restrict symbols,
                  const huff_ext_t   * __restrict extras,
                  uint16_t                      n) {
  uint_fast16_t l, i, prev_code, prev_sym_idx;
  uint_fast16_t count[HUFF_MAX_CODE_LENGTH   + 1] = {0};
  uint_fast16_t code[HUFF_MAX_CODE_LENGTH    + 1];
  uint_fast16_t sym_idx[HUFF_MAX_CODE_LENGTH + 1];

  (void)symbols; /* symbols auto-generated when NULL */

  prev_code    = 0;
  prev_sym_idx = 0;

  /* mark fast table as invalid */
  for (i = 0; i < (1U << HUFF_FAST_TABLE_BITS); i++)
    table->fast[i].len = 0;

  for (i = 0; i < n; i++)
    count[lengths[i]]++;

  count[0] = code[0] = sym_idx[0] = 0;

  for (l = 1; l <= HUFF_MAX_CODE_LENGTH; l++) {
    code[l]             = (prev_code   + count[l - 1]) << 1;
    sym_idx[l]          = prev_sym_idx + count[l - 1];
    table->sentinels[l] = code[l]      + count[l];
    table->offsets[l]   = sym_idx[l]   - code[l];

    prev_code    = code[l];
    prev_sym_idx = sym_idx[l];
  }

  table->extras = extras;
  table->offset = 0;

  /* fast table and symbol array */
  for (i = 0; i < n; i++) {
    if ((l = lengths[i])) {
      table->syms[sym_idx[l]++] = i;

      if (l <= HUFF_FAST_TABLE_BITS) {
        huff_ext_t ext;
        uint8_t    padlen, code8, index;
        uint16_t   pad;

        code8  = huff_rev8((uint8_t)code[l]++, l);
        padlen = HUFF_FAST_TABLE_BITS - l;

        for (pad = 0; pad < (1U << padlen); pad++) {
          index = (uint8_t)(code8 | (pad << l));
          ext   = extras[i];

          table->fast[index].sym   = i;
          table->fast[index].len   = (uint8_t)l;

          table->fast[index].value = (uint32_t)ext.base;
          table->fast[index].total = l + (uint8_t)ext.bits;
          table->fast[index].mask  = (1U<<ext.bits)-1;
        }
      }
    }
  }

  /* reverse bits only for slow path entries */
  for (i = 0; i < (1U << HUFF_FAST_TABLE_BITS); i++)
    if (!table->fast[i].len)
      table->fast[i].rev = huff_rev8full((uint8_t)i);

  return true;
}

HUFF_INLINE
bool
huff_init_lsb_extof(huff_table_ext_t   * __restrict table,
                    const uint8_t      * __restrict lengths,
                    const uint16_t     * __restrict symbols,
                    const huff_ext_t   * __restrict extras,
                    int                           offset,
                    uint16_t                      n) {
  uint_fast16_t l, i, prev_code, prev_sym_idx;
  uint_fast16_t count[HUFF_MAX_CODE_LENGTH   + 1] = {0};
  uint_fast16_t code[HUFF_MAX_CODE_LENGTH    + 1];
  uint_fast16_t sym_idx[HUFF_MAX_CODE_LENGTH + 1];

  (void)symbols; /* symbols auto-generated when NULL */

  prev_code    = 0;
  prev_sym_idx = 0;

  /* mark fast table as invalid */
  for (i = 0; i < (1U << HUFF_FAST_TABLE_BITS); i++)
    table->fast[i].len = 0;

  for (i = 0; i < n; i++)
    count[lengths[i]]++;

  count[0] = code[0] = sym_idx[0] = 0;

  for (l = 1; l <= HUFF_MAX_CODE_LENGTH; l++) {
    code[l]             = (prev_code   + count[l - 1]) << 1;
    sym_idx[l]          = prev_sym_idx + count[l - 1];
    table->sentinels[l] = code[l]      + count[l];
    table->offsets[l]   = sym_idx[l]   - code[l];

    prev_code    = code[l];
    prev_sym_idx = sym_idx[l];
  }

  table->extras = extras;
  table->offset = offset;

  /* fast table and symbol array */
  for (i = 0; i < n; i++) {
    if ((l = lengths[i])) {
      table->syms[sym_idx[l]++] = i;

      if (l <= HUFF_FAST_TABLE_BITS) {
        huff_ext_t ext;
        uint8_t    padlen, code8, index;
        uint16_t   pad;

        code8  = huff_rev8((uint8_t)code[l]++, l);
        padlen = HUFF_FAST_TABLE_BITS - l;

        for (pad = 0; pad < (1U << padlen); pad++) {
          index                  = (uint8_t)(code8 | (pad << l));
          table->fast[index].sym = i;
          table->fast[index].len = (uint8_t)l;
          
          if ((int)i >= offset) {
            ext                      = extras[i - offset];
            table->fast[index].value = (uint32_t)ext.base;
            table->fast[index].total = l + (uint8_t)ext.bits;
            table->fast[index].mask  = (1U<<ext.bits)-1;
          } else {
            table->fast[index].value = 0;
            table->fast[index].mask  = 0;
            table->fast[index].total = l;
          }
        }
      }
    }
  }

  /* reverse bits only for slow path entries */
  for (i = 0; i < (1U << HUFF_FAST_TABLE_BITS); i++)
    if (!table->fast[i].len)
      table->fast[i].rev = huff_rev8full((uint8_t)i);

  return true;
}

HUFF_INLINE
unsigned
huff_decode_lsb_ext(const huff_table_ext_t * __restrict table,
                    bitstream_t                         bitstream,
                    uint8_t                * __restrict used) {
  huff_fast_entry_ext_t fe;
  huff_ext_t            ext;
  uint16_t              l, code, bits, sym;
  uint8_t               bits8;

  /* align bits so LSB is always in the first position */
  bits8 = (uint8_t)bitstream;
  fe    = table->fast[bits8];

  if (likely(fe.len)) {
    *used = fe.total;
    return fe.value + (fe.mask & (unsigned)(bitstream >> fe.len));
  }

  bits = (uint16_t)(bitstream >> 8);
  code = fe.rev; /* huff_rev8full(bits8); */

  for (l = 9; l <= HUFF_MAX_CODE_LENGTH; l++) {
    code = (code << 1) | (bits & 1);
    if (code < table->sentinels[l]) {
      sym    = table->syms[(uint16_t)(table->offsets[l] + code)];
      ext    = table->extras[sym];
      *used  = l + (uint8_t)ext.bits;
      return (unsigned)(ext.base + (ext.mask & (unsigned)(bitstream >> l)));
    }
    bits >>= 1;
  }

  *used = 0;
  return -1;
}

HUFF_INLINE
uint_fast16_t
huff_decode_lsb_extof(const huff_table_ext_t * __restrict table,
                      bitstream_t                         bitstream,
                      uint8_t                * __restrict used,
                      unsigned               * __restrict value,
                      int                                 offset) {
  huff_fast_entry_ext_t fe;
  huff_ext_t            ext;
  uint16_t              l, code, bits, sym;
  uint8_t               bits8;

  /* align bits so LSB is always in the first position */
  bits8 = (uint8_t)bitstream;
  fe    = table->fast[bits8];

  if (likely(fe.len)) {
    *used  = fe.total;
    *value = fe.value + (fe.mask & (unsigned)(bitstream >> fe.len));
    return fe.sym;
  }

  bits = (uint16_t)(bitstream >> 8);
  code = fe.rev; /* huff_rev8full(bits8); */

  /* slow path */
  for (l = 9; l <= HUFF_MAX_CODE_LENGTH; l++) {
    code = (code << 1) | (bits & 1);
    if (code < table->sentinels[l]) {
      sym = table->syms[(uint16_t)(table->offsets[l] + code)];
      if (likely(sym >= offset)) {
        ext    = table->extras[sym - offset];
        *value = (unsigned)(ext.base + (ext.mask & (unsigned)(bitstream >> l)));
        *used  = l + (uint8_t)ext.bits;
      } else {
        *used  = (uint8_t)l;
        *value = 0;
      }
      return sym;
    }
    bits >>= 1;
  }

  *used = 0;
  return -1;
}

#ifdef __cplusplus
}
#endif
#endif /* huff_lsb_h */
