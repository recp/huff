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

#include "../include/huff.h"

/* 15: DEFLATE, 16: JPEG */
#define MAX_CODE_LENGTH 16  /* Maximum length of Huffman codes */
#define FAST_TABLE_BITS  8  /* Number of bits used for fast lookup */

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
