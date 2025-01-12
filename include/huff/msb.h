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

/**
 * TODO:
 */

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
//HUFF_EXPORT
//uint_fast16_t
//huff_decode_msb(const huff_table_t * __restrict table,
//                bitstream_t                     bitstream,
//                uint8_t                         bit_length,
//                uint8_t            * __restrict used_bits);

#ifdef __cplusplus
}
#endif
#endif /* huff_msb_h */
