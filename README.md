# 🗜️ Huffman Coding (In Progress)

This library is designed to make **Huffman coding** easy to use, while providing 
optimized and reusable utilities. It aims to simplify the integration of Huffman 
coding into projects, improve maintainability, and ensure robust testing.

Huffman coding is a cornerstone of many compression algorithms and formats, widely used in:
	•	🖼️ Images: JPEG, PNG, WebP, GIF.
	•	🎥 Video: MPEG, H.264, HEVC (H.265).
	•	🎵 Audio: MP3, AAC, Opus.
	•	📦 Compression: ZIP, DEFLATE, Brotli.
	•	🔗 Networking: HTTP/2, HTTP/3 (QUIC) header compression.

This library seeks to serve as a shared foundation for these and other 
applications, offering standardized, high-performance tools for encoding and 
decoding Huffman streams.

## ⚡ Features

- **Efficient Huffman Table Initialization**:
  - Supports sequential or custom symbol tables.
  - Optimized for both **LSB-first** and **MSB-first** bitstreams.
  
- **Convenient LSB and MSB Decoding**:
  - **LSB-first decoding**: Used in formats like DEFLATE (e.g., ZIP, PNG).
  - **MSB-first decoding**: Used in formats like JPEG.

- **Fast Lookup Table**:
  - Precomputes values for short codes to speed up decoding.
  - Highly optimized for common use cases with fallback for longer codes.

- **Flexible APIs**:
  - Dedicated functions for LSB-first and MSB-first streams (`huff_decode_lsb()` and `huff_decode_msb()`).
  - Utilities for bitstream manipulation and bit reversal (`huff_rev_bits()`).

## 🔧 Usage

### Initializing a Huffman Table
The library provides two initialization functions:
- `huff_init_lsb()` for LSB-first bitstreams.
- `huff_init_msb()` for MSB-first bitstreams.

```c
// Example: Initializing a Huffman table for LSB-first bitstreams (e.g., DEFLATE)
huff_table_t table;
uint8_t      lengths[] = {3, 3, 3, 3}; // Bit lengths for each symbol
uint16_t     symbols[] = {0, 1, 2, 3}; // Corresponding symbols
huff_init_lsb(&table, lengths, symbols, 4);
```

## Decoding a Symbol

```c
// LSB-first decoding
bitstream_t   bitstream = 0b10110010; // Example LSB-first bitstream
uint8_t       bit_length = 8; // Number of valid bits
uint8_t       used_bits;
uint_fast16_t symbol = huff_decode_lsb(&table, bitstream, bit_length, &used_bits);

// MSB-first decoding
bitstream_t bitstream_msb = 0b10110010; // Example MSB-first bitstream
symbol = huff_decode_msb(&table, bitstream_msb, bit_length, &used_bits);
```
