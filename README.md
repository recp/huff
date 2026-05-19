# 🗜️ Huffman Coding (In Progress)

This library is designed to make **Huffman coding** easy to use, while providing 
optimized and reusable utilities. It aims to simplify the integration of Huffman 
coding into projects, improve maintainability, and ensure robust testing.

This library seeks to serve as a shared foundation for these and other 
applications, offering standardized, high-performance tools for encoding and 
decoding Huffman streams.

🚨 Don't use this in production until tests are ready

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
```

## TODO

- [x] lsb
- [ ] sub tables?
- [ ] msb
- [ ] tests
- [ ] build
- [ ] documentation


