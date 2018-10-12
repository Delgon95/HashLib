#ifndef CRC_H_
#define CRC_H_

#include <chrono>       // std::chrono::system_clock / duration
#include <numeric>      // std::accumulate
#include <type_traits>  // std::enable_if

namespace hash {

enum CrcChunks {
  BYTE_BY_BYTE,
  CHUNKS_1x32b,
  CHUNKS_2x32b,
  CHUNKS_4x32b,
  CHUNKS_8x32b
};

struct OptionsCrc {
  constexpr OptionsCrc(uint64_t _polynomial, uint64_t _initial_crc,
                       uint64_t _xor_output, bool _reverse_data,
                       bool _reverse_out, CrcChunks _chunks = CHUNKS_4x32b);
  static constexpr OptionsCrc Crc16();
  static constexpr OptionsCrc Crc16_CCITT();
  static constexpr OptionsCrc Crc32();
  static constexpr OptionsCrc Crc64();
  static constexpr OptionsCrc Crc64_ISO();
  const uint64_t polynomial = 0;
  const uint64_t initial_crc = 0;
  const uint64_t xor_output = 0;
  const bool reverse_data = false;
  const bool reverse_out = false;
  CrcChunks chunks = CHUNKS_4x32b;
};

constexpr OptionsCrc::OptionsCrc(uint64_t _polynomial, uint64_t _initial_crc,
                                 uint64_t _xor_output, bool _reverse_data,
                                 bool _reverse_out, CrcChunks _chunks)
    : polynomial(_polynomial),
      initial_crc(_initial_crc),
      xor_output(_xor_output),
      reverse_data(_reverse_data),
      reverse_out(_reverse_out),
      chunks(_chunks) {}

template <typename T>
class Crc {
 public:
  Crc() = delete;
  Crc(const OptionsCrc& options = OptionsCrc::Crc32());
  ~Crc() = default;

  // Consume specified number of elements from given array of any type.
  template <typename Y>
  void Consume(const Y* data, size_t size);
  // Use not default processing method.
  template <typename Y>
  void Consume(const Y* data, size_t size, CrcChunks chunks);

  // Optimize CRC calculation by selecting processing method with most
  // performance. Calculate performance for 128 packages with 8kB of data.
  void Optimize(uint64_t buffer_size = 8 * 1024 - 1, uint64_t repeats = 128);

  // Retrieve CRC value.
  T crc() const noexcept;
  // Reset CRC value to initial one.
  void reset() noexcept;

 private:
  void GenerateLookupTable() noexcept;
  T CalculateTableValue(uint8_t value) const noexcept;
  T ReverseBits(T value, size_t bits = sizeof(T) * 8) const noexcept;
  template <typename Y>
  void Consume_byte_by_byte(const Y* data, size_t size);
  template <typename Y>
  void Consume_1x32b(const Y* data, size_t size);
  template <typename Y>
  void Consume_2x32b(const Y* data, size_t size);
  template <typename Y>
  void Consume_4x32b(const Y* data, size_t size);
  template <typename Y>
  void Consume_8x32b(const Y* data, size_t size);

  T crc_;
  const T initial_crc_;
  const T xor_output_;
  const T polynomial_;
  const bool reverse_data_;
  const bool reverse_out_;
  T lookup_table_[32][256];
  CrcChunks chunks_;
};

template <typename T>
Crc<T>::Crc(const OptionsCrc& options)
    : crc_((options.reverse_data)
               ? ReverseBits(static_cast<T>(options.initial_crc))
               : static_cast<T>(options.initial_crc)),
      initial_crc_(static_cast<T>(options.initial_crc)),
      xor_output_(static_cast<T>(options.xor_output)),
      polynomial_(static_cast<T>(options.polynomial)),
      reverse_data_(options.reverse_data),
      reverse_out_(options.reverse_out),
      chunks_(options.chunks) {
  GenerateLookupTable();
}

// Generates a lookup table for the checksums of all 8-bit values.
// Values can be genrated using reversed bit ordering depending on the
// standard.
template <typename T>
void Crc<T>::GenerateLookupTable() noexcept {
  for (size_t i = 0; i < 256; ++i) {
    lookup_table_[0][i] = CalculateTableValue(static_cast<uint8_t>(i));
  }
  // Precompute additional values, to allow computation with more 64b
  // chunks of data at the same time.
  if (reverse_data_) {
    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 1; j < 32; ++j) {
        lookup_table_[j][i] = (lookup_table_[j - 1][i] >> 8) ^
                              lookup_table_[0][lookup_table_[j - 1][i] & 0xFF];
      }
    }
  } else {
    constexpr static uint8_t shift = (sizeof(T) * 8) - 8;
    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 1; j < 32; ++j) {
        lookup_table_[j][i] =
            (lookup_table_[j - 1][i] << 8) ^
            lookup_table_[0][(lookup_table_[j - 1][i] >> shift) & 0xFF];
      }
    }
  }
}

template <typename T>
T Crc<T>::CalculateTableValue(uint8_t value) const noexcept {
  constexpr static uint8_t bits = sizeof(T) * 8;
  constexpr static T high_bit = static_cast<T>(1) << (bits - 1);
  T result = static_cast<T>(value);
  if (reverse_data_) {
    result = static_cast<T>(ReverseBits(value, 8));
  }
  result <<= (bits - 8);

  for (size_t i = 0; i < 8; ++i) {
    if (result & high_bit) {
      result <<= 1;
      result ^= polynomial_;
    } else {
      result <<= 1;
    }
  }

  if (reverse_data_) {
    return ReverseBits(result);
  }
  return result;
}

// Reverse bits in the 'value' parameter.
// If bits == sizeof(T) - reverse all bits
// 1000000011000011
// 1100001100000001 - all bits
// If bits < sizeof(T) only Reflect in range of specified bit size.
// 0000000011111001
// 0000000010011111 - only first 8 bits
template <typename T>
T Crc<T>::ReverseBits(T value, size_t bits) const noexcept {
  T reverse = 0;
  for (size_t i = 0; i < bits; ++i) {
    if (value & 1) {
      reverse |= (static_cast<T>(1) << ((bits - 1) - i));
    }
    value >>= 1;
  }
  return reverse;
}

// Swap endianess of a given type.
// Use builtin funcftions if possible.
template <typename T, std::enable_if_t<sizeof(T) == sizeof(uint64_t), int> = 0>
T SwapEndianess(T val) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_bswap64(val);
#elif defined(_MSC_VER)
  return _byteswap_uint64(val);
#else
  return val;
#endif
}

template <typename T, std::enable_if_t<sizeof(T) == sizeof(uint32_t), int> = 0>
T SwapEndianess(T val) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_bswap32(val);
#elif defined(_MSC_VER)
  return _byteswap_ulong(val);
#else
  return val;
#endif
}

template <typename T, std::enable_if_t<sizeof(T) == sizeof(uint16_t), int> = 0>
T SwapEndianess(T val) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_bswap16(val);
#elif defined(_MSC_VER)
  return _byteswap_ushort(val);
#else
  return val;
#endif
}

template <typename T, std::enable_if_t<sizeof(T) == sizeof(uint8_t), int> = 0>
T SwapEndianess(T val) {
  return val;
}

namespace {

class Timer {
 public:
  Timer() : start_(std::chrono::system_clock::now()) {}
  double elapsed() {
    return std::chrono::duration<double>(std::chrono::system_clock::now() -
                                         start_)
        .count();
  }

 private:
  std::chrono::system_clock::time_point start_;
};

}  // namespace

template <typename T>
void Crc<T>::Optimize(uint64_t buffer_size, uint64_t repeats) {
  double best_time = 1e16;
  CrcChunks best_type = CHUNKS_4x32b;
  // Delete later.
  char* buffer = static_cast<char*>(std::malloc(buffer_size * sizeof(*buffer)));
  Timer t;
  for (size_t i = 0; i < repeats; ++i) {
    Consume_1x32b(buffer, buffer_size);
  }
  double time = t.elapsed();
  if (time < best_time) {
    best_time = time;
    best_type = CHUNKS_1x32b;
  }
  t = Timer();
  for (size_t i = 0; i < repeats; ++i) {
    Consume_2x32b(buffer, buffer_size);
  }
  time = t.elapsed();
  if (time < best_time) {
    best_time = time;
    best_type = CHUNKS_2x32b;
  }
  t = Timer();
  for (size_t i = 0; i < repeats; ++i) {
    Consume_4x32b(buffer, buffer_size);
  }
  time = t.elapsed();
  if (time < best_time) {
    best_time = time;
    best_type = CHUNKS_4x32b;
  }
  t = Timer();
  for (size_t i = 0; i < repeats; ++i) {
    Consume_8x32b(buffer, buffer_size);
  }
  time = t.elapsed();
  if (time < best_time) {
    best_time = time;
    best_type = CHUNKS_8x32b;
  }
  t = Timer();
  for (size_t i = 0; i < repeats; ++i) {
    Consume_byte_by_byte(buffer, buffer_size);
  }
  time = t.elapsed();
  if (time < best_time) {
    best_type = BYTE_BY_BYTE;
  }
  reset();
  chunks_ = best_type;
  delete buffer;
}

template <typename T>
template <typename Y>
void Crc<T>::Consume(const Y* data, size_t size) {
  Consume(data, size, chunks_);
}

template <typename T>
template <typename Y>
void Crc<T>::Consume(const Y* data, size_t size, CrcChunks chunks_) {
  switch (chunks_) {
    case CHUNKS_8x32b:
      return Consume_8x32b(data, size);
    case CHUNKS_4x32b:
      return Consume_4x32b(data, size);
    case CHUNKS_2x32b:
      return Consume_2x32b(data, size);
    case CHUNKS_1x32b:
      return Consume_1x32b(data, size);
    case BYTE_BY_BYTE:
      return Consume_byte_by_byte(data, size);
  }
}

template <typename T>
template <typename Y>
void Crc<T>::Consume_byte_by_byte(const Y* data, size_t size) {
  // Cast provided data to match template type.
  const auto* casted_data_8 = reinterpret_cast<const uint8_t*>(data);
  const size_t bytes_left = size * sizeof(Y);
  if (reverse_data_) {
    crc_ = std::accumulate(casted_data_8, casted_data_8 + bytes_left, crc_,
                           [&](T crc, uint8_t value) {
                             return (crc >> 8) ^
                                    lookup_table_[0][(crc ^ value) & 0xFF];
                           });
  } else {
    static constexpr uint32_t shift = (sizeof(T) * 8) - 8;
    crc_ = std::accumulate(
        casted_data_8, casted_data_8 + bytes_left, crc_,
        [&](T crc, uint8_t value) {
          return (crc << 8) ^ lookup_table_[0][((crc >> shift) ^ value) & 0xFF];
        });
  }
}

// Notes:
// To obtain maximum performance code duplication was required. Adding any kind
// of generalized function to handle specific number of 32b words result in
// significant performance decrease.
// 1. Adding 'crc_ ^= ' after every 2 words provide significant increase in
// performance with MVSC compiler.
// 2. Switching from data++ to ++data porvided slight increase in performance.
// 2.1 Require uglu --data before processing.
// 3. Unrolling for. With optimizations turned on provide much faster code.
template <typename T>
template <typename Y>
void Crc<T>::Consume_1x32b(const Y* data, size_t size) {
  // Cast provided data to match template type.
  const auto* casted_data_32 = reinterpret_cast<const uint32_t*>(data);
  size_t bytes_left = size * sizeof(Y);
  static constexpr uint8_t unroll = 16;
  static constexpr uint32_t bytes_at_once = sizeof(uint32_t) * unroll;
  --casted_data_32;
  if (reverse_data_) {
    while (bytes_left >= bytes_at_once) {
      for (size_t i = 0; i < unroll; ++i) {
        const uint32_t word_1 = *++casted_data_32 ^ static_cast<uint32_t>(crc_);
        crc_ = lookup_table_[0][(word_1 >> 24) & 0xFF] ^
               lookup_table_[1][(word_1 >> 16) & 0xFF] ^
               lookup_table_[2][(word_1 >> 8) & 0xFF] ^
               lookup_table_[3][word_1 & 0xFF];
      }
      bytes_left -= bytes_at_once;
    }
  } else {
    while (bytes_left >= bytes_at_once) {
      for (size_t i = 0; i < unroll; ++i) {
        const uint32_t word_1 =
            *++casted_data_32 ^ static_cast<uint32_t>(SwapEndianess(crc_));
        crc_ = lookup_table_[0][(word_1 >> 24) & 0xFF] ^
               lookup_table_[1][(word_1 >> 16) & 0xFF] ^
               lookup_table_[2][(word_1 >> 8) & 0xFF] ^
               lookup_table_[3][word_1 & 0xFF];
      }
      bytes_left -= bytes_at_once;
    }
  }
  // Consume last bytes if any.
  const auto* casted_data_8 =
      reinterpret_cast<const uint8_t*>(++casted_data_32);
  Consume_byte_by_byte(casted_data_8, bytes_left);
}

template <typename T>
template <typename Y>
void Crc<T>::Consume_2x32b(const Y* data, size_t size) {
  // Cast provided data to match template type.
  const auto* casted_data_32 = reinterpret_cast<const uint32_t*>(data);
  size_t bytes_left = size * sizeof(Y);
  static constexpr uint8_t unroll = 8;
  static constexpr uint32_t bytes_at_once = sizeof(uint32_t) * 2 * unroll;
  --casted_data_32;
  if (reverse_data_) {
    while (bytes_left >= bytes_at_once) {
      for (size_t i = 0; i < unroll; ++i) {
        const uint32_t word_1 = *++casted_data_32 ^ static_cast<uint32_t>(crc_);
        const uint32_t word_2 =
            *++casted_data_32 ^
            static_cast<uint32_t>(static_cast<uint64_t>(crc_) >> 32);

        crc_ = lookup_table_[0][(word_2 >> 24) & 0xFF] ^
               lookup_table_[1][(word_2 >> 16) & 0xFF] ^
               lookup_table_[2][(word_2 >> 8) & 0xFF] ^
               lookup_table_[3][word_2 & 0xFF] ^
               lookup_table_[4][(word_1 >> 24) & 0xFF] ^
               lookup_table_[5][(word_1 >> 16) & 0xFF] ^
               lookup_table_[6][(word_1 >> 8) & 0xFF] ^
               lookup_table_[7][word_1 & 0xFF];
      }
      bytes_left -= bytes_at_once;
    }
  } else {
    while (bytes_left >= bytes_at_once) {
      for (size_t i = 0; i < unroll; ++i) {
        const T swapped = SwapEndianess(crc_);
        const uint32_t word_1 =
            *++casted_data_32 ^ static_cast<uint32_t>(swapped);
        const uint32_t word_2 =
            *++casted_data_32 ^
            static_cast<uint32_t>(static_cast<uint64_t>(swapped) >> 32);

        crc_ = lookup_table_[0][(word_2 >> 24) & 0xFF] ^
               lookup_table_[1][(word_2 >> 16) & 0xFF] ^
               lookup_table_[2][(word_2 >> 8) & 0xFF] ^
               lookup_table_[3][word_2 & 0xFF] ^
               lookup_table_[4][(word_1 >> 24) & 0xFF] ^
               lookup_table_[5][(word_1 >> 16) & 0xFF] ^
               lookup_table_[6][(word_1 >> 8) & 0xFF] ^
               lookup_table_[7][word_1 & 0xFF];
      }
      bytes_left -= bytes_at_once;
    }
  }
  // Consume last bytes if any.
  const auto* casted_data_8 =
      reinterpret_cast<const uint8_t*>(++casted_data_32);
  Consume_byte_by_byte(casted_data_8, bytes_left);
}

template <typename T>
template <typename Y>
void Crc<T>::Consume_4x32b(const Y* data, size_t size) {
  // Cast provided data to match template type.
  const auto* casted_data_32 = reinterpret_cast<const uint32_t*>(data);
  size_t bytes_left = size * sizeof(Y);
  static constexpr uint8_t unroll = 4;
  static constexpr uint32_t bytes_at_once = sizeof(uint32_t) * 4 * unroll;
  --casted_data_32;
  if (reverse_data_) {
    while (bytes_left >= bytes_at_once) {
      for (size_t i = 0; i < unroll; ++i) {
        const uint32_t word_1 = *++casted_data_32 ^ static_cast<uint32_t>(crc_);
        const uint32_t word_2 =
            *++casted_data_32 ^
            static_cast<uint32_t>(static_cast<uint64_t>(crc_) >> 32);
        const uint32_t word_3 = *++casted_data_32;
        const uint32_t word_4 = *++casted_data_32;

        crc_ = lookup_table_[0][(word_4 >> 24) & 0xFF] ^
               lookup_table_[1][(word_4 >> 16) & 0xFF] ^
               lookup_table_[2][(word_4 >> 8) & 0xFF] ^
               lookup_table_[3][word_4 & 0xFF] ^
               lookup_table_[4][(word_3 >> 24) & 0xFF] ^
               lookup_table_[5][(word_3 >> 16) & 0xFF] ^
               lookup_table_[6][(word_3 >> 8) & 0xFF] ^
               lookup_table_[7][word_3 & 0xFF];
        crc_ ^= lookup_table_[8][(word_2 >> 24) & 0xFF] ^
                lookup_table_[9][(word_2 >> 16) & 0xFF] ^
                lookup_table_[10][(word_2 >> 8) & 0xFF] ^
                lookup_table_[11][word_2 & 0xFF] ^
                lookup_table_[12][(word_1 >> 24) & 0xFF] ^
                lookup_table_[13][(word_1 >> 16) & 0xFF] ^
                lookup_table_[14][(word_1 >> 8) & 0xFF] ^
                lookup_table_[15][word_1 & 0xFF];
      }
      bytes_left -= bytes_at_once;
    }
  } else {
    while (bytes_left >= bytes_at_once) {
      for (size_t i = 0; i < unroll; ++i) {
        const T swapped = SwapEndianess(crc_);
        const uint32_t word_1 =
            *++casted_data_32 ^ static_cast<uint32_t>(swapped);
        const uint32_t word_2 =
            *++casted_data_32 ^
            static_cast<uint32_t>(static_cast<uint64_t>(swapped) >> 32);
        const uint32_t word_3 = *++casted_data_32;
        const uint32_t word_4 = *++casted_data_32;

        crc_ = lookup_table_[0][(word_4 >> 24) & 0xFF] ^
               lookup_table_[1][(word_4 >> 16) & 0xFF] ^
               lookup_table_[2][(word_4 >> 8) & 0xFF] ^
               lookup_table_[3][word_4 & 0xFF] ^
               lookup_table_[4][(word_3 >> 24) & 0xFF] ^
               lookup_table_[5][(word_3 >> 16) & 0xFF] ^
               lookup_table_[6][(word_3 >> 8) & 0xFF] ^
               lookup_table_[7][word_3 & 0xFF];
        crc_ ^= lookup_table_[8][(word_2 >> 24) & 0xFF] ^
                lookup_table_[9][(word_2 >> 16) & 0xFF] ^
                lookup_table_[10][(word_2 >> 8) & 0xFF] ^
                lookup_table_[11][word_2 & 0xFF] ^
                lookup_table_[12][(word_1 >> 24) & 0xFF] ^
                lookup_table_[13][(word_1 >> 16) & 0xFF] ^
                lookup_table_[14][(word_1 >> 8) & 0xFF] ^
                lookup_table_[15][word_1 & 0xFF];
      }
      bytes_left -= bytes_at_once;
    }
  }
  // Consume last bytes if any.
  const auto* casted_data_8 =
      reinterpret_cast<const uint8_t*>(++casted_data_32);
  Consume_byte_by_byte(casted_data_8, bytes_left);
}

template <typename T>
template <typename Y>
void Crc<T>::Consume_8x32b(const Y* data, size_t size) {
  // Cast provided data to match template type.
  const auto* casted_data_32 = reinterpret_cast<const uint32_t*>(data);
  size_t bytes_left = size * sizeof(Y);
  static constexpr uint8_t unroll = 2;
  static constexpr uint32_t bytes_at_once = sizeof(uint32_t) * 8 * unroll;
  --casted_data_32;
  if (reverse_data_) {
    while (bytes_left >= bytes_at_once) {
      for (size_t i = 0; i < unroll; ++i) {
        const uint32_t word_1 = *++casted_data_32 ^ static_cast<uint32_t>(crc_);
        const uint32_t word_2 =
            *++casted_data_32 ^
            static_cast<uint32_t>(static_cast<uint64_t>(crc_) >> 32);
        const uint32_t word_3 = *++casted_data_32;
        const uint32_t word_4 = *++casted_data_32;
        const uint32_t word_5 = *++casted_data_32;
        const uint32_t word_6 = *++casted_data_32;
        const uint32_t word_7 = *++casted_data_32;
        const uint32_t word_8 = *++casted_data_32;
        crc_ = lookup_table_[0][(word_8 >> 24) & 0xFF] ^
               lookup_table_[1][(word_8 >> 16) & 0xFF] ^
               lookup_table_[2][(word_8 >> 8) & 0xFF] ^
               lookup_table_[3][(word_8)&0xFF] ^
               lookup_table_[4][(word_7 >> 24) & 0xFF] ^
               lookup_table_[5][(word_7 >> 16) & 0xFF] ^
               lookup_table_[6][(word_7 >> 8) & 0xFF] ^
               lookup_table_[7][word_7 & 0xFF];
        crc_ ^= lookup_table_[8][(word_6 >> 24) & 0xFF] ^
                lookup_table_[9][(word_6 >> 16) & 0xFF] ^
                lookup_table_[10][(word_6 >> 8) & 0xFF] ^
                lookup_table_[11][(word_6)&0xFF] ^
                lookup_table_[12][(word_5 >> 24) & 0xFF] ^
                lookup_table_[13][(word_5 >> 16) & 0xFF] ^
                lookup_table_[14][(word_5 >> 8) & 0xFF] ^
                lookup_table_[15][word_5 & 0xFF];
        crc_ ^= lookup_table_[16][(word_4 >> 24) & 0xFF] ^
                lookup_table_[17][(word_4 >> 16) & 0xFF] ^
                lookup_table_[18][(word_4 >> 8) & 0xFF] ^
                lookup_table_[19][(word_4)&0xFF] ^
                lookup_table_[20][(word_3 >> 24) & 0xFF] ^
                lookup_table_[21][(word_3 >> 16) & 0xFF] ^
                lookup_table_[22][(word_3 >> 8) & 0xFF] ^
                lookup_table_[23][word_3 & 0xFF];
        crc_ ^= lookup_table_[24][(word_2 >> 24) & 0xFF] ^
                lookup_table_[25][(word_2 >> 16) & 0xFF] ^
                lookup_table_[26][(word_2 >> 8) & 0xFF] ^
                lookup_table_[27][(word_2)&0xFF] ^
                lookup_table_[28][(word_1 >> 24) & 0xFF] ^
                lookup_table_[29][(word_1 >> 16) & 0xFF] ^
                lookup_table_[30][(word_1 >> 8) & 0xFF] ^
                lookup_table_[31][word_1 & 0xFF];
      }
      bytes_left -= bytes_at_once;
    }
  } else {
    while (bytes_left >= bytes_at_once) {
      for (size_t i = 0; i < unroll; ++i) {
        const T swapped = SwapEndianess(crc_);
        const uint32_t word_1 =
            *++casted_data_32 ^ static_cast<uint32_t>(swapped);
        const uint32_t word_2 =
            *++casted_data_32 ^
            static_cast<uint32_t>(static_cast<uint64_t>(swapped) >> 32);
        const uint32_t word_3 = *++casted_data_32;
        const uint32_t word_4 = *++casted_data_32;
        const uint32_t word_5 = *++casted_data_32;
        const uint32_t word_6 = *++casted_data_32;
        const uint32_t word_7 = *++casted_data_32;
        const uint32_t word_8 = *++casted_data_32;

        crc_ = lookup_table_[0][(word_8 >> 24) & 0xFF] ^
               lookup_table_[1][(word_8 >> 16) & 0xFF] ^
               lookup_table_[2][(word_8 >> 8) & 0xFF] ^
               lookup_table_[3][(word_8)&0xFF] ^
               lookup_table_[4][(word_7 >> 24) & 0xFF] ^
               lookup_table_[5][(word_7 >> 16) & 0xFF] ^
               lookup_table_[6][(word_7 >> 8) & 0xFF] ^
               lookup_table_[7][word_7 & 0xFF];
        crc_ ^= lookup_table_[8][(word_6 >> 24) & 0xFF] ^
                lookup_table_[9][(word_6 >> 16) & 0xFF] ^
                lookup_table_[10][(word_6 >> 8) & 0xFF] ^
                lookup_table_[11][(word_6)&0xFF] ^
                lookup_table_[12][(word_5 >> 24) & 0xFF] ^
                lookup_table_[13][(word_5 >> 16) & 0xFF] ^
                lookup_table_[14][(word_5 >> 8) & 0xFF] ^
                lookup_table_[15][word_5 & 0xFF];
        crc_ ^= lookup_table_[16][(word_4 >> 24) & 0xFF] ^
                lookup_table_[17][(word_4 >> 16) & 0xFF] ^
                lookup_table_[18][(word_4 >> 8) & 0xFF] ^
                lookup_table_[19][(word_4)&0xFF] ^
                lookup_table_[20][(word_3 >> 24) & 0xFF] ^
                lookup_table_[21][(word_3 >> 16) & 0xFF] ^
                lookup_table_[22][(word_3 >> 8) & 0xFF] ^
                lookup_table_[23][word_3 & 0xFF];
        crc_ ^= lookup_table_[24][(word_2 >> 24) & 0xFF] ^
                lookup_table_[25][(word_2 >> 16) & 0xFF] ^
                lookup_table_[26][(word_2 >> 8) & 0xFF] ^
                lookup_table_[27][(word_2)&0xFF] ^
                lookup_table_[28][(word_1 >> 24) & 0xFF] ^
                lookup_table_[29][(word_1 >> 16) & 0xFF] ^
                lookup_table_[30][(word_1 >> 8) & 0xFF] ^
                lookup_table_[31][word_1 & 0xFF];
      }
      bytes_left -= bytes_at_once;
    }
  }
  // Consume last bytes if any.
  const auto* casted_data_8 =
      reinterpret_cast<const uint8_t*>(++casted_data_32);
  Consume_byte_by_byte(casted_data_8, bytes_left);
}

template <typename T>
T Crc<T>::crc() const noexcept {
  return (reverse_out_ ^ reverse_data_) ? ReverseBits(crc_) ^ xor_output_
                                        : crc_ ^ xor_output_;
}

template <typename T>
void Crc<T>::reset() noexcept {
  crc_ = (reverse_data_) ? ReverseBits(initial_crc_) : initial_crc_;
}

// Some of the most popoular CRC options.
OptionsCrc constexpr OptionsCrc::Crc16() {
  return OptionsCrc(static_cast<uint64_t>(0x8005),
                    static_cast<uint64_t>(0x0000),
                    static_cast<uint64_t>(0x0000), true, true);
}

OptionsCrc constexpr OptionsCrc::Crc16_CCITT() {
  return OptionsCrc(static_cast<uint64_t>(0x1021),
                    static_cast<uint64_t>(0xFFFF),
                    static_cast<uint64_t>(0x0000), false, false);
}

OptionsCrc constexpr OptionsCrc::Crc32() {
  return OptionsCrc(static_cast<uint64_t>(0x4C11DB7),
                    static_cast<uint64_t>(0xFFFFFFFF),
                    static_cast<uint64_t>(0xFFFFFFFF), true, true);
}

OptionsCrc constexpr OptionsCrc::Crc64() {
  return OptionsCrc(static_cast<uint64_t>(0x42F0E1EBA9EA3693),
                    static_cast<uint64_t>(0xFFFFFFFFFFFFFFFF),
                    static_cast<uint64_t>(0xFFFFFFFFFFFFFFFF), true, true);
}

OptionsCrc constexpr OptionsCrc::Crc64_ISO() {
  return OptionsCrc(static_cast<uint64_t>(0x000000000000001B),
                    static_cast<uint64_t>(0x0000000000000000),
                    static_cast<uint64_t>(0x0000000000000000), true, true);
}

// Create CRC class with CRC16 parameters.
static Crc<uint16_t> NewCrc16(const OptionsCrc& options = OptionsCrc::Crc16()) {
  return Crc<uint16_t>(options);
}

// Create CRC class with CRC32 parameters.
static Crc<uint32_t> NewCrc32(const OptionsCrc& options = OptionsCrc::Crc32()) {
  return Crc<uint32_t>(options);
}

// Create CRC class with CRC64 parameters.
static Crc<uint64_t> NewCrc64(const OptionsCrc& options = OptionsCrc::Crc64()) {
  return Crc<uint64_t>(options);
}

}  // namespace hash

#endif  // CRC_H_
