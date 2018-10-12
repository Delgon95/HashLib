#include <chrono>
#include <cstdlib>
#include <fstream>

#include "crc.h"

class Timer {
 public:
  Timer() : start_(std::chrono::high_resolution_clock::now()) {}
  double elapsed() {
    return std::chrono::duration<double>(
               std::chrono::high_resolution_clock::now() - start_)
        .count();
  }

 private:
  std::chrono::high_resolution_clock::time_point start_;
};

int main(int argc, char* argv[]) {
  auto crc16 = hash::NewCrc16();
  crc16.Optimize();
  auto crc16_ccitt = hash::NewCrc16(hash::OptionsCrc::Crc16_CCITT());
  crc16_ccitt.Optimize();
  auto crc32 = hash::NewCrc32();
  crc32.Optimize();
  auto crc64 = hash::NewCrc64();
  crc64.Optimize();
  auto crc64_iso = hash::NewCrc64(hash::OptionsCrc::Crc64_ISO());
  crc64_iso.Optimize();

  auto crc32_1 = hash::NewCrc32();
  auto crc32_2 = hash::NewCrc32();
  auto crc32_4 = hash::NewCrc32();
  auto crc32_8 = hash::NewCrc32();
  auto crc32_bbb = hash::NewCrc32();

  if (argc == 2) {
    std::basic_ifstream<char> file(argv[1],
                                   std::ios_base::in | std::ios_base::binary);
    const std::size_t buf_size = 1024 * 1024;
    char* buffer = static_cast<char*>(std::malloc(buf_size * sizeof(*buffer)));
    file.seekg(0, std::ios_base::end);
    const std::size_t size = file.tellg();
    file.seekg(0, std::ios_base::beg);

    double read_time = 0;
    double crc16_time = 0;
    double crc16_ccitt_time = 0;
    double crc32_time = 0;
    double crc64_time = 0;
    double crc64_iso_time = 0;

    double crc32_1_time = 0;
    double crc32_2_time = 0;
    double crc32_4_time = 0;
    double crc32_8_time = 0;
    double crc32_bbb_time = 0;

    Timer t;
    for (std::size_t i = 0; i < size; i += buf_size) {
      t = Timer();
      file.read(buffer, buf_size);
      read_time += t.elapsed();

      t = Timer();
      crc32_1.Consume(buffer, file.gcount(), hash::CHUNKS_1x32b);
      crc32_1_time += t.elapsed();
      t = Timer();
      crc32_2.Consume(buffer, file.gcount(), hash::CHUNKS_2x32b);
      crc32_2_time += t.elapsed();
      t = Timer();
      crc32_4.Consume(buffer, file.gcount(), hash::CHUNKS_4x32b);
      crc32_4_time += t.elapsed();
      t = Timer();
      crc32_8.Consume(buffer, file.gcount(), hash::CHUNKS_8x32b);
      crc32_8_time += t.elapsed();
      t = Timer();
      crc32_bbb.Consume(buffer, file.gcount(), hash::BYTE_BY_BYTE);
      crc32_bbb_time += t.elapsed();
      t = Timer();
      crc16.Consume(buffer, file.gcount());
      crc16_time += t.elapsed();
      t = Timer();
      crc16_ccitt.Consume(buffer, file.gcount());
      crc16_ccitt_time += t.elapsed();
      t = Timer();
      crc32.Consume(buffer, file.gcount());
      crc32_time += t.elapsed();
      t = Timer();
      crc64.Consume(buffer, file.gcount());
      crc64_time += t.elapsed();
      t = Timer();
      crc64_iso.Consume(buffer, file.gcount());
      crc64_iso_time += t.elapsed();
    }

    const double total_mb = static_cast<double>(size) / (1024.0 * 1024.0);
    printf("Read %lu B (%.3lf MiB) in %.3lfs (%lf MiB/s)\n", size, total_mb,
           read_time, static_cast<double>(total_mb) / read_time);

    printf("CRC32_bbb:        in %.6lfs (%.3lf MiB/s)\n", crc32_bbb_time,
           total_mb / crc32_bbb_time);
    printf("CRC32_1x32b:      in %.6lfs (%.3lf MiB/s)\n", crc32_1_time,
           total_mb / crc32_1_time);
    printf("CRC32_2x32b:      in %.6lfs (%.3lf MiB/s)\n", crc32_2_time,
           total_mb / crc32_2_time);
    printf("CRC32_4x32b:      in %.6lfs (%.3lf MiB/s)\n", crc32_4_time,
           total_mb / crc32_4_time);
    printf("CRC32_8x32b:      in %.6lfs (%.3lf MiB/s)\n", crc32_8_time,
           total_mb / crc32_8_time);

    printf("CRC16:            in %.6lfs (%.3lf MiB/s)\n", crc16_time,
           total_mb / crc16_time);
    printf("CRC16-CCITT:      in %.6lfs (%.3lf MiB/s)\n", crc16_ccitt_time,
           total_mb / crc16_ccitt_time);
    printf("CRC32:            in %.6lfs (%.3lf MiB/s)\n", crc32_time,
           total_mb / crc32_time);
    printf("CRC64:            in %.6lfs (%.3lf MiB/s)\n", crc64_time,
           total_mb / crc64_time);
    printf("CRC64-ISO:        in %.6lfs (%.3lf MiB/s)\n\n", crc64_iso_time,
           total_mb / crc64_iso_time);

    file.close();
    delete buffer;
  } else {
    crc16.Consume("1234567890", 10);
    crc16_ccitt.Consume("1234567890", 10);
    crc32.Consume("1234567890", 10);
    crc64.Consume("1234567890", 10);
    crc64_iso.Consume("1234567890", 10);
  }
  printf("CRC16:                            %.4X\n", crc16.crc());
  printf("CRC16-CCITT:                      %.4X\n", crc16_ccitt.crc());
  printf("CRC32:                        %.8X\n", crc32.crc());
  printf("CRC64:                %.16lX\n", crc64.crc());
  printf("CRC64-ISO:            %.16lX\n", crc64_iso.crc());

  return 0;
}
