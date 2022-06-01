#include "util.hpp"
#include <cstdarg>

namespace mog {
namespace util {

std::vector<int> range(int n) {
  std::vector<int> ret;
  for (int i = 0; i < n; ++i) ret.push_back(i);
  return ret;
}

//==================================================================================================
// String Utilities
//==================================================================================================
/**
 * @brief Safer alternative to sprintf.
 *
 * @param fmt
 * @param ...
 * @return char const*
 */
char const* format(char const* fmt, ...) {
  static std::string buffer;  // should have continuous memory space
  std::size_t buff_size = 256;
  va_list marker;

  while (true) {
    buffer.resize(buff_size);
    va_start(marker, fmt);
    int n = vsnprintf((char*)buffer.c_str(), buff_size, fmt, marker);
    va_end(marker);
    if (n >= 0) break;
    buff_size *= 2;  // double the buffer size
  }
  return buffer.c_str();
}
}  // namespace util
}  // namespace mog
