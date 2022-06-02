#include "bitmap.hpp"

namespace mog {
namespace data {

std::ostream& operator<<(std::ostream& os, Bitmap const& vm) {
  os << "Bitmap({";
  int i = 0;
  for (auto x : vm.to_vector()) {
    if (i++ > 0) os << ",";
    os << x;
  }
  return os << "})";
}

}  // namespace data
}  // namespace mog