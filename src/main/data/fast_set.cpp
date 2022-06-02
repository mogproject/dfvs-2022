#include "fast_set.hpp"

namespace mog {
namespace data {

/**
 * @brief Fixed-size, constant-time integer set representation.
 */
fast_set::fast_set(std::size_t size) : generation_(0) { data_.resize(size, -1); }

/**
 * @brief Clears the set
 *
 * O(1)
 */
void fast_set::clear() {
  ++generation_;
  if (generation_ == 0) {
    auto sz = data_.size();
    data_.clear();
    data_.resize(sz, -1);
  }
}

/**
 * @brief Resizes the capacity of the set.
 *
 * O(size)
 */
void fast_set::resize(std::size_t size) {
  generation_ = 0;
  data_.clear();
  data_.resize(size, -1);
}

/**
 * @brief Inserts one element to the set.
 *
 * O(1)
 */
void fast_set::set(int x) {
  if (x < 0 || static_cast<int>(data_.size()) <= x) throw std::invalid_argument("fast_set::set(): out of range");
  data_[x] = generation_;
}

/**
 * @brief Removes one element from the set.
 *
 * O(1)
 */
void fast_set::reset(int x) {
  if (x < 0 || static_cast<int>(data_.size()) <= x) throw std::invalid_argument("fast_set::reset(): out of range");
  data_[x] = -1;
}

/**
 * @brief Checks if the given element is in the set.
 *
 * O(1)
 */
bool fast_set::get(int x) const { return 0 <= x && x < static_cast<int>(data_.size()) && data_[x] == generation_; }

}  // namespace data
}  // namespace mog