#pragma once

#include <stdexcept>
#include <vector>

namespace mog {
namespace data {
/**
 * @brief Fixed-size, constant-time integer set representation.
 */
class fast_set {
 public:
  fast_set(std::size_t size = 0);

  /**
   * @brief Clears the set
   *
   * O(1)
   */
  void clear();

  /**
   * @brief Resizes the capacity of the set.
   *
   * O(size)
   */
  void resize(std::size_t size);

  /**
   * @brief Inserts one element to the set.
   *
   * O(1)
   */
  void set(int x);

  /**
   * @brief Removes one element from the set.
   *
   * O(1)
   */
  void reset(int x);

  /**
   * @brief Checks if the given element is in the set.
   *
   * O(1)
   */
  bool get(int x) const;

 private:
  int generation_;
  std::vector<int> data_;
};
}  // namespace data
}  // namespace mog