#pragma once
#undef _GLIBCXX_DEBUG                 // disable run-time bound checking, etc
#pragma GCC optimize("Ofast,inline")  // Ofast = O3,fast-math,allow-store-data-races,no-protect-parens

#pragma GCC target("bmi,bmi2,lzcnt,popcnt")                       // bit manipulation
#pragma GCC target("movbe")                                       // byte swap
#pragma GCC target("aes,pclmul,rdrnd")                            // encryption
#pragma GCC target("avx,avx2,f16c,fma,sse3,ssse3,sse4.1,sse4.2")  // SIMD

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "../util/util.hpp"

namespace mog {
namespace data {

typedef unsigned long long L;

/**
 * Bitmap of variable length.
 */
class Bitmap {
 private:
  static int const B = sizeof(L) * 8;
  std::size_t n_;
  std::vector<L> data_;

  inline void verify_argument_(int x, char const* label) const {
    if (0 <= x && x < static_cast<int>(size())) return;
    throw std::invalid_argument(util::format("`%s`: x (%d) must be between 0 and %d", label, x, size() - 1));
  }

 public:
  Bitmap(std::size_t n = 0) : n_(n) { data_.resize((n + (B - 1)) / B); }
  /**
   * Constructs a singleton.
   */
  Bitmap(std::size_t n, int x) : Bitmap(n) { *this |= x; }
  /**
   * Constructs from a list.
   */
  Bitmap(std::size_t n, std::vector<int> const& xs) : Bitmap(n) {
    for (int x : xs) *this |= x;
  }
  /**
   * Copy constructor.
   */
  Bitmap(Bitmap const& other) {
    this->n_ = other.n_;
    this->data_ = other.data_;
  }

  void operator=(Bitmap const& other) {
    n_ = other.n_;
    data_ = other.data_;
  }

  inline void clear() {
    for (std::size_t i = 0; i < data_.size(); ++i) data_[i] = 0;
  }

  inline bool empty() const {
    for (auto d : data_) {
      if (d) return false;
    }
    return true;
  }

  inline int front() const {
    int ret = -1;
    for (std::size_t i = 0; i < data_.size(); ++i) {
      if (data_[i]) {
        ret = i * B + __builtin_ctzll(data_[i]);
        break;
      }
    }
    return ret < static_cast<int>(size()) ? ret : -1;
  }

  inline int pop_front() {
    for (std::size_t i = 0; i < data_.size(); ++i)
      if (data_[i]) {
        int offset = __builtin_ctzll(data_[i]);
        int ret = i * B + offset;
        if (ret < static_cast<int>(size())) {
          data_[i] ^= 1ULL << offset;
          return ret;
        }
      }
    return -1;
  }

  //--------------------------------------------------------
  //    Operators
  //--------------------------------------------------------

  /**
   * negation (not)
   */
  Bitmap operator~() const {
    Bitmap ret(size());
    for (std::size_t i = 0; i < data_.size(); ++i) ret.data_[i] = ~data_[i];
    // trim extra bits
    if (size() % B) ret.data_.back() &= (1ULL << (size() % B)) - 1;
    return ret;
  }

  /**
   * set/union (or)
   */
  Bitmap& operator|=(int x) {
    verify_argument_(x, "|=");

    data_[x / B] |= 1ULL << (x % B);
    return *this;
  }

  Bitmap& operator|=(Bitmap const& rhs) {
    if (size() != rhs.size()) throw std::invalid_argument("inconsistent size");

    for (std::size_t i = 0; i < data_.size(); ++i) data_[i] |= rhs.data_[i];
    return *this;
  }

  /**
   * exclusive or (xor)
   */
  Bitmap& operator^=(int x) {
    verify_argument_(x, "^=");

    data_[x / B] ^= 1ULL << (x % B);
    return *this;
  }

  Bitmap& operator^=(Bitmap const& rhs) {
    if (size() != rhs.size()) throw std::invalid_argument("inconsistent size");

    for (std::size_t i = 0; i < data_.size(); ++i) data_[i] ^= rhs.data_[i];
    return *this;
  }

  /**
   * intersection (and)
   */
  Bitmap operator&=(int x) {
    verify_argument_(x, "&=");

    for (int i = 0; i < static_cast<int>(data_.size()); ++i) {
      if (i == x / B) {
        data_[x / B] &= 1ULL << (x % B);
      } else {
        data_[i] = 0;
      }
    }

    return *this;
  }

  Bitmap& operator&=(Bitmap const& rhs) {
    if (size() != rhs.size()) throw std::invalid_argument("inconsistent size");

    for (std::size_t i = 0; i < data_.size(); ++i) data_[i] &= rhs.data_[i];
    return *this;
  }

  /**
   * reset/set minus
   */
  Bitmap& operator-=(int x) {
    verify_argument_(x, "-=");

    data_[x / B] &= ~(1ULL << (x % B));
    return *this;
  }

  Bitmap& operator-=(Bitmap const& rhs) {
    if (size() != rhs.size()) throw std::invalid_argument("inconsistent size");

    for (std::size_t i = 0; i < data_.size(); ++i) data_[i] &= ~rhs.data_[i];
    return *this;
  }

  // clang-format off
  friend Bitmap operator|(Bitmap const& lhs, int x) { Bitmap ret(lhs); ret |= x; return ret; }
  friend Bitmap operator|(Bitmap const& lhs, Bitmap const& rhs) { Bitmap ret(lhs); ret |= rhs; return ret; }
  friend Bitmap operator^(Bitmap const& lhs, int x) { Bitmap ret(lhs); ret ^= x; return ret; }
  friend Bitmap operator^(Bitmap const& lhs, Bitmap const& rhs) { Bitmap ret(lhs); ret ^= rhs; return ret; }
  friend Bitmap operator&(Bitmap const& lhs, int x) { Bitmap ret(lhs); ret &= x; return ret; }
  friend Bitmap operator&(Bitmap const& lhs, Bitmap const& rhs) { Bitmap ret(lhs); ret &= rhs; return ret; }
  friend Bitmap operator-(Bitmap const& lhs, int x) { Bitmap ret(lhs); ret -= x; return ret; }
  friend Bitmap operator-(Bitmap const& lhs, Bitmap const& rhs) { Bitmap ret(lhs); ret -= rhs; return ret; }
  // clang-format on

  /**
   * get
   */
  bool operator[](int x) const {
    verify_argument_(x, "[]");
    return (data_[x / B] >> (x % B)) & 1ULL;
  }

  std::vector<int> to_vector() const {
    std::vector<int> ret;
    for (std::size_t i = 0; i < data_.size(); ++i) {
      if (data_[i]) {
        auto x = data_[i];
        while (x) {
          ret.push_back(i * B + __builtin_ctzll(x));
          x &= x - 1;
        }
      }
    }
    return ret;
  }

  /**
   * equality
   */
  friend inline bool operator==(Bitmap const& lhs, Bitmap const& rhs) {
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.data_.size(); ++i) {
      if (lhs.data_[i] != rhs.data_[i]) return false;
    }
    return true;
  }

  friend inline bool operator!=(Bitmap const& lhs, Bitmap const& rhs) { return !(lhs == rhs); }

  inline size_t size() const { return n_; }

  /**
   * @brief Returns the number of 1-bits.
   */
  std::size_t count() const {
    std::size_t ret = 0;
    for (std::size_t i = 0; i < data_.size(); ++i) ret += __builtin_popcountll(data_[i]);
    return ret;
  }

  inline bool subset(Bitmap const& rhs) const { return (*this & rhs) == *this; }

  inline bool superset(Bitmap const& rhs) const { return rhs.subset(*this); }

  std::string to_string() const {
    std::stringstream ss;
    for (int i = data_.size() - 1; i >= 0; --i) { ss << std::setfill('0') << std::setw(B / 4) << std::hex << data_[i]; }
    return ss.str();
  }

  void resize(std::size_t new_size) {
    if (n_ > new_size) throw std::invalid_argument("cannot shrink the data");
    if (n_ == new_size) return;  // do nothing
    n_ = new_size;
    data_.resize((n_ + (B - 1)) / B);
  }
};

std::ostream& operator<<(std::ostream& os, Bitmap const& vm);

}  // namespace data
}  // namespace mog