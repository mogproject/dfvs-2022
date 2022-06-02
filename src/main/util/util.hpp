#pragma once

#include <algorithm>
#include <cstring>
#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <unordered_set>
#include <vector>
#include <random>

#include "type_traits.hpp"

#define ANSI_RED "\x1b[31m"
#define ANSI_GREEN "\x1b[32m"
#define ANSI_YELLOW "\x1b[33m"
#define ANSI_BLUE "\x1b[34m"
#define ANSI_MAGENTA "\x1b[35m"
#define ANSI_CYAN "\x1b[36m"
#define ANSI_RESET "\x1b[0m"

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#if TRACE_ON
// #define TRACE(format, ...) printf("%s(%d) " format "\n", __FUNCTION__, __LINE__, __VA_ARGS__)
#define TRACE(format, ...) \
  printf("%s%s(%d):%s " format "\n", ANSI_BLUE, __FILENAME__, __LINE__, ANSI_RESET, __VA_ARGS__);
#define TRACE0(s) printf("%s%s(%d):%s %s\n", ANSI_BLUE, __FILENAME__, __LINE__, ANSI_RESET, s);
#else
#define TRACE(format, ...)
#define TRACE0(s)
#endif

namespace mog {
namespace util {

//==================================================================================================
// STL Extensions (assuming SFINAE)
//==================================================================================================
/** Implementation specialized for map and set. */
template <typename T, typename U>
typename std::enable_if<is_map<T>::value || is_set<T>::value, bool>::type contains(T const& col, U const& x) {
  return col.find(x) != col.end();
}

/** Otherwise, use std::find to search for the element */
template <typename T, typename U>
typename std::enable_if<!(is_map<T>::value || is_set<T>::value), bool>::type contains(T const& col, U const& x) {
  return std::find(col.begin(), col.end(), x) != col.end();
}

template <typename T, typename U>
inline void erase_first(T& col, U const& x) {
  col.erase(std::find(col.begin(), col.end(), x));
}

template <typename T>
void extend(std::vector<T>& xs, std::vector<T> const& ys) {
  xs.insert(xs.end(), ys.begin(), ys.end());
}

template <typename T>
std::vector<T> flatten(std::vector<std::vector<T>> const& xss) {
  std::vector<T> ret;
  for (auto const& xs : xss) extend(ret, xs);
  return ret;
}

/** Erase an element at the given index. */
template <typename T>
void erase(std::vector<T>& xs, int index) {
  xs.erase(xs.begin() + index);
}

template <typename K, typename V>
V get_or_else(std::map<K, V> const& m, K const& k, V const& default_value) {
  auto it = m.find(k);
  return m.find(k) == m.end() ? default_value : it->second;
}

template <typename T>
int max_size_element_index(std::vector<T> const& xs) {
  if (xs.empty()) throw std::invalid_argument("max_size_element_index(): xs cannot be empty");

  int ret = 0;
  for (int i = 1; i < static_cast<int>(xs.size()); ++i) {
    if (xs[i].size() > xs[ret].size()) ret = i;
  }
  return ret;
}

template <typename T>
T const& max_size_element(std::vector<T> const& xs) {
  return xs[max_size_element_index(xs)];
}

/**
 * @brief Returns a vector x such that x[i]=i for every i<n.
 * 
 * @param n size of the vector
 * @return std::vector<int> 
 */
std::vector<int> range(int n);

/**
 * @brief Returns the inverse map of the given vector, assuming no duplicates.
 *
 * @tparam T element type
 * @param xs vector with no duplicates
 * @return map from item to index
 */
template <typename T>
std::unordered_map<T, int> inverse(std::vector<T> const& xs) {
  std::unordered_map<T, int> ret;
  for (int i = 0; i < static_cast<int>(xs.size()); ++i) ret[xs[i]] = i;
  return ret;
}

template <typename T, typename P>
int count_if(std::vector<T> const& xs, P const& pred) {
  return std::count_if(xs.begin(), xs.end(), pred);
}

template <typename T, typename P>
int all_of(std::vector<T> const& xs, P const& pred) {
  return std::all_of(xs.begin(), xs.end(), pred);
}

template <typename T>
std::vector<T> sorted(std::vector<T> xs) {
  std::sort(xs.begin(), xs.end());
  return xs;
}

/**
 * @brief Randomly sample k elements from the given vector
 * 
 * @param xs given vector
 * @param k number of samples
 * @param gen pseudorandom number generator
 * @return randomly sampled elements
 */
template <typename T>
std::vector<T> sample(std::vector<T> const& xs, std::size_t k, std::default_random_engine &gen) {
  if (xs.size() <= k) return xs;  // k is large enough so we can sample everything
  std::vector<T> ret = xs;
  std::shuffle(ret.begin(), ret.end(), gen);
  ret.resize(k);
  return ret;
}

template <typename ForwardIter>
std::string to_string(ForwardIter begin, ForwardIter end) {
  std::stringstream ss;
  ss << "[";
  for (auto it = begin; it != end; ++it) {
    if (it != begin) ss << ", ";
    ss << (*it);
    // ss << (*it)->to_string();
  }
  ss << "]";
  return ss.str();
}

template <typename T>
std::string to_string(std::list<T> const& xs) {
  return to_string(xs.begin(), xs.end());
}
template <typename T>
std::string to_string(std::vector<T> const& xs) {
  return to_string(xs.begin(), xs.end());
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
char const* format(char const* fmt, ...);

}  // namespace util
}  // namespace mog

//==================================================================================================
// I/O Support
//==================================================================================================
// printers for STL types
// template <typename T>
// std::ostream& operator<<(std::ostream& stream, std::set<T> const& s);
// template <typename T>
// std::ostream& operator<<(std::ostream& stream, std::vector<T> const& s);
template <typename T>
std::ostream& operator<<(std::ostream& stream, std::vector<T> const& v);

template <typename A, typename B>
std::ostream& operator<<(std::ostream& stream, std::pair<A, B> const& p) {
  return stream << "(" << p.first << ", " << p.second << ")";
}
template <typename A, typename B>
std::ostream& operator<<(std::ostream& stream, std::map<A, B> const& p) {
  stream << "{";
  for (auto& x : p) stream << x.first << ":" << x.second << ", ";
  return stream << "}";
}
template <typename T>
std::ostream& operator<<(std::ostream& stream, std::vector<T> const& v) {
  stream << "[";
  for (auto i = v.begin(); i != v.end(); ++i) stream << ((i == v.begin()) ? "" : ", ") << *i;
  return stream << "]";
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, std::list<T> const& v) {
  stream << "[";
  for (auto i = v.begin(); i != v.end(); ++i) stream << ((i == v.begin()) ? "" : ", ") << *i;
  return stream << "]";
}
template <typename T>
std::ostream& operator<<(std::ostream& stream, std::set<T> const& s) {
  stream << "{";
  for (auto i = s.begin(); i != s.end(); ++i) stream << ((i == s.begin()) ? "" : ", ") << *i;
  return stream << "}";
}
template <typename T>
std::ostream& operator<<(std::ostream& stream, std::unordered_set<T> const& s) {
  stream << "{";
  for (auto i = s.begin(); i != s.end(); ++i) stream << ((i == s.begin()) ? "" : ", ") << *i;
  return stream << "}";
}
template <typename T>
std::ostream& operator<<(std::ostream& stream, std::queue<T>& q) {
  std::vector<T> v;
  for (; !q.empty(); q.pop()) v.push_back(q.front());
  for (auto& x : v) q.push(x);
  return stream << v;
}
template <typename T>
std::ostream& operator<<(std::ostream& stream, std::stack<T>& s) {
  std::vector<T> v;
  for (; !s.empty(); s.pop()) v.push_back(s.top());
  for (auto it = v.rbegin(); it != v.rend(); ++it) s.push(*it);
  return stream << v;
}
template <typename T>
std::ostream& operator<<(std::ostream& stream, std::priority_queue<T>& q) {
  std::vector<T> v;
  for (; !q.empty(); q.pop()) v.push_back(q.top());
  for (auto& x : v) q.push(x);
  return stream << v;
}

// #define TRC (std::cout << __FUNCTION__ << "(" << __LINE__ << ") ")
#if BUILD_DEV
#define PRINT(x)                                                                                                \
  {                                                                                                             \
    std::cerr << ANSI_BLUE << __FILENAME__ << ANSI_RESET << ":" << ANSI_YELLOW << __LINE__ << ANSI_RESET << ":" \
              << ANSI_GREEN << x << ANSI_RESET << std::endl;                                                  \
  }
#else
#define PRINT(x)

/* #define PRINT(x) \
   { std::cerr << __FILENAME__ << ":" << __LINE__ << ":" << x << std::endl; } */
#endif

#define cstr(x) (mog::util::to_string(x).c_str())
