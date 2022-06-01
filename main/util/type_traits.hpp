#pragma once

#include <map>
#include <set>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

/**
 * @brief Check if the given type parameter is std::map or std::unordered_map.
 *
 * @tparam T type parameter to check
 */
template <typename T>
struct is_map : std::false_type {};

template <typename K, typename V>
struct is_map<std::map<K, V>> : std::true_type {};

template <typename K, typename V>
struct is_map<std::unordered_map<K, V>> : std::true_type {};

/**
 * @brief Check if the given type parameter is std::set or std::unordered_set.
 *
 * @tparam T type parameter to check
 */
template <typename T>
struct is_set : std::false_type {};

template <typename U>
struct is_set<std::set<U>> : std::true_type {};

template <typename U>
struct is_set<std::unordered_set<U>> : std::true_type {};
