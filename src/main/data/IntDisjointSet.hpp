#pragma once

#include <unordered_map>
#include <stdexcept>
#include <vector>

namespace mog {
namespace data {
/**
 * Disjoint set for integers (i.e. Union-find data structure).
 */
class IntDisjointSet {
 public:
  IntDisjointSet(std::size_t size) : data(size, -1) {}

  bool Union(int x, int y) {
    if ((x = root(x)) == (y = root(y))) return false;
    if (data[y] < data[x]) std::swap(x, y);
    data[x] += data[y];
    data[y] = x;
    return true;
  }
  int root(int x) { return data[x] < 0 ? x : data[x] = root(data[x]); }
  int rank(int x) { return -data[root(x)]; }

  std::vector<std::vector<int>> partition() {
    std::vector<std::vector<int>> ret;
    std::unordered_map<int, int> mapping;
    for (int i = 0; i < static_cast<int>(data.size()); ++i) {
      if (rank(i) == 1) continue;  // exclude singletons

      auto r = root(i);
      int index = 0;
      auto it = mapping.find(r);
      if (it == mapping.end()) {
        index = mapping.size();
        mapping[r] = index;
        ret.push_back({i});
      } else {
        ret[it->second].push_back(i);
      }
    }
    return ret;
  }

 private:
  std::vector<int> data; /* stores parent or count */
};

}  // namespace data
}  // namespace mog