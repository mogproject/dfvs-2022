#include <queue>
#include <stack>

#include "../data/IntDisjointSet.hpp"
#include "clique.hpp"

namespace mog {
namespace algorithm {

// key:
// pair of
// 1. the number of common neighbors (the larger, the more priority) and
// 2. the negative of the sum of the degrees of u and v

static inline long long encode_edge(int u, int v) { return ((long long)std::min(u, v) << 32) | std::max(u, v); }

static inline std::pair<int, int> decode_edge(long long e) {
  int u = e >> 32;
  int v = e - (e >> 32 << 32);
  return {u, v};
}

/**
 * @brief Clique partitioning (not necessarily optimal)
 * Based on Tseng and Siewiorek 1986
 *
 * @tparam Graph grpah type
 * @param g input graph
 * @return std::vector<std::vector<int>> list of parts excluding isolates
 */

template <typename Graph>
std::vector<std::vector<int>> clique_partition(Graph& g) {
  if (g.capacity() != g.number_of_nodes()) throw std::invalid_argument("vertices must be compact");

  typedef std::pair<long long, long long> Key;
  int n = g.number_of_nodes();
  data::fast_set fset(n);
  data::IntDisjointSet dset(n);

  std::map<long long, long long> current;
  std::priority_queue<Key> q;

  auto f = [&](int i, bool remove_symmetry) {
    // printf("Update: %d\n", i);
    fset.clear();
    for (int j : g.neighbors(i)) fset.set(j);  // i's neighbors

    for (int j : g.neighbors(i)) {
      if (remove_symmetry && i >= j) continue;
      // consider edge ij (i<j)
      int cnt = 0;
      for (int w : g.neighbors(j)) {
        if (fset.get(w)) ++cnt;  // w is a common neighbor of ij
      }
      long long e = encode_edge(i, j);
      long long data = ((long long)cnt << 32) + (1LL << 32) - g.degree(i) - g.degree(j);
      current[e] = data;
      // printf("Push: edge=(%d-%d), common=%d, degsum=%d\n", i, j, cnt, g.degree(i) + g.degree(j));
      q.push({data, e});
    }
  };

  // initial edge info
  for (int i = 0; i < n; ++i) f(i, true);

  while (!q.empty()) {
    auto key = q.top();
    q.pop();

    if (key.first != current[key.second]) continue;  // outdated entry

    // contract edge and remove edges to uncommon neighbors
    auto p = decode_edge(key.second);
    int u = p.first;
    int v = p.second;
    dset.Union(u, v);
    // printf("Union: edge=(%d-%d), common=%d, degsum=%d\n", u, v, int(key.first >> 32), -(int)(key.first & 0xffffffff));

    std::vector<int> affected;  // list of affected vertices

    fset.clear();
    auto nbrs = g.neighbors(u);
    for (int j : nbrs) {
      g.remove_edge(u, j);
      current[encode_edge(u, j)] = -1;  // invalidate

      if (j == v) continue;
      affected.push_back(j);
      fset.set(j);  // u's neighbors
    }

    for (int w : g.neighbors(v)) {
      if (w == u) continue;
      if (fset.get(w)) continue;  // common neighbors; keep them as is

      affected.push_back(w);
      g.remove_edge(v, w);
      current[encode_edge(v, w)] = -1;  // invalidate
    }

    // register new info
    for (auto i : affected) f(i, false);
  }
  return dset.partition();
}

template <typename Graph>
std::vector<std::vector<int>> clique_partition(Graph const& g) {
  auto gg = g;
  return clique_partition(gg);
}

// explicit instantiation
template std::vector<std::vector<int>> clique_partition(data::graph::CLGraph& g);
template std::vector<std::vector<int>> clique_partition(data::graph::CLGraph const& g);
}  // namespace algorithm
}  // namespace mog