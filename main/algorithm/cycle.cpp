#include <cmath>
#include <queue>
#include <stack>

#include "../data/IntDisjointSet.hpp"
#include "../data/fast_set.hpp"
#include "cycle.hpp"

namespace mog {
namespace algorithm {

/**
 * @brief Finding a short cycle by sampling.
 *
 * @tparam Graph graph type
 * @param g input graph
 * @param sample number of samples
 * @return std::vector<int> vertices that form a cycle
 */
template <typename Graph>
std::vector<int> short_cycle(Graph const& g, int sample, int start) {
  static int const P = 1000000007;  // some prime number greater than n

  int n = g.number_of_nodes();
  if (n == 0) return {};

  int cap = g.capacity();
  auto vtcs = g.vertices();
  if (sample < 0) sample = static_cast<int>(std::ceil(std::sqrt(n)));

  std::vector<long long> dist(cap);
  std::vector<int> prev(cap);

  std::vector<int> ret;
  int best = n;
  int y = start % n;
  for (int t = 1; t <= sample; ++t) {
    // bfs
    std::queue<int> q;

    int s = vtcs[y];  // starting vertex
    dist[s] = cap * t;
    q.push(s);

    while (!q.empty()) {
      bool updated = false;
      auto u = q.front();
      q.pop();

      if (dist[u] % cap + 1 >= best) continue;  // cannot improve the best

      for (auto w : g.in_neighbors(u)) {
        if (dist[w] < cap * t) {
          // not visited
          dist[w] = dist[u] + 1;
          prev[w] = u;  // for backtracking
          q.push(w);
        } else if (w == s) {
          // found a cycle
          best = dist[u] % cap + 1;
          ret.clear();
          ret.push_back(s);
          for (int x = u; x != s; x = prev[x]) ret.push_back(x);
          updated = true;
          break;
        }
      }
      if (updated) break;
    }

    y = (y + P) % n;
  }
  return ret;
}

template <typename Graph>
std::vector<int> find_cycle(Graph const& g, int v) {
  if (g.has_edge(v, v)) return {v};  // found loop

  std::vector<int> prev(g.capacity(), -1);

  // bfs
  std::queue<int> q;
  q.push(v);

  while (!q.empty()) {
    auto u = q.front();
    q.pop();

    for (auto w : g.in_neighbors(u)) {
      if (w == v) {
        // found a cycle
        std::vector<int> ret = {v};
        for (int x = u; x != v; x = prev[x]) ret.push_back(x);
        return ret;
      } else if (prev[w] == -1) {
        // not visited
        prev[w] = u;  // for backtracking
        q.push(w);
      }
    }
  }

  return {};
}

// explicit instantiation
template std::vector<int> short_cycle(data::graph::CLDiGraph const& g, int sample, int start);

template std::vector<int> find_cycle(data::graph::CLDiGraph const& g, int v);
}  // namespace algorithm
}  // namespace mog