#include <queue>
#include <stack>

#include "../data/bitmap.hpp"
#include "../data/fast_set.hpp"
#include "../data/graph/graph.hpp"
#include "../util/profiler.hpp"
#include "../util/util.hpp"

namespace mog {
namespace algorithm {

template <typename Graph>
std::vector<std::vector<int>> components(Graph const &g) {
  std::vector<std::vector<int>> ret;
  data::fast_set visited(g.number_of_nodes());

  for (int v = 0; v < g.number_of_nodes(); ++v) {
    if (visited.get(v)) continue;
    ret.push_back(bfs(g, v, &visited));  // reuse `visited`
  }
  return ret;
}

template <typename DiGraph>
class StateTarjanSCC {
 public:
  DiGraph const &g;
  bool include_singletons;
  bool keep_result_set;
  std::vector<int> dfs_num;
  std::vector<int> dfs_low;
  std::stack<int> active;  // keep track of the current tree
  std::vector<std::vector<int>> result;

  StateTarjanSCC(DiGraph const &g, bool include_singletons, bool keep_result_set)
      : g(g),
        include_singletons(include_singletons),
        keep_result_set(keep_result_set),
        dfs_num(g.capacity(), -1),
        dfs_low(g.capacity()),
        dfs_counter_(0),
        is_active_(g.capacity()) {}

  void reset_dfs_num(int v) { dfs_num[v] = dfs_low[v] = dfs_counter_++; }

  inline bool is_active(int v) const { return is_active_[v]; }

  inline bool has_visited(int v) const { return dfs_num[v] != -1; }

  void push_active(int v) {
    active.push(v);
    is_active_ |= v;
  }

  int pop_active() {
    int ret = active.top();
    active.pop();
    is_active_ -= ret;  // reset flag
    return ret;
  }

 private:
  int dfs_counter_;
  data::Bitmap is_active_;  // bit on if the vertex is in `active`
};

template <typename DiGraph>
static void tarjan_scc(int u, StateTarjanSCC<DiGraph> &state) {
  state.reset_dfs_num(u);
  state.push_active(u);

  for (auto v : state.g.out_neighbors(u)) {
    if (!state.has_visited(v)) tarjan_scc(v, state);  // recursion
    if (state.is_active(v)) state.dfs_low[u] = std::min(state.dfs_low[u], state.dfs_low[v]);
  }

  if (state.dfs_low[u] == state.dfs_num[u]) {  // found SCC
    std::vector<int> scc;
    int v;
    do {
      v = state.pop_active();
      scc.push_back(v);
      if (!state.keep_result_set) state.dfs_low[v] = state.dfs_num[u];  // for scc_ids()
    } while (v != u);
    if (state.keep_result_set && (scc.size() > 1 || state.include_singletons)) state.result.push_back(scc);
  }
}

#define USE_SCC_IMPL 2
// TODO: find a faster implementation without recursion

/**
 * @brief Computes the strongly connected components (scc) by Tarjan's algorithm.
 *
 * @tparam DiGraph
 * @param g Input graph
 * @return list of scc's
 */
template <typename DiGraph>
std::vector<std::vector<int>> strongly_connected_components(DiGraph const &g, bool include_singletons, util::Profiler *prof) {
  PROF(if (prof) prof->start_timer("strongly_connected_comp()"));
  // int n = g.capacity();

#if USE_SCC_IMPL == 1
  // e: 183 us, h: 1219 us
  StateTarjanSCC<DiGraph> state(g, include_singletons, true);

  for (auto v : g.vertices()) {
    if (state.dfs_num[v] == -1) tarjan_scc(v, state);
  }
  PROF(if (prof) prof->stop_timer("strongly_connected_comp()"));
  return state.result;
#elif USE_SCC_IMPL == 2
  // e: 369 us, h: 1924 us
  StateTarjanSCC<DiGraph> state(g, include_singletons, true);
  std::stack<std::pair<int, int>> control;

  for (auto v : g.vertices()) {
    if (state.has_visited(v)) continue;  // already visited

    control.push({v, -2});

    while (!control.empty()) {
      auto p = control.top();
      control.pop();
      auto u = p.first;

      switch (p.second) {
        case -2:
          if (state.has_visited(u)) continue;
          state.reset_dfs_num(u);
          state.push_active(u);
          control.push({u, -1});

          for (auto w : g.out_neighbors(u)) {
            if (!state.has_visited(w)) {
              control.push({u, w});
              control.push({w, -2});  // recursion
            } else {
              if (state.is_active(w)) state.dfs_low[u] = std::min(state.dfs_low[u], state.dfs_low[w]);
            }
          }
          break;
        case -1:
          // check if there is an SCC
          if (state.dfs_low[u] == state.dfs_num[u]) {
            std::vector<int> scc;
            int w;
            do {
              w = state.pop_active();
              scc.push_back(w);
            } while (w != u);
            if (state.keep_result_set && (scc.size() > 1 || state.include_singletons)) state.result.push_back(scc);
          }
          break;
        default:
          // update dfs low
          if (state.is_active(p.second)) state.dfs_low[u] = std::min(state.dfs_low[u], state.dfs_low[p.second]);
      }
    }
  }
  PROF(if (prof) prof->stop_timer("strongly_connected_comp()"));
  return state.result;
#elif USE_SCC_IMPL == 3
  // e: 367 us, h: 1777 us
  std::vector<std::vector<int>> ret;
  std::vector<int> ts(n);  // vertex-id -> tree-number
  std::vector<int> ps(n);  // vertex-id -> parent
  std::vector<int> pre;    // dfs-index -> vertex-id
  std::vector<int> ws(n);  // dfs-low
  data::Bitmap visited(n);

  int tree_number = 0;

  for (auto v : g.vertices()) {
    if (visited[v]) continue;

    std::stack<int> st;
    st.push(v);

    while (!st.empty()) {
      auto u = st.top();
      st.pop();
      if (visited[u]) continue;

      visited |= u;
      ts[u] = tree_number;
      ws[u] = pre.size();  // dfs-index
      pre.push_back(u);

      for (auto w : g.out_neighbors(u)) {
        if (!visited[w]) {
          st.push(w);
          ps[w] = u;  // update parent
        }
      }
    }
    ++tree_number;
  }

  int nn = pre.size();

  // post-order traversal
  for (int i = nn - 1; i >= 0; --i) {
    int v = pre[i];
    for (auto u : g.out_neighbors(v)) {
      if (ts[v] == ts[u]) ws[v] = std::min(ws[v], ws[u]);
    }
  }

  // std::cout << "pre: " << pre << std::endl;
  // std::cout << "ps: " << ps << std::endl;
  // std::cout << "ws: " << ws << std::endl;

  // pre-order traversal
  int rep_idx = -1;
  std::stack<int> rep;  // stores the dfs-low of each representative
  for (int i = 0; i < nn; ++i) {
    int v = pre[i];
    if (i != ws[v]) {
      ws[v] = ws[ps[v]];  // inherit parent's dfs-low
      while (ws[v] != rep.top()) {
        rep.pop();
        --rep_idx;
      }
      ret[rep_idx].push_back(v);
    } else {  // representative
      rep.push(ws[v]);
      rep_idx = ret.size();
      ret.push_back({v});
    }
  }
  PROF(if (prof) prof->stop_timer("strongly_connected_comp()"));
  return ret;
#endif
}

template <typename DiGraph>
std::vector<int> scc_ids(DiGraph const &g) {
  std::vector<int> ret(g.capacity(), -1);
  auto sccs = strongly_connected_components(g, true, nullptr);
  for (int i = 0; i < static_cast<int>(sccs.size()); ++i) {
    for (auto v: sccs[i]) ret[v] = i;
  }
  return ret;
}

template <typename DiGraph>
bool is_acyclic(DiGraph const &g) {
  return strongly_connected_components<DiGraph>(g, false, nullptr).empty();
}

template <typename DiGraph>
static int number_of_reachable_nodes(DiGraph const &g, int v, int to_avoid, bool out_going) {
  // Assumption: v != to_avoid
  data::fast_set visited(g.capacity());
  std::queue<int> q;
  q.push(v);
  visited.set(v);
  int ret = 1;

  while (!q.empty()) {
    auto u = q.front();
    q.pop();
    for (auto w : out_going ? g.out_neighbors(u) : g.in_neighbors(u)) {
      if (w == to_avoid) continue;
      if (!visited.get(w)) {
        ++ret;
        q.push(w);
        visited.set(w);
      }
    }
  }
  return ret;
}

/**
 * @brief Tests if the given graph with an avoiding vertex is strongly connected,
 * using Kosaraju’s algorithm.
 *
 * @param g digraph
 * @param to_avoid avoiding vertex
 * @return true if g is strongly connected; false otherwise
 */
template <typename DiGraph>
static bool is_strongly_connected(DiGraph const &g, int to_avoid) {
  auto vtcs = g.vertices();
  int n = vtcs.size();

  if (n == 0) throw std::runtime_error("empty graph");
  if (n == 1) return true;  // empty graph is by definition connected

  auto v = vtcs[0] == to_avoid ? vtcs[1] : vtcs[0];
  return number_of_reachable_nodes(g, v, to_avoid, true) == n - 1 && number_of_reachable_nodes(g, v, to_avoid, false) == n - 1;
}

template <typename DiGraph>
std::vector<int> strong_articulation_points(DiGraph const &g) {
  // Assumptions:
  //   - g is strongly connected
  //   - g has is no loops
  std::vector<int> ret;

  // TODO: implement linear algorithm
  for (auto v : g.vertices()) {
    if (!is_strongly_connected(g, v)) ret.push_back(v);
  }
  return ret;
}

// explicit instantiation
template std::vector<std::vector<int>> strongly_connected_components(data::graph::CLDiGraph const &g,
                                                                     bool include_singletons, util::Profiler *prof);
template std::vector<int> scc_ids(data::graph::CLDiGraph const &g);
template bool is_acyclic(data::graph::CLDiGraph const &g);
template std::vector<int> strong_articulation_points(data::graph::CLDiGraph const &g);

}  // namespace algorithm
}  // namespace mog