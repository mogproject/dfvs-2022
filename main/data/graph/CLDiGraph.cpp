#include <algorithm>
#include <stdexcept>

#include "CLDiGraph.hpp"

namespace mog {
namespace data {
namespace graph {

//------------------------------------------------------------
//   Enumeration
//------------------------------------------------------------
// does not include the loop
std::vector<int> const& CLDiGraph::in_neighbors(int v) const {
  assert_vertex_existence(v, "in_neighbors");
  CLEAN(v);
  return in_edges_[v];
}

// does not include the loop
std::vector<int> const& CLDiGraph::out_neighbors(int v) const {
  assert_vertex_existence(v, "out_neighbors");
  CLEAN(v);
  return out_edges_[v];
}

/**
 * @brief Returns out-edges + in-edges.
 *
 * @param v
 * @return out-edges followed by in-edges
 * The same vertex may appear twice if there is a bidirected edge.
 */
std::vector<int> CLDiGraph::neighbors(int v) const {
  assert_vertex_existence(v, "neighbors");
  CLEAN(v);
  auto ret = out_edges_[v];
  ret.reserve(ret.size() + std::distance(in_edges_[v].begin(), in_edges_[v].end()));
  ret.insert(ret.end(), in_edges_[v].begin(), in_edges_[v].end());
  return ret;
}

/**
 * @brief Returns the intersection of out-neighbors and in-neighbors.
 *
 * O(deg^+(v) + deg^-(v))
 *
 * @param v
 * @return intersection of out-neighbors and in-neighbors
 */
std::vector<int> CLDiGraph::strong_neighbors(int v) const {
  assert_vertex_existence(v, "strong_neighbors");
  CLEAN(v);
  std::vector<int> ret;

  fset_.clear();
  for (int u : out_edges_[v]) fset_.set(u);
  for (int w : in_edges_[v]) {
    if (fset_.get(w)) ret.push_back(w);  // found an intersecting neighbor
  }
  return ret;
}

std::vector<int> CLDiGraph::in_neighbors_sorted(int v) const {
  assert_vertex_existence(v, "in_neighbors_sorted");
  CLEAN(v);
  std::vector<int> ret(in_edges_[v].begin(), in_edges_[v].end());
  std::sort(ret.begin(), ret.end());
  return ret;
}

std::vector<int> CLDiGraph::out_neighbors_sorted(int v) const {
  assert_vertex_existence(v, "out_neighbors_sorted");
  CLEAN(v);
  std::vector<int> ret(out_edges_[v].begin(), out_edges_[v].end());
  std::sort(ret.begin(), ret.end());
  return ret;
}

/**
 * @brief Enumerates all edges.
 */
std::vector<std::pair<int, int>> CLDiGraph::edges(bool use_in_neighbors) const {
  std::vector<std::pair<int, int>> ret;
  for (int i : vertices()) {
    if (use_in_neighbors) {
      for (int j : in_neighbors(i)) ret.push_back({j, i});
    } else {
      for (int j : out_neighbors(i)) ret.push_back({i, j});
    }
  }
  return ret;
}

//------------------------------------------------------------
//   Modification
//------------------------------------------------------------
void CLDiGraph::add_vertex(int v) {
  BaseDiGraph::add_vertex(v);

  in_edges_[v].clear();
  out_edges_[v].clear();
  has_loop_ -= v;
#if CLDIGRAPH_OPT == 1
  is_dirty_ -= v;
  removed_in_edges_[v].clear();
  removed_out_edges_[v].clear();
  clean_all();  // make sure that there is no edge between v and another
#endif
}

// O(deg(u))
void CLDiGraph::add_edge(int u, int v) {
  assert_vertex_existence(u);
  assert_vertex_existence(v);

  // needs to check if the edge exists to maintain data consistency
  if (has_edge(u, v)) throw std::invalid_argument("uv: already exists");
  add_edge_unsafe(u, v);
}

// O(1)
void CLDiGraph::add_edge_unsafe(int u, int v) {
  if (u == v) {
    has_loop_ |= u;
  } else {
    out_edges_[u].push_back(v);
    in_edges_[v].push_back(u);
  }
}

// O(deg(u) + deg(v))
void CLDiGraph::remove_edge(int u, int v) {
  assert_vertex_existence(u, "remove_edge");
  assert_vertex_existence(v, "remove_edge");

  if (u == v) {
    has_loop_ -= u;
  } else {
#if CLDIGRAPH_OPT == 1
    removed_out_edges_[u].push_back(v);
    removed_in_edges_[v].push_back(u);
    is_dirty_ |= u;
    is_dirty_ |= v;
#else
    util::erase_first(out_edges_[u], v);
    util::erase_first(in_edges_[v], u);
#endif
  }
}

//============================================================
//   Vertex Deletion
//============================================================
// O(|N^2[v]|)
void CLDiGraph::remove_vertex(int v) {
  assert_vertex_existence(v, "remove_vertex");

  // for (auto u : out_edges_[v]) {
  for (auto u : out_neighbors(v)) {
#if CLDIGRAPH_OPT == 1
    if (v != u) {
      removed_in_edges_[u].push_back(v);
      is_dirty_ |= u;
    }
#else
    if (v != u) util::erase_first(in_edges_[u], v);
#endif
  }
  // for (auto w : in_edges_[v]) {
  for (auto w : in_neighbors(v)) {
#if CLDIGRAPH_OPT == 1
    if (v != w) {
      removed_out_edges_[w].push_back(v);
      is_dirty_ |= w;
    }
#else
    if (v != w) util::erase_first(out_edges_[w], v);
#endif
  }
  BaseDiGraph::remove_vertex(v);
}

// CLDiGraph::Edges CLDiGraph::remove_vertex_with_trace(int v) {
//   CLDiGraph::Edges removed_edges;

//   // trace all the edges to remove
//   CLEAN(v);
//   for (auto u : out_edges_[v]) removed_edges.push_back({v, u});
//   for (auto w : in_edges_[v]) removed_edges.push_back({w, v});
//   if (has_loop_[v]) removed_edges.push_back({v, v});

//   remove_vertex(v);
//   return removed_edges;
// }

//============================================================
//   Edge Contraction
//============================================================
// O(deg(v) |N^2[v]|)
void CLDiGraph::contract_edge(int u, int v, bool keep_head) {
  assert_vertex_existence(u, "contract_edge");
  assert_vertex_existence(v, "contract_edge");
  if (u == v) return;  // do nothing

  int to_keep = keep_head ? u : v;
  int to_remove = keep_head ? v : u;

  for (auto x : in_neighbors(to_remove)) {
    if (x == to_keep) {
      if (!keep_head) has_loop_ |= to_keep;
    } else if (x != to_remove && !has_edge(x, to_keep)) {
      add_edge_unsafe(x, to_keep);
    }
  }
  for (auto y : out_neighbors(to_remove)) {
    if (y == to_keep) {
      if (keep_head) has_loop_ |= to_keep;
    } else if (y != to_remove && !has_edge(to_keep, y)) {
      add_edge_unsafe(to_keep, y);
    }
  }
  if (has_loop_[to_remove]) has_loop_ |= to_keep;
  remove_vertex(to_remove);
}

// std::pair<CLDiGraph::Edges, CLDiGraph::Edges> CLDiGraph::contract_edge_with_trace(int u, int v, bool keep_head) {
//   assert_vertex_existence(u);
//   assert_vertex_existence(v);
//   if (u == v) return {};  // do nothing

//   int to_keep = keep_head ? u : v;
//   int to_remove = keep_head ? v : u;

//   CLDiGraph::Edges added_edges;

//   for (auto x : in_neighbors(to_remove)) {
//     if (x == to_keep) {
//       if (!keep_head) has_loop_ |= to_keep;
//     } else if (x != to_remove && !has_edge(x, to_keep)) {
//       added_edges.push_back({x, to_keep});
//       add_edge_unsafe(x, to_keep);
//     }
//   }
//   for (auto y : out_neighbors(to_remove)) {
//     if (y == to_keep) {
//       if (keep_head) has_loop_ |= to_keep;
//     } else if (y != to_remove && !has_edge(to_keep, y)) {
//       added_edges.push_back({to_keep, y});
//       add_edge_unsafe(to_keep, y);
//     }
//   }
//   if (has_loop_[to_remove]) has_loop_ |= to_keep;
//   auto removed_edges = remove_vertex_with_trace(to_remove);

//   return {removed_edges, added_edges};
// }

//============================================================
//   Vertex Ignoration
//============================================================
// O(|N^2[v]|)
void CLDiGraph::ignore_vertex(int v) {
  assert_vertex_existence(v, "ignore_vertex");
  CLEAN(v);
  for (auto u : in_edges_[v]) {
    for (auto w : out_edges_[v]) {
      if (!has_edge(u, w)) add_edge_unsafe(u, w);
    }
  }
  remove_vertex(v);
}

// std::pair<CLDiGraph::Edges, CLDiGraph::Edges> CLDiGraph::ignore_vertex_with_trace(int v) {
//   assert_vertex_existence(v);

//   CLDiGraph::Edges added_edges;
//   CLEAN(v);
//   for (auto u : in_edges_[v]) {
//     for (auto w : out_edges_[v]) {
//       if (!has_edge(u, w)) {
//         added_edges.push_back({u, w});
//         add_edge_unsafe(u, w);
//       }
//     }
//   }
//   auto removed_edges = remove_vertex_with_trace(v);
//   return {removed_edges, added_edges};
// }

//============================================================
//   Lazy Edge Deletion
//============================================================
#if CLDIGRAPH_OPT == 1
void CLDiGraph::clean_all() const {
  for (auto v : is_dirty_.to_vector()) clean_(v);
}

void CLDiGraph::clean_(int v) const {
  if (is_dirty_[v]) {
    if (!removed_in_edges_[v].empty()) {
      // new implementation
      fset_.clear();
      for (auto x : removed_in_edges_[v]) fset_.set(x);

      int i = 0;
      while (i < static_cast<int>(in_edges_[v].size())) {
        auto u = in_edges_[v][i];
        if (fset_.get(u)) {
          // this vertex has been removed
          in_edges_[v][i] = in_edges_[v].back();
          in_edges_[v].pop_back();
        } else {
          ++i;
        }
      }

      // old implementation
      // int j = 0;
      // for (int i = 0; i < in_edges_[v].size(); ++i) {
      //   auto u = in_edges_[v][i];
      //   if (!util::contains(removed_in_edges_[v], u)) {
      //     if (i != j) in_edges_[v][j] = in_edges_[v][i];
      //     ++j;
      //   }
      // }
      // in_edges_[v].resize(j);
      removed_in_edges_[v].clear();
    }

    if (!removed_out_edges_[v].empty()) {
      // new implementation
      fset_.clear();
      for (auto x : removed_out_edges_[v]) fset_.set(x);

      int i = 0;
      while (i < static_cast<int>(out_edges_[v].size())) {
        auto u = out_edges_[v][i];
        if (fset_.get(u)) {
          // this vertex has been removed
          out_edges_[v][i] = out_edges_[v].back();
          out_edges_[v].pop_back();
        } else {
          ++i;
        }
      }

      // old implementation
      // int j = 0;
      // for (int i = 0; i < out_edges_[v].size(); ++i) {
      //   auto u = out_edges_[v][i];
      //   if (!util::contains(removed_out_edges_[v], u)) {
      //     if (i != j) out_edges_[v][j] = out_edges_[v][i];
      //     ++j;
      //   }
      // }
      // out_edges_[v].resize(j);
      removed_out_edges_[v].clear();
    }

    is_dirty_ -= v;
  }
}
#endif

void CLDiGraph::clear() {
  for (auto v : vertices()) remove_vertex(v);
}

void CLDiGraph::resize(int n) {
  if (n <= capacity()) return;  // do nothing

  clean_all();

  vertices_.resize(n);
  has_loop_.resize(n);
  fset_.resize(n);
  in_edges_.resize(n);
  out_edges_.resize(n);
  is_dirty_.resize(n);
  removed_in_edges_.resize(n);
  removed_out_edges_.resize(n);
}

//------------------------------------------------------------
//   Queries
//------------------------------------------------------------

// O(min{deg(u), deg(v)})
bool CLDiGraph::has_edge(int u, int v) const {
  assert_vertex_existence(u);
  assert_vertex_existence(v);

  if (u == v) return has_loop_[u];

  if (out_degree(u) < in_degree(v)) {  // cleans here
    return util::contains(out_edges_[u], v);
  } else {
    return util::contains(in_edges_[v], u);
  }
}

// O(1)
int CLDiGraph::in_degree(int v) const {
  assert_vertex_existence(v);
  CLEAN(v);
  return in_edges_[v].size() + (has_loop_[v] ? 1 : 0);
}

// O(1)
int CLDiGraph::out_degree(int v) const {
  assert_vertex_existence(v);
  CLEAN(v);
  return out_edges_[v].size() + (has_loop_[v] ? 1 : 0);
}

/**
 * @brief Checks if the given vertex is a strong vertex.
 * A vertex is strong if all its neighbors are strong neighbors.
 *
 * @param v vertex
 * @return true if v is a strong vertex
 */
bool CLDiGraph::is_strong_vertex(int v) const { return out_degree(v) == strong_degree(v); }

// O(deg+(v) + deg-(v))
// TODO: cache values
int CLDiGraph::strong_degree(int v) const {
  assert_vertex_existence(v);
  CLEAN(v);

  int ret = 0;
  fset_.clear();
  for (int u : out_edges_[v]) fset_.set(u);
  for (int w : in_edges_[v]) {
    if (fset_.get(w)) ++ret;  // found an intersecting neighbor
  }
  return ret;
}

//------------------------------------------------------------
//   Conversions
//------------------------------------------------------------
/**
 * @brief Creates a strong graph of this digraph. Loops are ignored.
 *
 * @return CLGraph undirected graph
 */
CLGraph CLDiGraph::strong_graph() const {
  CLGraph g(vertices_);  // g's edges <- h's strong edges

  for (auto v : vertices()) {
    for (auto u : strong_neighbors(v)) {
      if (v < u) g.add_edge(v, u);
    }
  }
  return g;
}

CLDiGraph CLDiGraph::weak_graph() const {
  CLDiGraph g(vertices_);

  clean_all();

  for (auto v : vertices()) {
    fset_.clear();
    for (int w : in_edges_[v]) fset_.set(w);
    for (int u : out_edges_[v]) {
      if (!fset_.get(u)) g.add_edge_unsafe(v, u);  // vu is a weak edge
    }
  }
  return g;
}

CLDiGraph CLDiGraph::reversal_graph() const {
  CLDiGraph g(vertices_);

  g.has_loop_ = has_loop_;

  for (auto v : vertices()) {
    for (auto u : out_neighbors(v)) g.add_edge_unsafe(u, v);
  }
  return g;
}

CLDiGraph CLDiGraph::compact() const {
  auto inv = util::inverse(vertices());
  CLDiGraph g(number_of_nodes());
  for (auto v : vertices()) {
    int i = inv[v];
    for (auto u : out_neighbors(v)) { g.add_edge_unsafe(i, inv[u]); }
  }
  return g;
}

CLGraph CLDiGraph::to_undirected() const {
  CLGraph g(vertices_);

  for (auto v : vertices()) {
    fset_.clear();
    for (auto u : in_neighbors(v)) fset_.set(u);
    for (auto u : out_neighbors(v)) {
      if (fset_.get(u)) {
        if (v < u) g.add_edge(v, u);  // strong edge
      } else {
        g.add_edge(v, u);  // weak edge
      }
    }
  }
  return g;
}

//------------------------------------------------------------
//   Debugging
//------------------------------------------------------------

void CLDiGraph::check_consistency() const {
  auto m = number_of_edges();
  long long out_degree_sum = 0, in_degree_sum = 0;
  for (auto i : vertices()) {
    fset_.clear();
    for (auto j : out_neighbors(i)) {
      if (!has_vertex(j)) throw std::invalid_argument(util::format("invalid out-neighbor: %d->%d", i, j));
      if (fset_.get(j)) throw std::invalid_argument(util::format("duplicate out-neighbor: %d->%d", i, j));
      if (!util::contains(in_neighbors(j), i)) {
        throw std::invalid_argument(util::format("edge inconsistency (%d->%d): %d not an in-neighbor of %d", i, j, j, i));
      }
      fset_.set(j);
      ++out_degree_sum;
    }
    fset_.clear();
    for (auto j : in_neighbors(i)) {
      if (!has_vertex(j)) throw std::invalid_argument(util::format("invalid in-neighbor: %d->%d", i, j));
      if (fset_.get(j)) throw std::invalid_argument(util::format("duplicate in-neighbor: %d->%d", j, i));
      if (!util::contains(out_neighbors(j), i)) {
        throw std::invalid_argument(util::format("edge inconsistency (%d->%d): %d not an out-neighbor of %d", j, i, i, j));
      }
      fset_.set(j);
      ++in_degree_sum;
    }
  }
  if (out_degree_sum != in_degree_sum || m != out_degree_sum) {
    throw std::invalid_argument(util::format("degree inconsistency: m=%lld, out_degree_sum=%lld, in_degree_sum=%lld", m,
                                             out_degree_sum, in_degree_sum));
  }
}

}  // namespace graph
}  // namespace data
}  // namespace mog