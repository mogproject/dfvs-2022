#include <algorithm>
#include <stdexcept>

#include "../../util/util.hpp"
#include "CLGraph.hpp"

namespace mog {
namespace data {
namespace graph {

//------------------------------------------------------------
//   Enumeration
//------------------------------------------------------------
std::vector<int> const& CLGraph::neighbors(int v) const {
  assert_vertex_existence(v, "neighbors");
  clean_(v);

  // debug
  // for (auto u: edges_[v]) {
  //   if (v == u) throw std::runtime_error("loop found");
  //   if (!has_vertex(u)) throw std::runtime_error("invalid vertex");
  // }


  return edges_[v];
}
std::vector<int> const& CLGraph::in_neighbors(int v) const { return neighbors(v); }
std::vector<int> const& CLGraph::out_neighbors(int v) const { return neighbors(v); }

std::vector<int> CLGraph::neighbors_sorted(int v) const {
  assert_vertex_existence(v, "neighbors_sorted");
  clean_(v);

  std::vector<int> ret(edges_[v].begin(), edges_[v].end());
  std::sort(ret.begin(), ret.end());
  return ret;
}

std::vector<std::vector<int>> const& CLGraph::get_adjacency_list() const {
  clean_all();
  return edges_;
}

/**
 * @brief Enumerates all edges.
 */
std::vector<std::pair<int, int>> CLGraph::edges() const {
  std::vector<std::pair<int, int>> ret;
  for (int i : vertices()) {
    for (int j : neighbors(i)) {
      if (i < j) ret.push_back({i, j});
    }
  }
  return ret;
}

//------------------------------------------------------------
//   Modification
//------------------------------------------------------------
void CLGraph::add_vertex(int v) {
  BaseGraph::add_vertex(v);

  edges_[v].clear();
  is_dirty_ -= v;
  removed_edges_[v].clear();
  clean_all();  // make sure that there is no edge between v and another
}

// O(1)
void CLGraph::add_edge(int u, int v) {
  assert_vertex_existence(u, "add_edge (u)");
  assert_vertex_existence(v, "add_edge (v)");
  if (u == v) throw std::invalid_argument("loops are not allowed");

  // never checks if the given edge already exists
  edges_[u].push_back(v);
  edges_[v].push_back(u);
}

void CLGraph::remove_edge(int u, int v) {
  assert_vertex_existence(u, "remove_edge (u)");
  assert_vertex_existence(v, "remove_edge (v)");
  if (u == v) throw std::invalid_argument("loops are not allowed");

  removed_edges_[u].push_back(v);
  removed_edges_[v].push_back(u);
  is_dirty_ |= u;
  is_dirty_ |= v;
  // util::erase_first(edges_[u], v);
  // util::erase_first(edges_[v], u);
}

void CLGraph::remove_vertex(int v) {
  assert_vertex_existence(v, "remove_vertex");

  for (auto u : edges_[v]) {
    removed_edges_[u].push_back(v);
    is_dirty_ |= u;
    // util::erase_first(edges_[u], v);
  }
  BaseGraph::remove_vertex(v);
}

//============================================================
//   Lazy Edge Deletion
//============================================================
void CLGraph::clean_all() const {
  for (auto v : is_dirty_.to_vector()) clean_(v);
}

void CLGraph::clean_(int v) const {
  if (is_dirty_[v]) {
    if (!removed_edges_[v].empty()) {
      fset_.clear();
      for (auto x : removed_edges_[v]) fset_.set(x);

      int i = 0;
      while (i < static_cast<int>(edges_[v].size())) {
        auto u = edges_[v][i];
        if (fset_.get(u)) {
          // this vertex has been removed
          edges_[v][i] = edges_[v].back();
          edges_[v].pop_back();
        } else {
          ++i;
        }
      }
      removed_edges_[v].clear();
    }

    is_dirty_ -= v;
  }
}

void CLGraph::resize(int n) {
  if (n <= capacity()) return;  // do nothing

  clean_all();

  vertices_.resize(n);
  fset_.resize(n);
  edges_.resize(n);
  is_dirty_.resize(n);
  removed_edges_.resize(n);
}
//------------------------------------------------------------
//   Queries
//------------------------------------------------------------

// O(deg(u) + deg(v))
bool CLGraph::has_edge(int u, int v) const {
  if (u == v) return false;

  if (degree(u) < degree(v)) {  // cleans here
    return util::contains(edges_[u], v);
  } else {
    return util::contains(edges_[v], u);
  }
}

int CLGraph::degree(int v) const {
  clean_(v);
  return edges_[v].size();
}

//------------------------------------------------------------
//   Conversions
//------------------------------------------------------------
CLGraph CLGraph::compact() const {
  auto inv = util::inverse(vertices());
  CLGraph g(number_of_nodes());
  for (auto v : vertices()) {
    int i = inv[v];
    for (auto u : neighbors(v)) {
      if (v < u) g.add_edge(i, inv[u]);
    }
  }
  return g;
}
}  // namespace graph
}  // namespace data
}  // namespace mog