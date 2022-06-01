#pragma once
#include <unordered_set>
#include <vector>

#include "../fast_set.hpp"
#include "BaseGraph.hpp"

namespace mog {
namespace data {
namespace graph {

/**
 * @brief Simple graph representation with adjacency lists (L for List). Vertices are labeled by 0-indexed integers.
 *
 * Space efficiency: O(n+m)
 */
class CLGraph : public BaseGraph {
 public:
  typedef int Vertex;
  typedef std::vector<Vertex> Vertices;
  typedef std::pair<Vertex, Vertex> Edge;
  typedef std::vector<Edge> Edges;

  CLGraph(std::size_t n = 0, std::vector<std::pair<int, int>> const& edges = {})
      : BaseGraph(n), edges_(n), is_dirty_(n), removed_edges_(n), fset_(n) {
    for (auto& p : edges) add_edge(p.first, p.second);
  }

  CLGraph(Bitmap const& vertices, std::vector<std::pair<int, int>> const& edges = {})
      : BaseGraph(vertices),
        edges_(vertices.size()),
        is_dirty_(vertices.size()),
        removed_edges_(vertices.size()),
        fset_(vertices.size()) {
    for (auto& p : edges) add_edge(p.first, p.second);
  }

  // Copy constructor
  CLGraph(CLGraph const& g)
      : BaseGraph(g), edges_(g.edges_), is_dirty_(g.is_dirty_), removed_edges_(g.removed_edges_), fset_(g.capacity()) {}

  //------------------------------------------------------------
  //   Properties
  //------------------------------------------------------------
  static constexpr bool const directed = false;
  static constexpr bool const use_matrix = false;

  //------------------------------------------------------------
  //   Enumeration
  //------------------------------------------------------------
  std::vector<int> const& in_neighbors(int v) const;   // alias of neighbors()
  std::vector<int> const& out_neighbors(int v) const;  // alias of neighbors()
  std::vector<int> neighbors_sorted(int v) const;
  std::vector<int> const& neighbors(int v) const;
  std::vector<std::vector<int>> const& get_adjacency_list() const;
  std::vector<std::pair<int, int>> edges() const;

  //------------------------------------------------------------
  //   Modification
  //------------------------------------------------------------
  void add_vertex(int v);
  void add_edge(int u, int v);
  void remove_edge(int u, int v);
  void remove_vertex(int v);

  void clean_all() const;
  void clean_(int v) const;

  void resize(int n);

  //------------------------------------------------------------
  //   Queries
  //------------------------------------------------------------
  bool has_edge(int u, int v) const;
  int degree(int v) const;

  //------------------------------------------------------------
  //   Conversions
  //------------------------------------------------------------
  CLGraph compact() const;

 private:
  mutable std::vector<std::vector<int>> edges_;
  mutable Bitmap is_dirty_;
  mutable std::vector<std::vector<int>> removed_edges_;
  mutable fast_set fset_;  // for general purpose; reset before use
};
}  // namespace graph
}  // namespace data
}  // namespace mog
