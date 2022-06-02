#pragma once
#include <unordered_set>
#include <vector>

#include "../../util/util.hpp"
#include "../fast_set.hpp"
#include "BaseDiGraph.hpp"
#include "CLGraph.hpp"

#define CLDIGRAPH_OPT 1  // optimize for edge removal

#if CLDIGRAPH_OPT == 1
#define CLEAN(v) clean_(v)
#else
#define CLEAN(v)
#endif

namespace mog {
namespace data {
namespace graph {

/**
 * @brief Simple digraph representation with adjacency lists. Vertices are labeled by 0-indexed integers.
 *
 * Space efficiency: O(n+m)
 */
class CLDiGraph : public BaseDiGraph {
 public:
  typedef int Vertex;
  typedef std::vector<Vertex> Vertices;
  typedef std::pair<Vertex, Vertex> Edge;
  typedef std::vector<Edge> Edges;

  CLDiGraph(std::size_t n = 0, std::vector<std::pair<int, int>> const& edges = {})
      : BaseDiGraph(n),
        has_loop_(n),
        fset_(n),
        in_edges_(n),
        out_edges_(n)
#if CLDIGRAPH_OPT == 1
        ,
        is_dirty_(n),
        removed_in_edges_(n),
        removed_out_edges_(n)
#endif
  {
    for (auto& p : edges) add_edge_unsafe(p.first, p.second);
  }

  CLDiGraph(Bitmap const& vertices, std::vector<std::pair<int, int>> const& edges = {})
      : BaseDiGraph(vertices),
        has_loop_(vertices.size()),
        fset_(vertices.size()),
        in_edges_(vertices.size()),
        out_edges_(vertices.size())
#if CLDIGRAPH_OPT == 1
        ,
        is_dirty_(vertices.size()),
        removed_in_edges_(vertices.size()),
        removed_out_edges_(vertices.size())
#endif
  {
    for (auto& p : edges) add_edge_unsafe(p.first, p.second);
  }

  // Copy constructor
  CLDiGraph(CLDiGraph const& g)
      : BaseDiGraph(g),
        has_loop_(g.has_loop_),
        fset_(g.capacity()),
        in_edges_(g.in_edges_),
        out_edges_(g.out_edges_)
#if CLDIGRAPH_OPT == 1
        ,
        is_dirty_(g.is_dirty_),
        removed_in_edges_(g.removed_in_edges_),
        removed_out_edges_(g.removed_out_edges_)
#endif
  {
  }

  // TODO: had error when I defined this...
  //  Assignment operator
  //   CLDiGraph& operator=(CLDiGraph const& g) {
  //     if (this != &g) {
  // #if CLDIGRAPH_OPT == 1
  //       is_dirty_ = g.is_dirty_;
  //       removed_in_edges_ = g.removed_in_edges_;
  //       removed_out_edges_ = g.removed_out_edges_;
  // #endif
  //       fset_.resize(g.capacity());
  //       in_edges_ = g.in_edges_;
  //       out_edges_ = g.out_edges_;
  //       has_loop_ = g.has_loop_;
  //     }
  //     return *this;
  //   }

  //------------------------------------------------------------
  //   Properties
  //------------------------------------------------------------
  static constexpr bool const directed = true;
  static constexpr bool const use_matrix = false;

  //------------------------------------------------------------
  //   Enumeration
  //------------------------------------------------------------
  std::vector<int> const& in_neighbors(int v) const;
  std::vector<int> const& out_neighbors(int v) const;
  std::vector<int> neighbors(int v) const;
  std::vector<int> in_neighbors_sorted(int v) const;
  std::vector<int> out_neighbors_sorted(int v) const;
  std::vector<int> strong_neighbors(int v) const;
  std::vector<std::pair<int, int>> edges(bool use_in_neighbors = false) const;

  //------------------------------------------------------------
  //   Modification
  //------------------------------------------------------------
  void add_vertex(int v);
  void add_edge(int u, int v);
  void add_edge_unsafe(int u, int v);
  void remove_edge(int u, int v);
  void remove_vertex(int v);
  // Edges remove_vertex_with_trace(int v);
  void contract_edge(int u, int v, bool keep_head = true);
  // std::pair<Edges, Edges> contract_edge_with_trace(int u, int v, bool keep_head = true);  // (removed edges, added edges)
  void ignore_vertex(int v);
  // std::pair<Edges, Edges> ignore_vertex_with_trace(int v);  // (removed edges, added edges)

  void clean_all() const;
  void clear();

  void resize(int n);

  //------------------------------------------------------------
  //   Queries
  //------------------------------------------------------------
  bool has_edge(int u, int v) const;
  int in_degree(int v) const;
  int out_degree(int v) const;
  bool is_strong_vertex(int v) const;
  int strong_degree(int v) const;

  //------------------------------------------------------------
  //   Conversions
  //------------------------------------------------------------
  CLGraph strong_graph() const;
  CLDiGraph weak_graph() const;
  CLDiGraph reversal_graph() const;
  CLDiGraph compact() const;
  CLGraph to_undirected() const;

  //------------------------------------------------------------
  //   Debugging
  //------------------------------------------------------------
  void check_consistency() const;
  void clean_(int v) const;

 private:
  Bitmap has_loop_;
  mutable fast_set fset_;  // for general purpose; reset before use

#if CLDIGRAPH_OPT == 1
  mutable std::vector<std::vector<int>> in_edges_;
  mutable std::vector<std::vector<int>> out_edges_;
  mutable Bitmap is_dirty_;
  mutable std::vector<std::vector<int>> removed_in_edges_;
  mutable std::vector<std::vector<int>> removed_out_edges_;

#else
  std::vector<std::vector<int>> in_edges_;
  std::vector<std::vector<int>> out_edges_;
#endif
};
}  // namespace graph
}  // namespace data
}  // namespace mog
