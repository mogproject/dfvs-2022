#pragma once
#include <unordered_set>
#include <vector>

#include "../../util/util.hpp"
#include "../bitmap.hpp"

namespace mog {
namespace data {
namespace graph {

class BaseGraph {
 public:
  BaseGraph(std::size_t n = 0) : vertices_(n) {
    vertices_ = ~vertices_;  // [0,n) are valid vertices
  }

  BaseGraph(Bitmap const& vertices) : vertices_(vertices) {}

  BaseGraph(BaseGraph const& g) : vertices_(g.vertices_) {}

  inline void assert_vertex_existence(int v, std::string const& where = "") const {
    if (!has_vertex(v)) throw std::invalid_argument(util::format("%s: invalid vertex (%d)", where.c_str(), v));
  }

  //------------------------------------------------------------
  //   Properties
  //------------------------------------------------------------

  int number_of_nodes() const { return vertices_.count(); }

  long long number_of_edges() const {
    long long ret = 0;
    for (auto v : vertices_.to_vector()) ret += degree(v);
    return ret / 2;
  }

  bool empty() const { return vertices_.empty(); }

  int capacity() const { return vertices_.size(); }

  bool is_compact() const { return number_of_nodes() == capacity(); }

  // O(n)
  int number_of_isolates() const {
    return util::count_if(vertices(), [&](int v) { return is_isolate(v); });
  }

  // O(n) w/ early return
  bool has_any_edge() const {
    return !util::all_of(vertices(), [&](int v) { return is_isolate(v); });
  }

  //------------------------------------------------------------
  //   Enumeration
  //------------------------------------------------------------

  std::vector<int> vertices() const { return vertices_.to_vector(); }

  //------------------------------------------------------------
  //   Modification
  //------------------------------------------------------------
  void add_vertex(int v) {
    if (v < 0 || capacity() <= v) throw std::invalid_argument(util::format("v=%d: out of range", v));
    if (has_vertex(v)) throw std::invalid_argument(util::format("%d: vertex already exists", v));

    vertices_ |= v;
  }

  void remove_vertex(int v) { vertices_ -= v; }

  //------------------------------------------------------------
  //   Queries
  //------------------------------------------------------------

  // O(1)
  bool has_vertex(int v) const { return 0 <= v && v < static_cast<int>(vertices_.size()) && vertices_[v]; }

  virtual int degree(int v) const = 0;

  bool is_isolate(int v) const { return degree(v) == 0; };

  virtual bool has_edge(int u, int v) const = 0;

  //------------------------------------------------------------
  //   Debugging
  //------------------------------------------------------------
  std::string to_matrix_string() const {
    if (capacity() > 100) return "too large";

    std::stringstream ss;
    for (int i = 0; i < capacity(); ++i) {
      for (int j = 0; j < capacity(); ++j) {
        // 0: non-edge, 1: edge, 2: removed vertex
        ss << ((vertices_[i] && vertices_[j]) ? has_edge(i, j) ? "1" : "0" : "2");
      }
      ss << std::endl;
    }
    return ss.str();
  }

 protected:
  Bitmap vertices_;
};
}  // namespace graph
}  // namespace data
}  // namespace mog