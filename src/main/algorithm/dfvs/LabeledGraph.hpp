#pragma once
#include "../../data/graph/CLDiGraph.hpp"
#include "../../data/graph/CLGraph.hpp"
#include "../component.hpp"

#define VERIFY_LABELED_GRAPH(x)
#define VERIFY_LABELED_GRAPH2(x)

namespace mog {
namespace algorithm {
namespace dfvs {

class LabeledGraph {
 public:
  struct History {
    data::graph::CLDiGraph::Vertices removed_vertices;  // stores vertex labels, not ids
    data::graph::CLDiGraph::Edges removed_edges;
    data::graph::CLDiGraph::Edges added_edges;
    int num_compactions;  // number of compactions

    void clear() {
      removed_vertices.clear();
      removed_edges.clear();
      added_edges.clear();
      num_compactions = 0;
    }
  };

  // compaction brings a new generation
  struct Labels {
    std::vector<int> labels;        // vertex-id -> label
    std::unordered_map<int, int> labels_inv;  // label -> vertex-id

    Labels(int n) : labels(util::range(n)), labels_inv(util::inverse(util::range(n))) {}

    Labels(std::vector<int> const& labels) : labels(labels), labels_inv(util::inverse(labels)) {}

    int get_index(int label) const { return labels_inv.at(label); }

    bool has_label(int label) const { return util::contains(labels_inv, label); }

    int get_label(int index) const { return labels[index]; }
  };

  data::graph::CLDiGraph g;
  data::graph::CLGraph sg;    // strong graph
  data::graph::CLDiGraph wg;  // weak graph

  //  private:
  std::vector<Labels> labels_;
  std::vector<History> history_;
  data::fast_set fset_;  // general-purpose

 public:
  LabeledGraph(data::graph::CLDiGraph const& g)
      : g(g), sg(g.strong_graph()), wg(g.weak_graph()), labels_({Labels(g.capacity())}), history_(1), fset_(g.capacity()) {
    // check_properties("construct 1");
  }

  LabeledGraph(std::vector<int> const& labels)
      : g(labels.size()), sg(labels.size()), wg(labels.size()), labels_({Labels(labels)}), history_(1), fset_(labels.size()) {
    // check_properties("construct 2");
  }

  /**
   * @brief Returns the index of the given label.
   *
   * @param label vertex label
   * @return index
   *
   * Precondition: label must be valid
   */
  int get_index(int label) const { return labels_.back().get_index(label); }

  bool has_label(int label) const { return labels_.back().has_label(label); }

  int get_label(int index) const { return labels_.back().get_label(index); }

  std::vector<int> get_labels() const { return labels_.back().labels; }

  /**
   * @brief Addes an edge to the graph assuming it does not exist currently.
   *
   * @param u predecessor
   * @param v successor
   */
  void add_edge(int u, int v, bool record_history) {
    if (record_history) add_edge_addition_to_history(u, v);
    g.add_edge_unsafe(u, v);
    if (u == v) return;
    if (g.has_edge(v, u)) {
      sg.add_edge(u, v);
      wg.remove_edge(v, u);
      // weak_edge_scores_.erase(encode_edge(v, u));
    } else {
      wg.clean_(u);
      wg.clean_(v);
      wg.add_edge_unsafe(u, v);
    }

    // update edge scores
    // for (auto w: wg.in_neighbors(u)) score_edge(w, u);
    // for (auto w: wg.out_neighbors(u)) score_edge(u, w);
    // for (auto w: wg.in_neighbors(v)) score_edge(w, v);
    // for (auto w: wg.out_neighbors(v)) score_edge(v, w);

    VERIFY_LABELED_GRAPH2(check_properties("add_edge"));
  }

  void add_edge_removal_to_history(int u, int v) {
    if (u == v) return;  // do not record loops
    auto& hist = history_.back();
    auto e = std::make_pair(get_label(u), get_label(v));
    if (util::contains(hist.added_edges, e)) {
      util::erase_first(hist.added_edges, e);
    } else {
      hist.removed_edges.push_back(e);
    }
  }

  void add_edge_addition_to_history(int u, int v) {
    if (u == v) return;  // do not record loops
    auto& hist = history_.back();
    auto e = std::make_pair(get_label(u), get_label(v));
    if (util::contains(hist.removed_edges, e)) {
      util::erase_first(hist.removed_edges, e);
    } else {
      hist.added_edges.push_back(e);
    }
  }

  /**
   * @brief Remove a vertex and its incident edges with logging.
   *
   * @param v vertex to remove
   */
  void remove_vertex(int v, bool record_history) {
    if (record_history) {
      history_.back().removed_vertices.push_back(get_label(v));
      for (auto u : g.out_neighbors(v)) add_edge_removal_to_history(v, u);
      for (auto w : g.in_neighbors(v)) add_edge_removal_to_history(w, v);
      // if (g.has_edge(v, v)) add_edge_removal_to_history(v, v);
    }

    g.remove_vertex(v);
    sg.remove_vertex(v);

    // update edge scores
    // std::vector<std::pair<int, int>> affected_edges;
    // for (auto u: wg.out_neighbors(v)) {
    //   for (auto w: wg.out_neighbors(u)) affected_edges.push_back({u, w});
    //   for (auto w: wg.in_neighbors(u)) {
    //     if (w != v) affected_edges.push_back({w, u});
    //   }
    //   weak_edge_scores_.erase(encode_edge(v, u));
    // }
    // for (auto u: wg.in_neighbors(v)) {
    //   for (auto w: wg.out_neighbors(u)) {
    //     if (w != v) affected_edges.push_back({u, w});
    //   }
    //   for (auto w: wg.in_neighbors(u)) affected_edges.push_back({w, u});
    //   weak_edge_scores_.erase(encode_edge(u, v));
    // }

    wg.remove_vertex(v);

    VERIFY_LABELED_GRAPH2(check_properties("remove_vertex"));
    // for (auto &e: affected_edges) score_edge(e.first, e.second);
  }

  void remove_edge(int u, int v, bool record_history) {
    if (record_history) add_edge_removal_to_history(u, v);

    g.remove_edge(u, v);
    if (u == v) return;

    if (g.has_edge(v, u)) {
      // uv was a strong edge
      sg.remove_edge(u, v);
      wg.add_edge_unsafe(v, u);
    } else {
      wg.remove_edge(u, v);
      // weak_edge_scores_.erase(encode_edge(u, v));
    }
    VERIFY_LABELED_GRAPH2(check_properties("remove_edge"));

    // update edge scores
    // for (auto w: wg.in_neighbors(u)) score_edge(w, u);
    // for (auto w: wg.out_neighbors(u)) score_edge(u, w);
    // for (auto w: wg.in_neighbors(v)) score_edge(w, v);
    // for (auto w: wg.out_neighbors(v)) score_edge(v, w);
  }

  void contract_edge(int u, int v, bool keep_head) {
    if (u == v) return;  // do nothing

    int to_keep = keep_head ? u : v;
    int to_remove = keep_head ? v : u;

    for (auto x : g.in_neighbors(to_remove)) {
      if (x == to_keep) {
        if (!keep_head && !g.has_edge(x, x)) g.add_edge_unsafe(x, x);
      } else if (!g.has_edge(x, to_keep)) {
        // PRINT("add_edge: " << (get_label(x)+1) << " " << (get_label(to_keep)+1));
        add_edge(x, to_keep, true);
      }
    }
    for (auto y : g.out_neighbors(to_remove)) {
      if (y == to_keep) {
        if (keep_head && !g.has_edge(y, y)) g.add_edge_unsafe(y, y);
      } else if (!g.has_edge(to_keep, y)) {
        // PRINT("add_edge: " << (get_label(to_keep)+1) << " " << (get_label(y)+1));
        add_edge(to_keep, y, true);
      }
    }
    if (g.has_edge(to_remove, to_remove) && !g.has_edge(to_keep, to_keep)) g.add_edge_unsafe(to_keep, to_keep);

    // PRINT("remove_vertex: " << (get_label(to_remove)+1));
    remove_vertex(to_remove, true);
    VERIFY_LABELED_GRAPH2(check_properties("contract_edge"));
  }

  void ignore_vertex(int v, bool record_history) {
    for (auto u : g.in_neighbors(v)) {
      for (auto w : g.out_neighbors(v)) {
        if (!g.has_edge(u, w)) add_edge(u, w, record_history);
      }
    }
    remove_vertex(v, record_history);
    VERIFY_LABELED_GRAPH2(check_properties("ignore_vertex"));
  }

  void remove_vertices(std::vector<int> vs, bool record_history) {
    sg.clean_all();  // why necessary?

    for (auto v : vs) remove_vertex(v, record_history);
    VERIFY_LABELED_GRAPH2(check_properties("remove_vertices"));
  }

  std::vector<int> ignore_and_reduce(int v, bool record_history) {
    // (1) Find the strong (bidirected) neighbors of v.
    auto snbrs = sg.neighbors(v);

    // (2) Remove the strong neighbors.
    remove_vertices(snbrs, record_history);

    // (3) Ignore this vertex. Then, there should not be a loop.
    ignore_vertex(v, record_history);

    // (4) Return the strong neighbors to be included in solution.
    VERIFY_LABELED_GRAPH2(check_properties("ignore_and_reduce"));
    return snbrs;
  }

  std::vector<int> ignore_and_reduce2(int u, int v, bool record_history) {
    // Assumption: uv is a weak edge
    std::vector<int> ret;

    // (1) Find directed WEAK triangles including u and v.
    fset_.clear();
    for (auto x : wg.in_neighbors(u)) fset_.set(x);
    for (auto x : wg.out_neighbors(v)) {
      if (fset_.get(x)) {
        // x-u-v forms a triangle
        ret.push_back(x);
      }
    }

    // (2) Find STRONG neighbors of u and v.
    fset_.clear();
    for (auto x : sg.neighbors(u)) {
      fset_.set(x);
      ret.push_back(x);
    }
    for (auto x : sg.neighbors(v)) {
      if (!fset_.get(x)) ret.push_back(x);  // avoid duplicate entries
    }

    // (3) Remove the reducible neighbors.
    remove_vertices(ret, record_history);

    // (4) Ignore two vertices. There should be no loops.
    for (auto x : g.in_neighbors(u)) {
      for (auto y : g.out_neighbors(u)) {
        if (!g.has_edge(x, y)) add_edge(x, y, record_history);  // cycles x-u-y-...-x
      }
      for (auto y : g.out_neighbors(v)) {
        if (!g.has_edge(x, y)) add_edge(x, y, record_history);  // cycles x-u-v-y-...-x
      }
    }
    for (auto x : g.in_neighbors(v)) {
      for (auto y : g.out_neighbors(v)) {
        if (!g.has_edge(x, y)) add_edge(x, y, record_history);  // cycles x-v-y-...-x
      }
    }
    remove_vertex(u, record_history);
    remove_vertex(v, record_history);

    // (5) Return the removed vertices to be included in solution.
    VERIFY_LABELED_GRAPH2(check_properties("ignore_and_reduce2"));
    return ret;
  }

  void compact() {
    VERIFY_LABELED_GRAPH(check_properties("compact start"));
    int n = g.number_of_nodes();
    int nn = g.capacity();
    if (n == nn) return;  // already compact

    // data::fast_set fset(nn);
    // for (auto v : g.vertices()) fset.set(v);

    // create a new set of labels

    // std::vector<int> new_labels = g.vertices().to_vector();
    // std::vector<int> removed_labels;
    // std::map<int, int> new_labels_inv;

    std::vector<int> new_label_entries;
    for (auto v : g.vertices()) new_label_entries.push_back(get_label(v));
    auto new_labels = Labels(new_label_entries);
    // for (int i = 0; i < static_cast<int>(labels_.size()); ++i) {
    //   if (i < nn && fset.get(i)) {
    //     new_labels_inv[get_label(i)] = new_labels.size();
    //     new_labels.push_back(get_label(i));
    //   } else {
    //     new_labels_inv[get_label(i)] = n + static_cast<int>(removed_labels.size());
    //     removed_labels.push_back(get_label(i));
    //   }
    // }

    // create new graph instances
    data::graph::CLDiGraph new_g(n);
    data::graph::CLGraph new_sg(n);
    data::graph::CLDiGraph new_wg(n);

    for (auto v : g.vertices()) {
      int i = new_labels.get_index(get_label(v));
      for (auto u : g.out_neighbors(v)) {
        new_g.add_edge_unsafe(i, new_labels.get_index(get_label(u)));  //
      }
      for (auto u : sg.neighbors(v)) {
        if (v < u) {
          if (!new_sg.has_vertex(i) || !new_sg.has_vertex(new_labels.get_index(get_label(u)))) {
            PRINT("n=" << n << ", nn=" << nn << ", v=" << v << ", u=" << u);
            PRINT("vertices=" << sg.vertices());
            PRINT(get_labels());
            PRINT(new_labels.labels);
          }
          new_sg.add_edge(i, new_labels.get_index(get_label(u)));
        }
      }
      for (auto u : wg.out_neighbors(v)) {
        new_wg.add_edge_unsafe(i, new_labels.get_index(get_label(u)));  //
      }
    }

    // replace objects
    g = new_g;
    sg = new_sg;
    wg = new_wg;
    labels_.push_back(new_labels);  // add a new layer of labels
    ++history_.back().num_compactions;

    // score edges
    // score_edge_all();
    // debug
    VERIFY_LABELED_GRAPH(check_properties("compact end"));
  }

  void clear() { remove_vertices(g.vertices(), true); }

  void commit() { history_.push_back({}); }

  /**
   * @brief Rolls back to the last commit point.
   *
   * Label mapping may be different from the original one.
   */
  void rollback() {
    VERIFY_LABELED_GRAPH(check_properties("rollback start"));
    auto& h = history_.back();

    // restore labels
    for (int t = 0; t < h.num_compactions; ++t) {
      auto& pre = labels_[labels_.size() - 2];
      auto& cur = labels_.back();
      for (int i = 0; i < static_cast<int>(cur.labels.size()); ++i) {
        // relocate pre.labels[i]
        auto v = cur.labels[i];
        auto j = pre.labels_inv[v];  // v's original position
        if (i == j) continue;

        auto u = pre.get_label(i);
        pre.labels[i] = v;
        pre.labels[j] = u;
        pre.labels_inv[v] = i;
        pre.labels_inv[u] = j;
      }
      labels_.pop_back();
    }

    // printf("rollback(): removed_vertices: %s\n", cstr(h.removed_vertices));

    // restore removed vertices
    int n = labels_.back().labels.size();
    if (n > g.capacity()) {  // increase the capacity if needed
      g.resize(n);
      sg.resize(n);
      wg.resize(n);
    }

    // PRINT("labels: " << get_labels());
    // PRINT("removed_vertices: " << h.removed_vertices);
    // PRINT("vertices: " << g.vertices());

    for (auto v : h.removed_vertices) {
      g.add_vertex(get_index(v));
      sg.add_vertex(get_index(v));
      wg.add_vertex(get_index(v));
    }

    // restore removed edges
    for (auto& e : h.removed_edges) {
      auto u = get_index(e.first);
      auto v = get_index(e.second);
      if (g.has_edge(u, v)) {
        PRINT("removed edges: " << h.removed_edges);
        throw std::runtime_error("edge exists");
      }
      add_edge(u, v, false);
    }
    // check_properties("rollback middle");

    // remove added edges
    wg.clean_all();  // necessary here

    for (auto& e : h.added_edges) {
      auto u = get_index(e.first);
      auto v = get_index(e.second);
      if (!g.has_edge(u, v)) {
        PRINT("added edges: " << h.added_edges);
        throw std::runtime_error("edge does not exist");
      }
      remove_edge(u, v, false);
    }

    // edge scoring
    // score_edge_all();

    history_.pop_back();
    if (history_.empty()) history_.push_back({});

    // debug
    VERIFY_LABELED_GRAPH(check_properties("rollback end"));

    // debug
    // for (auto v: g.vertices()) {
    //   for (auto u: g.out_neighbors(v)) {
    //     if (!g.has_vertex(u)) {
    //       PRINT("nn=" << g.capacity() << ", n=" << g.number_of_nodes() << ", v=" << v << ", u=" << u);
    //       PRINT("vs=" << g.vertices());
    //       PRINT("N[v]=" << g.out_neighbors(v));
    //       throw std::runtime_error("error 9001");
    //     }
    //   }
    //   for (auto u: g.in_neighbors(v)) {
    //     if (!g.has_vertex(u)) throw std::runtime_error("error 9002");
    //   }
    //   for (auto u: sg.neighbors(v)) {
    //     if (!sg.has_vertex(u)) throw std::runtime_error("error 9003");
    //   }
    // }
  }

  //------------------------------------------------------------
  //   Edge Ordering
  //------------------------------------------------------------

  // mutable std::priority_queue<std::pair<double, std::pair<int, int>>> edge_ordering_;
  // mutable std::unordered_map<long long, double> weak_edge_scores_;
  
  // long long encode_edge(int u, int v) const { return (static_cast<long long>(u) << 32) + v; }

  // void score_edge(int u, int v) {
  //   double score = wg.degree(u) + wg.degree(v);
  //   weak_edge_scores_[encode_edge(u, v)] = score;
  //   edge_ordering_.push({score, {u, v}});
  // }

  // void score_edge_all() {
  //   weak_edge_scores_.clear();
  //   for (auto u: wg.vertices()) {
  //     for (auto v: wg.out_neighbors(u)) score_edge(u, v);
  //   }
  // }

  // std::pair<int, int> get_best_edge() const {
  //   while (!edge_ordering_.empty()) {
  //     auto &p = edge_ordering_.top();
  //     auto enc = encode_edge(p.second.first, p.second.second);
  //     if (util::contains(weak_edge_scores_, enc) && weak_edge_scores_.at(enc) == p.first) return p.second;
  //     edge_ordering_.pop();  // outdated info
  //   }
  //   throw std::invalid_argument("no edge candidate");
  //   return {-1, -1};
  // }

  //------------------------------------------------------------
  //   Debugging
  //------------------------------------------------------------

  void check_properties(std::string const& where) const {
    // base graph
    if (util::sorted(g.edges(false)) != util::sorted(g.edges(true))) {
      PRINT("error in: " << where);
      PRINT("g.edges(false): " << util::sorted(g.edges(false)));
      PRINT("g.edges(true): " << util::sorted(g.edges(true)));
      throw std::runtime_error("inconsistent adjacency list: g");
    }

    // weak graph
    if (util::sorted(wg.edges(false)) != util::sorted(wg.edges(true))) {
      PRINT("error in: " << where);
      PRINT("wg.edges(false): " << util::sorted(wg.edges(false)));
      PRINT("wg.edges(true): " << util::sorted(wg.edges(true)));
      throw std::runtime_error("inconsistent adjacency list: wg");
    }

    // strong graph
    for (auto v : sg.vertices()) {
      for (auto u : sg.neighbors(v)) {
        if (!sg.has_edge(u, v)) {
          PRINT("error in: " << where);
          throw std::runtime_error("not symmetric");
        }
        if (v == u) {
          PRINT("error in: " << where);
          PRINT("v=" << v << ", u=" << u);
          throw std::runtime_error("found loop");
        }
        if (!sg.has_vertex(u)) {
          PRINT("error in: " << where);
          PRINT("v=" << v << ", u=" << u);
          throw std::runtime_error("invalid vertex");
        }
      }
    }

    auto strong_edges1 = g.strong_graph().edges();
    auto strong_edges2 = sg.edges();
    auto weak_edges1 = g.weak_graph().edges();
    auto weak_edges2 = wg.edges();

    std::sort(strong_edges1.begin(), strong_edges1.end());
    std::sort(strong_edges2.begin(), strong_edges2.end());
    std::sort(weak_edges1.begin(), weak_edges1.end());
    std::sort(weak_edges2.begin(), weak_edges2.end());

    if (strong_edges1 != strong_edges2) {
      PRINT("error in: " << where);
      PRINT("m1=" << strong_edges1.size() << ", m2=" << strong_edges2.size());
      PRINT("g.strong_graph().edges()=" << strong_edges1);
      PRINT("sg              .edges()=" << strong_edges2);
      throw std::runtime_error("broken strong graph");
    }
    if (weak_edges1 != weak_edges2) {
      PRINT("error in: " << where);
      PRINT("m1=" << weak_edges1.size() << ", m2=" << weak_edges2.size());
      PRINT("g.weak_graph().edges()=" << weak_edges1);
      PRINT("wg            .edges()=" << weak_edges2);
      throw std::runtime_error("broken weak graph");
    }
  }
};

}  // namespace dfvs
}  // namespace algorithm
}  // namespace mog
