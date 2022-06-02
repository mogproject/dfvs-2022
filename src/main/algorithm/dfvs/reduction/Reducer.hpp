#pragma once

#include <cassert>

#include "../Configuration.hpp"
#include "../LabeledGraph.hpp"

// #define TRACE_REDUCER(x) PRINT(x)
#define TRACE_REDUCER(x)

namespace mog {
namespace algorithm {
namespace dfvs {
class Reducer {
 private:
  Configuration conf_;
  util::Profiler* const prof_;
  std::vector<int> part_sol_;  // partial solution (list of vertex labels, not ids)
  int compaction_min_ = 20;

 public:
  Reducer(Configuration const& conf) : conf_(conf), prof_(conf.prof) {}

  std::vector<int> reduce(std::vector<LabeledGraph>& gs, int index, bool enforce_compaction = false) {
    part_sol_.clear();

    // the index of the currently handling subproblem (-1: given graph)
    auto offset = static_cast<int>(gs.size());
    index = index < 0 ? offset - 1 : index;
    auto cur_index = index;

    // gs[index].check_properties("reduce(): 1");

    // first decomposition
    decompose(gs, cur_index);
    /** @note now, gs[index], gs[offset], gs[offset+1]... are SCCs */
    // gs[index].check_properties("reduce(): 2");

    // main loop
    while (cur_index < static_cast<int>(gs.size())) {
      auto& h = gs[cur_index];
      bool updated = false;

      //------------------------------------------------------------
      //   (1) Compress one SCC
      //------------------------------------------------------------
      updated |= compress(h);

      if (h.g.empty()) goto erase_this;  // empty after compression
      /** @note now, h.g must be loop-less and have at least 3 vertices */

      if (conf_.reduce.complete_reduction) {
        //------------------------------------------------------------
        //   (2) Complete SCC reduction (special case of CORE reduction)
        //------------------------------------------------------------
        if (reduce_complete(h)) goto erase_this;

        //------------------------------------------------------------
        //   (3) n=4 reduction (assuming not complete)
        //------------------------------------------------------------
        if (conf_.reduce.n4_reduction && reduce_n4(h)) goto erase_this;
      }
      // gs[index].check_properties("reduce(): 3");

      //------------------------------------------------------------
      //   (4) PIE reduction
      //------------------------------------------------------------
      if (conf_.reduce.pie_reduction) {
        bool ret_reduce_pie = reduce_pi(h);
        PROF(if (prof_) prof_->count("reduce_pie()", ret_reduce_pie ? 1 : 0));
        updated |= ret_reduce_pie;
      }
      // gs[index].check_properties("reduce(): 4");

      //------------------------------------------------------------
      //   (5) CORE reduction
      //------------------------------------------------------------
      if (conf_.reduce.core_reduction) {
        bool ret_reduce_core = reduce_core(h);
        PROF(if (prof_) prof_->count("reduce_core()", ret_reduce_core ? 1 : 0));
        updated |= ret_reduce_core;
      }
      // gs[index].check_properties("reduce(): 5");

      //------------------------------------------------------------
      //   (6) DOME reduction
      //------------------------------------------------------------
      if (conf_.reduce.dome_reduction) {
        bool ret_reduce_dome = reduce_dome(h);
        PROF(if (prof_) prof_->count("reduce_dome()", ret_reduce_dome ? 1 : 0));
        updated |= ret_reduce_dome;
      }
      // gs[index].check_properties("reduce(): 6");

      //------------------------------------------------------------
      //   (7) Decompose into SCCs
      //------------------------------------------------------------
      if (!updated) goto move_next;  // done with this subproblem

      updated |= decompose(gs, cur_index);
      if (!updated) goto move_next;                  // decomposition did not happen
      if (gs[cur_index].g.empty()) goto erase_this;  // SCCs are all singletons

      continue;

    erase_this:
      if (cur_index == index) {
        gs[cur_index].clear();
        cur_index = offset;
      } else {
        util::erase(gs, cur_index);
      }
      continue;

    move_next:
      cur_index = cur_index == index ? offset : cur_index + 1;
    }

    // gs[index].check_properties("reduce(): 7");

    // PRINT(part_sol_);
    if (enforce_compaction || should_compact(gs[index])) {
      PROF(if (prof_) prof_->count("compact()"));
      PROF(if (prof_) prof_->start_timer("compact()"));
      gs[index].compact();
      PROF(if (prof_) prof_->stop_timer("compact()"));
    }
    return part_sol_;
  }

  void disable_compaction() { conf_.reduce.compact_threshold = 0; }

  bool should_compact(LabeledGraph const& h) const {
    if (h.g.empty()) return false;
    int nn = h.g.capacity();
    if (nn <= compaction_min_) return false;  // small enough capacity
    int n = h.g.number_of_nodes();
    return n <= nn * conf_.reduce.compact_threshold;
  }

 private:
  /**
   * @brief Naive but good enough compression (i.e. degree-based reductions).
   *
   * @param h
   * @return true
   * @return false
   */
  bool compress(LabeledGraph& h) {
    PROF(if (prof_) prof_->start_timer("compress()"));
    bool updated = true, ret = false;
    while (updated) {
      updated = false;
      for (auto v : h.g.vertices()) { updated |= reduce_loop(h, v) || reduce_deg_0(h, v) || reduce_deg_1(h, v); }
      ret |= updated;
    }
    PROF(if (prof_) prof_->stop_timer("compress()"));
    return ret;
  }

  /**
   * @brief Loop reduction
   *
   * @param h labeled graph
   * @return true reduced
   * @return false did not reduce
   */
  bool reduce_loop(LabeledGraph& h, int v) {
    if (h.g.has_edge(v, v)) {  // found a loop
      // h.g.remove_edge(v, v);  // I don't think this is necesary
      h.remove_vertex(v, true);
      part_sol_.push_back(h.get_label(v));
      TRACE_REDUCER("reduce loop: " << v);

      // PROF(if (prof_) prof_->count("reduce_loop()=true"));
      return true;
    }
    return false;
  }

  /**
   * @brief Degree-0 reduction.
   *
   * @param h labeled digraph
   * @return true reduced
   * @return false did not reduce
   */
  bool reduce_deg_0(LabeledGraph& h, int v) {
    if (h.g.in_degree(v) == 0 || h.g.out_degree(v) == 0) {
      TRACE_REDUCER(util::format("reduce_deg_0(): remove vertex (%d)", h.get_label(v)));
      h.remove_vertex(v, true);
      // PROF(if (prof_) prof_->count("reduce_deg_0()=true"));
      return true;
    }
    return false;
  }

  /**
   * @brief Degree-1 reduction.
   *
   * @param h labeled digraph
   * @return true reduced
   * @return false did not reduce
   */
  bool reduce_deg_1(LabeledGraph& h, int v) {
    bool ret = h.g.in_degree(v) == 1 || h.g.out_degree(v) == 1;
    if (h.g.in_degree(v) == 1) {
      TRACE_REDUCER(util::format("reduce_deg_1(): contract edge (%d, %d) keep-head", h.get_label(h.g.in_neighbors(v).front()), h.get_label(v)));
      h.contract_edge(h.g.in_neighbors(v).front(), v, true);
    } else if (h.g.out_degree(v) == 1) {
      TRACE_REDUCER(util::format("reduce_deg_1(): contract edge (%d, %d) keep-tail", h.get_label(v), h.get_label(h.g.out_neighbors(v).front())));
      h.contract_edge(v, h.g.out_neighbors(v).front(), false);
    }
    // if (ret) PROF(if (prof_) prof_->count("reduce_deg_1()=true"));
    return ret;
  }

  bool reduce_pi(LabeledGraph& h) {
    PROF(if (prof_) prof_->start_timer("reduce_pi()"));
    bool ret = false;

    auto& wg = h.wg;

    auto scc_id = scc_ids(wg);
    for (int v : wg.vertices()) {
      for (int u : wg.out_neighbors(v)) {
        if (scc_id[v] == -1 || scc_id[v] != scc_id[u]) {
          TRACE_REDUCER(util::format("reduce_pi(): remove edge (%d, %d)", h.get_label(v), h.get_label(u)));

          // h.check_properties("pi before");
          h.remove_edge(v, u, true);
          // h.check_properties("pi");
        }
      }
    }

    PROF(if (prof_) prof_->stop_timer("reduce_pi()"));
    return ret;
  }

  bool reduce_complete(LabeledGraph& h) {
    PROF(if (prof_) prof_->start_timer("reduce_complete()"));
    int n = h.g.number_of_nodes();
    if (h.g.number_of_edges() != n * (n - 1)) {
      PROF(if (prof_) prof_->stop_timer("reduce_complete()"));
      return false;
    }

    // pick any (n-1) vertices
    int i = 0;
    for (auto v : h.g.vertices()) {
      part_sol_.push_back(h.get_label(v));
      if (++i == n - 1) break;
    }
    PROF(if (prof_) prof_->stop_timer("reduce_complete()"));
    return true;
  }

  bool reduce_n4(LabeledGraph& h) {
    PROF(if (prof_) prof_->start_timer("reduce_n4()"));
    if (h.g.number_of_nodes() != 4) {
      PROF(if (prof_) prof_->stop_timer("reduce_n4()"));
      return false;
    }

    auto vs = h.g.vertices();

    // there must be a pair with no cycle as it is not complete
    for (int i = 0; i < 3; ++i) {
      for (int j = i + 1; j < 4; ++j) {
        if (!h.sg.has_edge(vs[i], vs[j])) {
          for (int k = 0; k < 4; ++k) {
            if (k != i && k != j) part_sol_.push_back(h.get_label(vs[k]));
          }
          PROF(if (prof_) prof_->stop_timer("reduce_n4()"));
          return true;
        }
      }
    }
    throw std::runtime_error("never happens");
    return false;
  }

  bool reduce_core(LabeledGraph& h) {
    PROF(if (prof_) prof_->start_timer("reduce_core()"));
    data::fast_set invalid(h.g.capacity());
    data::fast_set v_nbrs(h.g.capacity());
    std::vector<std::pair<int, int>> s_vtcs;

    auto& sg = h.sg;

    for (auto v : h.g.vertices()) {
      if (h.g.out_degree(v) == sg.degree(v)) s_vtcs.push_back({sg.degree(v), v});
    }
    std::sort(s_vtcs.begin(), s_vtcs.end());  // sort strong vertices by degree

    for (auto& p : s_vtcs) {
      auto v = p.second;
      if (invalid.get(v)) continue;  // already invalidated

      bool v_in_clique = true;
      v_nbrs.clear();
      for (auto u : sg.neighbors(v)) v_nbrs.set(u);  // v_nbrs = N(v)

      for (auto u : sg.neighbors(v)) {
        // count |N(u) intersect N(v)|
        int cnt = 0;
        for (auto w : sg.neighbors(u)) {
          if (v_nbrs.get(w)) ++cnt;
        }
        if (cnt != sg.degree(v) - 1) {
          v_in_clique = false;
          break;
        }
      }

      if (v_in_clique) {
        // include N(v) in solution
        for (auto u : sg.neighbors(v)) part_sol_.push_back(h.get_label(u));

        // remove N[v] from the graph
        h.remove_vertices(sg.neighbors(v), true);
        h.remove_vertex(v, true);

        PROF(if (prof_) prof_->stop_timer("reduce_core()"));
        return true;  // immediately return
      } else {
        // invalidate N(v)
        for (auto u : sg.neighbors(v)) invalid.set(u);
      }
    }
    PROF(if (prof_) prof_->stop_timer("reduce_core()"));
    return false;
  }

  bool reduce_dome(LabeledGraph& h) {
    PROF(if (prof_) prof_->start_timer("reduce_dome()"));
    bool ret = false;
    data::fast_set fset(h.g.capacity());

    for (auto u : h.g.vertices()) {
      for (auto v : h.wg.out_neighbors(u)) {
        // for every weak edge uv
        bool can_reduce = true;
        fset.clear();
        for (auto w : h.g.in_neighbors(v)) fset.set(w);  // fset := N^-(v)
        for (auto w : h.wg.in_neighbors(u)) {
          if (!fset.get(w)) {
            can_reduce = false;  // N^-weak(u) not in N^-(v)
            break;
          }
        }
        if (can_reduce) {
          // N^-weak(u) subseteq N^-(v)
          ret = true;
          h.remove_edge(u, v, true);
          TRACE_REDUCER(util::format("reduce_dome(): remove edge (%d, %d)", h.get_label(u), h.get_label(v)));
          // h.check_properties("dm");
          continue;
        }

        can_reduce = true;
        fset.clear();
        for (auto w : h.g.out_neighbors(v)) fset.set(w);  // fset := N^+(v)
        for (auto w : h.wg.out_neighbors(u)) {
          if (!fset.get(w)) {
            can_reduce = false;  // N^+weak(u) not in N^+(v)
            break;
          }
        }
        if (can_reduce) {
          // N^+weak(u) subseteq N^+(v)
          ret = true;
          h.remove_edge(u, v, true);
          TRACE_REDUCER(util::format("reduce_dome(): remove edge (%d, %d)", h.get_label(u), h.get_label(v)));
          // h.check_properties("dm2");
        }
      }
    }
    PROF(if (prof_) prof_->stop_timer("reduce_dome()"));
    return ret;
  }

  /**
   * @brief Decomposes a graph into a set of SCCs.
   * `gs[index]` will be replaced by the induced subgraph on the largest SCC.
   * Other SCCs are added to the back of `gs`.
   *
   * @param gs list of labeled digraphs
   * @param index index to process
   * @return true graph was updated
   * @return false graph was not updated
   */
  bool decompose(std::vector<LabeledGraph>& gs, int index) {
    PROF(if (prof_) prof_->start_timer("decompose()"));
    auto sccs = strongly_connected_components(gs[index].g, false, prof_);  // ignore singletons
    if (sccs.empty()) {                                                    // all done
      gs[index].clear();

      PROF(if (prof_) prof_->stop_timer("decompose()"));
      return true;
    }
    if (static_cast<int>(sccs[0].size()) == gs[index].g.number_of_nodes()) {
      PROF(if (prof_) prof_->stop_timer("decompose()"));
      return false;  // no change
    }

    int back_index = gs.size();

    // find the largest SCC
    int largest_scc = util::max_size_element_index(sccs);

    for (int i = 0; i < static_cast<int>(sccs.size()); ++i) TRACE("Decomposition[%d]: %s", i, cstr(sccs[i]));
    TRACE("Largest SCC: %d", largest_scc);

    // create subgraphs
    std::vector<int> scc_number(gs[index].g.capacity(), -1);  // vertex-index -> scc-id
    for (int i = 0; i < static_cast<int>(sccs.size()); ++i) {
      for (auto j : sccs[i]) scc_number[j] = i;
      if (i == largest_scc) continue;  // reuse `this` for the largest SCC

      std::vector<int> labels;
      for (auto j : sccs[i]) {
        assert(gs[index].g.has_vertex(j));
        labels.push_back(gs[index].get_label(j));
      }
      gs.push_back(LabeledGraph(labels));  // create new graph instance
    }

    // add edges
    std::vector<int> to_remove;
    for (auto i : gs[index].g.vertices()) {
      int scc_index = scc_number[i];
      if (scc_index == largest_scc) continue;
      to_remove.push_back(i);
      if (scc_index == -1) continue;  // singleton

      for (int j : gs[index].g.out_neighbors(i)) {
        if (scc_index == scc_number[j]) {
          auto offset = scc_index < largest_scc ? scc_index : scc_index - 1;
          auto& k = gs[back_index + offset];
          k.add_edge(k.get_index(gs[index].get_label(i)), k.get_index(gs[index].get_label(j)), false);
        }
      }
    }

    // update this graph with logging
    gs[index].remove_vertices(to_remove, true);

    // gs[index].check_properties("decompose() end");
    // for (int i = back_index; i < static_cast<int>(gs.size()); ++i) gs[i].check_properties("decompose() end+");

    PROF(if (prof_) prof_->stop_timer("decompose()"));
    return true;
  }
};

}  // namespace dfvs
}  // namespace algorithm
}  // namespace mog