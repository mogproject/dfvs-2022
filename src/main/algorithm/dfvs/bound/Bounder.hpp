#pragma once

#include <cassert>
#include <random>

#include "../../../data/fast_set.hpp"
#include "../../../data/graph/CLGraph.hpp"
#include "../../clique.hpp"
#include "../../cycle.hpp"
#include "../../vertex_cover.hpp"
#include "../Configuration.hpp"
#include "../LabeledGraph.hpp"
#include "../reduction/Reducer.hpp"

namespace mog {
namespace algorithm {
namespace dfvs {

class Bounder {
 public:
  Configuration conf_;
 private:
  Reducer reducer_;
  util::Profiler* const prof_;
  std::default_random_engine rand_gen_;
  std::uniform_real_distribution<> real_dist_;

 public:
  Bounder(Configuration const& conf)
      : conf_(conf), reducer_(conf), prof_(conf.prof), rand_gen_(conf.seed), real_dist_(0, 1) {}

  int lower_bound(LabeledGraph const& h) {
    int ret = lower_bound_const();
    if (conf_.bound.exec_rate < 1.0 && real_dist_(rand_gen_) >= conf_.bound.exec_rate) return ret;

    switch (conf_.bound.strategy) {
      case BoundStrategy::BoundConstant: ret = lower_bound_const(); break;
      case BoundStrategy::BoundPacking: ret = lower_bound_clique_and_cycle(h, false); break;
      case BoundStrategy::BoundPackingFull: ret = lower_bound_clique_and_cycle(h, true); break;
      case BoundStrategy::BoundVertexCoverFull: ret = lower_bound_vc_full(h); break;
      case BoundStrategy::BoundVertexCoverOnly:
       {
        // ret = lower_bound_vc(h, 1);
        // ret = lower_bound_vc(h, 10);
        // int ret3 = lower_bound_vc_full(h);
        ret = lower_bound_vc_full(h);
        // int ret4 = lower_bound_cycle_packing(h.g, 10);
        // PROF(if (prof_) prof_->count("ret2-ret (+ is good)", ret2 - ret));
        // PROF(if (prof_) prof_->count("ret3-ret (+ is good)", ret3 - ret));
        // PROF(if (prof_) prof_->count("ret4-ret (+ is good)", ret4 - ret));
        // ret = std::max(ret, std::max(ret2, ret3));

        // ret = lower_bound_cycle_packing(h.g, 1);
        break;
      }
      default: throw std::invalid_argument("unknown strategy");
    }
    return ret;
  }

  constexpr int lower_bound_const() const { return 2; }

  int lower_bound_clique_and_cycle(LabeledGraph const& hh, bool full_iteration) {
    int ret = 0;

    // PRINT("V:" << g.vertices());
    // PRINT("E:" << g.edges());

    std::vector<LabeledGraph> gs = {hh};
    gs[0].compact();

    while (!gs.empty()) {
      auto& h = gs.back();
      if (h.g.empty()) {
        gs.pop_back();
        continue;
      }

      if (full_iteration) {
        //------------------------------------------------------------
        // Costly version
        //------------------------------------------------------------
        h.compact();

        // (1) find a large (not neccesarily the largest) clique
        PROF(if (prof_) prof_->start_timer("find_large_clique()"));
        auto clique = find_large_clique(h.sg);
        PROF(if (prof_) prof_->stop_timer("find_large_clique()"));

        // PRINT(clique);

        if (clique.empty()) {
          // (2) find a short cycle
          PROF(if (prof_) prof_->start_timer("short_cycle()"));
          auto cycle = short_cycle(h.g);
          PROF(if (prof_) prof_->stop_timer("short_cycle()"));
          if (cycle.empty()) {
            // no cycle found
          } else {
            for (auto c : cycle) h.remove_vertex(c, false);
            ret += 1;
          }
        } else {
          for (auto c : clique) h.remove_vertex(c, false);
          ret += clique.size() - 1;
        }

        // (3) reduction
        ret += reducer_.reduce(gs, -1).size();
      } else {
        //------------------------------------------------------------
        // Cheap version
        //------------------------------------------------------------
        // find clique partition
        PROF(if (prof_) prof_->start_timer("clique_partition()"));
        auto parts = clique_partition(h.sg);
        PROF(if (prof_) prof_->stop_timer("clique_partition()"));

        for (auto& cs : parts) {
          ret += cs.size() - 1;
          // PROF(if (prof_) prof_->count("clique size", cs.size()));
          for (auto x : cs) h.remove_vertex(x, false);
        }

        while (!h.g.empty() && true) {
          PROF(if (prof_) prof_->start_timer("short_cycle()"));
          auto cycle = short_cycle(h.g);
          PROF(if (prof_) prof_->stop_timer("short_cycle()"));
          if (cycle.empty()) {
            break;
          } else {
            for (auto x : cycle) h.remove_vertex(x, false);
            ++ret;
          }
        }
        // at this point, g is very likely acyclic
        gs.pop_back();
        // ret += reducer_.reduce_lb(gs, -1);
      }
    }
    // PRINT("ret:" << ret);
    return ret;
  }

  int lower_bound_vc_full(LabeledGraph const& h) {
    PROF(if (prof_) prof_->start_timer("lb: lower_bound_vc_full()"));

    int ret = 0;
    std::vector<LabeledGraph> gs = {h};  // create a copy of the graph

    while (!gs.empty()) {
      // process one SCC
      int index = gs.size() - 1;

      if (gs[index].g.empty()) {
        gs.pop_back();
        continue;
      }

      if (gs[index].sg.has_any_edge()) {
        // (1) Find the exact vertex cover number of the strong graph.
        PROF(if (prof_) prof_->start_timer("lb: solve_vertex_cover_num"));
        ret += solve_vertex_cover_number(gs[index].sg);
        PROF(if (prof_) prof_->stop_timer("lb: solve_vertex_cover_num"));

        // (2) Remove all strong vertices, i.e. the vertices incident to any strong edge
        std::vector<int> to_remove;
        for (auto x: gs[index].sg.vertices()) {
          if (!gs[index].sg.is_isolate(x)) to_remove.push_back(x);
        }
        for (auto x: to_remove) gs[index].remove_vertex(x, false);
      } else {
        // (3) Find one short cycle
        auto cycle = short_cycle(gs[index].wg);
        if (cycle.empty()) {
          // done with this subgraph
          gs[index].clear();
          continue;
        }
        // PRINT("cycle: " << cycle);

        // debug
        // for (int i = 0; i < cycle.size(); ++i) {
        //   int u = cycle[i];
        //   int v = cycle[(i+1) % cycle.size()];
        //   if (!gs[index].sg.is_isolate(u)) throw std::runtime_error("u is not weak");
        //   if (!gs[index].sg.is_isolate(v)) throw std::runtime_error("v is not weak");
        //   if (!gs[index].wg.has_edge(u, v)) throw std::runtime_error("not a cycle");
        // }

        // (4) Remove the cycle
        for (auto x: cycle) gs[index].remove_vertex(x, false);
        ++ret;  // count up the lower-bound
      }

      // (5) Reduce and repeat
      ret += reducer_.reduce(gs, index).size();
    }

    PROF(if (prof_) prof_->stop_timer("lb: lower_bound_vc_full()"));
    return ret;
  }

  int lower_bound_vc(LabeledGraph const& h, int num_cycle_attempts) {
    PROF(if (prof_) prof_->start_timer("lb: lower_bound_vc()"));

    // Find the exact vertex cover number of the strong graph.
    PROF(if (prof_) prof_->start_timer("lb: solve_vertex_cover()"));
    int ret = solve_vertex_cover(h.sg).size();
    PROF(if (prof_) prof_->stop_timer("lb: solve_vertex_cover()"));

    if (h.sg.number_of_isolates() != h.g.number_of_isolates()) {
      // Then, find short cycles in the remaining graph.
      PROF(if (prof_) prof_->start_timer("lb: vc prep"));
      auto gg = h.wg;
      std::vector<int> to_remove;
      for (auto v : h.sg.vertices()) {
         // remove all strong vertices (i.e. the vertices with any strong edge)
        if (!h.sg.is_isolate(v)) to_remove.push_back(v);
      }
      for (auto x: to_remove) gg.remove_vertex(x);
      PROF(if (prof_) prof_->stop_timer("lb: vc prep"));

      PROF(if (prof_) prof_->start_timer("lb: vc cycle"));
      ret += lower_bound_cycle_packing(gg, num_cycle_attempts);
      PROF(if (prof_) prof_->stop_timer("lb: vc cycle"));
    }

    // experimental
    // LB2: VC of the entire graph -> not advantageus at all...
    // int ret2 = solve_vertex_cover(h.g).size();

    PROF(if (prof_) prof_->stop_timer("lb: lower_bound_vc()"));
    return ret;
  }

 private:
  int lower_bound_cycle_packing(data::graph::CLDiGraph const& g, int num_attempts) {
    PROF(if (prof_) prof_->start_timer("lb: cycle packing"));
    int n = g.number_of_nodes();

    int best = 0;
    if (n == 0) return best;

    for (int t = 0; t < num_attempts; ++t) {
      int ret = 0;
      auto gg = g;
      while (!gg.empty()) {
        PROF(if (prof_) prof_->start_timer("lb: short_cycle()"));
        auto cycle = short_cycle(gg, -1, t * 331 % n);
        PROF(if (prof_) prof_->stop_timer("lb: short_cycle()"));
        if (cycle.empty()) break;  // done
        for (auto x: cycle) gg.remove_vertex(x);
        ++ret;
      }
      best = std::max(best, ret);
    }
    PROF(if (prof_) prof_->stop_timer("lb: cycle packing"));
    return best;
  }
  
  std::vector<int> find_large_clique(data::graph::CLGraph const& g) {
    // find clique partition
    auto parts = clique_partition(g);
    if (parts.empty()) return {};

    // pick the largest one
    std::vector<int> ret;
    for (auto i : util::max_size_element(parts)) {
      ret.push_back(i);  // no need to convert labels
    }
    return ret;
  }
};

}  // namespace dfvs
}  // namespace algorithm
}  // namespace mog