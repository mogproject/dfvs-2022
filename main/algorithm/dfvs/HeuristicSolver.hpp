#pragma once

#include <memory>
#include <random>

#include "../../util/profiler.hpp"
#include "../component.hpp"
#include "../cycle.hpp"
#include "../vertex_cover.hpp"
#include "Configuration.hpp"
#include "LabeledGraph.hpp"
#include "branch/EdgeBrancher.hpp"
#include "branch/VertexBrancher.hpp"
#include "reduction/Reducer.hpp"

namespace mog {
namespace algorithm {
namespace dfvs {

enum HeuristicStrategy { VertexInclude, VertexExclude, EdgeInclude, EdgeExclude };
enum GuessResult { GuessUpdated, GuessNotUpdated, GuessTimeout };

/**
 * @brief Non-exact DFVS Solver
 *
 */
class HeuristicSolver {
 public:
  static constexpr int const INF = 1000000000;  // 1e9

 private:
  Configuration conf_;
  util::Profiler* const prof_;
  std::default_random_engine rand_gen_;
  std::uniform_real_distribution<> real_dist_;
  Reducer reducer_;
  int best_obj_;                // best objective value
  std::vector<int> best_cert_;  // best certificate
  std::chrono::system_clock::time_point guess_time_limit_;

 public:
  HeuristicSolver(Configuration conf)
      : conf_(conf),           //
        prof_(conf.prof),      //
        rand_gen_(conf.seed),  //
        real_dist_(0, 1),      //
        reducer_(conf)         //
  {
    //
  }

  /**
   * @brief Guesses good solution to each root problem.
   *
   * @param roots
   * @return std::vector<std::vector<int>>
   */
  std::vector<std::pair<std::shared_ptr<Brancher>, std::vector<int>>> guess(std::vector<LabeledGraph> const& roots) {
    // Assumption: every LabeledGraph has already been reduced and is compact.
    std::vector<std::pair<std::shared_ptr<Brancher>, std::vector<int>>> ret(roots.size());

    for (int index = 0; index < static_cast<int>(roots.size()); ++index) {
      // process each problem root
      auto start_time = std::chrono::system_clock::now();
      guess_time_limit_ = start_time + std::chrono::seconds(conf_.heuristic.guess_time_limit_sec);

      auto n = roots[index].g.number_of_nodes();
      auto m = roots[index].g.number_of_edges();

      bool is_large = n >= conf_.heuristic.num_guesses_large_thres_n && m >= conf_.heuristic.num_guesses_large_thres_m;

      if (n >= 100) {
        PRINT("Guessing root problem: index=" << index << " (n=" << n << ", m=" << m << ", m_s=" << roots[index].sg.number_of_edges()
                                              << ", m_w=" << roots[index].wg.number_of_edges() << ")");
      }

      uint32_t seed = conf_.seed;
      best_obj_ = INF;

      std::shared_ptr<Brancher> best_brancher;
      best_brancher = std::make_shared<VertexBrancher>(VertexBrancher(conf_, seed, true, 0.7));
      best_cert_.clear();
      // auto &h = roots[index];
      // PRINT("h    : index=" << 0 << " (n=" << h.g.number_of_nodes()
      //                                  << ", m=" << h.g.number_of_edges()
      //                                  << ", m_s=" << h.sg.number_of_edges()
      //                                  << ", m_w=" << h.wg.number_of_edges() << ")");

      for (int t = 0; t < conf_.heuristic.num_guesses; ++t) {
        // set up branchers
        std::vector<std::shared_ptr<Brancher>> branchers;
        branchers.push_back(std::make_shared<VertexBrancher>(VertexBrancher(conf_, seed, true, 0.7)));
        if (!is_large) branchers.push_back(std::make_shared<VertexBrancher>(VertexBrancher(conf_, seed, true, 0.5)));
        if (!is_large) branchers.push_back(std::make_shared<VertexBrancher>(VertexBrancher(conf_, seed, true, 0.9)));
        // if (!is_large) branchers.push_back(std::make_shared<VertexBrancher>(VertexBrancher(conf_, seed, true, 1.0)));
        if (!is_large) branchers.push_back(std::make_shared<VertexBrancher>(VertexBrancher(conf_, seed, false, 0.5)));
        // if (!is_large) branchers.push_back(std::make_shared<VertexBrancher>(VertexBrancher(conf_, seed, false, 0.3)));
        // if (!is_large) branchers.push_back(std::make_shared<VertexBrancher>(VertexBrancher(conf_, seed, false, 0.1)));
        if (!is_large) branchers.push_back(std::make_shared<EdgeBrancher>(EdgeBrancher(conf_, seed, true, 1.0)));
        if (!is_large) branchers.push_back(std::make_shared<EdgeBrancher>(EdgeBrancher(conf_, seed, true, 0.9)));
        // branchers.push_back(std::make_shared<EdgeBrancher>(EdgeBrancher(conf_, seed, false, 0.0)));
        // branchers.push_back(std::make_shared<EdgeBrancher>(EdgeBrancher(conf_, seed, false, 0.1)));
        ++seed;

        GuessResult guess_ret = GuessResult::GuessTimeout;
        for (auto& brancher : branchers) {
          // PRINT("Guessing: index=" << index << ",brancher=" << brancher->to_string());

          auto guess_ret = guess_once(roots[index], brancher);
          if (guess_ret == GuessResult::GuessTimeout) {
            PRINT("[index=" << index << ",brancher=" << brancher->to_string() << "] Time limit exceeded");
            break;
          } else if (guess_ret == GuessResult::GuessUpdated) {
            PRINT("[index=" << index << ",brancher=" << brancher->to_string() << "] Updated: obj=" << best_obj_);
            best_brancher = brancher;
          }
        }
        if (guess_ret == GuessResult::GuessTimeout) break;
      }

      if (best_obj_ == INF) {
        // time out during hte first guess
        ret[index] = {best_brancher, {-1}};
      } else {
        ret[index] = {best_brancher, best_cert_};
      }
    }
    return ret;
  }

  bool time_limit_exceeded() {
    return std::chrono::system_clock::now() > guess_time_limit_;
  }

  GuessResult guess_once(LabeledGraph const& h, std::shared_ptr<Brancher> brancher) {
    //------------------------------------------------------------
    // (1) Set up
    //------------------------------------------------------------
    std::vector<int> part_sol;  // partial solution
    std::vector<LabeledGraph> gs = {h};

    //------------------------------------------------------------
    // (2) Main loop
    //------------------------------------------------------------
    while (!gs.empty()) {
      // remove empty subgraphs
      if (gs.back().g.empty()) {
        gs.pop_back();
        continue;
      }

      // process one (non-empty) SCC
      int index = gs.size() - 1;

      // process until there are no weak edges
      while (gs[index].wg.has_any_edge()) {
        if (time_limit_exceeded()) return GuessTimeout;

        // (2-1) Pick next target
        brancher->find_target(gs[index], false, true);

        // (2-2) Branch on the next target
        brancher->branch_left(gs[index], part_sol, false);

        // (2-3) Reduce the resulting graph
        PROF(if (prof_) prof_->start_timer("guess: reduce"));
        auto sol_by_reduction = reducer_.reduce(gs, index);  // returns labels, not ids
        util::extend(part_sol, sol_by_reduction);
        PROF(if (prof_) prof_->stop_timer("guess: reduce"));
      }

      if (time_limit_exceeded()) return GuessResult::GuessTimeout;

      // (2-4) Find a vertex cover of the strong graph.
      // PRINT("partial solution size before vc: " << part_sol.size());
      // PRINT("guess: vc (n=" << gs[index].sg.number_of_nodes() << ", m=" << gs[index].sg.number_of_edges() << ")");
      PROF(if (prof_) prof_->start_timer("guess: vc"));
      auto vc = solve_vertex_cover(gs[index].sg, conf_.heuristic.vc_time_limit);
      PROF(if (prof_) prof_->stop_timer("guess: vc"));

      if (vc == std::vector<int>({-1})) {
        // VC solver could not find a solution; find cliques instead?
        // TODO
      } else {
        for (auto x : vc) part_sol.push_back(gs[index].get_label(x));
      }
      // PRINT("partial solution size after vc : " << part_sol.size());

      gs[index].clear();
    }

    //------------------------------------------------------------
    // (3) Local Search
    //------------------------------------------------------------
    if (time_limit_exceeded()) return GuessResult::GuessTimeout;
  
    part_sol = local_search(h, part_sol);

    if (time_limit_exceeded()) return GuessResult::GuessTimeout;

    // update result
    if (best_obj_ > static_cast<int>(part_sol.size())) {
      best_obj_ = part_sol.size();
      best_cert_ = part_sol;
      return GuessResult::GuessUpdated;
    }
    return GuessResult::GuessNotUpdated;
  }

  std::vector<int> local_search(LabeledGraph const& root, std::vector<int> current_solution) {
    PROF(if (prof_) prof_->start_timer("guess: local_search()"));
    // auto before_size = current_solution.size();

    std::vector<int> vtcs;  // vertex ids
    for (auto x : current_solution) vtcs.push_back(root.get_index(x));

    auto restore_vertex = [&](data::graph::CLDiGraph& g, int u) {
      g.add_vertex(u);
      for (auto w : root.g.in_neighbors(u)) {
        if (g.has_vertex(w)) g.add_edge_unsafe(w, u);
      }
      for (auto w : root.g.out_neighbors(u)) {
        if (g.has_vertex(w)) g.add_edge_unsafe(u, w);
      }
    };

    bool updated = true;
    while (updated) {
      if (time_limit_exceeded()) return {};
      updated = false;
      auto g = root.g;
      int s = current_solution.size();

      // remove all vertices in the current solution
      for (int i = 0; i < s; ++i) g.remove_vertex(vtcs[i]);

      //------------------------------------------------------------
      // 1-out
      //------------------------------------------------------------

      auto samples = util::sample(util::range(s), conf_.heuristic.local_search_num_samples_for_one_out, rand_gen_);
      int k = samples.size();

      for (int i = 0; i < k; ++i) {
        if (time_limit_exceeded()) return {};

        // restore i-th vertex
        auto ii = samples[i];
        auto v = vtcs[ii];
        restore_vertex(g, v);

        // check if S \ v is a solution
        if (strongly_connected_components(g, false).empty()) {
          updated = true;
          util::erase(current_solution, ii);
          util::erase(vtcs, ii);
          break;
        }

        // remove i-th vertex
        g.remove_vertex(v);
      }

      // at this point, all solution vertices must be removed
      //------------------------------------------------------------
      // 2-out, 1-in
      //------------------------------------------------------------
      // for (int i = 0; i < s; ++i) {
      //   if (g.has_vertex(vtcs[i])) throw std::runtime_error("sanity check");
      // }

      // pick two vertices to restore
      s = current_solution.size();  // refresh s (important!)
      auto samples2 = util::sample(util::range(s), conf_.heuristic.local_search_num_samples_for_two_out, rand_gen_);
      int k2 = samples2.size();
      data::fast_set fset(root.g.capacity());
      for (int i = 0; i < k2 - 1; ++i) {
        if (time_limit_exceeded()) return {};

        // restore i
        int ii = samples2[i];
        auto u = vtcs[ii];
        restore_vertex(g, u);

        for (int j = i + 1; j < k2; ++j) {
          // restore j
          int jj = samples2[j];
          auto v = vtcs[jj];
          restore_vertex(g, v);

          // find cycles passing u and v, respectively
          fset.clear();
          std::vector<int> candidates;
          for (auto w : find_cycle(g, u)) fset.set(w);
          for (auto w : find_cycle(g, v)) {
            if (fset.get(w)) candidates.push_back(w);
          }

          // check the intersection of those cycles
          for (auto w : candidates) {
            g.remove_vertex(w);

            // check if S \ {v,u} union {w} is a solution
            if (strongly_connected_components(g, false).empty()) {
              updated = true;
              // PRINT(s << " " << current_solution.size() << " " << ii << " " << jj);
              util::erase(current_solution, std::max(ii, jj));  // remove the larger index first
              util::erase(current_solution, std::min(ii, jj));
              current_solution.push_back(root.get_label(w));
              util::erase(vtcs, std::max(ii, jj));
              util::erase(vtcs, std::min(ii, jj));
              vtcs.push_back(w);
              // PRINT("found 2-out 1-in!");
              break;
            }

            restore_vertex(g, w);
          }
          if (updated) break;
          g.remove_vertex(v);
        }
        if (updated) break;
        g.remove_vertex(u);
      }
    }
    PROF(if (prof_) prof_->stop_timer("guess: local_search()"));
    // PRINT("local search: before=" << before_size << ", after=" << current_solution.size());
    return current_solution;
  }
};
}  // namespace dfvs
}  // namespace algorithm
}  // namespace mog