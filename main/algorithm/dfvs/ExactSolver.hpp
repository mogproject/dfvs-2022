#pragma once

#include <cassert>
#include <random>
#include <memory>

#include "../../util/profiler.hpp"
#include "../vertex_cover.hpp"
#include "Configuration.hpp"
#include "HeuristicSolver.hpp"
#include "LabeledGraph.hpp"
#include "bound/Bounder.hpp"
#include "branch/Brancher.hpp"
#include "reduction/Reducer.hpp"

// #include "../../readwrite/c_pace.hpp"

#define INF_SET std::vector<int>({-1})

namespace mog {
namespace algorithm {
namespace dfvs {

/**
 * @brief DFVS Solver.
 *
 */
class ExactSolver {
 public:
  static constexpr int const INF = 1000000000;  // 1e9
  static int cardinality(std::vector<int> const& xs) { return !xs.empty() && xs[0] < 0 ? INF : xs.size(); }

 private:
  util::Profiler* const prof_;
  Configuration conf_;
  Reducer reducer_;
  Bounder bounder_;

  int best_obj_;                // the currently best objective value
  std::vector<int> best_cert_;  // the currently best certificate
  std::vector<int> part_sol_;   // partial solution

 public:
  ExactSolver(Configuration const& conf)
      : prof_(conf.prof),  //
        conf_(conf),       //
        reducer_(conf),    //
        bounder_(conf)    //
  {}

  /**
   * @brief Checks if the given vertex set is a solution (not necessarily optimal)
   *        to the given Directed Feedback Vertex Set instance.
   *
   * @param g directed graph
   * @param vs vertices
   * @return true if the vertex set is a feedback set
   */
  bool is_solution(data::graph::CLDiGraph g, std::vector<int> vs) const {
    for (auto v : vs) {
      if (!g.has_vertex(v)) return false;  // invalid index or duplicate
      g.remove_vertex(v);
    } 
    return algorithm::strongly_connected_components(g, false).empty();
  }

  std::vector<int> solve(data::graph::CLDiGraph const& g) {
    // PRINT("lb=" << bounder_.lower_bound(LabeledGraph(g)));
    // return {};

    //------------------------------------------------------------
    // (1) Reduce the original instance.
    //------------------------------------------------------------
    PROF(if (prof_) prof_->start_timer("initial reduction"));
    std::vector<LabeledGraph> roots = {LabeledGraph(g)};
    auto root_solution = reducer_.reduce(roots, 0, true);  // enforce compaction
    PROF(if (prof_) prof_->stop_timer("initial reduction"));

    // now, all subgraphs are compact
    PRINT("Initially-reduced solution size: " << root_solution.size());
    // PRINT("reduced: " << root_solution);

    //------------------------------------------------------------
    // (2) Find reducible components.
    //------------------------------------------------------------
    std::vector<LabeledGraph> reduced_roots;
    for (auto& h : roots) {
      if (!h.wg.has_any_edge()) {
        // reducible
        if (h.sg.number_of_nodes() >= 100) {
          PRINT("exact: vc (n=" << h.sg.number_of_nodes() << ", m=" << h.sg.number_of_edges() << ")");
        }
        PROF(if (prof_) prof_->start_timer("vc: reducible"));
        auto vc = solve_vertex_cover(h.sg);
        PROF(if (prof_) prof_->stop_timer("vc: reducible"));

        for (auto x : vc) root_solution.push_back(h.get_label(x));  // extend the solution
      } else {
        reduced_roots.push_back(h);
      }
    }

    // if the entire graph is reducible, we are done
    if (reduced_roots.empty()) {
      PRINT("Reducible instance");
      return root_solution;
    }
    PRINT("# of irreducible instances: " << reduced_roots.size());

    //------------------------------------------------------------
    // (3) Guess good solutions.
    //------------------------------------------------------------
    std::vector<std::pair<std::shared_ptr<Brancher>, std::vector<int>>> guess_result(reduced_roots.size(), {std::make_shared<VertexBrancher>(VertexBrancher(conf_, conf_.seed, true, 0.7)), INF_SET});

    if (conf_.heuristic.enabled) {
      HeuristicSolver hsolver(conf_);

      PROF(if (prof_) prof_->start_timer("guess()"));
      guess_result = hsolver.guess(reduced_roots);
      PROF(if (prof_) prof_->stop_timer("guess()"));

      int guess_solution_size = root_solution.size();
      for (auto& r : guess_result) guess_solution_size += r.second.size();
      PRINT("Guessed solution size: " << guess_solution_size);
    }

    // finally, we are ready to search for exact solution

    //------------------------------------------------------------
    // (4) Run branch-and-reduce.
    //------------------------------------------------------------
    for (int index = 0; index < static_cast<int>(reduced_roots.size()); ++index) {
      auto brancher = guess_result[index].first;
      PRINT("Solving root problem: index=" << index << " (n=" << reduced_roots[index].g.number_of_nodes()
                                           << ", m=" << reduced_roots[index].g.number_of_edges()
                                           << ", m_s=" << reduced_roots[index].sg.number_of_edges()
                                           << ", m_w=" << reduced_roots[index].wg.number_of_edges() << ")"
                                           << " brancher=" << brancher->to_string());

      
      int best_card = cardinality(guess_result[index].second);

      // leaders must be included in the solution
      std::vector<int> ldr_solution;
      if (best_card <= 1000) {
        ldr_solution = find_leaders(reduced_roots[index], guess_result[index].second);
        PRINT("Number of leaders: " << ldr_solution.size());
      }

      // non-recursion
      brancher->set_guess_solution(guess_result[index].second);
      
      auto bnr_solution = solve_nonrec(reduced_roots, index, best_card, brancher, ldr_solution);

      if (cardinality(bnr_solution) == INF) {
        PRINT("Using guessed solution: index=" << index);
        util::extend(root_solution, guess_result[index].second);  // guess solution is actually optimal
      } else {
        PRINT("Solution found: index=" << index << ", size=" << bnr_solution.size());
        util::extend(root_solution, bnr_solution);
      }
    }

    return root_solution;
  }

  // returns the labels of leader vertices
  std::vector<int> find_leaders(LabeledGraph &h, std::vector<int> const& guess_solution) {
    std::vector<int> ret;
    if (guess_solution == INF_SET) return ret;

    int ub = guess_solution.size();

    for (auto x: guess_solution) {
      h.commit();
      auto red = h.ignore_and_reduce(h.get_index(x), true);
      int lb = bounder_.lower_bound(h) + red.size();
      // PRINT("find_leaders: x=" << x << ": ub=" << ub << ", lb=" << lb);
      if (lb >= ub) ret.push_back(x);
      h.rollback();
    }
    return ret;
  }

 private:
  /** Information for a blue (branch) node */
  struct BlueInfo {
    int num_done;          // number of finished children
    int index;             // target graph
    int ub;                // upper-bound
    std::vector<int> sol;  // best solution so far
    int u_label;           // branch on edge u-v
    int v_label;
    bool include_first;

    BlueInfo(int index, int ub)
        : num_done(0), index(index), ub(ub), sol(INF_SET), u_label(-1), v_label(-1), include_first(true) {}

    bool has_solution() const { return sol != INF_SET; }

    std::string to_string() const {
      return std::string(
          util::format("BlueInfo(num_done=%d, index=%d, ub=%d, sol=%s, u_label=%d, v_label=%d, include_first=%d)",  //
                       num_done, index, ub, cstr(sol), u_label, v_label, include_first));
    }
  };

  /** Information for a red (reduction) node */
  struct RedInfo {
    int num_done;                  // number of finished subgraphs
    std::vector<int> target_idx;   // target indices
    std::vector<int> lower_bound;  // lower-bound for each subgraph
    int lower_bound_sum;           // sum of the lower-bounds
    std::vector<int> part_sol;     // partial solution
    int offset;                    // original number of labeled graphs when this node was initialized

    RedInfo(int index, int offset) : num_done(0), target_idx({index}), lower_bound_sum(0), offset(offset) {}

    bool has_solution() const { return part_sol != INF_SET; }
    int get_sol_size() const { return cardinality(part_sol); }

    std::string to_string() const {
      return std::string(util::format(
          "RedInfo(num_done=%d, target_idx=%s, lower_bound=%s, lower_bound_sum=%d, part_sol=%s, offset=%d)",  //
          num_done, cstr(target_idx), cstr(lower_bound), lower_bound_sum, cstr(part_sol), offset));
    }
  };

  std::vector<int> solve_nonrec(std::vector<LabeledGraph>& gs, int init_index, int init_ub, std::shared_ptr<Brancher> brancher, std::vector<int> const& leaders) {
    std::vector<BlueInfo> blue_stack = {BlueInfo(init_index, init_ub)};
    std::vector<RedInfo> red_stack = {RedInfo(init_index, gs.size())};

    // PRINT("init_ub=" << init_ub << ", leaders=" << leaders);

    // must remove leaders
    for (auto x: leaders) gs[init_index].remove_vertex(gs[init_index].get_index(x), false);
    red_stack.back().part_sol = leaders;
    gs[init_index].commit();

    std::vector<int> q;
    q.push_back(2);  // blue back
    q.push_back(3);  // red init

    while (!q.empty()) {
      PROF(if (prof_) prof_->count("solve iterations"));
      int op = q.back();
      q.pop_back();

      // debug print
      // PRINT("op[" << op << "], q=" << q);
      // for (int i = 0; i < static_cast<int>(gs.size()); ++i) {
      //   // PRINT("  gs[" << i << "]: labels=" << gs[i].get_labels() << ", edges=" << gs[i].g.edges());
      // }
      // for (int i = 0; i < static_cast<int>(blue_stack.size()); ++i) {
      //   PRINT("  b [" << i << "]: " << blue_stack[i].to_string());
      // }
      // for (int i = 0; i < static_cast<int>(red_stack.size()); ++i) {
      //   PRINT("  r [" << i << "]: " << red_stack[i].to_string());
      // }
      // if (blue_stack.size() != red_stack.size() && blue_stack.size() != red_stack.size() + 1) throw std::runtime_error("stack inconsistency");

      switch (op) {
        case 0: solve_blue_init(gs, q, blue_stack, brancher); break;
        case 1: solve_blue_forward(gs, q, blue_stack, red_stack, brancher); break;
        case 2: solve_blue_backward(gs, q, blue_stack, red_stack); break;
        case 3: solve_red_init(gs, q, blue_stack, red_stack); break;
        case 4: solve_red_forward(gs, q, blue_stack, red_stack); break;
        case 5: solve_red_backward(gs, q, blue_stack, red_stack); break;
        default: throw std::runtime_error("unknown op code"); ;
      }
    }
    // PRINT("done");
    // for (int i = 0; i < static_cast<int>(gs.size()); ++i) {
    //   PRINT("  gs[" << i << "]: labels=" << gs[i].get_labels() << ", edges=" << gs[i].g.edges());
    // }
    // for (int i = 0; i < static_cast<int>(blue_stack.size()); ++i) {
    //   PRINT("  b [" << i << "]: " << blue_stack[i].to_string());
    // }
    // for (int i = 0; i < static_cast<int>(red_stack.size()); ++i) {
    //   PRINT("  r [" << i << "]: " << red_stack[i].to_string());
    // }
    return blue_stack[0].sol;
  }

  /**
   * @brief Blue Init: New Blue -> Filled Blue
   *
   * @param gs labeled graphs
   * @param q operation stack
   * @param blue_stack blue stack
   */
  void solve_blue_init(std::vector<LabeledGraph>& gs, std::vector<int>& q, std::vector<BlueInfo>& blue_stack, std::shared_ptr<Brancher> brancher) {
    PROF(if (prof_) prof_->start_timer("solve blue init"));
    auto& bn = blue_stack.back();

    PROF(if (prof_) prof_->start_timer("branch target"));
    // auto branch_target = brancher_.branch_target2(gs[bn.index]);
    brancher->find_target(gs[bn.index], true, false);
    PROF(if (prof_) prof_->stop_timer("branch target"));

    /** @note stores label instead of index because indices may change after `rollback()` */
    bn.u_label = brancher->get_u();
    bn.v_label = brancher->get_v();
    bn.include_first = brancher->get_include();

    q.push_back(2);
    q.push_back(1);  // to the right child
    q.push_back(2);
    q.push_back(1);  // to the left child
    PROF(if (prof_) prof_->stop_timer("solve blue init"));
  }

  /**
   * @brief Blue Forward: Filled Blue -> New Red
   *
   * @param gs labeled graphs
   * @param q operation stack
   * @param blue_stack blue stack
   * @param red_stack red stack
   */
  void solve_blue_forward(std::vector<LabeledGraph>& gs, std::vector<int>& q, std::vector<BlueInfo>& blue_stack,
                          std::vector<RedInfo>& red_stack, std::shared_ptr<Brancher> brancher) {
    PROF(if (prof_) prof_->start_timer("solve blue forward"));
    auto& bn = blue_stack.back();
    red_stack.push_back(RedInfo(bn.index, gs.size()));  // new red node
    auto& rn = red_stack.back();

    gs[bn.index].commit();  // take snapshot

    brancher->set_u(bn.u_label);
    brancher->set_v(bn.v_label);
    brancher->set_include(bn.include_first);

    if (bn.num_done == 0) {
      brancher->branch_left(gs[bn.index], rn.part_sol, true);
    } else {
      brancher->branch_right(gs[bn.index], rn.part_sol, true);
    }

    if (rn.get_sol_size() >= bn.ub) {
      // cannot be improved and do not proceed
      rn.part_sol = INF_SET;
      PROF(if (prof_) prof_->stop_timer("solve blue forward"));
      return;
    }

    // gs[bn.index].check_properties(util::format("solve_blue_forward:after %d, %s", bn.num_done, brancher->to_string().c_str()));
    // if ((bn.num_done == 0) == bn.include_first) {
    //   // include
    //   PROF(if (prof_) prof_->start_timer("branch graft process"));
    //   gs[bn.index].add_edge(v, u, true);  // add reverse edge to make strong
    //   PROF(if (prof_) prof_->stop_timer("branch graft process"));
    // } else {
    //   // exclude
    //   PROF(if (prof_) prof_->start_timer("branch exclude process"));
    //   auto exclude_sol = gs[bn.index].ignore_and_reduce2(u, v, true);  // add solution to red node if exists
    //   if (static_cast<int>(exclude_sol.size()) >= bn.ub) {
        
    //     rn.part_sol = INF_SET;
        
    //   }
    //   for (auto x : exclude_sol) rn.part_sol.push_back(gs[bn.index].get_label(x));
      // PROF(if (prof_) prof_->stop_timer("branch exclude process"));
    // }
    q.push_back(3);
    PROF(if (prof_) prof_->stop_timer("solve blue forward"));
  }
  /**
   * @brief Blue Backword: Finished Red -> Blue
   *
   * @param gs labeled graphs
   * @param q operation stack
   * @param blue_stack blue stack
   * @param red_stack red stack
   */
  void solve_blue_backward(std::vector<LabeledGraph>& gs, std::vector<int>& q, std::vector<BlueInfo>& blue_stack,
                           std::vector<RedInfo>& red_stack) {
    PROF(if (prof_) prof_->start_timer("solve blue backward"));
    auto& bn = blue_stack.back();
    auto& rn = red_stack.back();

    // collect results from red child
    if (rn.has_solution()) {
      // if (bn.ub <= static_cast<int>(rn.part_sol.size())) {
      //   //throw std::runtime_error("upper-bound inconsistency");
      //   PRINT(bn.ub);
      //   PRINT(rn.part_sol.size());
      //   PRINT(rn.lower_bound);
      // } else {
      bn.sol = red_stack.back().part_sol;
      bn.ub = bn.sol.size();
      // }
    }

    // increment done
    ++bn.num_done;

    // shrink graphs and remove red node
    while (rn.offset < static_cast<int>(gs.size())) gs.pop_back();
    red_stack.pop_back();

    // rollback
    PROF(if (prof_) prof_->start_timer("branch rollback"));
    gs[bn.index].rollback();
    PROF(if (prof_) prof_->stop_timer("branch rollback"));
    PROF(if (prof_) prof_->stop_timer("solve blue backward"));
  }

  /**
   * @brief Red Init: New Red -> Filled Red
   *
   * @param gs labeled graphs
   * @param q operation stack
   * @param blue_stack blue stack
   * @param red_stack red stack
   */
  void solve_red_init(std::vector<LabeledGraph>& gs, std::vector<int>& q, std::vector<BlueInfo>& blue_stack,
                      std::vector<RedInfo>& red_stack) {
    PROF(if (prof_) prof_->start_timer("solve red init"));
    auto& bn = blue_stack.back();
    auto& rn = red_stack.back();
    int index = rn.target_idx[0];

    // gs[index].check_properties("solve_red_init");
    // reduction
    PROF(if (prof_) prof_->start_timer("reduce()"));
    auto red_sol = reducer_.reduce(gs, index);  // reduce returns labels
    PROF(if (prof_) prof_->stop_timer("reduce()"));

    if (rn.get_sol_size() + static_cast<int>(red_sol.size()) >= bn.ub) {
      rn.part_sol = INF_SET;
      PROF(if (prof_) prof_->count("branch cut: after reduction"));
      PROF(if (prof_) prof_->stop_timer("solve red init"));
      return;
    }
    util::extend(rn.part_sol, red_sol);

    // list of (graph size, index) examined at this level
    int process_sign = conf_.process_small_subgraph_first ? 1 : -1;
    std::vector<std::pair<int, int>> target = {{process_sign * gs[index].g.number_of_nodes(), index}};
    for (int j = rn.offset; j < static_cast<int>(gs.size()); ++j) {
      int n = gs[j].g.number_of_nodes();
      if (n > 0) target.push_back({process_sign * n, j});  // filter non-empty graphs
    }
    std::sort(target.begin(), target.end());  // order subgraphs

    // check for reducibility
    rn.target_idx.clear();  // reset target idx
    for (int i = 0; i < static_cast<int>(target.size()); ++i) {
      bool reduced = false;
      int idx = target[i].second;

      if (conf_.use_vc_solver) {
        auto& sg = gs[idx].sg;
        if (sg.number_of_edges() * 2 == gs[idx].g.number_of_edges()) {
          //------------------------------------------------------------
          // Strong graph -> reduce to Vertex Cover
          //------------------------------------------------------------
          PROF(if (prof_) prof_->start_timer("solve_vertex_cover()"));
          auto vc_sol = solve_vertex_cover(sg);
          PROF(if (prof_) prof_->stop_timer("solve_vertex_cover()"));
          for (auto x : vc_sol) rn.part_sol.push_back(gs[idx].get_label(x));
          reduced = true;
        }
      }

      // compute lower-bounds
      if (!reduced) {
        PROF(if (prof_) prof_->start_timer("lower_bound()"));
        int lb = bounder_.lower_bound(gs[idx]);
        PROF(if (prof_) prof_->stop_timer("lower_bound()"));

        rn.target_idx.push_back(idx);
        rn.lower_bound.push_back(lb);
        rn.lower_bound_sum += lb;
      }

      // bound check
      if (rn.get_sol_size() + rn.lower_bound_sum >= bn.ub) {
        // PRINT("idx: " << idx);
        // // PRINT("labels: " << gs[idx].get_labels());
        // // PRINT("g edges: " << gs[idx].g.edges());
        // // PRINT("sg edges: " << gs[idx].sg.edges());
        // // PRINT("wg edges: " << gs[idx].wg.edges());
        // PRINT("lb: " << rn.lower_bound);
        // PRINT("lb sum: " << rn.lower_bound_sum);
        // PRINT("part sol: " << rn.part_sol);
        // PRINT("ub: " << bn.ub);

        // readwrite::save_pace_2022("debug_out2.in", gs[idx].g, true);

        rn.part_sol = INF_SET;
        PROF(if (prof_) prof_->count("branch cut: lower-bound"));
        PROF(if (prof_) prof_->stop_timer("solve red init"));
        return;
      }
    }

    if (!rn.target_idx.empty()) {
      q.push_back(5);
      q.push_back(4);
    }
    PROF(if (prof_) prof_->stop_timer("solve red init"));
  }

  /**
   * @brief Red Forward: Filled Red -> New Blue
   *
   * @param gs labeled graphs
   * @param q operation stack
   * @param blue_stack blue stack; will be pushed
   * @param red_stack red stack
   */
  void solve_red_forward(std::vector<LabeledGraph>& gs, std::vector<int>& q, std::vector<BlueInfo>& blue_stack,
                         std::vector<RedInfo>& red_stack) {
    PROF(if (prof_) prof_->start_timer("solve red forward"));
    auto& rn = red_stack.back();
    auto ub = blue_stack.back().ub - rn.get_sol_size() - rn.lower_bound_sum + rn.lower_bound[rn.num_done];
    auto idx = rn.target_idx[rn.num_done];

    // if (ub <= 0) throw std::runtime_error("solve_red_forward: nonpositive upper-bound");
    blue_stack.push_back(BlueInfo(idx, ub));

    q.push_back(0);
    PROF(if (prof_) prof_->stop_timer("solve red forward"));
  }

  /**
   * @brief Red Backward: Finished Blue -> Red
   *
   * @param gs labeled graphs
   * @param q operation stack
   * @param blue_stack blue stack; will be popped
   * @param red_stack red stack
   */
  void solve_red_backward(std::vector<LabeledGraph>& gs, std::vector<int>& q, std::vector<BlueInfo>& blue_stack,
                          std::vector<RedInfo>& red_stack) {
    PROF(if (prof_) prof_->start_timer("solve red backward"));
    auto& bn = blue_stack.back();
    auto& rn = red_stack.back();
    int i = rn.num_done;
    ++rn.num_done;

    if (!bn.has_solution()) {
      // child has no better solution
      rn.part_sol = INF_SET;
    } else {
      // update lower-bound sum
      rn.lower_bound_sum -= rn.lower_bound[i];
      util::extend(rn.part_sol, bn.sol);

      if (rn.num_done == static_cast<int>(rn.target_idx.size())) {
        // all done
      } else {
        // move to next subgraph
        q.push_back(5);
        q.push_back(4);
      }
    }

    blue_stack.pop_back();
    PROF(if (prof_) prof_->stop_timer("solve red backward"));
  }
};

}  // namespace dfvs
}  // namespace algorithm
}  // namespace mog