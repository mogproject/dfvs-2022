#pragma once

#include "Brancher.hpp"

namespace mog {
namespace algorithm {
namespace dfvs {

class VertexBrancher : public Brancher {
 private:
  uint32_t seed_;
  bool include_first_;
  double include_rate_;
  std::unordered_map<int, std::pair<int, bool>> record_;
  std::pair<int, bool> next_target_;

 public:
  VertexBrancher(Configuration const& conf, uint32_t seed, bool include_first, double include_rate)
      : Brancher(conf, seed),          //
        seed_(seed),                    //
        include_first_(include_first),  //
        include_rate_(include_rate)     //
  {
    //
  }

  void clear_record() { record_.clear(); }

  int get_u() const { return next_target_.first; }
  int get_v() const { return next_target_.first; }
  bool get_include() const { return next_target_.second; }

  void set_u(int u) { next_target_.first = u; }
  void set_v(int v) { next_target_.first = v; }
  void set_include(bool should_include) { next_target_.second = should_include; }

  void branch_left(LabeledGraph& h, std::vector<int>& part_sol, bool record_history) {
    auto v = h.get_index(next_target_.first);
    auto should_include = next_target_.second;

    if (should_include) {
      branch_include(h, part_sol, v, record_history);
    } else {
      branch_exclude(h, part_sol, v, record_history);
    }
  }

  void branch_right(LabeledGraph& h, std::vector<int>& part_sol, bool record_history) {
    auto v = h.get_index(next_target_.first);
    auto should_include = next_target_.second;

    if (should_include) {
      branch_exclude(h, part_sol, v, record_history);
    } else {
      branch_include(h, part_sol, v, record_history);
    }
  }

  void branch_include(LabeledGraph& h, std::vector<int>& part_sol, int v, bool record_history) {
    PROF(if (prof_) prof_->start_timer("branch: vtx include"));
    part_sol.push_back(h.get_label(v));
    h.remove_vertex(v, record_history);
    PROF(if (prof_) prof_->stop_timer("branch: vtx include"));
  }

  void branch_exclude(LabeledGraph& h, std::vector<int>& part_sol, int v, bool record_history) {
    PROF(if (prof_) prof_->start_timer("branch: vtx exclude"));
    auto sol_by_exclude = h.ignore_and_reduce(v, record_history);  // returns ids
    for (auto x : sol_by_exclude) part_sol.push_back(h.get_label(x));
    PROF(if (prof_) prof_->stop_timer("branch: vtx exclude"));
  }

  std::string to_string() const {
    return util::format("VertexBrancher(%u,%s,%.1f,record=%lu)", seed_, include_first_ ? "true" : "false", include_rate_, record_.size());
  }

  // O(n)
  void find_target(LabeledGraph& h, bool use_record, bool keep_record) {
    PROF(if (prof_) prof_->start_timer("branch: vtx find_target()"));

    int bestv = -1;
    double best_prior = 1e9;
    double best_score = -1e300;
    bool should_include = rand() < include_rate_;

    // computes weak degrees plus random tie-breaking
    for (auto v : h.wg.vertices()) {
      double score = 0;

      int v_label = h.get_label(v);
      if (is_in_guess_solution(v_label)) score += 100000;

      if (use_record) {
        auto& r = record_[h.get_label(v)];
        score += r.first;
        if (r.first > 0 && score < best_prior) {  // priority-based decision
          best_prior = score;
          bestv = v;
          should_include = r.second;
          continue;
        }
      }

      score += h.wg.degree(v) + rand();
      if (!include_first_) score *= -1;
      // s += std::max(h.wg.in_degree(v), h.wg.out_degree(v)) * 0.8;  // favor unbalanced vertex
      if (score > best_score) {
        best_score = score;
        bestv = v;
      }
    }

    if (bestv < 0) {
      throw std::invalid_argument("no branch target");
    }

    auto v_label = h.get_label(bestv);
    next_target_ = {v_label, should_include};

    if (keep_record) {
      int priority = record_.size() + 1;
      record_[v_label] = {priority, should_include};
    }

    PROF(if (prof_) prof_->stop_timer("branch: vtx find_target()"));
  }

 private:
};

}  // namespace dfvs
}  // namespace algorithm
}  // namespace mog