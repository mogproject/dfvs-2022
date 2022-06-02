#pragma once

#include "Brancher.hpp"

namespace mog {
namespace algorithm {
namespace dfvs {

class EdgeBrancher : public Brancher {
 private:
  uint32_t seed_;
  bool include_first_;
  double include_rate_;
  std::unordered_map<long long, std::pair<int, bool>> record_;
  std::unordered_set<int> has_record_;
  std::pair<std::pair<int, int>, bool> next_target_;

 public:
  EdgeBrancher(Configuration const& conf, uint32_t seed, bool include_first, double include_rate)
      : Brancher(conf, seed),          //
        seed_(seed),                    //
        include_first_(include_first),  //
        include_rate_(include_rate)     //
  {
    //
  }

  void clear_record() {
    record_.clear();
    has_record_.clear();
  }

  int get_u() const { return next_target_.first.first; }
  int get_v() const { return next_target_.first.second; }
  bool get_include() const { return next_target_.second; }

  void set_u(int u) { next_target_.first.first = u; }
  void set_v(int v) { next_target_.first.second = v; }
  void set_include(bool should_include) { next_target_.second = should_include; }

  void branch_left(LabeledGraph& h, std::vector<int>& part_sol, bool record_history) {
    auto u = h.get_index(next_target_.first.first);
    auto v = h.get_index(next_target_.first.second);
    auto should_include = next_target_.second;

    if (should_include) {
      branch_include(h, part_sol, u, v, record_history);
    } else {
      branch_exclude(h, part_sol, u, v, record_history);
    }
  }

  void branch_right(LabeledGraph& h, std::vector<int>& part_sol, bool record_history) {
    auto u = h.get_index(next_target_.first.first);
    auto v = h.get_index(next_target_.first.second);
    auto should_include = next_target_.second;

    if (should_include) {
      branch_exclude(h, part_sol, u, v, record_history);
    } else {
      branch_include(h, part_sol, u, v, record_history);
    }
  }

  void branch_include(LabeledGraph& h, std::vector<int>& part_sol, int u, int v, bool record_history) {
    PROF(if (prof_) prof_->start_timer("branch: edge include"));
    h.add_edge(v, u, record_history);  // add the reversal edge
    PROF(if (prof_) prof_->stop_timer("branch: edge include"));
  }

  void branch_exclude(LabeledGraph& h, std::vector<int>& part_sol, int u, int v, bool record_history) {
    // PRINT("branch exclude: u=" << u << ", v=" << v);
    PROF(if (prof_) prof_->start_timer("branch: edge exclude"));
    auto sol_by_exclude = h.ignore_and_reduce2(u, v, record_history);  // returns ids
    for (auto x : sol_by_exclude) part_sol.push_back(h.get_label(x));
    PROF(if (prof_) prof_->stop_timer("branch: edge exclude"));
  }

  std::string to_string() const {
    return util::format("EdgeBrancher(%u,%s,%.1f,record=%lu)", seed_, include_first_ ? "true" : "false", include_rate_, record_.size());
  }

  // O(n) (want to avoid O(m))
  void find_target(LabeledGraph& h, bool use_record, bool keep_record) {
    PROF(if (prof_) prof_->start_timer("branch: edge find_target()"));
    int bestu = -1, bestv = -1;
    int best_prior = 1000000000;
    double best_score = -1e300;
    bool should_include = rand() < include_rate_;

    // computes weak degrees plus random tie-breaking

    // fix one vertex first
    for (auto u : h.wg.vertices()) {
      if (h.wg.out_degree(u) == 0) continue;

      auto score = h.wg.degree(u) + rand();
      // s += std::max(h.wg.in_degree(v), h.wg.out_degree(v)) * 0.8;  // favor unbalanced vertex
      if (!include_first_) score *= -1;

      if (use_record) {
        if (util::contains(has_record_, h.get_label(u))) score += 1e6;
      }

      if (score > best_score) {
        best_score = score;
        bestu = u;
      }
    }

    best_score = -1e300;
    for (auto v : h.wg.out_neighbors(bestu)) {
      if (use_record) {
        auto& r = record_[encode(h.get_label(bestu), h.get_label(v))];
        auto priority = r.first;
        if (priority > 0 && priority < best_prior) {  // priority-based decision
          best_prior = priority;
          bestv = v;
          should_include = r.second;
          continue;
        }
      }

      auto score = h.wg.degree(v) + rand();
      if (!include_first_) score *= -1;

      if (score > best_score) {
        best_score = score;
        bestv = v;
      }
    }

    if (bestu < 0 || bestv < 0) throw std::invalid_argument("no branch target");

    int u_label = h.get_label(bestu);
    int v_label = h.get_label(bestv);
    next_target_ = {{u_label, v_label}, should_include};

    if (keep_record) {
      int priority = record_.size() + 1;
      record_[encode(u_label, v_label)] = {priority, should_include};
      has_record_.insert(u_label);
    }
    PROF(if (prof_) prof_->stop_timer("branch: edge find_target()"));
  }

 private:
  inline long long encode(int u, int v) const { return (static_cast<long long>(u) << 32) | v; }

};

}  // namespace dfvs
}  // namespace algorithm
}  // namespace mog