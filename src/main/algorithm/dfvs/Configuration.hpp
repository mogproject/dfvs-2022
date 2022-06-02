#pragma once

#include "../../util/profiler.hpp"

namespace mog {
namespace algorithm {
namespace dfvs {

enum BoundStrategy { BoundConstant, BoundPacking, BoundPackingFull, BoundVertexCoverFull, BoundVertexCoverOnly };
enum BranchOrder { BranchOrderHeavy, BranchOrderLight, BranchOrderRandom };
enum SampleStrategy { SampleLeastHeavy, SampleSeparator };

class Configuration {
 private:
  util::Profiler profiler_;

 public:
  /** Pointer to the profiler. */
  util::Profiler *prof = nullptr;

  uint32_t seed;
  bool use_vc_solver;  // use external Vertex Cover solver
  double vc_threshold;
  bool process_small_subgraph_first;

  static constexpr int const reduce_offset = 10;
  static constexpr int const bound_offset = 16;
  static constexpr int const branch_offset = 24;

  struct {
    bool enabled = true;
    int guess_time_limit_sec = 300;  // guess time limit in seconds
    int num_guesses = 50;
    int num_guesses_large_thres_n = 1600; // 2000;
    int num_guesses_large_thres_m = 10000; // 30000;
    int vc_time_limit = 1;  // in minutes
    int local_search_num_samples_for_one_out = 2000;
    int local_search_num_samples_for_two_out = 100;
  } heuristic;

  struct {
    bool complete_reduction;
    bool n4_reduction;
    bool pie_reduction;
    bool core_reduction;
    bool dome_reduction;
    double compact_threshold;
  } reduce;

  struct {
    BoundStrategy strategy;
    double exec_rate;
  } bound;

  struct {
    BranchOrder order;
    bool include_first;
    double alpha;  // weight scale for strong edges
    SampleStrategy sample_strategy;
    double sample_rate;
    int sample_min;
    int sample_max;
    bool sample_ignore;  // ignore or remove
    bool sample_reduce;  // perform reduction or just decompose into SCCs
  } branch;

  Configuration(uint64_t params, uint32_t seed) : profiler_(50), seed(seed) {
    int offset = 0;
    if (get_bool(params, offset)) prof = &profiler_;

    use_vc_solver = get_bool(params, offset);
    vc_threshold = get_negative_power(params, offset);
    process_small_subgraph_first = get_bool(params, offset);
    heuristic.enabled = get_bool(params, offset);

    offset = reduce_offset;
    reduce.complete_reduction = true;
    reduce.n4_reduction = true;
    reduce.pie_reduction = get_bool(params, offset);   // 1 bit
    reduce.core_reduction = get_bool(params, offset);  // 1 bit
    reduce.dome_reduction = get_bool(params, offset);  // 1 bit
    reduce.compact_threshold = 0.5;

    offset = bound_offset;
    switch (get_int(params, offset, 2)) {  // 2 bits
      case 0: bound.strategy = BoundStrategy::BoundConstant; break;
      case 1: bound.strategy = BoundStrategy::BoundPacking; break;
      case 2: bound.strategy = BoundStrategy::BoundPackingFull; break;
      case 3: bound.strategy = BoundStrategy::BoundVertexCoverFull; break;
      default: throw std::invalid_argument("unknown bound strategy");
    }
    bound.exec_rate = get_negative_power(params, offset);  // 4 bits

    offset = branch_offset;
    switch (get_int(params, offset, 2)) {  // 2 bits
      case 0: branch.order = BranchOrder::BranchOrderHeavy; break;
      case 1: branch.order = BranchOrder::BranchOrderLight; break;
      case 2: branch.order = BranchOrder::BranchOrderRandom; break;
      default: throw std::invalid_argument("unknown branch order");
    }
    branch.include_first = get_bool(params, offset);  // 1 bit
    branch.alpha = get_int(params, offset, 3) - 2;    // 3 bits

    switch (get_int(params, offset, 2)) {  // 2 bits
      case 0: branch.sample_strategy = SampleStrategy::SampleLeastHeavy; break;
      case 1: branch.sample_strategy = SampleStrategy::SampleSeparator; break;
      default: throw std::invalid_argument("unknown sample strategy");
    }

    branch.sample_rate = get_negative_power(params, offset);  // 4 bits
    branch.sample_min = get_int(params, offset, 10);          // 10 bits
    branch.sample_max = get_int(params, offset, 10);          // 10 bits
    branch.sample_ignore = get_bool(params, offset);          // 1 bit
    branch.sample_reduce = get_bool(params, offset);          // 1 bit
  }

  static constexpr uint64_t default_params() {
    uint64_t params = 0;
    params |= 1ULL << 1;
    params |= 6ULL << 2;                    // almost strong 1/64
    params |= 1ULL << 6;                    // process small subgraph first
    params |= 1ULL << 7;                    // guess
    params |= 1ULL << (reduce_offset);      // PIE
    params |= 1ULL << (reduce_offset + 1);  // CORE
    params |= 1ULL << (reduce_offset + 2);  // DOME

    // params |= 1ULL << (bound_offset);      // bound strategy: pack
    params |= 3ULL << (bound_offset);      // bound strategy: VC
    params |= 0ULL << (bound_offset + 2);  // bound exec rate: 1/1

    params |= 0ULL << (branch_offset);       // branch strategy
    params |= 1ULL << (branch_offset + 2);   // include first
    params |= 0ULL << (branch_offset + 3);   // branch alpha
    params |= 0ULL << (branch_offset + 6);   // sample strategy
    params |= 6ULL << (branch_offset + 8);   // sample rate 1/64
    params |= 1ULL << (branch_offset + 12);  // sample min
    params |= 1ULL << (branch_offset + 22);  // sample max (1: no sampling)
    params |= 1ULL << (branch_offset + 32);  // ignore (1) or remove(0)
    params |= 1ULL << (branch_offset + 33);  // reduce after sampling enabled (1)
    return params;
  }

  std::string to_string() const {
    std::stringstream ss;
    ss << "[Configuration]\n";
    ss << util::format("Profiler              : %s\n", prof ? "on" : "off");
    ss << util::format("Seed                  : %u\n", seed);
    ss << util::format("Vertex Cover threshold: %.3f\n", vc_threshold);
    ss << util::format("Small subgraph first  : %s\n", process_small_subgraph_first ? "yes" : "no");
    ss << util::format("Guess by heuristics   : %s\n", heuristic.enabled ? "on" : "off");
    ss << util::format("Complete reduction    : %s\n", reduce.complete_reduction ? "on" : "off");
    ss << util::format("N4 reduction          : %s\n", reduce.n4_reduction ? "on" : "off");
    ss << util::format("PIE reduction         : %s\n", reduce.pie_reduction ? "on" : "off");
    ss << util::format("CORE reduction        : %s\n", reduce.core_reduction ? "on" : "off");
    ss << util::format("DOME reduction        : %s\n", reduce.dome_reduction ? "on" : "off");
    ss << util::format("Compact threshold     : %.3f\n", reduce.compact_threshold);
    ss << util::format("Lower-bound strategy  : %d\n", bound.strategy);
    ss << util::format("Lower-bound rate      : %.3f\n", bound.exec_rate);
    ss << util::format("Branch strategy       : %d\n", branch.order);
    ss << util::format("Include first         : %s\n", branch.include_first ? "yes" : "no");
    ss << util::format("Branch alpha value    : %.3f\n", branch.alpha);
    ss << util::format("Branch sample strategy: %d\n", branch.sample_strategy);
    ss << util::format("Branch sample rate    : %.3f\n", branch.sample_rate);
    ss << util::format("Branch sample minimum : %d\n", branch.sample_min);
    ss << util::format("Branch sample maximum : %d\n", branch.sample_max);
    ss << util::format("Branch sample ignore  : %s\n", branch.sample_ignore ? "yes" : "no");
    ss << util::format("Branch sample reduce  : %s", branch.sample_reduce ? "yes" : "no");
    return ss.str();
  }

 private:
  int get_int(uint64_t params, int &offset, int size) {
    int ret = (params >> offset) & ((1ULL << size) - 1);
    offset += size;
    return ret;
  }

  bool get_bool(uint64_t params, int &offset) { return get_int(params, offset, 1) == 1; }

  /**
   * @brief Returns 2^(-x), where x is a 4-bit integer.
   */
  double get_negative_power(uint64_t params, int &offset, int size = 4) {
    return 1.0 / (1 << get_int(params, offset, size));
  }
};

}  // namespace dfvs
}  // namespace algorithm
}  // namespace mog