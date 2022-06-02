#pragma once

#include <cassert>
#include <random>

#include "../../component.hpp"
#include "../Configuration.hpp"
#include "../LabeledGraph.hpp"
#include "../reduction/Reducer.hpp"

namespace mog {
namespace algorithm {
namespace dfvs {

class Brancher {
 protected:
  Configuration conf_;
  util::Profiler* const prof_;

 private:
  std::default_random_engine rand_gen_;
  std::uniform_real_distribution<> real_dist_;
  std::unordered_set<int> guess_solution_;

 public:
  Brancher(Configuration const& conf, uint32_t seed)
      : conf_(conf),       //
        prof_(conf.prof),  //
        rand_gen_(seed),   //
        real_dist_(0, 1)   //
  {
    //
  }
  double rand() { return real_dist_(rand_gen_); }
  virtual void clear_record() = 0;
  virtual void find_target(LabeledGraph& h, bool use_record, bool keep_record) = 0;
  virtual void branch_left(LabeledGraph& h, std::vector<int>& part_sol, bool record_history) = 0;
  virtual void branch_right(LabeledGraph& h, std::vector<int>& part_sol, bool record_history) = 0;
  virtual int get_u() const = 0;
  virtual int get_v() const = 0;
  virtual bool get_include() const = 0;
  virtual void set_u(int u) = 0;
  virtual void set_v(int v) = 0;
  virtual void set_include(bool should_include) = 0;
  virtual std::string to_string() const = 0;
  void set_guess_solution(std::vector<int> const& guess_solution) {
    guess_solution_ = std::unordered_set<int>(guess_solution.begin(), guess_solution.end());
  }
  bool is_in_guess_solution(int v) const { return util::contains(guess_solution_, v); }
};

}  // namespace dfvs
}  // namespace algorithm
}  // namespace mog