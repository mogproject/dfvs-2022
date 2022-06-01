#include "algorithm/dfvs/ExactSolver.hpp"
#include "readwrite/c_pace.hpp"

using namespace std;
using namespace mog::readwrite;
using namespace mog::algorithm::dfvs;

/**
 * See Configuration.hpp for the flags.
 */
int main(int argc, char* argv[]) {
  // set up configuration
  uint64_t option = Configuration::default_params();
  uint32_t seed = 12345;
  Configuration conf(option, seed);

  // load graph
  auto g = read_pace_2022(std::cin);

  // run algorithm
  auto solver = ExactSolver(conf);
  std::vector<int> ret;

  while (true) {
    ret = solver.solve(g);
    if (solver.is_solution(g, ret)) break;
  }

  // output result
  for (auto r : ret) printf("%d\n", r + 1);

  return 0;
}
