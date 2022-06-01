#pragma once
#include <vector>

#include "../data/fast_set.hpp"
#include "../data/graph/graph.hpp"

namespace mog {
namespace algorithm {

template <typename Graph>
std::vector<std::vector<int>> clique_partition(Graph &g);

template <typename Graph>
std::vector<std::vector<int>> clique_partition(Graph const&g);

}  // namespace algorithm
}  // namespace mog
