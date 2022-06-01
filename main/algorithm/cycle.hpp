#pragma once
#include <vector>

#include "../data/graph/graph.hpp"

namespace mog {
namespace algorithm {

template <typename Graph>
std::vector<int> short_cycle(Graph const& g, int sample = -1, int start = 0);

/**
 * @brief Finds a shortest cycle passing through the given vertex.
 * 
 * @tparam Graph 
 * @param g graph
 * @param v vertex
 * @return list of vertices in the cycle in order
 */
template <typename Graph>
std::vector<int> find_cycle(Graph const& g, int v);

}  // namespace algorithm
}  // namespace mog
