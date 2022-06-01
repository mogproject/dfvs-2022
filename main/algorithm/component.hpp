#pragma once
#include "../util/profiler.hpp"
#include <vector>

namespace mog {
namespace algorithm {

template <typename Graph>
std::vector<std::vector<int>> components(Graph const& g);

template <typename DiGraph>
std::vector<std::vector<int>> strongly_connected_components(DiGraph const& g, bool include_singletons = true,
                                                            util::Profiler* prof = nullptr);
template <typename DiGraph>
std::vector<int> scc_ids(DiGraph const& g);

template <typename DiGraph>
bool is_acyclic(DiGraph const& g);

template <typename DiGraph>
std::vector<int> strong_articulation_points(DiGraph const& g);

}  // namespace algorithm
}  // namespace mog
