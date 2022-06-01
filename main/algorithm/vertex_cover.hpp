#pragma once
#include <vector>

#include "../data/graph/graph.hpp"

namespace mog {
namespace algorithm {

std::vector<int> solve_vertex_cover(data::graph::CLGraph const& g, int time_limit_min=0);
std::vector<int> solve_vertex_cover(data::graph::CLDiGraph const& g, int time_limit_min=0);
int solve_vertex_cover_number(data::graph::CLGraph const& g, int time_limit_min=0);
int solve_vertex_cover_number(data::graph::CLDiGraph const& g, int time_limit_min=0);

}  // namespace algorithm
}  // namespace mog
