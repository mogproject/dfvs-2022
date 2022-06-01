#pragma once
#include "../data/graph/CLDiGraph.hpp"
#include "../data/graph/CLGraph.hpp"
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace mog {
namespace readwrite {

typedef std::pair<data::graph::CLGraph, std::pair<std::vector<std::string>, std::vector<int>>> PACE2016Result;

/**
 * @brief Loads a graph instance in the PACE 2016 format.
 * A graph may contain loops but cannot have multi-edges.
 *
 * @param path path to the input file
 * @return pair of [undirected graph] and [a pair of label mapping, and vertices with a loop]
 */
std::pair<data::graph::CLGraph, std::pair<std::vector<std::string>, std::vector<int>>> load_pace_2016(char const *path);
std::pair<data::graph::CLGraph, std::pair<std::vector<std::string>, std::vector<int>>> read_pace_2016(std::istream &is);
data::graph::CLDiGraph load_pace_2022(char const *path);
data::graph::CLDiGraph read_pace_2022(std::istream &is);

std::ostream &write_pace_2022(data::graph::CLDiGraph const &graph, std::ostream &os, bool compact = false);
std::ostream &write_tedder_2008(data::graph::CLGraph const &graph, std::ostream &os);
void save_pace_2022(char const *path, data::graph::CLDiGraph const &graph, bool compact = false);

}  // namespace readwrite
}  // namespace mog
