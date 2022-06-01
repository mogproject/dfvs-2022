#include <algorithm>

#include "../util/util.hpp"
#include "c_pace.hpp"
#include <cctype>
#include <fstream>
#include <map>
#include <sstream>

namespace mog {
namespace readwrite {

PACE2016Result read_pace_2016(std::istream &is) {
  // label map utilities
  std::map<std::string, size_t> label_map;
  auto f = [&](std::string const &s) {  // get index with update
    auto p = label_map.find(s);
    if (p == label_map.end()) {
      auto idx = label_map.size();
      label_map[s] = idx;
      return idx;
    }
    return p->second;
  };

  // parse edges
  std::vector<std::pair<int, int>> edges;
  std::vector<int> loops;

  for (std::string line; std::getline(is, line);) {
    if (line.empty()) continue;

    auto ss = std::stringstream(line);
    std::string u, v;
    ss >> u >> v;
    auto uu = f(u);
    auto vv = f(v);
    if (uu == vv) {
      loops.push_back(uu);
    } else {
      edges.push_back({uu, vv});
    }
  }

  // construct graph without loops
  int n = label_map.size();
  auto G = data::graph::CLGraph(n);
  for (auto &e : edges) G.add_edge(e.first, e.second);

  // construct label mapping
  std::vector<std::string> labels(n);
  for (auto &p : label_map) labels[p.second] = p.first;

  return {G, {labels, loops}};
}

data::graph::CLDiGraph read_pace_2022(std::istream &is) {
  int lineno = -1;
  int n = -1;
  long long m = -1;
  std::vector<std::pair<int, int>> edges;

  for (std::string line; std::getline(is, line);) {
    int i = 0;
    while (i < static_cast<int>(line.size()) && !std::isprint(line[i])) ++i;
    if (line[i] == '%') continue;  // ignore comments

    auto ss = std::stringstream(line);

    if (lineno == -1) {
      // header: number of vertices, number of edges, identifier (should be 0)
      int t;
      ss >> n >> m >> t;
    } else {
      // adjacency list (1-indexed)
      int v;
      while (ss >> v) { edges.push_back({lineno, v - 1}); }
    }
    ++lineno;
  }

  // construct digraph
  return data::graph::CLDiGraph(n, edges);
}

PACE2016Result load_pace_2016(char const *path) {
  std::ifstream f(path);
  if (f.fail()) throw std::invalid_argument(util::format("Failed to open file: %s", path));
  return read_pace_2016(f);
}

data::graph::CLDiGraph load_pace_2022(char const *path) {
  std::ifstream f(path);
  if (f.fail()) throw std::invalid_argument(util::format("Failed to open file: %s", path));
  return read_pace_2022(f);
}

std::ostream &write_pace_2022(data::graph::CLDiGraph const &graph, std::ostream &os, bool compact) {
  auto nn = graph.capacity();
  auto n = graph.number_of_nodes();

  std::vector<int> vtcs;
  std::unordered_map<int, int> mapping;

  if (compact) {
    nn = n;
    vtcs = graph.vertices();
    std::sort(vtcs.begin(), vtcs.end());
    for (int i = 0; i < n; ++i) mapping[vtcs[i]] = i;
  }

  auto m = graph.number_of_edges();
  os << nn << " " << m << " 0" << std::endl;

  for (int i = 0; i < nn; ++i) {
    auto ii = compact ? vtcs[i] : i;
    if (graph.has_vertex(ii)) {
      bool first = true;
      for (int j : graph.out_neighbors(ii)) {
        if (!first) os << " ";
        os << ((compact ? mapping[j] : j) + 1);
        first = false;
      }
    }
    os << std::endl;
  }

  return os;
}

std::ostream &write_tedder_2008(data::graph::CLGraph const &graph, std::ostream &os) {
  auto n = graph.number_of_nodes();
  for (int i = 0; i < n; ++i) {
    os << i << "->";
    bool first = true;
    for (int j : graph.neighbors(i)) {
      if (!first) os << ",";
      os << j;
      first = false;
    }
    os << std::endl;
  }
  return os;
}

void save_pace_2022(char const *path, data::graph::CLDiGraph const &graph, bool compact) {
  std::ofstream f(path);
  if (f.fail()) throw std::invalid_argument(util::format("Failed to open file: %s", path));
  write_pace_2022(graph, f, compact);
}

}  // namespace readwrite
}  // namespace mog
