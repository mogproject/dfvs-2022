#include "vertex_cover.hpp"

#include "../../pace-2019/lib/mis/exact_mis.h"
#include "../../pace-2019/extern/KaHIP/interface/kaHIP_interface.h"

#include "../../pace-2019/app/configuration_mis.h"  // must be after including kaHIP_interface.h

namespace mog {
namespace algorithm {

std::vector<int> solve_vertex_cover(data::graph::CLGraph const& g, int time_limit_min) {
  if (g.empty()) return {};

  std::vector<int> ret;
  MISConfig mis_config;
  configuration_mis cfg;
  cfg.standard(mis_config);
  // mis_config.ils_iterations = 0;
  if (time_limit_min > 0) mis_config.time_limit = time_limit_min;

  // PRINT("edges" << g.compact().edges());

  auto adj = g.is_compact() ? g.get_adjacency_list() : g.compact().get_adjacency_list();

  // if (static_cast<int>(adj.size()) != g.number_of_nodes()) throw std::runtime_error("wrong list size");
  // std::vector<std::vector<bool>> adjm(adj.size(), std::vector<bool>(adj.size()));
  // for (int i = 0; i < static_cast<int>(adj.size()); ++i) {
  //   for (int j: adj[i]) {
  //     if (adjm[i][j]) throw std::runtime_error("found duplicates");
  //     adjm[i][j] = true;
  //   }
  // }
  // for (int i = 0; i < static_cast<int>(adj.size()); ++i) {
  //   for (int j = 0; j < static_cast<int>(adj.size()); ++j) {
  //     if (i == j && adjm[i][j]) {
  //       // PRINT(i);
  //       // PRINT(adj[i]);
  //       throw std::runtime_error("found loop");
  //     }
  //     if (adjm[i][j] != adjm[j][i]) throw std::runtime_error("found asymmetric");
  //   }
  // }

  // PRINT("n=" << adj.size() << ", m=" << g.number_of_edges());

  auto mis = getExactMISCombined(adj, mis_config);
  for (int i = 0; i < static_cast<int>(mis.size()); ++i) {
    if (!mis[i]) ret.push_back(i);
  }

  // check if the result is really a vertex cover
  // if (static_cast<int>(mis.size()) != g.number_of_nodes()) {
  //   // throw std::runtime_error("mis result size error");
  //   return {-1};
  // }

  if (g.is_compact()) return ret;

  std::vector<int> vtcs = g.vertices();
  // auto vr = util::inverse(vtcs);
  // for (int i = 0; i < g.number_of_nodes(); ++i) {
  //   int u = vtcs[i];
  //   for (int v: g.out_neighbors(u)) {
  //     int j = vr[v];
  //     if (mis[i] && mis[j]) {
  //       // throw std::runtime_error("mis result error");
  //       return {-1};
  //     }
  //   }
  // }

  std::vector<int> converted_ret;
  for (auto i : ret) converted_ret.push_back(vtcs[i]);
  return converted_ret;
}

int solve_vertex_cover_number(data::graph::CLGraph const& g, int time_limit_min) {
  if (g.empty()) return 0;

  int ret = 0;
  MISConfig mis_config;
  configuration_mis cfg;
  cfg.standard(mis_config);

  if (time_limit_min > 0) mis_config.time_limit = time_limit_min;

  auto adj = g.is_compact() ? g.get_adjacency_list() : g.compact().get_adjacency_list();

  // if (static_cast<int>(adj.size()) != g.number_of_nodes()) throw std::runtime_error("wrong list size");
  // std::vector<std::vector<bool>> adjm(adj.size(), std::vector<bool>(adj.size()));
  // for (int i = 0; i < static_cast<int>(adj.size()); ++i) {
  //   for (int j: adj[i]) {
  //     if (adjm[i][j]) throw std::runtime_error("found duplicates");
  //     adjm[i][j] = true;
  //   }
  // }
  // for (int i = 0; i < static_cast<int>(adj.size()); ++i) {
  //   for (int j = 0; j < static_cast<int>(adj.size()); ++j) {
  //     if (i == j && adjm[i][j]) {
  //       // PRINT(i);
  //       // PRINT(adj[i]);
  //       throw std::runtime_error("found loop");
  //     }
  //     if (adjm[i][j] != adjm[j][i]) throw std::runtime_error("found asymmetric");
  //   }
  // }

  // PRINT("n=" << adj.size() << ", m=" << g.number_of_edges());

  auto mis = getExactMISCombined(adj, mis_config);
  for (int i = 0; i < static_cast<int>(mis.size()); ++i) {
    if (!mis[i]) ++ret;
  }

  // check if the result is really a vertex cover
  if (static_cast<int>(mis.size()) != g.number_of_nodes()) throw std::runtime_error("mis result size error");

  // std::vector<int> vvv = g.vertices();
  // auto vvvr = util::inverse(vvv);
  // for (int i = 0; i < g.number_of_nodes(); ++i) {
  //   int u = vvv[i];
  //   for (int v: g.out_neighbors(u)) {
  //     int j = vvvr[v];
  //     if (mis[i] && mis[j]) throw std::runtime_error("mis result error");
  //   }
  // }

  return ret;
}

std::vector<int> solve_vertex_cover(data::graph::CLDiGraph const& g, int time_limit_min) { return solve_vertex_cover(g.to_undirected(), time_limit_min); }
int solve_vertex_cover_number(data::graph::CLDiGraph const& g, int time_limit_min) { return solve_vertex_cover_number(g.to_undirected(), time_limit_min); }

}  // namespace algorithm
}  // namespace mog