// Copyright 2026 fastgap contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PYMATCHING2_FASTGAP_INTERIOR_GRAPH_H
#define PYMATCHING2_FASTGAP_INTERIOR_GRAPH_H

#include <functional>
#include <vector>

#include "pymatching/fastgap/types.h"
#include "pymatching/sparse_blossom/flooder/graph.h"

namespace pm {
namespace fastgap {

/// The detector graph G in internal integer units, read directly from the MatchingGraph that the
/// decoder floods, so distances use exactly the decoder's edges and weights. Interior edges are
/// stored in CSR form; the (at most one) boundary edge of each detector is stored separately with
/// its side L or R (CLAUDE.md §2.3).
struct InteriorGraph {
    size_t num_nodes = 0;
    std::vector<uint32_t> offsets;  // size num_nodes + 1
    std::vector<uint32_t> neighbors;
    std::vector<dist_int> weights;
    std::vector<dist_int> boundary_weight;  // DIST_INF if no boundary edge
    std::vector<uint8_t> boundary_side;     // SIDE_NONE / SIDE_L / SIDE_R
    double normalising_constant = 0;

    /// Reads the graph. `observable` is the (gauge-fixed) observable whose bit labels side L.
    /// Throws if an interior edge carries that observable (canonicalisation failed).
    static InteriorGraph from_matching_graph(const pm::MatchingGraph& graph, size_t observable);

    /// 64-bit fingerprint of edges, weights, sides and the normalising constant.
    uint64_t fingerprint() const;

    /// Multi-source Dijkstra over interior edges only (never through the boundary). `init` gives
    /// (node, starting distance) pairs. Nodes farther than `cap` are left at DIST_INF.
    std::vector<dist_int> dijkstra(const std::vector<std::pair<uint32_t, dist_int>>& init, dist_int cap = DIST_INF)
        const;

    /// h_X(x) = dist(x, X) for side X: shortest path from x ending with a boundary edge on side X.
    std::vector<dist_int> side_potential(uint8_t side) const;

    /// Exact interior distance from a to b with A*; `heuristic(x)` must be a consistent lower bound
    /// on dist(x, b) (e.g. a landmark bound, which is 1-Lipschitz). Uses `scratch` as workspace.
    struct AStarScratch {
        std::vector<dist_int> g;
        std::vector<uint32_t> touched;
    };
    dist_int astar(
        uint32_t a, uint32_t b, const std::function<dist_int(uint32_t)>& heuristic, AStarScratch& scratch) const;

    /// Median interior-or-boundary edge weight (internal units).
    dist_int median_edge_weight() const;
};

}  // namespace fastgap
}  // namespace pm

#endif  // PYMATCHING2_FASTGAP_INTERIOR_GRAPH_H
