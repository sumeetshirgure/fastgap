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

#include "pymatching/fastgap/interior_graph.h"

#include <algorithm>
#include <cstring>
#include <queue>

namespace pm {
namespace fastgap {

InteriorGraph InteriorGraph::from_matching_graph(const pm::MatchingGraph& graph, size_t observable) {
    InteriorGraph g;
    g.num_nodes = graph.nodes.size();
    g.normalising_constant = graph.normalising_constant;
    g.offsets.assign(g.num_nodes + 1, 0);
    g.boundary_weight.assign(g.num_nodes, DIST_INF);
    g.boundary_side.assign(g.num_nodes, SIDE_NONE);
    pm::obs_int mask = (pm::obs_int)1 << observable;
    const pm::DetectorNode* base = graph.nodes.data();
    for (size_t i = 0; i < g.num_nodes; i++) {
        const auto& node = graph.nodes[i];
        for (size_t j = 0; j < node.neighbors.size(); j++) {
            dist_int w = (dist_int)node.neighbor_weights[j];
            if (node.neighbors[j] == nullptr) {
                if (g.boundary_side[i] != SIDE_NONE)
                    throw std::logic_error("detector has two boundary edges in the matching graph");
                g.boundary_weight[i] = w;
                g.boundary_side[i] = (node.neighbor_observables[j] & mask) ? SIDE_L : SIDE_R;
            } else {
                if (node.neighbor_observables[j] & mask)
                    throw std::logic_error(
                        "interior edge carries the gauge-fixed observable; canonicalisation was not applied");
                g.neighbors.push_back((uint32_t)(node.neighbors[j] - base));
                g.weights.push_back(w);
            }
        }
        g.offsets[i + 1] = (uint32_t)g.neighbors.size();
    }
    return g;
}

uint64_t InteriorGraph::fingerprint() const {
    // FNV-1a over the raw arrays.
    uint64_t h = 1469598103934665603ULL;
    auto mix = [&](const void* data, size_t len) {
        auto p = (const uint8_t*)data;
        for (size_t i = 0; i < len; i++) {
            h ^= p[i];
            h *= 1099511628211ULL;
        }
    };
    mix(&num_nodes, sizeof(num_nodes));
    mix(offsets.data(), offsets.size() * sizeof(uint32_t));
    mix(neighbors.data(), neighbors.size() * sizeof(uint32_t));
    mix(weights.data(), weights.size() * sizeof(dist_int));
    mix(boundary_weight.data(), boundary_weight.size() * sizeof(dist_int));
    mix(boundary_side.data(), boundary_side.size());
    mix(&normalising_constant, sizeof(normalising_constant));
    return h;
}

std::vector<dist_int> InteriorGraph::dijkstra(
    const std::vector<std::pair<uint32_t, dist_int>>& init, dist_int cap) const {
    std::vector<dist_int> dist(num_nodes, DIST_INF);
    typedef std::pair<dist_int, uint32_t> item;
    std::priority_queue<item, std::vector<item>, std::greater<item>> pq;
    for (auto [v, d] : init) {
        if (v >= num_nodes)
            throw std::invalid_argument("dijkstra source out of range");
        if (d <= cap && d < dist[v]) {
            dist[v] = d;
            pq.push({d, v});
        }
    }
    while (!pq.empty()) {
        auto [d, x] = pq.top();
        pq.pop();
        if (d != dist[x])
            continue;
        for (uint32_t e = offsets[x]; e < offsets[x + 1]; e++) {
            dist_int nd = d + weights[e];
            uint32_t y = neighbors[e];
            if (nd <= cap && nd < dist[y]) {
                dist[y] = nd;
                pq.push({nd, y});
            }
        }
    }
    return dist;
}

std::vector<dist_int> InteriorGraph::side_potential(uint8_t side) const {
    std::vector<std::pair<uint32_t, dist_int>> init;
    for (uint32_t i = 0; i < num_nodes; i++) {
        if (boundary_side[i] == side)
            init.push_back({i, boundary_weight[i]});
    }
    return dijkstra(init);
}

dist_int InteriorGraph::astar(
    uint32_t a, uint32_t b, const std::function<dist_int(uint32_t)>& heuristic, AStarScratch& scratch) const {
    if (a == b)
        return 0;
    if (scratch.g.size() != num_nodes)
        scratch.g.assign(num_nodes, DIST_INF);
    typedef std::pair<dist_int, uint32_t> item;  // (g + h, node)
    std::priority_queue<item, std::vector<item>, std::greater<item>> pq;
    scratch.g[a] = 0;
    scratch.touched.push_back(a);
    pq.push({heuristic(a), a});
    dist_int result = DIST_INF;
    while (!pq.empty()) {
        auto [f, x] = pq.top();
        pq.pop();
        dist_int gx = scratch.g[x];
        if (x == b) {
            result = gx;
            break;
        }
        dist_int hx = heuristic(x);
        if (f != gx + hx)
            continue;  // stale entry
        for (uint32_t e = offsets[x]; e < offsets[x + 1]; e++) {
            uint32_t y = neighbors[e];
            dist_int ng = gx + weights[e];
            if (ng < scratch.g[y]) {
                if (is_inf(scratch.g[y]))
                    scratch.touched.push_back(y);
                scratch.g[y] = ng;
                pq.push({ng + heuristic(y), y});
            }
        }
    }
    for (auto t : scratch.touched)
        scratch.g[t] = DIST_INF;
    scratch.touched.clear();
    return result;
}

dist_int InteriorGraph::median_edge_weight() const {
    std::vector<dist_int> all(weights.begin(), weights.end());
    for (auto w : boundary_weight) {
        if (!is_inf(w))
            all.push_back(w);
    }
    if (all.empty())
        return 0;
    std::nth_element(all.begin(), all.begin() + all.size() / 2, all.end());
    return all[all.size() / 2];
}

}  // namespace fastgap
}  // namespace pm
