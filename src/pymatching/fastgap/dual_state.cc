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

#include "pymatching/fastgap/dual_state.h"

#include <algorithm>

#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"

namespace pm {
namespace fastgap {

void DualState::clear() {
    defects.clear();
    y_leaf.clear();
    rho.clear();
    anc_offset.assign(1, 0);
    anc_blossom.clear();
    anc_prefix.clear();
    y_blossom.clear();
    blossom_depth.clear();
    mate.clear();
    w_star = 0;
    obs_mask = 0;
}

namespace {

dist_int frozen_radius(const pm::GraphFillRegion* region, pm::cumulative_time_int now) {
    if (!region->radius.is_frozen())
        throw InvariantViolation("region is not frozen after sparse blossom finished");
    dist_int y = region->radius.get_distance_at_time(now);
    if (y != region->radius.y_intercept())
        throw InvariantViolation("frozen region radius differs from its y-intercept");
    if (y < 0)
        throw InvariantViolation("negative region radius");
    return y;
}

}  // namespace

void extract_dual_state(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    size_t observable,
    DualState& out,
    DualStateScratch& scratch) {
    if (!mwpm.flooder.negative_weight_detection_events.empty() || mwpm.flooder.negative_weight_sum != 0)
        throw std::invalid_argument("fastgap v1 does not support graphs with negative edge weights.");

    out.clear();
    auto& graph = mwpm.flooder.graph;
    if (scratch.defect_of_node.size() != graph.nodes.size())
        scratch.defect_of_node.assign(graph.nodes.size(), -1);
    scratch.blossom_ids.clear();
    auto release = [&]() {
        for (auto d : out.defects)
            scratch.defect_of_node[d] = -1;
    };

    // Validate and index the defects before flooding, so that a bad shot never leaves the Mwpm
    // half-processed.
    for (auto d : detection_events) {
        if (d >= graph.nodes.size()) {
            release();
            throw std::invalid_argument(
                "The detection event with index " + std::to_string(d) +
                " does not correspond to a node in the graph, which only has " + std::to_string(graph.nodes.size()) +
                " nodes.");
        }
        if (!(d + 1 > graph.is_user_graph_boundary_node.size() || !graph.is_user_graph_boundary_node[d]))
            continue;  // upstream skips detection events on user-declared boundary nodes
        if (scratch.defect_of_node[d] != -1) {
            release();
            throw std::invalid_argument("duplicate detection event D" + std::to_string(d));
        }
        scratch.defect_of_node[d] = (int32_t)out.defects.size();
        out.defects.push_back((uint32_t)d);
    }
    try {
        process_timeline_until_completion(mwpm, detection_events);
    } catch (...) {
        release();
        throw;
    }
    pm::cumulative_time_int now = mwpm.flooder.queue.cur_time;

    // Walk each defect's leaf -> blossom_parent chain before anything is shattered.
    dist_int total = 0;
    for (size_t i = 0; i < out.defects.size(); i++) {
        uint32_t d = out.defects[i];
        const pm::DetectorNode& node = graph.nodes[d];
        const pm::GraphFillRegion* leaf = node.region_that_arrived;
        // The source node of a detection event stays owned by its own leaf region for the whole
        // decode; the checks below fail loudly if that ever stops being true.
        if (leaf == nullptr || !leaf->blossom_children.empty() || node.reached_from_source != &node)
            throw InvariantViolation("defect node is not owned by its own leaf region after flooding");
        dist_int y = frozen_radius(leaf, now);
        out.y_leaf.push_back(y);
        total += y;

        scratch.chain.clear();
        for (auto r = leaf->blossom_parent; r != nullptr; r = r->blossom_parent)
            scratch.chain.push_back(r);
        const pm::GraphFillRegion* top = scratch.chain.empty() ? leaf : scratch.chain.back();
        if (node.region_that_arrived_top != top)
            throw InvariantViolation("region_that_arrived_top disagrees with the blossom_parent chain");

        dist_int rho = y;
        for (size_t k = scratch.chain.size(); k-- > 0;) {
            const pm::GraphFillRegion* b = scratch.chain[k];
            auto [it, inserted] = scratch.blossom_ids.emplace(b, (int32_t)out.y_blossom.size());
            int32_t depth = (int32_t)(scratch.chain.size() - 1 - k);
            if (inserted) {
                dist_int yb = frozen_radius(b, now);
                out.y_blossom.push_back(yb);
                out.blossom_depth.push_back(depth);
                total += yb;
            } else if (out.blossom_depth[it->second] != depth) {
                throw InvariantViolation("inconsistent blossom depth");
            }
            rho += out.y_blossom[it->second];
            out.anc_blossom.push_back(it->second);
            out.anc_prefix.push_back(rho - y);
        }
        out.anc_offset.push_back((uint32_t)out.anc_blossom.size());
        out.rho.push_back(rho);
    }
    out.w_star = total;

    // Shatter with the upstream routine to get the matching.
    mwpm.flooder.match_edges.clear();
    shatter_blossoms_for_all_detection_events_and_extract_match_edges(mwpm, detection_events);

    size_t n = out.defects.size();
    out.mate.assign(n, INT32_MIN);
    pm::obs_int bit = (pm::obs_int)1 << observable;
    const pm::DetectorNode* base = graph.nodes.data();
    for (const auto& e : mwpm.flooder.match_edges) {
        int32_t i = scratch.defect_of_node[e.loc_from - base];
        if (i < 0 || out.mate[i] != INT32_MIN)
            throw InvariantViolation("match edge does not start at an unmatched defect");
        out.obs_mask ^= e.obs_mask;
        if (e.loc_to == nullptr) {
            out.mate[i] = (e.obs_mask & bit) ? MATE_L : MATE_R;
        } else {
            int32_t j = scratch.defect_of_node[e.loc_to - base];
            if (j < 0 || out.mate[j] != INT32_MIN)
                throw InvariantViolation("match edge does not end at an unmatched defect");
            if (e.obs_mask & bit)
                throw InvariantViolation("defect-defect relay crosses the gauge-fixed observable");
            out.mate[i] = j;
            out.mate[j] = i;
        }
    }
    release();
    for (size_t i = 0; i < n; i++) {
        if (out.mate[i] == INT32_MIN)
            throw InvariantViolation("defect left unmatched after shattering");
    }
}

}  // namespace fastgap
}  // namespace pm
