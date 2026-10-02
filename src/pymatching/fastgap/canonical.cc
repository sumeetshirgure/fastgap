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

#include "pymatching/fastgap/canonical.h"

#include <cmath>
#include <deque>
#include <map>

namespace pm {
namespace fastgap {

namespace {

/// One graphlike piece of an error instruction (the targets between two separators).
struct RawComponent {
    std::vector<uint64_t> dets;
    std::vector<uint64_t> observables;  // in order, possibly repeated
    bool obs_k = false;                 // parity of observable k in this component
};

std::vector<RawComponent> split_components(const stim::DemInstruction& instruction, size_t observable) {
    std::vector<RawComponent> out(1);
    for (const auto& t : instruction.target_data) {
        if (t.is_separator()) {
            out.emplace_back();
        } else if (t.is_relative_detector_id()) {
            out.back().dets.push_back(t.val());
        } else if (t.is_observable_id()) {
            out.back().observables.push_back(t.val());
            if (t.val() == observable)
                out.back().obs_k ^= true;
        }
    }
    return out;
}

std::string edge_name(uint64_t u, uint64_t v) {
    if (v == UINT64_MAX)
        return "(D" + std::to_string(u) + ", boundary)";
    return "(D" + std::to_string(u) + ", D" + std::to_string(v) + ")";
}

/// Shared parallel-edge check (CLAUDE.md §2.3 step 4). Keys are (u, v) with u < v, or (u, MAX)
/// for a boundary edge at u.
struct ParallelEdgeChecker {
    std::map<std::pair<uint64_t, uint64_t>, bool> seen;
    void add(uint64_t u, uint64_t v, bool obs_prime) {
        if (v != UINT64_MAX && u > v)
            std::swap(u, v);
        auto [it, inserted] = seen.emplace(std::make_pair(u, v), obs_prime);
        if (!inserted && it->second != obs_prime) {
            if (v == UINT64_MAX) {
                throw ParallelEdgeConflictError(
                    "Detector D" + std::to_string(u) +
                    " has boundary edges on both sides L and R of the gauge-fixed observable. This is a local "
                    "weight-2 logical (PyMatching would merge both into one boundary edge, so the code has "
                    "distance 2), which fastgap does not support.");
            }
            throw ParallelEdgeConflictError(
                "Parallel error mechanisms on edge " + edge_name(u, v) +
                " differ in the gauge-fixed observable. This is a local weight-2 logical (the code has "
                "distance 2), which fastgap does not support.");
        }
    }
};

}  // namespace

size_t GraphSpec::baseline_obs_detector() const {
    if (kind == FROM_REFERENCE_EDGES)
        return ref_obs_detector;
    return num_detectors;
}

GraphSpec canonicalize_dem(const stim::DetectorErrorModel& dem, size_t observable) {
    GraphSpec spec;
    spec.kind = GraphSpec::FROM_DEM;
    stim::DetectorErrorModel flat = dem.flattened();
    spec.num_detectors = flat.count_detectors();
    spec.num_observables = flat.count_observables();
    spec.observable = observable;
    if (observable >= spec.num_observables && spec.num_observables > 0)
        throw std::invalid_argument(
            "observable index " + std::to_string(observable) + " is out of range: the DEM has " +
            std::to_string(spec.num_observables) + " observables.");
    if (spec.num_observables == 0)
        throw std::invalid_argument("The DEM has no observables, so there is no logical gap to compute.");
    if (spec.num_observables > 64)
        throw std::invalid_argument("fastgap v1 supports at most 64 observables.");

    // Pass 1: collect the interior edges with their obs_k parity (step 1).
    size_t n = spec.num_detectors;
    std::vector<std::vector<std::pair<uint64_t, uint8_t>>> adj(n);
    for (const auto& instruction : flat.instructions) {
        if (instruction.type != stim::DemInstructionType::DEM_ERROR)
            continue;
        double p = instruction.arg_data[0];
        if (p > 0.5)
            throw std::invalid_argument(
                "fastgap v1 does not support negative edge weights, but the DEM has an error with probability " +
                std::to_string(p) + " > 0.5.");
        if (!(p > 0))
            continue;  // PyMatching ignores these
        for (auto& c : split_components(instruction, observable)) {
            if (c.dets.size() == 2) {
                if (c.dets[0] == c.dets[1])
                    throw std::invalid_argument(
                        "Error component with a repeated detector D" + std::to_string(c.dets[0]) +
                        " (a self-loop) is not supported.");
                adj[c.dets[0]].push_back({c.dets[1], c.obs_k});
                adj[c.dets[1]].push_back({c.dets[0], c.obs_k});
            } else if (c.dets.size() > 2) {
                spec.num_hyperedges_ignored++;
            }
        }
    }

    // BFS 2-colouring so that obs_k(u, v) = s(u) XOR s(v) on every interior edge (steps 1 and 2).
    spec.gauge.assign(n, 0);
    std::vector<uint8_t> visited(n, 0);
    std::deque<uint64_t> queue;
    for (uint64_t root = 0; root < n; root++) {
        if (visited[root])
            continue;
        visited[root] = 1;
        queue.push_back(root);
        while (!queue.empty()) {
            uint64_t x = queue.front();
            queue.pop_front();
            for (auto [y, bit] : adj[x]) {
                uint8_t want = spec.gauge[x] ^ bit;
                if (!visited[y]) {
                    visited[y] = 1;
                    spec.gauge[y] = want;
                    queue.push_back(y);
                } else if (spec.gauge[y] != want) {
                    throw BoundarylessLogicalError(
                        "Observable L" + std::to_string(observable) +
                        " can be flipped by a closed loop of interior errors that never touches a boundary (an "
                        "odd-parity cycle through D" +
                        std::to_string(x) + " and D" + std::to_string(y) +
                        "). The code has a logical that is not a boundary-to-boundary chain (e.g. a toroidal "
                        "code), which fastgap v1 does not support.");
                }
            }
        }
    }

    // Pass 2: relabel (step 3), check parallel edges (step 4), and write the canonical DEM.
    ParallelEdgeChecker checker;
    std::vector<stim::DemTarget> targets;
    for (const auto& instruction : flat.instructions) {
        if (instruction.type != stim::DemInstructionType::DEM_ERROR) {
            spec.canonical_dem.append_dem_instruction(instruction);
            continue;
        }
        double p = instruction.arg_data[0];
        auto components = split_components(instruction, observable);
        targets.clear();
        for (size_t ci = 0; ci < components.size(); ci++) {
            auto& c = components[ci];
            if (ci > 0)
                targets.push_back(stim::DemTarget::separator());
            for (auto d : c.dets)
                targets.push_back(stim::DemTarget::relative_detector_id(d));
            bool graphlike = c.dets.size() == 1 || c.dets.size() == 2;
            bool obs_prime = c.obs_k;
            if (graphlike) {
                // Interior edges get obs'_k = 0; a boundary edge at u gets obs_k(u, b) XOR s(u).
                for (auto d : c.dets)
                    obs_prime ^= (bool)spec.gauge[d];
                if (p > 0) {
                    if (c.dets.size() == 2) {
                        if (obs_prime)
                            throw std::logic_error("gauge fixing left observable k on an interior edge");
                        checker.add(c.dets[0], c.dets[1], false);
                    } else {
                        checker.add(c.dets[0], UINT64_MAX, obs_prime);
                    }
                }
            }
            // Components with 0 or > 2 detectors are dropped by PyMatching's decoder; their labels
            // are left untouched (a detector-free component has no gauge shift anyway).
            for (auto o : c.observables) {
                if (o != observable)
                    targets.push_back(stim::DemTarget::observable_id(o));
            }
            if (obs_prime)
                targets.push_back(stim::DemTarget::observable_id(observable));
        }
        spec.canonical_dem.append_error_instruction(p, targets, instruction.tag);
    }

    if (spec.canonical_dem.count_detectors() != spec.num_detectors)
        throw std::logic_error("canonical DEM changed the number of detectors");
    if (spec.canonical_dem.count_observables() < spec.num_observables)
        spec.canonical_dem.append_logical_observable_instruction(
            stim::DemTarget::observable_id(spec.num_observables - 1), "");
    return spec;
}

GraphSpec spec_from_reference_edges(size_t num_detectors, size_t obs_detector, std::vector<RefEdge> edges) {
    GraphSpec spec;
    spec.kind = GraphSpec::FROM_REFERENCE_EDGES;
    spec.num_detectors = num_detectors;
    spec.num_observables = 1;
    spec.observable = 0;
    spec.gauge.assign(num_detectors, 0);
    spec.ref_obs_detector = obs_detector;
    if (obs_detector >= num_detectors)
        throw std::invalid_argument("obs_detector must be a detector of the reference graph.");
    ParallelEdgeChecker checker;
    for (auto& e : edges) {
        if (e.u == SIZE_MAX && e.v == SIZE_MAX)
            continue;
        if (e.u >= num_detectors || (e.v != SIZE_MAX && e.v >= num_detectors))
            throw std::invalid_argument("reference edge refers to a detector outside the graph");
        if (e.u == obs_detector || e.v == obs_detector)
            throw std::invalid_argument(
                "reference edges must already be converted to the L/R model (no edge may touch obs_detector)");
        if (e.weight < 0)
            throw std::invalid_argument("fastgap v1 does not support negative edge weights.");
        if (e.v == SIZE_MAX) {
            checker.add(e.u, UINT64_MAX, e.side_l);
        } else {
            checker.add(e.u, e.v, false);
        }
    }
    spec.ref_edges = std::move(edges);
    return spec;
}

stim::DetectorErrorModel make_baseline_dem(const GraphSpec& spec) {
    if (spec.kind != GraphSpec::FROM_DEM)
        throw std::invalid_argument("make_baseline_dem requires a DEM-derived GraphSpec");
    size_t obs_det = spec.baseline_obs_detector();
    stim::DetectorErrorModel out;
    std::vector<stim::DemTarget> targets;
    for (const auto& instruction : spec.canonical_dem.instructions) {
        if (instruction.type == stim::DemInstructionType::DEM_DETECTOR) {
            out.append_dem_instruction(instruction);
            continue;
        }
        if (instruction.type != stim::DemInstructionType::DEM_ERROR)
            continue;  // logical_observable declarations are dropped with the observables
        targets.clear();
        bool first = true;
        for (auto& c : split_components(instruction, spec.observable)) {
            if (!first)
                targets.push_back(stim::DemTarget::separator());
            first = false;
            for (auto d : c.dets)
                targets.push_back(stim::DemTarget::relative_detector_id(d));
            if (c.obs_k)
                targets.push_back(stim::DemTarget::relative_detector_id(obs_det));
        }
        out.append_error_instruction(instruction.arg_data[0], targets, instruction.tag);
    }
    std::vector<double> no_coords;
    out.append_detector_instruction(no_coords, stim::DemTarget::relative_detector_id(obs_det), "");
    return out;
}

pm::UserGraph make_decode_user_graph(const GraphSpec& spec) {
    if (spec.kind == GraphSpec::FROM_DEM)
        return pm::detector_error_model_to_user_graph(spec.canonical_dem, false, pm::NUM_DISTINCT_WEIGHTS);
    pm::UserGraph g(spec.num_detectors, 1);
    for (auto& e : spec.ref_edges) {
        if (e.u == SIZE_MAX) {
            // Keep the detector-free L-R edge where the reference graph has it (a boundary edge of
            // the observable detector, which is never a defect here), so both graphs share the
            // same heaviest edge and therefore the same internal units.
            g.add_or_merge_boundary_edge(spec.ref_obs_detector, {}, e.weight, e.error_probability, pm::INDEPENDENT);
            continue;
        }
        if (e.v == SIZE_MAX) {
            std::vector<size_t> obs;
            if (e.side_l)
                obs.push_back(0);
            g.add_or_merge_boundary_edge(e.u, obs, e.weight, e.error_probability, pm::INDEPENDENT);
        } else {
            g.add_or_merge_edge(e.u, e.v, {}, e.weight, e.error_probability, pm::INDEPENDENT);
        }
    }
    g.loaded_from_dem_without_correlations = true;
    return g;
}

pm::UserGraph make_baseline_user_graph(const GraphSpec& spec) {
    if (spec.kind == GraphSpec::FROM_DEM)
        return pm::detector_error_model_to_user_graph(make_baseline_dem(spec), false, pm::NUM_DISTINCT_WEIGHTS);
    size_t obs_det = spec.baseline_obs_detector();
    pm::UserGraph g(spec.num_detectors, 0);
    for (auto& e : spec.ref_edges) {
        if (e.u == SIZE_MAX) {
            g.add_or_merge_boundary_edge(obs_det, {}, e.weight, e.error_probability, pm::INDEPENDENT);
        } else if (e.v == SIZE_MAX) {
            if (e.side_l) {
                g.add_or_merge_edge(e.u, obs_det, {}, e.weight, e.error_probability, pm::INDEPENDENT);
            } else {
                g.add_or_merge_boundary_edge(e.u, {}, e.weight, e.error_probability, pm::INDEPENDENT);
            }
        } else {
            g.add_or_merge_edge(e.u, e.v, {}, e.weight, e.error_probability, pm::INDEPENDENT);
        }
    }
    g.loaded_from_dem_without_correlations = true;
    return g;
}

uint8_t gauge_shift(const GraphSpec& spec, const std::vector<uint64_t>& detection_events) {
    uint8_t s = 0;
    for (auto d : detection_events) {
        if (d < spec.gauge.size())
            s ^= spec.gauge[d];
    }
    return s;
}

}  // namespace fastgap
}  // namespace pm
