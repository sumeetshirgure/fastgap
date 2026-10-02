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

#ifndef PYMATCHING2_FASTGAP_CANONICAL_H
#define PYMATCHING2_FASTGAP_CANONICAL_H

#include <cstdint>
#include <string>
#include <vector>

#include "pymatching/sparse_blossom/driver/user_graph.h"
#include "stim.h"

namespace pm {
namespace fastgap {

/// Raised when observable k can be flipped by a closed loop of interior edges (CLAUDE.md §2.3 step 2).
struct BoundarylessLogicalError : std::invalid_argument {
    explicit BoundarylessLogicalError(const std::string& what) : std::invalid_argument(what) {
    }
};

/// Raised when two raw mechanisms on the same detectors differ in obs'_k (CLAUDE.md §2.3 step 4).
struct ParallelEdgeConflictError : std::invalid_argument {
    explicit ParallelEdgeConflictError(const std::string& what) : std::invalid_argument(what) {
    }
};

/// One edge of an externally supplied reference graph (CLAUDE.md §7.3), already in the L/R model.
/// `u == SIZE_MAX && v == SIZE_MAX` is a detector-free L–R mechanism (an edge from the observable
/// detector straight to the boundary).
struct RefEdge {
    size_t u;
    size_t v;  // SIZE_MAX means a boundary edge at u
    bool side_l;  // only meaningful for boundary edges
    double weight;
    double error_probability;
};

/// The canonicalised description of one (DEM, observable) pair. Every consumer (decode graph,
/// index, exact baseline) is built from this one object, so they all see identical edges,
/// weights and internal units.
struct GraphSpec {
    enum Kind : uint8_t { FROM_DEM = 0, FROM_REFERENCE_EDGES = 1 };
    Kind kind = FROM_DEM;
    size_t num_detectors = 0;
    /// Number of observables of the decode graph (the size of the prediction vector).
    size_t num_observables = 0;
    /// The observable index k whose gap is computed.
    size_t observable = 0;
    /// Gauge bits s(x): obs'_k(E) = obs_k(E) XOR (XOR_{x in D} s(x)).
    std::vector<uint8_t> gauge;

    /// FROM_DEM: the canonical DEM (observable k only on boundary edges).
    stim::DetectorErrorModel canonical_dem;
    /// FROM_REFERENCE_EDGES: the reference edges, plus the detector that played the role of the
    /// observable in the reference DEM. Detection events on that detector are ignored.
    std::vector<RefEdge> ref_edges;
    size_t ref_obs_detector = SIZE_MAX;

    /// Number of raw components with more than two detectors (PyMatching drops these).
    size_t num_hyperedges_ignored = 0;

    /// Index of the extra detector that represents observable k in the exact baseline graph.
    size_t baseline_obs_detector() const;
};

/// Canonicalise observable `observable` of `dem` onto the boundary (CLAUDE.md §2.3).
/// Works on the raw (flattened) error list, not on PyMatching's merged graph.
GraphSpec canonicalize_dem(const stim::DetectorErrorModel& dem, size_t observable);

/// Build a GraphSpec from the edges PyMatching keeps for a reference DEM (CLAUDE.md §7.3).
/// Runs the parallel-edge check. Gauge bits are all zero.
GraphSpec spec_from_reference_edges(size_t num_detectors, size_t obs_detector, std::vector<RefEdge> edges);

/// The DEM used by the exact two-decode baseline: observable k becomes detector
/// `spec.baseline_obs_detector()`; all observables are dropped. FROM_DEM only.
stim::DetectorErrorModel make_baseline_dem(const GraphSpec& spec);

/// The UserGraph fastgap decodes on (built exactly as upstream builds it from a DEM).
pm::UserGraph make_decode_user_graph(const GraphSpec& spec);

/// The UserGraph of the exact two-decode baseline.
pm::UserGraph make_baseline_user_graph(const GraphSpec& spec);

/// The gauge shift XOR_{x in D} s(x) for one shot.
uint8_t gauge_shift(const GraphSpec& spec, const std::vector<uint64_t>& detection_events);

}  // namespace fastgap
}  // namespace pm

#endif  // PYMATCHING2_FASTGAP_CANONICAL_H
