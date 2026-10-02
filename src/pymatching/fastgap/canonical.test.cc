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

#include <gtest/gtest.h>

#include "pymatching/fastgap/fastgap_test_util.test.h"
#include "pymatching/fastgap/interior_graph.h"

using namespace pm::fastgap;

namespace {

/// Checks that the canonical DEM has observable k only on boundary components.
void expect_observable_only_on_boundary(const GraphSpec& spec) {
    spec.canonical_dem.iter_flatten_error_instructions([&](const stim::DemInstruction& e) {
        size_t dets = 0;
        bool obs = false;
        for (auto& t : e.target_data) {
            if (t.is_separator()) {
                EXPECT_FALSE(dets == 2 && obs);
                dets = 0;
                obs = false;
            } else if (t.is_relative_detector_id()) {
                dets++;
            } else if (t.is_observable_id() && t.val() == spec.observable) {
                obs = !obs;
            }
        }
        EXPECT_FALSE(dets == 2 && obs);
    });
}

}  // namespace

TEST(FastgapCanonical, GaugeMovesObservableToBoundary) {
    std::mt19937_64 rng(5);
    for (int trial = 0; trial < 20; trial++) {
        auto dem = testutil::random_grid_dem(rng, 4 + trial % 3, 3 + trial % 2);
        auto spec = canonicalize_dem(dem, 0);
        expect_observable_only_on_boundary(spec);
        EXPECT_EQ(spec.num_detectors, dem.count_detectors());
        EXPECT_EQ(spec.canonical_dem.count_detectors(), dem.count_detectors());

        // The canonical DEM decodes on exactly the same graph as the original, apart from labels.
        auto g1 = pm::detector_error_model_to_user_graph(dem, false, pm::NUM_DISTINCT_WEIGHTS);
        auto g2 = make_decode_user_graph(spec);
        auto m1 = g1.to_matching_graph(pm::NUM_DISTINCT_WEIGHTS);
        auto m2 = g2.to_matching_graph(pm::NUM_DISTINCT_WEIGHTS);
        ASSERT_EQ(m1.normalising_constant, m2.normalising_constant);
        ASSERT_EQ(m1.nodes.size(), m2.nodes.size());
        for (size_t i = 0; i < m1.nodes.size(); i++)
            ASSERT_EQ(m1.nodes[i].neighbor_weights, m2.nodes[i].neighbor_weights);
    }
}

TEST(FastgapCanonical, ShiftRelabelsClassesPerShot) {
    // obs'(E) = obs(E) XOR (XOR_{x in D} s(x)) for any error set E.
    std::mt19937_64 rng(11);
    auto dem = testutil::random_grid_dem(rng, 5, 4);
    auto spec = canonicalize_dem(dem, 0);
    std::vector<std::pair<std::vector<uint64_t>, uint8_t>> raw, canon;
    auto collect = [](const stim::DetectorErrorModel& m, std::vector<std::pair<std::vector<uint64_t>, uint8_t>>& out) {
        m.iter_flatten_error_instructions([&](const stim::DemInstruction& e) {
            std::vector<uint64_t> dets;
            uint8_t o = 0;
            for (auto& t : e.target_data) {
                if (t.is_relative_detector_id())
                    dets.push_back(t.val());
                if (t.is_observable_id() && t.val() == 0)
                    o ^= 1;
            }
            out.push_back({dets, o});
        });
    };
    collect(dem, raw);
    collect(spec.canonical_dem, canon);
    ASSERT_EQ(raw.size(), canon.size());
    for (int trial = 0; trial < 200; trial++) {
        std::set<uint64_t> syndrome;
        uint8_t o_raw = 0, o_canon = 0;
        for (size_t i = 0; i < raw.size(); i++) {
            if (rng() % 7 == 0) {
                for (auto d : raw[i].first) {
                    if (!syndrome.erase(d))
                        syndrome.insert(d);
                }
                o_raw ^= raw[i].second;
                o_canon ^= canon[i].second;
            }
        }
        std::vector<uint64_t> events(syndrome.begin(), syndrome.end());
        ASSERT_EQ(o_canon, o_raw ^ gauge_shift(spec, events));
    }
}

TEST(FastgapCanonical, BoundarylessLogicalRaises) {
    // A ring with odd observable parity and no boundary: a toroidal-style logical.
    stim::DetectorErrorModel dem(R"DEM(
        error(0.1) D0 D1 L0
        error(0.1) D1 D2
        error(0.1) D2 D3
        error(0.1) D3 D0
    )DEM");
    EXPECT_THROW(canonicalize_dem(dem, 0), BoundarylessLogicalError);
    // Even parity is fine.
    stim::DetectorErrorModel ok(R"DEM(
        error(0.1) D0 D1 L0
        error(0.1) D1 D2 L0
        error(0.1) D2 D3
        error(0.1) D3 D0
        error(0.1) D0
    )DEM");
    EXPECT_NO_THROW(canonicalize_dem(ok, 0));
}

TEST(FastgapCanonical, ToroidalSurfaceCodeRaises) {
    // A 2x2 periodic repetition-like lattice: every cycle around the torus flips L0.
    stim::DetectorErrorModel dem(R"DEM(
        error(0.1) D0 D1
        error(0.1) D1 D0 L0
        error(0.1) D2 D3
        error(0.1) D3 D2 L0
        error(0.1) D0 D2
        error(0.1) D1 D3
    )DEM");
    EXPECT_THROW(canonicalize_dem(dem, 0), std::invalid_argument);
}

TEST(FastgapCanonical, ParallelEdgeConflictRaises) {
    stim::DetectorErrorModel boundary_both_sides(R"DEM(
        error(0.1) D0 L0
        error(0.1) D0
        error(0.1) D0 D1
        error(0.1) D1
    )DEM");
    EXPECT_THROW(canonicalize_dem(boundary_both_sides, 0), ParallelEdgeConflictError);
    stim::DetectorErrorModel interior(R"DEM(
        error(0.1) D0 D1 L0
        error(0.1) D1 D0
        error(0.1) D0
        error(0.1) D1
    )DEM");
    // The two parallel mechanisms form an odd cycle D0-D1-D0: also a weight-2 logical.
    EXPECT_THROW(canonicalize_dem(interior, 0), std::invalid_argument);
}

TEST(FastgapCanonical, NegativeWeightsRaise) {
    stim::DetectorErrorModel dem(R"DEM(
        error(0.6) D0 L0
        error(0.1) D0 D1
        error(0.1) D1
    )DEM");
    EXPECT_THROW(canonicalize_dem(dem, 0), std::invalid_argument);
}

TEST(FastgapCanonical, HyperedgesAndDecomposition) {
    stim::DetectorErrorModel dem(R"DEM(
        error(0.1) D0 D1 D2
        error(0.1) D0 D1 L0 ^ D2
        error(0.1) D0
        error(0.1) D2
        error(0.1) D1 D2
    )DEM");
    auto spec = canonicalize_dem(dem, 0);
    EXPECT_EQ(spec.num_hyperedges_ignored, 1u);
    expect_observable_only_on_boundary(spec);
    // The baseline graph shares internal units with the decode graph.
    auto b = make_baseline_user_graph(spec).to_matching_graph(pm::NUM_DISTINCT_WEIGHTS);
    auto g = make_decode_user_graph(spec).to_matching_graph(pm::NUM_DISTINCT_WEIGHTS);
    EXPECT_EQ(b.normalising_constant, g.normalising_constant);
    EXPECT_EQ(b.nodes.size(), g.nodes.size() + 1);
}

TEST(FastgapCanonical, RepeatBlocksAndShifts) {
    auto dem = testutil::surface_code_dem(3, 3, 0.01, "rotated_memory_x");
    auto spec = canonicalize_dem(dem, 0);
    expect_observable_only_on_boundary(spec);
    EXPECT_EQ(spec.num_detectors, dem.count_detectors());
    auto g1 = pm::detector_error_model_to_user_graph(dem, false, pm::NUM_DISTINCT_WEIGHTS);
    auto g2 = make_decode_user_graph(spec);
    EXPECT_EQ(
        g1.to_matching_graph(pm::NUM_DISTINCT_WEIGHTS).normalising_constant,
        g2.to_matching_graph(pm::NUM_DISTINCT_WEIGHTS).normalising_constant);
}

TEST(FastgapCanonical, ReferenceEdgesParallelCheck) {
    std::vector<RefEdge> edges = {
        {0, SIZE_MAX, true, 1.0, 0.1},
        {0, 1, false, 1.0, 0.1},
        {1, SIZE_MAX, false, 1.0, 0.1},
    };
    EXPECT_NO_THROW(spec_from_reference_edges(3, 2, edges));
    edges.push_back({0, SIZE_MAX, false, 1.0, 0.1});
    EXPECT_THROW(spec_from_reference_edges(3, 2, edges), ParallelEdgeConflictError);
}
