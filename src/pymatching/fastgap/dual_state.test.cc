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

#include <gtest/gtest.h>

#include "pymatching/fastgap/fastgap_test_util.test.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"

using namespace pm::fastgap;

namespace {

/// Checks Sigma y == W* against upstream, the prediction against upstream on the original DEM,
/// relay tightness and all-pairs sigma >= 0 with brute-force Dijkstra.
size_t check_against_upstream(const stim::DetectorErrorModel& dem, size_t shots, uint64_t seed) {
    auto idx = GapIndex::create_from_dem(dem.str(), 0);
    pm::Mwpm ours = idx->make_mwpm();
    pm::Mwpm upstream = pm::detector_error_model_to_mwpm(dem, pm::NUM_DISTINCT_WEIGHTS);
    DualState ds;
    DualStateScratch scratch;
    std::mt19937_64 rng(seed);
    size_t blossom_shots = 0;
    for (size_t s = 0; s < shots; s++) {
        auto events = testutil::sample_dem(dem, rng);
        extract_dual_state(ours, events, 0, ds, scratch);
        auto ref = pm::decode_detection_events_for_up_to_64_observables(upstream, events, false);
        EXPECT_EQ(ds.w_star, ref.weight);
        pm::obs_int shift = gauge_shift(idx->spec, events);
        EXPECT_EQ(ds.obs_mask ^ shift, ref.obs_mask);
        blossom_shots += ds.num_blossoms() > 0;

        size_t n = ds.num_defects();
        dist_int sum = 0;
        for (auto y : ds.y_leaf)
            sum += y;
        for (auto y : ds.y_blossom)
            sum += y;
        EXPECT_EQ(sum, ds.w_star);
        for (size_t i = 0; i < n; i++) {
            auto row = idx->graph.dijkstra({{ds.defects[i], 0}});
            EXPECT_GE(idx->h_l[ds.defects[i]] - ds.rho[i], 0);
            EXPECT_GE(idx->h_r[ds.defects[i]] - ds.rho[i], 0);
            if (ds.mate[i] == MATE_L)
                EXPECT_EQ(idx->h_l[ds.defects[i]], ds.rho[i]);
            if (ds.mate[i] == MATE_R)
                EXPECT_EQ(idx->h_r[ds.defects[i]], ds.rho[i]);
            for (size_t j = 0; j < n; j++) {
                if (i == j)
                    continue;
                dist_int d = row[ds.defects[j]];
                dist_int sigma = d - ds.rho[i] - ds.rho[j] + 2 * ds.beta(i, j);
                EXPECT_GE(sigma, 0);
                EXPECT_EQ(ds.beta(i, j), ds.beta(j, i));
                if (ds.mate[i] == (int32_t)j) {
                    EXPECT_EQ(sigma, 0);
                    EXPECT_EQ(d, ds.r(i) + ds.r(j));
                }
            }
        }
    }
    return blossom_shots;
}

}  // namespace

TEST(FastgapDualState, TinyHandExample) {
    // D0 -- D1 -- D2 with boundaries at both ends.
    stim::DetectorErrorModel dem(R"DEM(
        error(0.1) D0 L0
        error(0.1) D0 D1
        error(0.1) D1 D2
        error(0.1) D2
    )DEM");
    auto idx = GapIndex::create_from_dem(dem.str(), 0);
    pm::Mwpm mwpm = idx->make_mwpm();
    DualState ds;
    DualStateScratch scratch;
    extract_dual_state(mwpm, {0, 1}, 0, ds, scratch);
    ASSERT_EQ(ds.num_defects(), 2u);
    EXPECT_EQ(ds.mate[0], 1);
    EXPECT_EQ(ds.mate[1], 0);
    EXPECT_EQ(ds.w_star, idx->graph.weights[idx->graph.offsets[0]]);
    extract_dual_state(mwpm, {0}, 0, ds, scratch);
    EXPECT_EQ(ds.mate[0], MATE_L);
    EXPECT_EQ(ds.obs_mask, 1u);
    extract_dual_state(mwpm, {}, 0, ds, scratch);
    EXPECT_EQ(ds.num_defects(), 0u);
    EXPECT_EQ(ds.w_star, 0);
}

TEST(FastgapDualState, RandomGridsAgainstUpstream) {
    std::mt19937_64 rng(17);
    size_t blossoms = 0;
    for (int trial = 0; trial < 12; trial++) {
        auto dem = testutil::random_grid_dem(rng, 5 + trial % 4, 4 + trial % 3, 0.5);
        blossoms += check_against_upstream(dem, 150, trial);
    }
    EXPECT_GT(blossoms, 0u);
}

TEST(FastgapDualState, SurfaceCodeAgainstUpstream) {
    for (auto task : {"rotated_memory_x", "rotated_memory_z"}) {
        auto dem = testutil::surface_code_dem(5, 5, 0.01, task);
        size_t blossoms = check_against_upstream(dem, 150, 99);
        EXPECT_GT(blossoms, 0u);
    }
}

TEST(FastgapDualState, OddComponentWithoutBoundaryPropagates) {
    stim::DetectorErrorModel dem(R"DEM(
        error(0.1) D0 L0
        error(0.1) D0 D1
        error(0.1) D2 D3
    )DEM");
    auto idx = GapIndex::create_from_dem(dem.str(), 0);
    pm::Mwpm mwpm = idx->make_mwpm();
    DualState ds;
    DualStateScratch scratch;
    EXPECT_THROW(extract_dual_state(mwpm, {2}, 0, ds, scratch), std::invalid_argument);
    // The Mwpm remains usable afterwards.
    extract_dual_state(mwpm, {2, 3}, 0, ds, scratch);
    EXPECT_EQ(ds.mate[0], 1);
}

TEST(FastgapDualState, BadShotsLeaveMwpmUsable) {
    stim::DetectorErrorModel dem(R"DEM(
        error(0.1) D0 L0
        error(0.1) D0 D1
        error(0.1) D1 D2
        error(0.1) D2
    )DEM");
    auto idx = GapIndex::create_from_dem(dem.str(), 0);
    pm::Mwpm mwpm = idx->make_mwpm();
    DualState ds;
    DualStateScratch scratch;
    EXPECT_THROW(extract_dual_state(mwpm, {1, 1}, 0, ds, scratch), std::invalid_argument);
    EXPECT_THROW(extract_dual_state(mwpm, {0, 7}, 0, ds, scratch), std::invalid_argument);
    extract_dual_state(mwpm, {0, 1}, 0, ds, scratch);
    EXPECT_EQ(ds.mate[0], 1);
    EXPECT_EQ(ds.w_star, idx->graph.weights[idx->graph.offsets[0]]);
}
