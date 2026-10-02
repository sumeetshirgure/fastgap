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

#include "pymatching/fastgap/exact.h"

#include <gtest/gtest.h>

#include "pymatching/fastgap/fastgap_test_util.test.h"

using namespace pm::fastgap;

TEST(FastgapExact, BaselineEqualsBruteForce) {
    std::mt19937_64 rng(2024);
    size_t compared = 0;
    for (int trial = 0; trial < 10; trial++) {
        auto dem = testutil::random_grid_dem(rng, 4 + trial % 3, 3 + trial % 3, 0.4, 0.02, 0.3);
        auto idx = testutil::make_index(dem);
        ExactBaseline baseline(idx);
        EXPECT_EQ(baseline.normalising_constant(), idx->graph.normalising_constant);
        for (int s = 0; s < 100; s++) {
            auto events = testutil::sample_dem(dem, rng);
            if (events.size() > 12)
                continue;
            auto a = baseline.decode(events);
            auto b = brute_force_gap(*idx, events);
            ASSERT_EQ(a.w_even, b.w_even);
            ASSERT_EQ(a.w_odd, b.w_odd);
            compared++;
        }
    }
    EXPECT_GT(compared, 500u);
}

TEST(FastgapExact, GaugeLeavesGapUnchanged) {
    // Same code with two different random gauges: identical gaps, predictions flipped by the
    // per-shot shift so that both agree in the original labels.
    std::mt19937_64 rng_a(7), rng_b(7);
    auto dem_a = testutil::random_grid_dem(rng_a, 5, 4, 0.4);
    // Re-label the observable by an explicit gauge on the original DEM: flip L0 on every
    // mechanism touching D3 (an odd number of times).
    std::string text = dem_a.str();
    stim::DetectorErrorModel dem_b;
    dem_a.iter_flatten_error_instructions([&](const stim::DemInstruction& e) {
        std::vector<stim::DemTarget> t;
        bool touches = false, obs = false;
        for (auto& x : e.target_data) {
            if (x.is_relative_detector_id()) {
                t.push_back(x);
                touches ^= x.val() == 3;
            } else if (x.is_observable_id()) {
                obs ^= true;
            }
        }
        if (obs ^ touches)
            t.push_back(stim::DemTarget::observable_id(0));
        dem_b.append_error_instruction(e.arg_data[0], t, "");
    });
    auto idx_a = testutil::make_index(dem_a);
    auto idx_b = testutil::make_index(dem_b);
    ExactBaseline ba(idx_a), bb(idx_b);
    std::mt19937_64 rng(1);
    for (int s = 0; s < 300; s++) {
        auto events = testutil::sample_dem(dem_a, rng);
        auto ga = ba.decode(events);
        auto gb = bb.decode(events);
        ASSERT_EQ(ga.gap, gb.gap);
        ASSERT_EQ(ga.w_star, gb.w_star);
        bool has3 = std::find(events.begin(), events.end(), 3u) != events.end();
        if (ga.gap > 0)
            ASSERT_EQ(ga.prediction, gb.prediction ^ (uint8_t)has3);
    }
}
