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

#include <chrono>

#include "pymatching/fastgap/gap_decoder.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"

namespace pm {
namespace fastgap {

namespace {

void finish(ExactGap& g, uint8_t shift) {
    if (is_inf(g.w_even) && is_inf(g.w_odd))
        throw std::invalid_argument("no perfect matching exists for this syndrome");
    bool odd_wins = g.w_odd < g.w_even;
    g.w_star = odd_wins ? g.w_odd : g.w_even;
    dist_int other = odd_wins ? g.w_even : g.w_odd;
    g.gap = is_inf(other) ? DIST_INF : other - g.w_star;
    g.prediction = (uint8_t)odd_wins ^ shift;
}

}  // namespace

ExactBaseline::ExactBaseline(std::shared_ptr<const GapIndex> index, size_t num_threads)
    : index_(std::move(index)), pool_(num_threads) {
    for (size_t t = 0; t < pool_.size(); t++) {
        pm::UserGraph g = make_baseline_user_graph(index_->spec);
        mwpms_.push_back(g.to_mwpm(pm::NUM_DISTINCT_WEIGHTS, false));
        if (mwpms_.back().flooder.graph.normalising_constant != index_->graph.normalising_constant)
            throw std::logic_error("exact baseline graph has different internal units from fastgap's graph");
    }
    scratch_.resize(pool_.size());
}

double ExactBaseline::normalising_constant() const {
    return mwpms_[0].flooder.graph.normalising_constant;
}

ExactGap ExactBaseline::decode_with(pm::Mwpm& mwpm, std::vector<uint64_t>& events) {
    auto t0 = std::chrono::steady_clock::now();
    ExactGap g;
    size_t obs_det = index_->spec.baseline_obs_detector();
    int64_t t_parity[2];
    for (int parity = 0; parity < 2; parity++) {
        auto tp = std::chrono::steady_clock::now();
        if (parity == 1)
            events.push_back(obs_det);
        try {
            auto r = pm::decode_detection_events_for_up_to_64_observables(mwpm, events, false);
            (parity ? g.w_odd : g.w_even) = r.weight;
        } catch (const std::invalid_argument&) {
            // No matching in this class (e.g. no path from the observable detector to anything).
        }
        if (parity == 1)
            events.pop_back();
        t_parity[parity] =
            std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - tp).count();
    }
    finish(g, gauge_shift(index_->spec, events));
    g.t_exact_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t0).count();
    bool odd_wins = g.w_odd < g.w_even;  // as in finish()
    g.t_mstar_ns = t_parity[odd_wins ? 1 : 0];
    g.t_comp_ns = t_parity[odd_wins ? 0 : 1];
    return g;
}

ExactGap ExactBaseline::decode(const std::vector<uint64_t>& detection_events) {
    auto& events = scratch_[0];
    events.clear();
    for (auto d : detection_events) {
        if (d != index_->spec.ref_obs_detector)
            events.push_back(d);
    }
    return decode_with(mwpms_[0], events);
}

void ExactBaseline::decode_batch(
    const uint8_t* shots, size_t num_shots, size_t bytes_per_shot, std::vector<ExactGap>& out) {
    out.resize(num_shots);
    size_t num_detectors = index_->num_detectors();
    size_t ignored = index_->spec.ref_obs_detector;
    pool_.run(num_shots, [&](size_t worker, size_t shot) {
        auto& events = scratch_[worker];
        unpack_shot(shots + shot * bytes_per_shot, bytes_per_shot, num_detectors, ignored, events);
        out[shot] = decode_with(mwpms_[worker], events);
    });
}

ExactGap brute_force_gap(const GapIndex& index, const std::vector<uint64_t>& detection_events) {
    std::vector<uint32_t> defects;
    for (auto d : detection_events) {
        if (d != index.spec.ref_obs_detector)
            defects.push_back((uint32_t)d);
    }
    size_t n = defects.size();
    if (n > 20)
        throw std::invalid_argument("brute_force_gap supports at most 20 defects");
    std::vector<std::vector<dist_int>> dist(n);
    for (size_t i = 0; i < n; i++) {
        auto row = index.graph.dijkstra({{defects[i], 0}});
        for (size_t j = 0; j < n; j++)
            dist[i].push_back(row[defects[j]]);
    }
    auto add = [](dist_int a, dist_int b) { return (is_inf(a) || is_inf(b)) ? DIST_INF : a + b; };
    size_t full = ((size_t)1 << n) - 1;
    // f[mask * 2 + p]: min weight matching the defects in mask, with L-edge parity p.
    std::vector<dist_int> f((full + 1) * 2, DIST_INF);
    f[0] = 0;
    for (size_t mask = 1; mask <= full; mask++) {
        size_t i = __builtin_ctzll(mask);
        size_t rest = mask & ~((size_t)1 << i);
        for (int p = 0; p < 2; p++) {
            dist_int best = DIST_INF;
            best = std::min(best, add(f[rest * 2 + (p ^ 1)], index.h_l[defects[i]]));
            best = std::min(best, add(f[rest * 2 + p], index.h_r[defects[i]]));
            for (size_t j = i + 1; j < n; j++) {
                if (rest & ((size_t)1 << j))
                    best = std::min(best, add(f[(rest & ~((size_t)1 << j)) * 2 + p], dist[i][j]));
            }
            f[mask * 2 + p] = best;
        }
    }
    ExactGap g;
    // A matching may also include one bare L-R chain, which flips the class.
    g.w_even = std::min(f[full * 2 + 0], add(f[full * 2 + 1], index.d_lr));
    g.w_odd = std::min(f[full * 2 + 1], add(f[full * 2 + 0], index.d_lr));
    finish(g, gauge_shift(index.spec, detection_events));
    return g;
}

}  // namespace fastgap
}  // namespace pm
