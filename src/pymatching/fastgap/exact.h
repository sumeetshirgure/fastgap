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

#ifndef PYMATCHING2_FASTGAP_EXACT_H
#define PYMATCHING2_FASTGAP_EXACT_H

#include <memory>
#include <vector>

#include "pymatching/fastgap/gap_index.h"
#include "pymatching/fastgap/thread_pool.h"

namespace pm {
namespace fastgap {

struct ExactGap {
    /// Predicted flip of observable k in the original (un-gauged) labels. On a tie (gap 0) the
    /// even class (canonical labels) is reported.
    uint8_t prediction = 0;
    dist_int w_star = 0;
    dist_int gap = 0;
    /// Minimum weights over matchings with an even / odd number of L edges (canonical labels).
    dist_int w_even = DIST_INF;
    dist_int w_odd = DIST_INF;
    int64_t t_exact_ns = 0;
    /// The two decodes timed separately: the one in the class of M* (it yields W* and the
    /// prediction, like an ordinary decode) and the one in the complementary class (W^c). The
    /// latter is what the baseline pays for the gap on top of an ordinary decode.
    int64_t t_mstar_ns = 0;
    int64_t t_comp_ns = 0;
};

/// The exact two-decode baseline (CLAUDE.md §5.2): PyMatching on the gauge-fixed graph with the
/// observable turned into an extra detector, decoded with that detector off and on. This is the
/// method of Gidney's cultivation and yoking samplers.
class ExactBaseline {
   public:
    explicit ExactBaseline(std::shared_ptr<const GapIndex> index, size_t num_threads = 1);
    ExactGap decode(const std::vector<uint64_t>& detection_events);
    void decode_batch(const uint8_t* shots, size_t num_shots, size_t bytes_per_shot, std::vector<ExactGap>& out);
    double normalising_constant() const;

   private:
    ExactGap decode_with(pm::Mwpm& mwpm, std::vector<uint64_t>& events);
    std::shared_ptr<const GapIndex> index_;
    ThreadPool pool_;
    std::vector<pm::Mwpm> mwpms_;
    std::vector<std::vector<uint64_t>> scratch_;
};

/// Brute-force exact gap by enumerating all matchings in the L/R model (subset DP over at most
/// 20 defects), using exact interior distances. Includes the bare L-R chain of weight d_LR.
ExactGap brute_force_gap(const GapIndex& index, const std::vector<uint64_t>& detection_events);

}  // namespace fastgap
}  // namespace pm

#endif  // PYMATCHING2_FASTGAP_EXACT_H
