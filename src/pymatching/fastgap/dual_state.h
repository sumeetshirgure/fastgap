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

#ifndef PYMATCHING2_FASTGAP_DUAL_STATE_H
#define PYMATCHING2_FASTGAP_DUAL_STATE_H

#include <unordered_map>
#include <vector>

#include "pymatching/fastgap/types.h"
#include "pymatching/sparse_blossom/matcher/mwpm.h"

namespace pm {
namespace fastgap {

/// The decoder's final dual state and matching (CLAUDE.md §2.4), read before shattering.
struct DualState {
    /// Detector id of each defect, in detection-event order (events on user boundary nodes skipped).
    std::vector<uint32_t> defects;
    std::vector<dist_int> y_leaf;
    /// rho_v = y_v + sum over blossoms S containing v of y_S.
    std::vector<dist_int> rho;
    /// Ancestor blossoms of each defect, ordered from the top-level blossom down to the leaf's
    /// immediate parent, so common ancestors of two defects form a common prefix. CSR layout.
    std::vector<uint32_t> anc_offset;
    std::vector<int32_t> anc_blossom;
    /// anc_prefix[anc_offset[i] + j] = sum of y over the first j + 1 ancestors of defect i.
    std::vector<dist_int> anc_prefix;
    std::vector<dist_int> y_blossom;
    std::vector<int32_t> blossom_depth;  // 0 = top level
    /// mate[i]: index of the defect matched to defect i, or MATE_L / MATE_R.
    std::vector<int32_t> mate;
    /// W* = sum of all y in internal units.
    dist_int w_star = 0;
    /// XOR of the observable masks of all match edges (canonical labels).
    pm::obs_int obs_mask = 0;

    size_t num_defects() const {
        return defects.size();
    }
    size_t num_blossoms() const {
        return y_blossom.size();
    }
    int32_t top_blossom(size_t i) const {
        return anc_offset[i] == anc_offset[i + 1] ? -1 : anc_blossom[anc_offset[i]];
    }
    /// beta(i, j): summed y_S over blossoms containing both defects.
    inline dist_int beta(size_t i, size_t j) const {
        uint32_t a = anc_offset[i], ae = anc_offset[i + 1];
        uint32_t b = anc_offset[j], be = anc_offset[j + 1];
        if (a == ae || b == be || anc_blossom[a] != anc_blossom[b])
            return 0;  // fast path: different (or no) top-level blossoms
        uint32_t k = 0;
        while (a + k + 1 < ae && b + k + 1 < be && anc_blossom[a + k + 1] == anc_blossom[b + k + 1])
            k++;
        return anc_prefix[a + k];
    }
    /// r(i) = rho_i - beta(i, mate(i)), so that dist(i, mate(i)) = r(i) + r(mate(i)).
    inline dist_int r(size_t i) const {
        int32_t m = mate[i];
        return m >= 0 ? rho[i] - beta(i, (size_t)m) : rho[i];
    }

    void clear();
};

/// Reusable per-thread scratch for `extract_dual_state`.
struct DualStateScratch {
    std::vector<int32_t> defect_of_node;
    std::unordered_map<const pm::GraphFillRegion*, int32_t> blossom_ids;
    std::vector<const pm::GraphFillRegion*> chain;
};

/// Runs sparse blossom to completion on `detection_events`, reads the dual state before
/// shattering, then shatters to obtain the matching. `observable` is the gauge-fixed observable
/// whose bit on a boundary match edge means side L. Leaves the Mwpm ready for the next decode.
void extract_dual_state(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    size_t observable,
    DualState& out,
    DualStateScratch& scratch);

}  // namespace fastgap
}  // namespace pm

#endif  // PYMATCHING2_FASTGAP_DUAL_STATE_H
