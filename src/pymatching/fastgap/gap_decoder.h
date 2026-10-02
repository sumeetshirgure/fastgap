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

#ifndef PYMATCHING2_FASTGAP_GAP_DECODER_H
#define PYMATCHING2_FASTGAP_GAP_DECODER_H

#include <memory>
#include <vector>

#include "pymatching/fastgap/dual_state.h"
#include "pymatching/fastgap/gap_index.h"
#include "pymatching/fastgap/thread_pool.h"

namespace pm {
namespace fastgap {

enum class Pricing : uint8_t {
    /// sigma_hat from tables and landmarks (CLAUDE.md §3.4): the number fastgap reports.
    INDEX = 0,
    /// Exact sigma from a per-shot Dijkstra: the slow reference Delta_sigma (CLAUDE.md §5.2).
    EXACT = 1,
};

struct DecodeOptions {
    /// Early-termination threshold tau in internal units (CLAUDE.md §3.6); < 0 means none.
    dist_int threshold = -1;
    /// Compute gap_ub (CLAUDE.md §3.5).
    bool upper_bound = false;
    Pricing pricing = Pricing::INDEX;
    /// Run the §5.1 debug invariants on every shot (always on in tests).
    bool check = CHECKS_DEFAULT_ON;
    /// Keep the walk in the result (single-shot API).
    bool keep_walk = false;
};

/// One step of the returned walk P_H: a hop from `from` to `to`. Endpoints are detector ids, or
/// -1 for L and -2 for R.
struct Hop {
    int64_t from;
    int64_t to;
    dist_int price;  // c(from, to) used by the search
    bool exact;      // the price used an exact distance
};

struct GapResult {
    pm::obs_int prediction = 0;  // original (un-gauged) observable labels
    dist_int w_star = 0;
    dist_int gap_lb = 0;
    dist_int gap_ub = DIST_INF;
    bool censored = false;
    bool walk_simple = false;
    bool all_hops_exact = false;
    bool blossom_pass = false;  // some hop-in/hop-out run through a blossom
    uint32_t num_defects = 0;
    uint32_t num_blossoms = 0;
    uint32_t num_hops = 0;
    dist_int hop_slack = -1;  // sum over hops of (dist - d_hat), only with upper_bound
    int64_t t_decode_ns = 0;
    int64_t t_gap_ns = 0;
    std::vector<Hop> walk;
};

/// Per-thread state: a private Mwpm plus scratch buffers reused across shots.
struct GapWorker {
    pm::Mwpm mwpm;
    DualState ds;
    DualStateScratch ds_scratch;
    std::vector<uint64_t> events;
    std::vector<dist_int> label;
    std::vector<uint8_t> settled;
    std::vector<int32_t> pred_state;
    std::vector<int32_t> pred_hop;
    std::vector<dist_int> row;  // exact pricing: distances from the current state to all detectors
    InteriorGraph::AStarScratch astar;
    // Per-shot pricing scratch (index pricing).
    std::vector<int32_t> node_to_defect;  // -1 for non-defects; reset after every shot
    std::vector<int32_t> act;             // active hop targets (defect indices), compacted
    std::vector<int32_t> pos;             // position of each defect in `act`, or -1
    std::vector<int32_t> gathered;        // saturated potentials of active targets, potential-major (k * n + p)
    std::vector<int32_t> gathered_rows;   // saturated potentials of all defects, defect-major
    std::vector<int32_t> row_lb;
    std::vector<dist_int> row_price;      // prices of the current row
    std::vector<uint8_t> row_is_exact;
    std::vector<int32_t> top;             // top-level blossom of each defect, -1 if none
};

class GapDecoder {
   public:
    GapDecoder(std::shared_ptr<const GapIndex> index, size_t num_threads = 1, bool pin_threads = false);

    /// Decode one shot given as detection-event indices.
    GapResult decode(const std::vector<uint64_t>& detection_events, const DecodeOptions& options);

    /// Decode `num_shots` bit-packed shots (`bytes_per_shot` bytes each, little-endian bit order as
    /// produced by stim). Shots are spread over the thread pool. `out` is resized to num_shots.
    void decode_batch(
        const uint8_t* shots,
        size_t num_shots,
        size_t bytes_per_shot,
        const DecodeOptions& options,
        std::vector<GapResult>& out);

    const GapIndex& index() const {
        return *index_;
    }
    size_t num_threads() const {
        return pool_.size();
    }

   private:
    /// `parallel_pricing` lets a single-shot decode fill the price matrix on the pool; it must be
    /// false inside a pool task (batch), since the pool is not re-entrant.
    void decode_with(
        GapWorker& w,
        const std::vector<uint64_t>& events,
        const DecodeOptions& o,
        GapResult& res,
        bool parallel_pricing);
    std::shared_ptr<const GapIndex> index_;
    ThreadPool pool_;
    std::vector<std::unique_ptr<GapWorker>> workers_;
};

/// Unpacks one bit-packed shot into detection events (dropping `ignored_detector` if given).
void unpack_shot(
    const uint8_t* shot, size_t bytes_per_shot, size_t num_detectors, size_t ignored_detector,
    std::vector<uint64_t>& out);

}  // namespace fastgap
}  // namespace pm

#endif  // PYMATCHING2_FASTGAP_GAP_DECODER_H
