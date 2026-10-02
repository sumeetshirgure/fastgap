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

#include "pymatching/fastgap/gap_decoder.h"

#include <algorithm>
#include <chrono>

namespace pm {
namespace fastgap {

namespace {

constexpr int32_t STATE_L = -1;
constexpr int32_t END_R = -2;
constexpr int32_t PRED_FREE_RELAY = -3;  // standing at u with mate(u) = L, reached for free

/// Rows of the price matrix are prefilled in parallel only when a shot is this large.
constexpr size_t MATRIX_MIN_DEFECTS = 48;

inline int64_t elapsed_ns(std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
}

/// Hop prices c(x, u) (CLAUDE.md §3.3, §3.4). Every price is a lower bound on the true reduced
/// cost, clamped at 0 (safe because sigma >= 0).
///
/// The search reads one row at a time ("standing at x") and only needs prices to *active* hop
/// targets: defects u whose mate is R or an unsettled state (once state m is settled, hopping
/// into mate(m) can no longer improve anything). The active targets are kept compacted, with
/// their saturated potentials in a potential-major block, so `begin_row(x)` computes the landmark
/// bound for all of them in one vectorisable pass, then scans x's exact local table once.
/// `price_at` (current row) and `hop` (any pair, read-only) return identical prices.
struct Pricer {
    const GapIndex& idx;
    const DualState& ds;
    GapWorker& w;
    Pricing mode;
    // Optional prefilled matrix of prices between defects (n x n), with exactness flags.
    const std::vector<dist_int>* matrix = nullptr;
    const std::vector<uint8_t>* matrix_exact = nullptr;
    int32_t row_state = INT32_MIN;
    size_t stride = 0;
    size_t n = 0;

    Pricer(const GapIndex& idx, const DualState& ds, GapWorker& w, Pricing mode)
        : idx(idx), ds(ds), w(w), mode(mode), stride(idx.potential_stride()), n(ds.num_defects()) {
        w.top.resize(n);
        for (size_t i = 0; i < n; i++)
            w.top[i] = ds.top_blossom(i);
        // Hops into a defect matched to L only lead back to L, so those are never targets.
        w.act.clear();
        w.pos.assign(n, -1);
        for (size_t i = 0; i < n; i++) {
            if (ds.mate[i] != MATE_L) {
                w.pos[i] = (int32_t)w.act.size();
                w.act.push_back((int32_t)i);
            }
        }
        if (mode != Pricing::INDEX)
            return;
        if (w.node_to_defect.size() != idx.graph.num_nodes)
            w.node_to_defect.assign(idx.graph.num_nodes, -1);
        w.gathered.resize(n * stride);
        w.gathered_rows.resize(n * stride);
        w.row_price.resize(n);
        w.row_lb.resize(n);
        w.row_is_exact.resize(n);
        for (size_t i = 0; i < n; i++) {
            uint32_t u = ds.defects[i];
            w.node_to_defect[u] = (int32_t)i;
            std::copy_n(idx.potentials32.data() + (size_t)u * stride, stride, w.gathered_rows.data() + i * stride);
        }
        for (size_t p = 0; p < w.act.size(); p++) {
            const int32_t* src = w.gathered_rows.data() + (size_t)w.act[p] * stride;
            for (size_t k = 0; k < stride; k++)
                w.gathered[k * n + p] = src[k];
        }
    }
    ~Pricer() {
        if (mode == Pricing::INDEX) {
            for (auto u : ds.defects)
                w.node_to_defect[u] = -1;
        }
    }
    Pricer(const Pricer&) = delete;

    size_t num_active() const {
        return w.act.size();
    }

    /// Removes u from the active targets (swap with the last one).
    void deactivate(int32_t u) {
        int32_t p = w.pos[u];
        if (p < 0)
            return;
        size_t last = w.act.size() - 1;
        int32_t moved = w.act[last];
        w.act[p] = moved;
        w.pos[moved] = p;
        w.pos[u] = -1;
        w.act.pop_back();
        if (mode == Pricing::INDEX) {
            for (size_t k = 0; k < stride; k++)
                w.gathered[k * n + p] = w.gathered[k * n + last];
        }
        row_state = INT32_MIN;  // row arrays are indexed by position
    }

    void begin_row(int32_t x) {
        if (x < 0 || row_state == x)
            return;
        row_state = x;
        if (mode == Pricing::EXACT) {
            w.row = idx.graph.dijkstra({{ds.defects[x], 0}});
            return;
        }
        if (matrix != nullptr)
            return;
        // Landmark bound for every active target at once (vectorisable: contiguous over p).
        size_t a = w.act.size();
        int32_t* lb32 = w.row_lb.data();
        std::fill(lb32, lb32 + a, 0);
        const int32_t* px = w.gathered_rows.data() + (size_t)x * stride;
        for (size_t k = 0; k < stride; k++) {
            const int32_t* col = w.gathered.data() + k * n;
            int32_t pk = px[k];
            for (size_t p = 0; p < a; p++) {
                int32_t diff = col[p] - pk;
                diff = diff < 0 ? -diff : diff;
                lb32[p] = diff > lb32[p] ? diff : lb32[p];
            }
        }
        dist_int* price = w.row_price.data();
        for (size_t p = 0; p < a; p++)
            price[p] = lb32[p];
        // Exact table hits override the bound.
        std::fill(w.row_is_exact.begin(), w.row_is_exact.begin() + a, 0);
        uint32_t u = ds.defects[x];
        for (uint64_t e = idx.table_offset[u]; e < idx.table_offset[u + 1]; e++) {
            int32_t j = w.node_to_defect[idx.table_ids[e]];
            if (j >= 0 && w.pos[j] >= 0) {
                price[w.pos[j]] = idx.table_dist[e];
                w.row_is_exact[w.pos[j]] = 1;
            }
        }
        dist_int floor = idx.table_radius + 1;
        dist_int rx = ds.rho[x];
        for (size_t p = 0; p < a; p++) {
            int32_t j = w.act[p];
            dist_int d = w.row_is_exact[p] ? price[p] : std::max(price[p], floor);
            if (is_inf(d)) {
                price[p] = DIST_INF;
                continue;
            }
            dist_int sv = d - rx - ds.rho[j];
            if (w.top[x] >= 0 && w.top[x] == w.top[j])
                sv += 2 * ds.beta(x, j);
            price[p] = sv > 0 ? sv : 0;
        }
    }

    /// Price of the hop from the current row x to the active target at position p (requires
    /// begin_row(x)).
    inline dist_int price_at(int32_t x, size_t p, bool& exact) const {
        int32_t j = w.act[p];
        if (matrix != nullptr)
            return from_matrix(x, j, exact);
        if (mode == Pricing::EXACT) {
            exact = true;
            return reduce(x, j, w.row[ds.defects[j]]);
        }
        exact = w.row_is_exact[p];
        return w.row_price[p];
    }

    inline dist_int beta(int32_t i, int32_t j) const {
        return (w.top[i] < 0 || w.top[i] != w.top[j]) ? 0 : ds.beta(i, j);
    }

    inline dist_int reduce(int32_t i, int32_t j, dist_int d) const {
        if (is_inf(d))
            return DIST_INF;
        dist_int s = d - ds.rho[i] - ds.rho[j] + 2 * beta(i, j);
        return s > 0 ? s : 0;
    }

    /// Price of any hop i -> j; read-only, so safe to call from several threads in INDEX mode.
    /// In EXACT mode it requires begin_row(i). Optionally returns the distance estimate.
    inline dist_int hop(int32_t i, int32_t j, bool& exact, dist_int* d_out = nullptr) const {
        if (matrix != nullptr && d_out == nullptr)
            return from_matrix(i, j, exact);
        dist_int d;
        if (mode == Pricing::EXACT) {
            exact = true;
            d = w.row[ds.defects[j]];
        } else {
            d = idx.d_hat(ds.defects[i], ds.defects[j], exact);
        }
        if (d_out)
            *d_out = d;
        return reduce(i, j, d);
    }

    inline dist_int from_matrix(int32_t i, int32_t j, bool& exact) const {
        size_t n = ds.num_defects();
        exact = (*matrix_exact)[(size_t)i * n + j];
        return (*matrix)[(size_t)i * n + j];
    }

    inline dist_int from_l(int32_t j) const {
        dist_int h = idx.h_l[ds.defects[j]];
        if (is_inf(h))
            return DIST_INF;
        dist_int s = h - ds.rho[j];
        return s > 0 ? s : 0;
    }

    inline dist_int to_r(int32_t i) const {
        if (i == STATE_L)
            return idx.d_lr;
        dist_int h = idx.h_r[ds.defects[i]];
        if (is_inf(h))
            return DIST_INF;
        dist_int s = h - ds.rho[i];
        return s > 0 ? s : 0;
    }
};

struct InternalHop {
    int32_t from;  // defect index or STATE_L
    int32_t to;    // defect index or END_R
    dist_int price;
    dist_int dhat;  // the distance estimate behind the price
    bool exact;
};

/// Debug invariants of CLAUDE.md §5.1: relays tight, sigma >= 0, d_hat <= dist.
void check_dual_invariants(const GapIndex& idx, const DualState& ds, GapWorker& w) {
    size_t n = ds.num_defects();
    for (size_t i = 0; i < n; i++) {
        uint32_t u = ds.defects[i];
        if (idx.h_l[u] - ds.rho[i] < 0 || idx.h_r[u] - ds.rho[i] < 0)
            throw InvariantViolation("boundary reduced cost sigma(u, X) < 0 at D" + std::to_string(u));
        int32_t m = ds.mate[i];
        if (m == MATE_L && idx.h_l[u] != ds.rho[i])
            throw InvariantViolation("relay to L is not tight at D" + std::to_string(u));
        if (m == MATE_R && idx.h_r[u] != ds.rho[i])
            throw InvariantViolation("relay to R is not tight at D" + std::to_string(u));
    }
    size_t pair_rows = std::min<size_t>(n, 64);
    std::vector<uint8_t> relay_checked(n, 0);
    for (size_t i = 0; i < n; i++) {
        int32_t m = ds.mate[i];
        bool need_row = i < pair_rows || (m >= 0 && (size_t)m > i);
        if (!need_row)
            continue;
        auto row = idx.graph.dijkstra({{ds.defects[i], 0}});
        for (size_t j = i + 1; j < n; j++) {
            bool is_mate = (int32_t)j == m;
            if (i >= pair_rows && !is_mate)
                continue;
            dist_int d = row[ds.defects[j]];
            dist_int reduced = d - ds.rho[i] - ds.rho[j] + 2 * ds.beta(i, j);
            if (!is_inf(d) && reduced < 0)
                throw InvariantViolation(
                    "dual infeasible: sigma(D" + std::to_string(ds.defects[i]) + ", D" + std::to_string(ds.defects[j]) +
                    ") = " + std::to_string(reduced) + " < 0");
            if (is_mate) {
                if (reduced != 0)
                    throw InvariantViolation("relay between matched defects is not tight");
                relay_checked[i] = relay_checked[j] = 1;
            }
            bool exact;
            dist_int dh = idx.d_hat(ds.defects[i], ds.defects[j], exact);
            if (dh > d || (exact && dh != d))
                throw InvariantViolation("d_hat overestimates the interior distance");
        }
    }
    for (size_t i = 0; i < n; i++) {
        if (ds.mate[i] >= 0 && !relay_checked[i])
            throw InvariantViolation("relay tightness was not checked");
    }
}

}  // namespace

void unpack_shot(
    const uint8_t* shot, size_t bytes_per_shot, size_t num_detectors, size_t ignored_detector,
    std::vector<uint64_t>& out) {
    out.clear();
    size_t limit = std::min(num_detectors, bytes_per_shot * 8);
    for (size_t b = 0; b < bytes_per_shot; b++) {
        uint8_t byte = shot[b];
        while (byte) {
            int bit = __builtin_ctz(byte);
            byte &= byte - 1;
            size_t d = b * 8 + bit;
            if (d >= limit)
                throw std::invalid_argument(
                    "shot has a detection event at index " + std::to_string(d) + " but the graph only has " +
                    std::to_string(num_detectors) + " detectors");
            if (d != ignored_detector)
                out.push_back(d);
        }
    }
}

GapDecoder::GapDecoder(std::shared_ptr<const GapIndex> index, size_t num_threads, bool pin_threads)
    : index_(std::move(index)), pool_(num_threads, pin_threads) {
    if (!index_->finalized)
        throw std::invalid_argument("GapIndex has no landmarks yet (call set_landmarks first)");
    for (size_t t = 0; t < pool_.size(); t++) {
        auto w = std::make_unique<GapWorker>();
        w->mwpm = index_->make_mwpm();
        workers_.push_back(std::move(w));
    }
}

GapResult GapDecoder::decode(const std::vector<uint64_t>& detection_events, const DecodeOptions& options) {
    GapResult res;
    auto& w = *workers_[0];
    w.events.clear();
    for (auto d : detection_events) {
        if (d != index_->spec.ref_obs_detector)
            w.events.push_back(d);
    }
    decode_with(w, w.events, options, res, true);
    return res;
}

void GapDecoder::decode_batch(
    const uint8_t* shots,
    size_t num_shots,
    size_t bytes_per_shot,
    const DecodeOptions& options,
    std::vector<GapResult>& out) {
    out.resize(num_shots);
    size_t num_detectors = index_->num_detectors();
    size_t ignored = index_->spec.ref_obs_detector;
    pool_.run(num_shots, [&](size_t worker, size_t shot) {
        auto& w = *workers_[worker];
        unpack_shot(shots + shot * bytes_per_shot, bytes_per_shot, num_detectors, ignored, w.events);
        out[shot] = GapResult();
        decode_with(w, w.events, options, out[shot], false);
    });
}

void GapDecoder::decode_with(
    GapWorker& w, const std::vector<uint64_t>& events, const DecodeOptions& o, GapResult& res, bool parallel_pricing) {
    const GapIndex& idx = *index_;
    auto t0 = std::chrono::steady_clock::now();

    // Stage 1: the ordinary, untruncated PyMatching decode, keeping the dual state.
    extract_dual_state(w.mwpm, events, idx.spec.observable, w.ds, w.ds_scratch);
    const DualState& ds = w.ds;
    auto t1 = std::chrono::steady_clock::now();

    pm::obs_int shift = (pm::obs_int)gauge_shift(idx.spec, events) << idx.spec.observable;
    res.prediction = ds.obs_mask ^ shift;
    res.w_star = ds.w_star;
    res.num_defects = (uint32_t)ds.num_defects();
    res.num_blossoms = (uint32_t)ds.num_blossoms();
    if (o.check)
        check_dual_invariants(idx, ds, w);

    // Stage 2: dense Dijkstra over the |D| + 2 states "standing at x" (CLAUDE.md §3.3).
    const int32_t n = (int32_t)ds.num_defects();
    Pricer pricer(idx, ds, w, o.pricing);
    std::vector<dist_int> matrix;
    std::vector<uint8_t> matrix_exact;
    if (parallel_pricing && o.pricing == Pricing::INDEX && pool_.size() > 1 && (size_t)n >= MATRIX_MIN_DEFECTS) {
        // Single-shot latency path: fill the upper triangle by row blocks on the persistent pool.
        matrix.assign((size_t)n * n, 0);
        matrix_exact.assign((size_t)n * n, 0);
        size_t blocks = std::min<size_t>(pool_.size() * 4, (size_t)n);
        pool_.run(blocks, [&](size_t, size_t b) {
            for (int32_t i = (int32_t)b; i < n; i += (int32_t)blocks) {
                for (int32_t j = i + 1; j < n; j++) {
                    bool exact;
                    dist_int c = pricer.hop(i, j, exact);
                    matrix[(size_t)i * n + j] = matrix[(size_t)j * n + i] = c;
                    matrix_exact[(size_t)i * n + j] = matrix_exact[(size_t)j * n + i] = exact;
                }
            }
        });
        pricer.matrix = &matrix;
        pricer.matrix_exact = &matrix_exact;
    }

    // State index n is "standing at L".
    w.label.assign(n + 1, DIST_INF);
    w.settled.assign(n + 1, 0);
    w.pred_state.assign(n + 1, INT32_MIN);
    w.pred_hop.assign(n + 1, -1);
    w.label[n] = 0;
    for (int32_t j = 0; j < n; j++) {
        if (ds.mate[j] == MATE_L) {
            w.label[j] = 0;
            w.pred_state[j] = PRED_FREE_RELAY;
        }
    }
    dist_int tau = o.threshold >= 0 ? o.threshold : DIST_INF;
    dist_int best = DIST_INF;
    int32_t best_from = INT32_MIN;
    int32_t best_hop = END_R;  // END_R: hop straight to R; otherwise the defect whose mate is R

    dist_int stop_label = DIST_INF;
    while (true) {
        int32_t x = -1;
        dist_int lx = DIST_INF;
        for (int32_t s = 0; s <= n; s++) {
            if (!w.settled[s] && w.label[s] < lx) {
                lx = w.label[s];
                x = s;
            }
        }
        if (x < 0 || lx >= best || lx >= tau) {
            stop_label = lx;
            break;
        }
        w.settled[x] = 1;
        int32_t xs = x == n ? STATE_L : x;
        // Move B: hop to R.
        dist_int c_r = pricer.to_r(xs);
        if (!is_inf(c_r) && lx + c_r < best) {
            best = lx + c_r;
            best_from = xs;
            best_hop = END_R;
        }
        // Move A: hop to defect u, then relay to mate(u). Settling x retires the target mate(x).
        if (xs >= 0 && ds.mate[xs] >= 0)
            pricer.deactivate(ds.mate[xs]);
        pricer.begin_row(xs);
        size_t a = pricer.num_active();
        for (size_t p = 0; p < a; p++) {
            int32_t u = w.act[p];
            if (u == xs)
                continue;
            bool exact;
            dist_int c = xs == STATE_L ? pricer.from_l(u) : pricer.price_at(xs, p, exact);
            if (is_inf(c))
                continue;
            dist_int cand = lx + c;
            int32_t m = ds.mate[u];
            if (m == MATE_R) {
                if (cand < best) {
                    best = cand;
                    best_from = xs;
                    best_hop = u;
                }
            } else if (cand < w.label[m]) {
                w.label[m] = cand;
                w.pred_state[m] = xs;
                w.pred_hop[m] = u;
            }
        }
    }
    if (best <= stop_label) {
        res.gap_lb = best;
        res.censored = false;
    } else {
        // Every unsettled label is >= tau, so the answer is >= tau (CLAUDE.md §3.6).
        res.gap_lb = tau;
        res.censored = true;
    }

    // Reconstruct the walk P_H (only when it finished).
    std::vector<InternalHop> hops;
    int32_t start_state = STATE_L;
    if (!res.censored && best_from != INT32_MIN) {
        auto make_hop = [&](int32_t a, int32_t b) {
            InternalHop h{a, b, 0, 0, true};
            pricer.begin_row(a);
            if (a == STATE_L && b == END_R) {
                h.price = h.dhat = idx.d_lr;
            } else if (a == STATE_L) {
                h.price = pricer.from_l(b);
                h.dhat = idx.h_l[ds.defects[b]];
            } else if (b == END_R) {
                h.price = pricer.to_r(a);
                h.dhat = idx.h_r[ds.defects[a]];
            } else {
                h.price = pricer.hop(a, b, h.exact, &h.dhat);
            }
            return h;
        };
        hops.push_back(make_hop(best_from, best_hop));
        int32_t s = best_from;
        while (s != STATE_L) {
            int32_t p = w.pred_state[s];
            if (p == PRED_FREE_RELAY) {
                start_state = s;
                break;
            }
            hops.push_back(make_hop(p, w.pred_hop[s]));
            s = p;
        }
        std::reverse(hops.begin(), hops.end());
    }

    if (!hops.empty()) {
        // Node sequence z_0 = L, ..., z_last = R, with the type of the edge entering each node.
        struct Step {
            int32_t node;  // defect index, STATE_L or END_R
            bool via_hop;
        };
        std::vector<Step> seq;
        seq.push_back({STATE_L, false});
        if (start_state != STATE_L)
            seq.push_back({start_state, false});
        for (auto& h : hops) {
            seq.push_back({h.to, true});
            if (h.to != END_R) {
                int32_t m = ds.mate[h.to];
                seq.push_back({m >= 0 ? m : END_R, false});
            }
        }
        res.num_hops = (uint32_t)hops.size();
        res.all_hops_exact = std::all_of(hops.begin(), hops.end(), [](const InternalHop& h) { return h.exact; });

        // Simple: every defect at most once.
        std::vector<int32_t> visited;
        for (auto& st : seq) {
            if (st.node >= 0)
                visited.push_back(st.node);
        }
        std::sort(visited.begin(), visited.end());
        res.walk_simple = std::adjacent_find(visited.begin(), visited.end()) == visited.end();

        // Blossom passes: runs inside a blossom S entered and left by hops (n_S of §3.2).
        if (ds.num_blossoms() > 0) {
            std::vector<int32_t> blossoms;
            for (auto v : visited) {
                for (uint32_t a = ds.anc_offset[v]; a < ds.anc_offset[v + 1]; a++)
                    blossoms.push_back(ds.anc_blossom[a]);
            }
            std::sort(blossoms.begin(), blossoms.end());
            blossoms.erase(std::unique(blossoms.begin(), blossoms.end()), blossoms.end());
            auto in_blossom = [&](int32_t v, int32_t b) {
                if (v < 0)
                    return false;
                for (uint32_t a = ds.anc_offset[v]; a < ds.anc_offset[v + 1]; a++) {
                    if (ds.anc_blossom[a] == b)
                        return true;
                }
                return false;
            };
            for (auto b : blossoms) {
                for (size_t p = 0; p < seq.size(); p++) {
                    if (!in_blossom(seq[p].node, b) || (p > 0 && in_blossom(seq[p - 1].node, b)))
                        continue;
                    size_t q = p;
                    while (q + 1 < seq.size() && in_blossom(seq[q + 1].node, b))
                        q++;
                    if (seq[p].via_hop && q + 1 < seq.size() && seq[q + 1].via_hop)
                        res.blossom_pass = true;
                }
            }
        }

        // Upper bound Delta_h on the same walk (CLAUDE.md §3.5).
        if (o.upper_bound) {
            dist_int ub = 0, slack = 0;
            for (auto& h : hops) {
                dist_int d;
                if (h.from == STATE_L && h.to == END_R) {
                    d = idx.d_lr;
                } else if (h.from == STATE_L) {
                    d = idx.h_l[ds.defects[h.to]];
                } else if (h.to == END_R) {
                    d = idx.h_r[ds.defects[h.from]];
                } else if (h.exact) {
                    d = h.dhat;
                } else {
                    d = idx.exact_distance(ds.defects[h.from], ds.defects[h.to], w.astar);
                }
                if (is_inf(d)) {
                    ub = DIST_INF;
                    break;
                }
                if (o.check && d < h.dhat)
                    throw InvariantViolation("d_hat exceeds the exact hop distance");
                slack += d - h.dhat;
                dist_int ra = h.from >= 0 ? ds.r(h.from) : 0;
                dist_int rb = h.to >= 0 ? ds.r(h.to) : 0;
                ub += d - ra - rb;
            }
            res.hop_slack = is_inf(ub) ? -1 : slack;
            res.gap_ub = (res.walk_simple && !is_inf(ub)) ? ub : DIST_INF;
            if (o.check && res.gap_ub < res.gap_lb)
                throw InvariantViolation("gap_ub < gap_lb");
        }

        if (o.keep_walk) {
            auto code = [&](int32_t v) -> int64_t {
                if (v == STATE_L)
                    return -1;
                if (v == END_R)
                    return -2;
                return ds.defects[v];
            };
            for (auto& h : hops)
                res.walk.push_back({code(h.from), code(h.to), h.price, h.exact});
        }
    }
    auto t2 = std::chrono::steady_clock::now();
    res.t_decode_ns = elapsed_ns(t0, t1);
    res.t_gap_ns = elapsed_ns(t1, t2);
}

}  // namespace fastgap
}  // namespace pm
