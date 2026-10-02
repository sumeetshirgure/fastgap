# CLAUDE.md — fastgap

This file tells Claude Code how to build and extend **fastgap**, a fork of PyMatching
(https://github.com/oscarhiggott/PyMatching) that returns, alongside every minimum-weight
perfect matching (MWPM) decode, a **certified lower bound on the complementary gap** computed
from the decoder's own final state. It never runs a second decode and never sweeps the empty
spacetime volume.

Read the whole file before writing code. Sections 2–4 are the specification. Section 5 is the
implementation plan, section 7 is the landmark-placement framework that later experiments
will extend, and section 9 lists the invariants that must never regress.

---

## 0. Context and non-goals

- The upstream baseline is PyMatching v2 (sparse blossom; Higgott & Gidney, *Quantum* 9, 1600,
  2025). The fork point used for this plan is upstream commit
  `6f63b2b9474ba0fa7e511fe52bffdce858a06984` (2026-05-21). If you fork a newer commit, re-check
  the file and function names in §5.1 before relying on them.
- The earlier SpecMatching work (a speculative, truncated-dual decoder) is **not** used here.
  At d ≤ 17 it gives no latency gain over PyMatching, and applications such as code
  concatenation cannot discard shots. fastgap therefore runs the ordinary, untruncated
  PyMatching decode and adds the gap bound after it. Do not port SpecMatching code.
- Non-goals for v1: correlated matching (`enable_correlations=True`), negative edge weights,
  codes whose observable can be flipped by a closed loop that never touches a boundary (see
  §2.3), and hardware/FPGA work.

---

## 1. Glossary (use these names in code, comments and output files)

| Name in docs | Code name | Meaning |
|---|---|---|
| detector graph `G` | `graph` | PyMatching's matching graph built from a Stim detector error model (DEM). Nodes are detectors; edges are error mechanisms with weight `w_e = ln((1-p_e)/p_e)`; one virtual boundary node `b`. |
| `dist(u,v)` | `dist` | Shortest-path length between detectors `u,v` in `G`, **not** passing through `b`. |
| `dist(u,X)` | `dist_to_side` | For `X ∈ {L,R}` (§2.3): length of the shortest path from `u` that ends with a boundary edge on side `X`. |
| `d_LR` | `d_lr` | `dist(L,R)`: weight of the lightest logical operator (a boundary-to-boundary chain). A matching may include one bare `L`–`R` chain, so `W^c ≤ W* + d_LR`. |
| defect | `defect` | A detector that fired in this shot (a "detection event"). Set `D`. |
| `M*`, `W*` | `mstar`, `w_star` | Minimum-weight matching and its weight. |
| mate | `mate[u]` | The partner of defect `u` in `M*`: another defect, `L`, or `R`. |
| `W^c` | `w_comp` | Minimum weight over matchings in the **other** logical class. |
| gap `Δ` | `gap` | `Δ = W^c − W*` (the complementary gap). |
| `y_v`, `y_S` | `y_leaf`, `y_blossom` | Dual variables: one per defect region, one per blossom region. |
| blossom `S` | `blossom` | An odd set of defects merged by sparse blossom, with `y_S ≥ 0`. Nested or disjoint. |
| total radius `ρ_v` | `rho[v]` | `y_v + Σ_{S ∋ v} y_S`. |
| shared mass `β(u,v)` | `beta(u,v)` | `Σ_{S ∋ u, S ∋ v} y_S` (0 if no common blossom). |
| reduced cost `σ(u,v)` | `sigma` | `dist(u,v) − ρ_u − ρ_v + 2β(u,v) ≥ 0`. Boundary version: `σ(u,X) = dist(u,X) − ρ_u` for `X ∈ {L,R}`. |
| hop / relay | `hop`, `relay` | Along an alternating path: an edge added (in `M^c` only) / an edge of `M*` removed. |
| landmark potential `h_ℓ` | `potential[l]` | `h_ℓ(x) = dist(ℓ, x)` for a landmark `ℓ` (a node or a set of nodes). Precomputed per DEM. |
| `LB(u,v)` | `lb` | `max_ℓ |h_ℓ(u) − h_ℓ(v)| ≤ dist(u,v)`. |
| `Δ_σ,H` | `gap_lb` | The number fastgap reports: a certified lower bound on `Δ`. |
| `Δ_h` | `gap_ub` | Optional upper bound on `Δ` (§3.5). |

Never call the reduced costs or the dual radii "Johnson potentials". The correct term is
**reduced costs** (LP-dual reduced costs, as in Edmonds' blossom algorithm).

---

## 2. Problem setup

### 2.1 Units

PyMatching discretises weights as `w_int = 2 * round(w * c)`:
- `c` comes from `UserGraph::get_edge_weight_normalising_constant`; for non-integral weights it is
  `(NUM_DISTINCT_WEIGHTS − 1) / max|w|` with `NUM_DISTINCT_WEIGHTS = 2^24`;
- the `×2` and the returned `2c` are in `iter_discretized_edges` (`driver/user_graph.h`);
- `to_matching_graph` (`driver/user_graph.cc`) stores the result as `normalising_constant = 2c`.

`c` depends on the **heaviest edge in the graph**. Two graphs that differ by even one edge can
therefore have different internal units.
All radii, distances and reduced costs in fastgap are **integers in these internal units**. Do all
comparisons in integers. Internal distances can reach about `2^25` per edge times the path length,
so:
- store potentials and table distances as `int64`, or as `uint32` **saturated** at a cap `C`;
- saturating is still a valid lower bound, because `|min(a,C) − min(b,C)| ≤ |a − b|`;
- never use 16-bit storage. Convert to user units only at the API boundary by dividing by
`normalising_constant`. Also expose a decibel value, `gap_db = 10*log10(e) * gap`, because the
cultivation and yoking papers bin gaps in dB.

### 2.2 One observable at a time

fastgap computes the gap for one logical observable index `k` per call (default 0). For `K`
observables, call it once per observable. All experiments planned so far use one observable.

### 2.3 Canonicalise the observable onto the boundary (once per DEM)

Stim DEMs can attach an observable flip to interior edges. fastgap first applies a **gauge change**
so that observable `k` appears only on boundary edges.

Work on the **raw DEM error list** (after `decompose_errors`), not on PyMatching's merged graph.
PyMatching's `UserGraph` does two things that would hide problems:
- it merges parallel edges (`merge_edge_or_boundary_edge`, `INDEPENDENT`) and keeps the *first*
  edge's observables;
- it silently drops errors with more than two detectors (`handle_dem_instruction`).

The steps are:

1. Consider the interior graph (all 2-detector edges, no boundary). Assign each detector a bit
   `s(x)` by BFS so that every interior edge `(u,v)` satisfies `obs_k(u,v) = s(u) XOR s(v)`.
2. If some interior cycle has odd observable parity, no such `s` exists. Raise a clear error: the
   code has a logical that never touches a boundary, which v1 does not support.
3. New observable labels: interior edges get `obs'_k = 0`; a boundary edge at `u` gets
   `obs'_k(u,b) = obs_k(u,b) XOR s(u)`.
4. **Parallel-edge check.** If two raw mechanisms share the same detector set but differ in
   `obs'_k`, raise an error. This is a local weight-2 logical: PyMatching's comment on merging
   says "the code has distance 2". In particular, a detector must not have boundary edges on both
   `L` and `R`, because PyMatching would merge them into one boundary edge and the `L`/`R` model
   below could no longer describe the decoder.

For any error set `E` with syndrome `D`, `obs'_k(E) = obs_k(E) XOR (XOR_{x ∈ D} s(x))`. The shift is
a constant per shot, so the two logical classes are relabelled but not changed, and `W*`, `W^c`,
`Δ` are unchanged. When fastgap reports a prediction it must XOR the shift back in.

**One graph for everything.** After canonicalisation, write the result out as a new DEM (the
*canonical DEM*). Build every consumer from that same canonical DEM, so they all see identical
edges, weights and internal units:
- fastgap's PyMatching decode;
- the landmark and table index;
- the exact two-decode baseline, where the observable becomes an extra detector; this is valid
  because the canonical DEM has the observable only on boundary edges, so no 3-detector errors
  appear.

Assert that the baseline graph and fastgap's graph have the same `normalising_constant`. For
externally supplied reference DEMs (the cultivation `gap_dem`), see §7.3.

After canonicalisation the single boundary `b` splits into two:
- `L`: boundary edges with `obs' = 1`;
- `R`: boundary edges with `obs' = 0`.

A matching's class is the parity of its number of `L`-edges. Interior distances no longer depend on
parity. Every experiment DEM we plan to use (surface-code memory, yoked patches, phenomenological
inner patches, the cultivation escape-stage matchable DEM) has boundaries and passes step 2. Add a
test for each.

### 2.4 Decoder state after the first decode

After sparse blossom finishes:

- every defect `v` has a leaf region with radius `y_v`;
- every blossom region `S` has radius `y_S` (PyMatching stores a blossom's *extra* radius, which is
  exactly `y_S`);
- the matching `M*` pairs each defect with another defect or with the boundary, along a path whose
  observable mask is known.

Standard LP facts that the algorithm relies on (Edmonds 1965; Higgott & Gidney 2025):

- **Feasibility:** `σ(u,v) ≥ 0` for every pair and `σ(u,X) ≥ 0` for `X ∈ {L,R}`.
- **Tightness:** `σ = 0` on every edge of `M*`, and on every edge inside a blossom's odd cycle.
- **Strong duality:** `W* = Σ_v y_v + Σ_S y_S`. (PyMatching itself computes the weight this way in
  `shatter_blossom_and_extract_matches`.)
- **Complementary slackness:** if `y_S > 0`, exactly one edge of `M*` crosses `S`, called the
  **base** of `S`.

---

## 3. The algorithm

### 3.1 Why a shortest path

`M* XOR M^c` (the edges in exactly one of the two matchings) decomposes into components. Each
defect has degree 0 or 2 in it (one hop, one relay), and `L`/`R` can have any degree. Components
are cycles, paths from a boundary back to itself, or `L`–`R` paths. After §2.3 only `L`–`R` paths
change the class. Flipping any component costs `Σ hops dist − Σ relays dist ≥ 0`, because `M*` is
optimal. For the optimal `M^c` the difference is therefore a single simple `L`–`R` alternating
path `P`, so

```
Δ = min over simple alternating L–R paths P of [ Σ_{hops} dist − Σ_{relays} dist ].
```

Hops and relays alternate along `P`. The two end edges touch `L` or `R` and may be of either
kind: a defect that `M*` matches to `L` can be released by a relay at the start of the path. A
path with no defects is the bare `L`–`R` chain of weight `d_LR`; this is why `W^c` is defined over
matchings that may include one such chain.

### 3.2 Reduced costs turn it into a non-negative shortest path

Relays are tight, `dist(a,b) = ρ_a + ρ_b − 2β(a,b)`, and every interior defect of `P` lies on one
hop and one relay. Moving `ρ` from each relay onto its neighbouring hops telescopes the sum:

```
cost(P) = Σ_{hops h} σ(h) + 2 Σ_S y_S · n_S(P),
```

Here `n_S(P) ≥ 0` counts the maximal runs of `P` inside blossom `S` that are entered **and** left
by hops. Runs that use the base contribute 0, and a simple path cannot use the base twice. So
every simple `P` satisfies `cost(P) ≥ Σ_hops σ(h)`, with every term non-negative. The proof is
in the talk slides (Blossoms 1–3) and the advisor report.

### 3.3 The search: "stand at x, hop to u, relay to mate(u)"

The search below is written with a generic hop price `c(x,u) ≥ 0`.
- With exact reduced costs, `c = max(σ, 0) = σ`, the answer is called `Δ_σ`.
- With the priced reduced costs `c = σ̂` of §3.4, the answer is `Δ_σ,H`, the number fastgap
  computes.

The walks are of the following kind:

- **State:** "standing at `x`", with `x ∈ D ∪ {L}`. Standing at a defect `x` means we arrived at
  `x` by a relay from `mate(x)` and must leave `x` by a hop.
- **Start:** standing at `L` with cost 0. Also standing at every defect `u` with `mate(u) = L`,
  with cost 0 (a free relay `L ⇒ u`).
- **Move A (hop to a defect):** from `x`, hop to a defect `u ≠ x`, `u ≠ mate(x)`, paying
  `c(x,u)`. Then relay to `m = mate(u)`:
  - if `m` is a defect, the new state is standing at `m`;
  - if `m = R`, the walk is finished;
  - if `m = L`, discard the move (it only returns to the start).
- **Move B (hop to `R`):** from `x`, pay `c(x,R)`. The walk is finished. From `x = L` this
  price is `d_LR = dist(L,R) = min_x (h_L(x) + h_R(x))`, the weight of the lightest logical,
  stored once in the index. It is the answer when a shot has no defects.
- **Answer:** the cost of the cheapest finished walk.

Boundary prices are `c(L,u) = dist(L,u) − ρ_u` and `c(x,R) = dist(x,R) − ρ_x`, always exact. The
graph is complete on `|D|+2` states with non-negative costs, so use **dense (array-based)
Dijkstra**: `O(|D|^2)` time, no heap, and cache-friendly rows. The matrix `P[x][u] = c(x,u)` is
symmetric. Row `x` is read when standing at `x`.

Why this is a lower bound: every simple alternating path is one of these walks, with hop prices
`≤` its true cost. The search also admits walks that reuse a matched pair (land on `u`, later
land on `mate(u)`), and those only lower the minimum. Hence `Δ_σ ≤ Δ` on every shot, blossoms
included.

### 3.4 Pricing a hop without visiting the volume

```
σ̂(u,v) = max( d̂(u,v) − ρ_u − ρ_v + 2β(u,v), 0 )
d̂(u,v) = dist(u,v)                  if v is in u's local table (radius R_tab), exact
        = max(LB(u,v), R_tab + 1)    otherwise (not in the table ⇒ dist > R_tab; integers)
LB(u,v) = max_ℓ |h_ℓ(u) − h_ℓ(v)|
```

The boundary prices `dist(u,L)` and `dist(u,R)` are always exact, because `L` and `R` are
mandatory landmarks (§7.1).

Since `d̂ ≤ dist`, it follows that `σ̂ ≤ σ` and therefore **`Δ_σ,H ≤ Δ_σ ≤ Δ`**. This is the number
fastgap reports as `gap_lb`. The clamp at 0 is safe because `σ ≥ 0`.

### 3.5 Optional upper bound on the same walk

Let `r(u) = ρ_u − β(u, mate(u))`, with `β = 0` when the mate is a boundary. Relays then split
exactly: `dist(u, mate(u)) = r(u) + r(mate(u))`. If the returned walk `P_H` is simple, meaning it
visits each defect at most once (no matched pair is used twice), then

```
Δ ≤ cost(P_H) = Σ_{hops (a,b) of P_H} ( dist(a,b) − r(a) − r(b) ) =: Δ_h,     with r(L) = r(R) = 0.
```

Compute `dist` exactly only for the hops of `P_H`: typically 1–5 targeted queries per shot (§5.4).
Report `gap_ub = Δ_h` when the walk is simple and `+inf` otherwise. The interval
`[gap_lb, gap_ub]` always contains `Δ`. For post-selection at threshold `τ`: keep if
`gap_lb ≥ τ`, discard if `gap_ub < τ`, and otherwise mark the shot uncertain. Log the uncertain
fraction.

### 3.6 Early termination

When a threshold `τ` is given, stop the Dijkstra as soon as the smallest unsettled label is
`≥ τ`, and report `gap_lb = τ` with the flag `censored = True`. Post-selection needs nothing more.

---

## 4. Public API (target)

Python, in a new package `fastgap` that ships in the same distribution as the forked `pymatching`:

```python
import fastgap

idx = fastgap.GapIndex.from_dem(
    dem,                        # stim.DetectorErrorModel (graphlike after decompose_errors=True)
    observable=0,
    strategy="auto",            # landmark strategy name or LandmarkStrategy object (§7)
    num_landmarks=16,
    table_radius=None,          # internal-unit radius of exact local tables; None = auto (§7.4)
)
idx.save(path); idx = fastgap.GapIndex.load(path)   # tables are expensive; cache per DEM

dec = fastgap.GapDecoder(idx, num_threads=1)
res = dec.decode(syndrome, threshold=None, upper_bound=False)
# res.prediction   : uint8 array of observable flips (original, un-gauged labels)
# res.weight       : W* in user units
# res.gap_lb       : Δ_σ,H in user units          res.gap_lb_db : same in dB
# res.gap_ub       : Δ_h or inf (if upper_bound)  res.censored  : bool (early stop)
# res.walk_simple, res.all_hops_exact, res.num_defects, res.num_blossoms
# res.t_decode_ns, res.t_gap_ns                    (per-stage timings)

batch = dec.decode_batch(shots_bit_packed, threshold=None, upper_bound=False)  # arrays of the above

# Exact baseline (two PyMatching decodes on the gauge-fixed DEM with the observable turned into
# an extra detector, the method used in Gidney's cultivation and yoking code):
exact = fastgap.exact_gap_batch(idx, shots_bit_packed)   # returns (prediction, w_star, gap)
```

C++ core under `src/pymatching/fastgap/` (namespace `pm::fastgap`). Keep the upstream
`pymatching` package importable and unmodified in behaviour, so upstream tests still pass and
future rebases stay easy. Do not install upstream `pymatching` in the same environment.

---

## 5. Implementation plan

Work in phases. Each phase ends with its tests green (§9) and a short entry in `NOTES.md`.

### 5.1 Phase 1: extract the dual state from PyMatching

Relevant upstream code:

- `src/pymatching/sparse_blossom/driver/mwpm_decoding.cc`: `process_timeline_until_completion`
  runs sparse blossom to termination. `shatter_blossoms_for_all_detection_events_and_extract_match_edges`
  then walks the region trees and **deletes** the regions (`flooder.region_arena.del`). Both are
  file-local free functions, not declared in `mwpm_decoding.h`. Add declarations there; this is
  the one small upstream edit. Keep it minimal so rebases stay easy.
- `src/pymatching/sparse_blossom/flooder/graph_fill_region.h`: `GraphFillRegion` has
  `radius` (`VaryingCT`), `blossom_parent`, `blossom_parent_top`, `blossom_children` (vector of
  `RegionEdge {region, edge}` around the odd cycle), and `match` (`{region, edge}`; `region ==
  nullptr` means matched to the boundary).
- `src/pymatching/sparse_blossom/matcher/mwpm.cc`: `create_detection_event` allocates one leaf
  region per defect. `shatter_blossom_and_extract_match_edges` produces the final
  `CompressedEdge {loc_from, loc_to, obs_mask}` list, where `loc_to == nullptr` means the
  boundary.

Add `pm::fastgap::extract_dual_state(Mwpm&, const std::vector<uint64_t>& detection_events, DualState&)`:

1. Record the leaf region of each defect when it is created. Add a hook in
   `Mwpm::create_detection_event`, or read `nodes[i].region_that_arrived` right after event
   creation. Do not assume `region_that_arrived` still points at the leaf after flooding
   without testing it.
2. Call `process_timeline_until_completion`.
3. **Before shattering**, walk each defect's leaf → `blossom_parent` chain up to the top.
   - Assert every region is frozen.
   - Read each value with `radius.get_distance_at_time(flooder.queue.cur_time)`. For frozen
     regions this equals `y_intercept()`.
   - Assign dense ids to blossom regions.
   - Store per defect: `y_leaf`, the ancestor list `[(blossom_id, y_S)]`, and `rho`.
   - Store per blossom: its `y_S` and its depth.
   - `β(u,v)` is the summed `y_S` of the common ancestors. With a per-defect ancestor list and
     blossom depths this is a short merge. Most pairs share no blossom, so add a fast path:
     `top(u) != top(v) → β = 0`.
4. Shatter with the existing `shatter_blossoms_for_all_detection_events_and_extract_match_edges`
   to get `mwpm.flooder.match_edges`. From these, build `mate[u]` (defect index, or `L`/`R` after
   gauge-fixing the edge's `obs_mask` bit `k` together with `s(u)`) and each relay's observable
   parity.
5. Take `W* = Σ y` (all leaf and blossom values) from the pre-shatter walk. This is exactly what
   upstream computes: it sums `y_intercept()` while shattering. After shattering to match edges
   the regions are gone and cannot be shattered again, so this is the only point where `W*` can be
   read. The prediction is the XOR of the match edges' `obs_mask`, plus the gauge shift.
6. Handle the upstream corner cases:
   - detection events on user-declared boundary nodes are skipped by upstream; do the same;
   - negative-weight graphs throw an error in v1;
   - an odd component with no boundary makes upstream throw; propagate that error.

Debug invariants (behind `FASTGAP_CHECK`, always on in tests):
- `Σ y == W*` in internal units;
- for every `M*` relay, `ρ_u + ρ_v − 2β(u,v) == dist(u,v)` (take `dist` from a brute-force
  Dijkstra in tests);
- `σ ≥ 0` for all defect pairs in small tests.

### 5.2 Phase 2: exact reduced-cost Dijkstra (no landmarks yet)

- Implement §3.3 with exact `dist` from a per-shot brute-force multi-source Dijkstra. This is
  slow, and it is the reference `Δ_σ`.
- Implement the exact baseline `Δ`: two PyMatching decodes on the gauge-fixed DEM with the
  observable converted into an extra detector toggled on or off. This mirrors
  `cultiv/_decoding/_pymatching_gap_sampler.py` and `_desaturation_sampler.py` in
  `Strilanc/magic-state-cultivation`, which require the observable to live on boundary edges;
  §2.3 guarantees that.
- Add a brute-force exact `Δ` (enumerate all matchings) for graphs with ≤ 10 defects.

### 5.3 Phase 3: index (tables and landmarks) and the priced search

- `GapIndex` stores:
  - the gauge bits `s(x)`;
  - boundary potentials `h_L`, `h_R`;
  - `K` landmark potentials (`int64`, or `uint32` saturated at a cap; see §2.1);
  - exact local tables: for each detector, a sorted list of `(neighbour id, dist)` for
    `dist ≤ R_tab`, excluding transit through `b`;
  - the strategy config and its hash.
  Build it by Dijkstra over the interior graph; each landmark is a multi-source Dijkstra
  (a set of nodes at distance 0).
- Pair lookup `d̂(u,v)`: search `v` in `u`'s table (binary search, or a small hash for large
  tables), and fall back to `LB`. Store potentials node-major (`K` values contiguous per node),
  so `LB` is two cache lines.
- Price matrix: fill the upper triangle, parallelised by row blocks with a persistent thread pool
  created once per `GapDecoder`. Never spawn threads per shot: that alone costs microseconds.
- Dense Dijkstra: single-threaded, over `|D|+2` states, with predecessor tracking for path output.

### 5.4 Phase 4: upper bound and diagnostics

- Exact distance for a hop: build the upstream `SearchFlooder` (`get_mwpm_with_search_graph`) and
  use its shortest-path routine between two detectors. Alternatively, write an A* on the flooder
  graph using `LB` as the heuristic; it is consistent because `LB` is 1-Lipschitz. Either way,
  exclude transit through the boundary, and verify against brute-force Dijkstra.
- Per-shot diagnostics, kept cheap:
  - `walk_simple`;
  - `all_hops_exact` (every hop of `P_H` came from a table);
  - number of hops;
  - `Σ (dist − d̂)` over the hops (only when `upper_bound=True`);
  - `num_blossoms`;
  - whether any hop was a hop-in/hop-out pass through a blossom.

### 5.5 Phase 5: Python bindings, CLI, benchmark

- pybind11 bindings following `user_graph.pybind.cc`. Batch entry points take bit-packed shots,
  like upstream `decode_batch`.
- CLI:
  - `python -m fastgap bench --dem X.dem --shots N --threads T` writes per-shot CSV rows:
    `shot, num_defects, w_star, gap_exact, gap_lb, gap_ub, t_decode_ns, t_gap_ns, t_exact_ns,
    walk_simple, all_hops_exact, censored`;
  - `python -m fastgap build-index --dem X.dem --strategy S --out X.idx`.
- Timing rules:
  - use `std::chrono::steady_clock` around each stage;
  - warm up with 1000 shots;
  - pin threads when the OS allows;
  - report median, p99 and p99.9;
  - report single-thread and `T`-thread numbers, plus core-seconds per shot;
  - report compute-only latency separately from end-to-end latency (syndrome arriving from
    outside the process, answer leaving it). Craig Gidney pointed out that the end-to-end number
    is the relevant one for real-time use.

---

## 6. Correctness rules (must always hold)

1. `gap_lb ≤ gap_exact` on **every** shot of every test DEM. Compare in **integer internal
   units**, after asserting that both come from graphs with identical edges and the same
   `normalising_constant`. A single violation is a bug. Stop and investigate; never loosen the
   test or add a float tolerance.
2. `gap_ub ≥ gap_exact` whenever `walk_simple`.
3. The prediction and `W*` are bit-identical to upstream PyMatching on the same DEM and shot
   (after undoing the gauge shift).
4. Every distance used in pricing is `≤` the true interior distance. Never use a heuristic that
   can overestimate, for example a closed-form lattice formula on a DEM with defects, hook edges
   or postselected detectors, unless it is proven to be a lower bound.
5. Integer arithmetic only inside the core. No floats in pricing or Dijkstra.

---

## 7. Landmark placement is DEM-dependent (framework for future experiments)

Landmark quality decides how tight `gap_lb` is, and its correctness never depends on the choice:
any set of exact potentials gives a valid lower bound. Different DEMs need different strategies,
so placement is a pluggable component, and every experiment records which strategy produced its
numbers.

### 7.1 Interface

```python
class LandmarkStrategy(Protocol):
    name: str
    def select(self, graph: InteriorGraph, coords: dict[int, list[float]] | None,
               meta: dict, k: int) -> list[LandmarkSource]: ...
# LandmarkSource = frozenset[int] of detector ids (distance 0 at every member), plus a label.
```

- `L` and `R` (the two boundary sides from §2.3) are **always** added by the framework, outside
  any strategy. Two reasons:
  - boundary hop prices must be exact;
  - the yoking distance-preservation argument needs `|h_L(u) − h_L(v)| ≤ d̂(u,v)`. Keeping the
    `L` potential in the landmark set guarantees that the bound "sees" the full `L`–`R`
    separation.
- A landmark may be a **set** of nodes (a face, a time slice). Distance-to-a-set is 1-Lipschitz,
  so `LB` stays valid.
- The strategy config (name, parameters, chosen sources) is serialised into the index file and
  hashed. Output CSVs carry the hash.

### 7.2 Built-in strategies

| Strategy | Use for | Sources |
|---|---|---|
| `boundary_only` | ablation baseline | nothing beyond `L`, `R` |
| `faces` | surface-code memory, yoked patches | the time faces (first-round and last-round detectors) and the two spatial faces not already covered by `L`/`R`, found from detector coordinates |
| `faces_corners` | default for surface-code memory | `faces` plus the 8 spacetime corners (nearest detector to each bounding-box corner) |
| `fps` | any DEM, needs no coordinates | farthest-point sampling in the graph metric: start from `L ∪ R`, repeatedly add the detector farthest from all current landmarks |
| `greedy_cover` | tuning on a new DEM | from a candidate pool (e.g. `faces_corners ∪ fps(4k)`), greedily add the landmark that most increases `Σ LB` over hop pairs sampled from real shots (the hops of `P_H` on a few thousand shots). This is the "maxcover" idea for ALT landmarks (Goldberg & Harrelson, SODA 2005). |
| `grafted` | magic-state-cultivation escape DEM | see §7.3 |
| `phenomenological` | inner patches for qLDPC concatenation | uniform `L1` lattice: `faces_corners` is usually already near-exact; confirm with the tightness report |

`auto` picks by inspecting the DEM:
- with `coords` and a regular lattice → `faces_corners`;
- with cultivation coordinate annotations (a 5th coordinate, or virtual pair nodes) → `grafted`;
- otherwise → `fps`.

### 7.3 The grafted cultivation DEM

The escape-stage DEM from Gidney's `DesaturationSampler` is a grafted code: a color-code region
made matchable by postselection and **virtual pair ("doublet") nodes**, embedded in a surface
code, plus a transition slab in time.

- The irregular region is about `d1^2` in space by `r1` rounds in time. Distances bend around it,
  so face landmarks alone are loose there.
- Place landmark sets on:
  - the virtual pair nodes (as one set, and optionally per colour);
  - the time slice where the code switches from grafted to matchable;
  - the first and last rounds;
  - the surface-code faces;
  - then fill the rest of the budget with `fps` restricted to the irregular region.
- Use a larger `table_radius` inside the irregular region if memory allows. Exact tables are the
  cheapest way to be tight where the geometry is complicated.
- Build fastgap's graph from **exactly the edges PyMatching keeps** when it loads Gidney's
  desaturation `gap_dem` (`pymatching.Matching.from_detector_error_model(gap_dem).edges()`):
  - edges to the observable detector become `L` boundary edges (`obs' = 1`);
  - ordinary boundary edges become `R` edges;
  - drop the observable-detector node;
  - keep weights exactly.

  Why not reconstruct from the DEM text: `_dem_with_obs_detector` keeps the `L0` targets, and
  PyMatching drops any 3-detector error it creates, so a DEM-level "inverse" would contain edges
  the reference decoder never saw.

  Add a test that fastgap's edge set, weights and `normalising_constant` match the reference
  graph exactly, and run the §2.3 parallel-edge check on it.
- Check how postselected detectors appear: clipped errors become synthetic boundary edges, which
  can sit in the middle of the patch. Potentials computed by Dijkstra on this graph remain valid
  automatically; only coordinate-based shortcuts would break.

### 7.4 Table radius and memory

- `table_radius` is chosen per DEM:
  - default: about 4 edge weights, i.e. `4 · median(edge weight)`;
  - for high-`p` concatenation DEMs, large enough that most overlapping-region pairs are exact.
- Report table memory (bytes per detector) and potential memory (`K × 4` bytes per detector) at
  build time.

### 7.5 Tightness report (run for every new DEM or strategy)

`python -m fastgap landmark-report --dem X.dem --idx X.idx --shots 10000` writes:
- the distribution of `LB/dist` on hop pairs sampled from real shots;
- the fraction of shots with `gap_lb == Δ_σ` (exact-distance reference) and with
  `gap_lb == gap_exact`;
- a histogram of `gap_exact − gap_lb` in dB, split into:
  - the `LB` part: `Δ_σ − Δ_σ,H`;
  - the alternation/blossom part: `Δ − Δ_σ`;
- the time per shot.

Future experiments add strategies by implementing `LandmarkStrategy` and registering them in
`fastgap/landmarks/__init__.py`. Never edit the core.

---

## 8. Repository conventions

- New C++ lives in `src/pymatching/fastgap/`, with `*.test.cc` beside each file (GoogleTest, like
  upstream).
- Register new C++ sources in `CMakeLists.txt`; `setup.py` delegates to CMake. The Bazel `BUILD`
  file globs `src/**/*.cc`, so check that new test files do not end up in the library target.
- New Python lives in `src/fastgap/`, with tests in `tests/fastgap/`.
- Keep every upstream test passing: `pytest tests`, plus the C++ test target.
- Update `README.md` with an "fastgap" section once Phase 5 lands.
- Benchmarks go in `benchmarks/fastgap/`, with the DEM generation script, the index config, and
  the raw CSV. Never commit only a summary.

---

## 9. Test plan (all must pass before any experiment uses fastgap)

1. **Tiny exact:** random small DEMs (≤ 10 defects) cover:
   - all-pairs `σ ≥ 0`;
   - relays tight;
   - `Σ y = W*`;
   - `gap_exact` (two-decode) equals the brute-force `Δ`;
   - `gap_lb ≤ Δ_σ ≤ Δ`;
   - `gap_ub ≥ Δ` when `walk_simple`.
2. **Gauge:** gauge-fixing leaves `Δ` unchanged, and flips the prediction exactly by the
   per-shot shift. It raises the documented error on a toroidal/boundaryless DEM.
3. **Blossom-heavy:** rotated surface-code memory (`stim.Circuit.generated`,
   `rotated_memory_x/z`), `d ∈ {3,5,7}`, `p ∈ {5e-3, 1e-2}`, 10^5 shots each. Check
   `gap_lb ≤ gap_exact` on every shot, and log how often blossoms occur.
4. **Scale:** `d ∈ {9,13,17}` at `p ∈ {5e-4, 1e-3}`, 10^5 shots, same soundness check plus
   the latency report.
5. **Upstream parity:** predictions and weights identical to upstream on all of the above.
6. **Index round-trip:** save and load produce bit-identical results.
