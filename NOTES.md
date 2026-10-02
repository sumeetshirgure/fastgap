# fastgap implementation notes

Fork point: upstream PyMatching `6f63b2b9474ba0fa7e511fe52bffdce858a06984`. The file and function
names in CLAUDE.md §5.1 were checked against this commit and are unchanged.

Environment used for every number below: Apple M-series laptop (arm64, 14 cores), macOS 15,
Apple clang 16, Python 3.12, stim 1.16. The timer is `std::chrono::steady_clock`, which ticks in
steps of ~41.7 ns on this machine, so sub-microsecond medians are quantised. macOS has no hard
thread-affinity API, so `--pin` is a no-op here; pinning only takes effect on Linux.

## Layout

| Path | What |
|---|---|
| `src/pymatching/fastgap/canonical.{h,cc}` | §2.3 gauge fixing on the raw flattened DEM, the parallel-edge check, canonical and baseline DEMs/graphs |
| `src/pymatching/fastgap/interior_graph.{h,cc}` | CSR copy of the decoder's own `MatchingGraph` (integer weights, sides L/R), Dijkstra, A* |
| `src/pymatching/fastgap/dual_state.{h,cc}` | §5.1 `extract_dual_state`: dual radii before shattering, then the upstream shatter for `mate` |
| `src/pymatching/fastgap/gap_index.{h,cc}` | §5.3 index: `h_L`, `h_R`, landmark potentials, exact local tables, save/load |
| `src/pymatching/fastgap/gap_decoder.{h,cc}` | §3.3–3.6 dense Dijkstra over hop/relay states, the upper bound, diagnostics, batch decoding |
| `src/pymatching/fastgap/exact.{h,cc}` | §5.2 two-decode baseline and the brute-force subset DP |
| `src/pymatching/fastgap/thread_pool.h` | persistent pool, one per decoder |
| `src/pymatching/fastgap/fastgap.pybind.cc` | `fastgap._cpp_fastgap` (internal integer units only) |
| `src/fastgap/` | Python API (`GapIndex`, `GapDecoder`, `exact_gap_batch`), landmark strategies, CLI |
| `tests/fastgap/` | §9 tests in Python; C++ tests are the `*.test.cc` files next to each source |
| `benchmarks/fastgap/` | DEM generator, index config, run script, raw CSVs and summaries |

The only upstream edit is the pair of declarations added to `driver/mwpm_decoding.h`. Build
plumbing changed too: `CMakeLists.txt` adds `FASTGAP_*` source lists, the `_cpp_fastgap` module
and `-DFASTGAP_CHECK` on the test target; `setup.py` builds both extensions from one CMake tree;
`BUILD` adds `-DFASTGAP_CHECK` to `cc_test`. Upstream tests still pass: 116 upstream C++ tests,
92 upstream Python tests (13 skipped by upstream for optional dependencies).

Note for macOS (case-insensitive file systems): `cmake -B build` collides with Bazel's `BUILD`
file. Use another directory name, e.g. `cmake -B cmake-build` (already git-ignored).

## Phase 1: dual state (§5.1)

- `extract_dual_state` validates and indexes the events first, then calls the upstream
  `process_timeline_until_completion`. Before shattering, it walks every defect's leaf →
  `blossom_parent` chain. It checks that each region is frozen, that
  `get_distance_at_time(now) == y_intercept()`, and that `region_that_arrived_top` agrees with
  the chain. Each defect stores its ancestors top-first with prefix sums, so `β(u,v)` is a
  common-prefix walk; when the top blossoms differ it returns 0 immediately.
- Leaf identification: the leaf is read from `nodes[d].region_that_arrived` *after* flooding,
  without hooking `create_detection_event`. The code checks on every shot that this region has
  no children and that `reached_from_source == &node`. The checks are cheap and have not fired on
  any of the several million shots decoded by the test suite: the source node of a detection
  event is never released by its own leaf region.
- `W* = Σ y` is read before shattering. Tests compare it bit-for-bit with upstream's weight on
  random grids (blossom-heavy) and on surface codes.
- Upstream corner cases: events on user-declared boundary nodes are skipped exactly as upstream
  does. Negative weights are rejected (p > 0.5 at canonicalisation; any negative-weight state in
  the Mwpm). "No perfect matching" propagates, and the Mwpm stays usable afterwards (tested). Bad
  shots (duplicate or out-of-range events) are rejected before flooding, so the Mwpm is never
  left half-processed.

## Phase 2: references (§5.2)

- `Pricing::EXACT` runs the same search with exact `dist`, computing one Dijkstra row per settled
  state. This gives `Δ_σ`.
- `ExactBaseline` runs two upstream decodes on the baseline graph, with the observable detector
  off and on. The baseline graph is built from the same canonical spec. Its normalising constant
  is asserted equal to the decode graph's, both when the index is built and when the baseline is
  constructed.
- `brute_force_gap`: a subset DP over ≤ 20 defects in the L/R model, including the bare L–R
  chain. On random small DEMs it equals the two-decode baseline in both `W_even` and `W_odd`.

## Phase 3: index and priced search (§5.3)

- Canonicalisation works on `dem.flattened()`, splitting components at separators, as
  PyMatching's own iterator does. The canonical DEM keeps every instruction in the same order and
  changes only the observable-k labels of graphlike components, so PyMatching builds a
  structurally identical graph. Merges and the first-edge observables are the same, and so are
  `M*` and the tie-breaking. Upstream PyMatching itself parses `str(dem)`, which is not an exact
  round trip of the probabilities, so fastgap also builds from `str(dem)`. The index stores that
  text and re-canonicalises it on load.
- The distance graph is read from the decoder's own `MatchingGraph`, so its edges, weights and
  units match the decoder's by construction. A 64-bit fingerprint of the graph is stored in the
  index file and checked on load.
- Potentials: exact int64 in the index and the file. Pricing uses a saturated int32 copy (cap
  2³¹−1, §2.1), which halves the work and lets the bound vectorise; saturation keeps the bound
  valid.
- Tables: a truncated Dijkstra from every detector, stored as CSR sorted by neighbour id. The
  default radius is 4 × the median edge weight (§7.4).
- Search: dense Dijkstra over `|D|+2` states, as in §3.3. Implementation details that matter for
  speed (prices are identical to the plain per-pair formula; tests compare against it):
  - Only *active* hop targets are priced. Once state `m` is settled, `mate(m)` is swap-removed
    from a compacted list, since hopping into it can no longer improve anything.
  - Each settled row is priced in one pass. The landmark bound is computed over a
    potential-major block of the active targets. The row's exact table is then scanned once
    through a node→defect map, instead of one binary search per pair.
  - With more than one thread, a single-shot `decode` fills the price matrix by row blocks on the
    persistent pool. Batch decoding spreads shots over the pool instead. The pool is created once
    per decoder, and a nested `run` falls back to serial.
- Speed history at d=17, p=1e-3 (88 defects per shot on average; median gap-stage time per shot):

  | step | time |
  |---|---|
  | first version: per-pair binary search plus int64 LB | 58.5 µs |
  | row scan of the table, gathered potentials | 37.6 µs |
  | int32 saturated potentials | 31.8 µs |
  | active-target compaction | 27.9 µs |

  For comparison, the decode is about 9.5 µs and the exact baseline about 600 µs. The remaining
  time splits as follows:
  - ~15 µs is the dense O(|D|²) skeleton itself (measured with no landmarks and no tables);
  - ~10 µs is table scans;
  - ~4 µs is the landmark bound.
- **Rejected idea, recorded so nobody retries it:** A* with the potential `h_R − ρ`. The intended
  consistency argument needs an *upper* bound on `h_R(u)` across the relay `u ⇒ m`, but tightness
  only gives `h_R(u) ≤ h_R(m) + ρ_u + ρ_m`, which leaves a slack of 2(ρ_u + ρ_m). Relays are free
  in reduced-cost space, so there is no consistent goal potential of this form. The §9 soundness
  test (`gap_lb ≤ gap_exact`) failed immediately on surface codes, and the change was reverted.

## Phase 4: upper bound and diagnostics (§5.4)

- Exact hop distances use an A* on the interior graph with the landmark bound as heuristic.
  The heuristic is consistent because a max of 1-Lipschitz potentials is 1-Lipschitz, saturated
  or not. It is verified against brute-force Dijkstra in `interior_graph.test.cc`. Table hits are
  reused, and L/R hops use `h_L`/`h_R`. These are the only targeted distance queries.
- Diagnostics per shot: `walk_simple`, `all_hops_exact`, `num_hops`, `hop_slack = Σ(dist − d̂)`
  (only with `upper_bound`), `num_blossoms`, and `blossom_pass` (a run inside a blossom entered
  and left by hops, i.e. `n_S > 0`).
- Debug invariants (`check=True`; always on in the C++ tests through `FASTGAP_CHECK`):
  - boundary reduced costs are ≥ 0;
  - relays to L and R are tight;
  - for every defect-defect relay, `ρ_u + ρ_v − 2β == dist` (brute-force Dijkstra);
  - all-pairs `σ ≥ 0` (all pairs when there are ≤ 64 defects; for larger shots, all pairs among
    the first 64);
  - `d̂ ≤ dist`, with equality on table hits;
  - `gap_ub ≥ gap_lb`.

## Phase 5: bindings, CLI, benchmark (§5.5)

- `python -m fastgap bench | build-index | landmark-report | serve`. The bench CSV has exactly
  the columns of §5.5 plus a trailing `strategy_hash` column, so every row carries the hash (§7.1).
- Timing:
  - compute-only per-stage timers in C++;
  - 1000 warm-up shots;
  - median, p99 and p99.9;
  - single-thread and T-thread wall time, and core-seconds per shot;
  - end-to-end latency (`--e2e`): a `python -m fastgap serve` subprocess reads bit-packed shots
    from a pipe and writes answers back. The bench measures the round trip per shot, so this
    number includes process-boundary I/O and the Python wrapper, as Craig Gidney asked for.

## Results

### Soundness (§9.3, §9.4, §9.5): 10⁵ shots per configuration

Every configuration below passed all of these checks:
- `gap_lb ≤ gap_exact` in internal integer units on every shot;
- `gap_ub ≥ gap_exact` on every simple walk;
- predictions and weights bit-identical to `pymatching.Matching.decode_batch` on the original DEM;
- the in-decoder invariants on every shot for d ≤ 7, and on the first 2000 shots for d ≥ 9.

The raw per-configuration stats are in `benchmarks/fastgap/results/soundness_1e5.jsonl`. Strategy
`auto` (= `faces_corners` here), K ≤ 16, default table radius. Times are the per-shot
compute-only stage timers, measured while 4 worker threads decoded concurrently.

| basis | d | p | blossom shots | blossom passes | walk simple | gap_lb == exact | t_decode (median) | t_gap (median / p99) | t_exact (median) |
|---|---|---|---|---|---|---|---|---|---|
| X | 3 | 0.005 | 1.4% | 0.01% | 100.00% | 100.0% | 0.12 µs | 0.17 / 1.2 µs | 1.0 µs |
| Z | 3 | 0.005 | 0.9% | 0.02% | 100.00% | 100.0% | 0.08 µs | 0.12 / 1.0 µs | 1.0 µs |
| X | 5 | 0.005 | 9.7% | 0.22% | 99.81% | 99.6% | 1.00 µs | 1.21 / 4.6 µs | 9.7 µs |
| Z | 5 | 0.005 | 8.2% | 0.22% | 99.88% | 99.6% | 1.00 µs | 1.08 / 4.2 µs | 8.9 µs |
| X | 7 | 0.005 | 29.1% | 0.93% | 98.95% | 95.4% | 4.04 µs | 5.17 / 13.8 µs | 32.7 µs |
| Z | 7 | 0.005 | 26.4% | 0.88% | 99.19% | 95.1% | 3.92 µs | 4.42 / 11.5 µs | 30.8 µs |
| X | 3 | 0.01 | 3.6% | 0.07% | 99.98% | 99.9% | 0.25 µs | 0.29 / 1.6 µs | 1.2 µs |
| Z | 3 | 0.01 | 3.0% | 0.04% | 100.00% | 100.0% | 0.21 µs | 0.25 / 1.3 µs | 1.2 µs |
| X | 5 | 0.01 | 24.5% | 0.65% | 99.12% | 98.7% | 2.83 µs | 2.50 / 7.5 µs | 13.2 µs |
| Z | 5 | 0.01 | 23.3% | 0.77% | 99.38% | 98.9% | 2.79 µs | 2.21 / 6.8 µs | 12.5 µs |
| X | 7 | 0.01 | 66.3% | 3.01% | 95.77% | 94.0% | 11.75 µs | 10.29 / 24.4 µs | 43.2 µs |
| Z | 7 | 0.01 | 64.9% | 2.96% | 96.76% | 94.9% | 11.33 µs | 9.00 / 21.5 µs | 41.3 µs |
| X | 9 | 0.0005 | 0.7% | 0.02% | 99.99% | 88.9% | 0.54 µs | 0.96 / 6.7 µs | 38.3 µs |
| Z | 9 | 0.0005 | 0.7% | 0.01% | 100.00% | 91.3% | 0.54 µs | 0.83 / 6.3 µs | 34.5 µs |
| X | 13 | 0.0005 | 1.8% | 0.04% | 99.97% | 45.3% | 2.00 µs | 6.58 / 19.6 µs | 163.0 µs |
| Z | 13 | 0.0005 | 1.9% | 0.06% | 99.97% | 53.5% | 1.88 µs | 5.83 / 16.8 µs | 150.2 µs |
| X | 17 | 0.0005 | 3.7% | 0.08% | 99.92% | 17.8% | 4.58 µs | 16.17 / 36.1 µs | 498.0 µs |
| Z | 17 | 0.0005 | 3.9% | 0.10% | 99.94% | 25.5% | 4.54 µs | 14.62 / 30.9 µs | 451.0 µs |
| X | 9 | 0.001 | 2.7% | 0.07% | 99.95% | 78.5% | 1.21 µs | 2.71 / 8.8 µs | 48.2 µs |
| Z | 9 | 0.001 | 2.9% | 0.08% | 99.96% | 82.0% | 1.21 µs | 2.33 / 8.3 µs | 43.9 µs |
| X | 13 | 0.001 | 7.7% | 0.20% | 99.83% | 40.0% | 4.08 µs | 11.54 / 23.3 µs | 207.2 µs |
| Z | 13 | 0.001 | 8.0% | 0.20% | 99.87% | 46.0% | 4.08 µs | 10.33 / 21.5 µs | 191.7 µs |
| X | 17 | 0.001 | 15.7% | 0.39% | 99.65% | 24.4% | 10.58 µs | 31.71 / 55.1 µs | 704.1 µs |
| Z | 17 | 0.001 | 16.7% | 0.37% | 99.72% | 31.0% | 9.88 µs | 27.42 / 47.2 µs | 628.6 µs |

Reading the table:
- Blossoms are common at d=7, p=1e-2 (65% of shots), but the bound almost never loses anything
  to them. The `Δ − Δ_σ` part of the shortfall is ~0 on every DEM measured.
- All of the looseness at large d comes from the landmark bound on long hops (`Δ_σ − Δ_σ,H`).

### Tightness vs index parameters (d=13, p=1e-3, rotated memory X, 5000 shots)

| strategy | K | table radius | gap_lb == exact | mean shortfall | p90 shortfall | gap stage (median, before speed-ups) | table memory |
|---|---|---|---|---|---|---|---|
| boundary_only | 0 | 4w | 16.2% | 16.9 dB | 37.7 dB | 12.2 µs | 3.9 MB |
| faces_corners | 12 | 4w | 40.0% | 5.9 dB | 18.1 dB | 13.0 µs | 3.9 MB |
| fps | 16 | 4w | 43.3% | 5.4 dB | 16.5 dB | 12.6 µs | 3.9 MB |
| faces_corners + fps | 32 | 4w | 47.8% | 4.8 dB | 14.8 dB | 14.0 µs | 3.9 MB |
| faces_corners + fps | 64 | 4w | 60.7% | 3.4 dB | 12.5 dB | 17.2 µs | 3.9 MB |
| greedy_cover | 16 | 4w | 58.2% | 2.8 dB | 10.1 dB | 14.2 µs | 3.9 MB |
| faces_corners + fps | 16 | **8w** | **92.7%** | **0.55 dB** | **0.0 dB** | 18.0 µs | 16.3 MB |
| faces_corners + fps | 16 | 0 | 12.8% | 10.4 dB | 25.5 dB | 7.3 µs | 0 |

(`w` = median edge weight. The landmark-report histogram puts the entire shortfall in the LB
part: `Δ_σ == Δ` on > 99.9% of shots.)

Take-aways for future experiments:
1. **The table radius is the main tightness lever at large distance.** The default of 4 edge
   weights (§7.4) is enough for d ≤ 7. For d ≥ 13, memory of about 16 MB at R = 8w buys exact
   bounds on more than 90% of shots.
2. **Among landmark-only options, `greedy_cover` is the best use of a fixed K.** Plain FPS and
   face/corner placements need about 4× more landmarks for the same tightness.
3. Post-selection decisions (`threshold`) are cheap regardless. Early termination stops at τ, and
   `classify(τ)` is never wrong, which is tested against the exact gap.

### Benchmarks

Every DEM in `benchmarks/fastgap/index_config.json` was benchmarked with 10⁴ shots, 1000
warm-up shots, T = 4, the `auto` strategy and the default table radius. Shots are sampled from
the DEM. The raw per-shot CSVs and summaries are in `benchmarks/fastgap/results/`.

| DEM | detectors | mean defects | t_decode med / p99 | t_gap med / p99 / p99.9 | t_exact med | end-to-end med / p99 | 1-thread µs/shot | 4-thread µs/shot | core-µs/shot (4T) | gap_lb == exact | LB/dist on hops (median) | violations |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| rmx_d5_p5e-3 | 120 | 8.4 | 0.96 / 5.8 µs | 1.12 / 4.4 / 8.1 µs | 9.2 µs | 14.5 / 25.7 µs | 2.79 | 0.77 | 3.06 | 99.6% | 1.000 | 0 |
| rmx_d7_p5e-3 | 336 | 25.3 | 3.83 / 14.8 µs | 4.92 / 12.3 / 16.9 µs | 31.8 µs | 23.8 / 43.5 µs | 10.11 | 2.67 | 10.66 | 95.4% | 1.000 | 0 |
| rmx_d5_p1e-3 | 120 | 1.8 | 0.12 / 1.3 µs | 0.12 / 1.5 / 2.5 µs | 4.3 µs | 12.5 / 20.2 µs | 0.55 | 0.19 | 0.75 | 99.9% | 1.000 | 0 |
| rmx_d9_p1e-3 | 720 | 12.1 | 1.21 / 4.2 µs | 2.67 / 8.7 / 12.4 µs | 47.0 µs | 17.0 / 29.4 µs | 4.58 | 1.26 | 5.02 | 78.0% | 0.997 | 0 |
| rmx_d13_p1e-3 | 2184 | 38.1 | 3.92 / 9.3 µs | 11.00 / 22.5 / 31.5 µs | 200.9 µs | 28.3 / 45.5 µs | 16.14 | 4.26 | 17.04 | 38.7% | 0.996 | 0 |
| rmx_d17_p1e-3 | 4896 | 87.8 | 9.29 / 16.8 µs | 29.54 / 50.3 / 60.0 µs | 677.2 µs | 51.3 / 78.4 µs | 40.89 | 10.78 | 43.11 | 24.2% | 1.000 | 0 |

Columns:
- `t_*` are compute-only stage timers measured inside a single-threaded `decode_batch`.
- End-to-end is the round trip of one shot through the `serve` pipe protocol: syndrome in from
  another process, answer out. It is dominated by the pipe and the Python wrapper, roughly 10 µs
  of fixed cost.
- µs/shot is the batch wall time divided by the number of shots.

## Known gaps / not done

- **The real cultivation escape-stage DEM has not been tested.** The DEM comes from Gidney's
  `magic-state-cultivation` repository, which is not available in this environment. What is in
  place:
  - the `grafted` strategy, tested with virtual-node annotations supplied through `meta`;
  - `GapIndex.from_gap_dem`, which builds fastgap's graph from exactly the edges PyMatching keeps
    and asserts that edges, weights and `normalising_constant` match the reference graph;
  - tests of `from_gap_dem` on synthetic gap DEMs that contain 3-detector errors PyMatching drops.

  Run `from_gap_dem` on the real `gap_dem` before using fastgap for cultivation, and check which
  coordinate the virtual pair nodes carry. The `grafted` default assumes a 5th coordinate.
- **Yoked patches** are approximated by two independent patches whose logicals are XORed into one
  observable. It passes canonicalisation and soundness, but it is not Gidney's actual yoked
  circuit.
- **Bazel** is not installed here. The `BUILD` change is a one-line `copts` addition, and the
  globs already exclude `*.test.cc`, `*.test.h` and `*.pybind.cc` from the library target.
- **CPU pinning** is implemented for Linux only.
