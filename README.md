# fastgap

This fork adds **fastgap**, a Python package shipped in the same distribution as `pymatching`.
It decodes each shot with the ordinary sparse-blossom decode and returns a
**certified lower bound on the complementary gap** `Δ = W^c − W*` (the weight difference between
the best matching in the other logical class and the minimum-weight matching). The bound comes
from the decoder's own final dual state. fastgap never runs a second decode and never sweeps the
empty spacetime volume.

Why: the complementary gap is the standard per-shot confidence signal for post-selection in magic-state cultivation
and decoder confidence in yoked surface code and other code concatenation schemes.
The usual way to compute it is to decode a second time, forced into the other logical class.
That second decode costs as much as the first, or more, and grows with the spacetime volume.
fastgap replaces it with a small search over the defects alone.
The search is sound by construction: `gap_lb ≤ Δ` on every shot, checked against the exact
two-decode baseline on 10⁵ shots per test configuration. Predictions and
`W*` stay bit-identical to upstream PyMatching.

```python
import stim, fastgap

circuit = stim.Circuit.generated("surface_code:rotated_memory_x", distance=7, rounds=7,
                                 after_clifford_depolarization=1e-3)
dem = circuit.detector_error_model(decompose_errors=True)

idx = fastgap.GapIndex.from_dem(dem, observable=0, strategy="auto", num_landmarks=16)
idx.save("rmx_d7.idx")                      # tables are the expensive part: cache per DEM
dec = fastgap.GapDecoder(idx, num_threads=4)  # 4 threads share each shot's gap search

shots = circuit.compile_detector_sampler().sample(10_000, bit_packed=True)
res = dec.decode_batch(shots, threshold=None, upper_bound=True)
res.prediction    # identical to pymatching.Matching.decode_batch
res.weight        # W*, identical to upstream
res.gap_lb        # certified lower bound on the gap (nats); res.gap_lb_db in dB
res.gap_ub        # upper bound when the returned walk is simple, else inf

prediction, w_star, gap = fastgap.exact_gap_batch(idx, shots)   # exact two-decode baseline
assert (res.gap_lb_int <= fastgap.exact_gap_batch(idx, shots).gap_int).all()
```

For post-selection at threshold `τ`, pass `threshold=τ`: the search stops as soon as the bound
reaches `τ` (`censored=True`), and `res.classify(τ)` sorts shots into keep, discard and uncertain.

How it works:

1. Observable `k` is gauge-fixed onto the boundary, which splits the boundary into sides `L` and `R`.
2. The radii `y_v`, `y_S` are read before shattering the blossoms.
3. A dense Dijkstra runs over alternating "hop-relay-hop-relay" walks from `L` to `R`,
    priced by reduced costs `σ(u,v) = dist(u,v) − ρ_u − ρ_v + 2β(u,v)`.
    Each hop distance is a lower bound taken from landmark potentials.
4. `num_threads = T` threads cooperate on *each* shot to lower its latency; shots are decoded one
   after another. The calling thread runs the Dijkstra, and the other `T − 1` threads price, in
   advance, the rows of the states it is about to settle. Results are bit-identical for every
   `T`. Use at most the size of the fastest core cluster (e.g. 4–5 on an Apple M5 Pro).

### Gap-computation latency

The baseline runs two PyMatching decodes with the observable as an extra detector, toggled off and
on. Both methods run one ordinary decode, which gives the prediction (min weight matching).
The baseline then pays a second, complementary-class decode for the gap, while fastgap pays its gap search,
which includes `gap_ub`. The table compares only those two costs:
speed-up = baseline complementary decode / fastgap gap stage.

Setup: rotated surface-code memory X, `d` rounds, uniform circuit noise `p`, 2×10⁴ shots per DEM,
compute-only per-shot timers (median / p99), Apple M5 Pro. There were 0 soundness violations.

| DEM | mean defects | baseline gap (complementary decode) | fastgap gap, T=1 | fastgap gap, T=4 | speed-up T=1 (med / p99) | speed-up T=4 (med / p99) |
|---|---|---|---|---|---|---|
| d=5, p=5e-3 | 8.4 | 7.17 / 14.6 µs | 0.75 / 2.5 µs | 0.75 / 2.7 µs | 9.6× / 5.7× | 9.6× / 5.5× |
| d=7, p=5e-3 | 25.3 | 24.67 / 44.5 µs | 2.75 / 7.1 µs | 2.62 / 5.6 µs | 9.0× / 6.3× | 9.4× / 7.9× |
| d=5, p=1e-3 | 1.8 | 3.67 / 10.0 µs | 0.12 / 0.9 µs | 0.12 / 0.9 µs | 29.3× / 10.9× | 29.3× / 10.9× |
| d=9, p=1e-3 | 12.0 | 38.50 / 62.7 µs | 1.50 / 5.9 µs | 1.54 / 5.7 µs | 25.7× / 10.7× | 25.0× / 11.0× |
| d=13, p=1e-3 | 38.3 | 158.71 / 232.8 µs | 6.71 / 14.8 µs | 5.50 / 13.0 µs | 23.7× / 15.8× | 28.9× / 17.9× |
| d=17, p=1e-3 | 87.8 | 535.88 / 808.5 µs | 18.88 / 34.2 µs | 11.92 / 21.9 µs | 28.4× / 23.6× | 45.0× / 37.0× |

The shared ordinary decode costs about the same in both methods (e.g. 7.8 µs for fastgap vs
8.7 µs for the baseline at d=17), so it is left out. Raw per-shot CSVs and summaries are in
`benchmarks/fastgap/results/`; more detail is in `NOTES.md`.

Landmark placement is pluggable (`fastgap.landmarks`). It changes only how tight the bound is,
never whether it is valid.

Command line:

```
python -m fastgap build-index --dem X.dem --strategy auto --out X.idx
python -m fastgap bench --dem X.dem --idx X.idx --shots 100000 --threads 4 --out X.csv --e2e
python -m fastgap landmark-report --dem X.dem --idx X.idx --shots 10000
```

Supported in v1: graphlike DEMs with boundaries, non-negative weights, and one observable per
index. Not supported: correlated matching, and codes with a logical that never touches a
boundary (toric codes raise `fastgap.BoundarylessLogicalError`). The design document is
`CLAUDE.md`. Implementation notes and measured results are in `NOTES.md`, and benchmarks are in
`benchmarks/fastgap/`.

## Attribution

Cooked with Claude Opus 5.5. Paper for `fastgap` is coming soon.

Please also cite the following as our method is based on it.

```
@article{Higgott2025sparseblossom,
  doi = {10.22331/q-2025-01-20-1600},
  url = {https://doi.org/10.22331/q-2025-01-20-1600},
  title = {Sparse {B}lossom: correcting a million errors per core second with minimum-weight matching},
  author = {Higgott, Oscar and Gidney, Craig},
  journal = {{Quantum}},
  issn = {2521-327X},
  publisher = {{Verein zur F{\"{o}}rderung des Open Access Publizierens in den Quantenwissenschaften}},
  volume = {9},
  pages = {1600},
  month = jan,
  year = {2025}
}
```
