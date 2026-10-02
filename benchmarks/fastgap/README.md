# fastgap benchmarks

Everything here can be regenerated from the repository root, with the package installed:

```
PYTHON=python benchmarks/fastgap/run_bench.sh
```

| File | What |
|---|---|
| `index_config.json` | DEM list (rotated surface-code memory, uniform circuit noise `p`, `d` rounds unless `rounds` is set), index parameters (strategy, K, table radius) and bench parameters (shots, warm-up, threads per shot, seed) |
| `gen_dems.py` | writes `dems/<name>.dem` from the config (git-ignored; deterministic) |
| `run_bench.sh` | for every DEM: `build-index`, `bench` and `landmark-report` |
| `results/<name>.csv` | raw per-shot rows: `shot, num_defects, w_star, gap_exact, gap_lb, gap_ub, t_decode_ns, t_gap_ns, t_exact_ns, walk_simple, all_hops_exact, censored, strategy_hash`, then `t_gap_ns_t<T>` for each thread count `T` in `bench.threads` (`t_gap_ns` is T = 1), then the exact baseline's two decodes timed separately: `t_exact_mstar_ns` (the decode in the class of M*) and `t_exact_comp_ns` (the complementary-class decode, i.e. the baseline's cost of the gap) |
| `results/<name>.summary.json` | per-shot latency percentiles only: compute-only (`latency_ns_by_threads`: decode, gap and total stages, at T = 1 and at each `T` in `bench.threads`) and end-to-end (`end_to_end_ns_by_threads`), the exact baseline's latency (`exact_ns`), its two decodes separately (`exact_mstar_decode_ns`, `exact_comp_decode_ns`), the gap-computation speed-up (`gap_speedup_vs_exact_by_threads`: baseline complementary-decode latency / fastgap gap-stage latency, at the median, p99 and p99.9), plus the soundness count |
| `results/<name>.index.json` | index build report: landmarks, memory, hash |
| `results/<name>.landmarks.json` | tightness report (CLAUDE.md §7.5) |
| `results/soundness_1e5.jsonl` | stats logged by the 10^5-shot §9.3/§9.4 test run (`FASTGAP_TEST_LOG=...`) |

`T` is the number of threads that cooperate on *one* shot's gap search (intra-shot parallelism).
Shots are always decoded one after another, so there are no throughput numbers.

Gaps and weights in the CSV are in nats. The exact gap comes from the two-decode baseline on the
same canonical graph.
