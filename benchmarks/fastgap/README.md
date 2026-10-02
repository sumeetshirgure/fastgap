# fastgap benchmarks

Everything here can be regenerated from the repository root, with the package installed:

```
PYTHON=python benchmarks/fastgap/run_bench.sh
```

| File | What |
|---|---|
| `index_config.json` | DEM list (rotated surface-code memory, uniform circuit noise `p`, `d` rounds unless `rounds` is set), index parameters (strategy, K, table radius) and bench parameters (shots, warm-up, threads, seed) |
| `gen_dems.py` | writes `dems/<name>.dem` from the config (git-ignored; deterministic) |
| `run_bench.sh` | for every DEM: `build-index`, `bench` and `landmark-report` |
| `results/<name>.csv` | raw per-shot rows: `shot, num_defects, w_star, gap_exact, gap_lb, gap_ub, t_decode_ns, t_gap_ns, t_exact_ns, walk_simple, all_hops_exact, censored, strategy_hash` |
| `results/<name>.summary.json` | latency percentiles (compute-only and end-to-end), throughput at 1 thread and at each thread count in `bench.threads` (`throughput.by_threads`), core-seconds per shot, soundness count |
| `results/<name>.index.json` | index build report: landmarks, memory, hash |
| `results/<name>.landmarks.json` | tightness report (CLAUDE.md §7.5) |
| `results/soundness_1e5.jsonl` | stats logged by the 10^5-shot §9.3/§9.4 test run (`FASTGAP_TEST_LOG=...`) |

Gaps and weights in the CSV are in nats. The exact gap comes from the two-decode baseline on the
same canonical graph.
