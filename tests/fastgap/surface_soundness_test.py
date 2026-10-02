"""CLAUDE.md §9.3 (blossom-heavy), §9.4 (scale) and §9.5 (upstream parity).

Each configuration decodes SHOTS (default 10^5) shots and checks, in integer internal units:
gap_lb <= gap_exact on every shot, gap_ub >= gap_exact when the walk is simple, and predictions
and weights bit-identical to upstream PyMatching on the original DEM.
"""

import json
import os

import numpy as np
import pymatching
import pytest

import fastgap
from fastgap_helpers import SHOTS, sample, surface_dem

LOG = os.environ.get("FASTGAP_TEST_LOG")


def _check(task, d, p, shots, threads=4, check_shots=None):
    """``check_shots``: how many shots also run the in-decoder invariants (all-pairs sigma >= 0,
    relays tight, d_hat <= dist); None = all of them."""
    dem = surface_dem(task, d, p)
    idx = fastgap.GapIndex.from_dem(dem)
    dets, _ = sample(dem, shots, seed=d * 1000 + int(p * 1e5))
    dec = fastgap.GapDecoder(idx, num_threads=threads)
    dec.parallel_grain = 1  # every shot takes the intra-shot parallel search
    res = dec.decode_batch(dets, upper_bound=True)
    single = fastgap.GapDecoder(idx, num_threads=1).decode_batch(dets, upper_bound=True)
    for field in ("gap_lb_int", "gap_ub_int", "prediction", "weight_int", "walk_simple", "all_hops_exact", "censored"):
        assert np.array_equal(getattr(res, field), getattr(single, field)), field
    n_check = len(dets) if check_shots is None else min(check_shots, len(dets))
    checked = fastgap.GapDecoder(idx, num_threads=threads, check=True).decode_batch(dets[:n_check], upper_bound=True)
    assert np.array_equal(checked.gap_lb_int, res.gap_lb_int[:n_check])
    exact = fastgap.exact_gap_batch(idx, dets, num_threads=threads)

    violations = np.flatnonzero(res.gap_lb_int > exact.gap_int)
    assert len(violations) == 0, f"gap_lb > gap_exact on shots {violations[:10]}"
    ub_bad = np.flatnonzero(res.walk_simple & (res.gap_ub_int < exact.gap_int))
    assert len(ub_bad) == 0, f"gap_ub < gap_exact on simple walks {ub_bad[:10]}"
    assert np.array_equal(res.weight_int, exact.w_star_int)

    matching = pymatching.Matching.from_detector_error_model(dem)
    pred, weights = matching.decode_batch(dets, bit_packed_shots=True, return_weights=True)
    assert np.array_equal(res.prediction, pred)
    assert np.array_equal(res.weight, weights)

    stats = {
        "task": task, "d": d, "p": p, "shots": shots,
        "blossom_shot_fraction": float(np.mean(res.num_blossoms > 0)),
        "blossom_pass_fraction": float(np.mean(res.blossom_pass)),
        "walk_simple_fraction": float(np.mean(res.walk_simple)),
        "gap_lb_equals_exact": float(np.mean(res.gap_lb_int == exact.gap_int)),
        "median_t_decode_ns": float(np.median(res.t_decode_ns)),
        "median_t_gap_ns": float(np.median(res.t_gap_ns)),
        "p99_t_gap_ns": float(np.percentile(res.t_gap_ns, 99)),
        "median_t_exact_ns": float(np.median(exact.t_exact_ns)),
        "strategy_hash": idx.strategy_hash,
    }
    print(json.dumps(stats))
    if LOG:
        with open(LOG, "a") as f:
            f.write(json.dumps(stats) + "\n")
    return stats


@pytest.mark.parametrize("task", ["rotated_memory_x", "rotated_memory_z"])
@pytest.mark.parametrize("d", [3, 5, 7])
@pytest.mark.parametrize("p", [5e-3, 1e-2])
def test_blossom_heavy(task, d, p):
    stats = _check(task, d, p, SHOTS)
    assert stats["blossom_shot_fraction"] > 0


@pytest.mark.parametrize("task", ["rotated_memory_x", "rotated_memory_z"])
@pytest.mark.parametrize("d", [9, 13, 17])
@pytest.mark.parametrize("p", [5e-4, 1e-3])
def test_scale(task, d, p):
    _check(task, d, p, SHOTS, check_shots=2000)
