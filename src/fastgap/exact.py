"""Exact references: the two-decode baseline and brute-force enumeration (CLAUDE.md §5.2)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from fastgap import _cpp_fastgap
from fastgap.decoder import pack_shots
from fastgap.index import GapIndex


@dataclass
class ExactBatch:
    """Exact gaps. ``prediction`` is observable k in the original labels; weights in user
    units; ``*_int`` in internal units (use these for soundness comparisons)."""

    prediction: np.ndarray
    w_star: np.ndarray
    gap: np.ndarray
    w_star_int: np.ndarray
    gap_int: np.ndarray
    t_exact_ns: np.ndarray

    def __iter__(self):
        # Allows ``prediction, w_star, gap = exact_gap_batch(...)``.
        return iter((self.prediction, self.w_star, self.gap))


def _baseline(idx: GapIndex, num_threads: int):
    core = idx._exact_cache.get(num_threads)
    if core is None:
        core = _cpp_fastgap.ExactBaselineCore(idx._core, int(num_threads))
        if core.normalising_constant != idx.normalising_constant:
            raise AssertionError("exact baseline graph and fastgap graph have different internal units")
        idx._exact_cache[num_threads] = core
    return core


def exact_gap_batch(idx: GapIndex, shots, *, bit_packed_shots: bool = True, num_threads: int = 1) -> ExactBatch:
    """Exact gap by two PyMatching decodes on the gauge-fixed graph, with the observable turned
    into an extra detector that is off / on (the method of Gidney's cultivation and yoking
    samplers). Built from the same canonical graph as fastgap, with the same internal units.

    Returns an :class:`ExactBatch`, which unpacks as ``(prediction, w_star, gap)``.
    """
    packed = pack_shots(shots, idx.num_detectors, bit_packed_shots)
    r = _baseline(idx, num_threads).decode_batch(packed, idx.num_detectors)
    return ExactBatch(
        prediction=r["prediction"],
        w_star=r["w_star"] / idx.normalising_constant,
        gap=idx.to_user_units(r["gap"]),
        w_star_int=r["w_star"],
        gap_int=np.where(r["gap"] < 0, np.inf, r["gap"].astype(np.float64)),
        t_exact_ns=r["t_exact_ns"],
    )


def brute_force_gap(idx: GapIndex, detection_events) -> dict:
    """Exact gap by enumerating all matchings (at most 20 defects); internal units."""
    events = np.asarray(detection_events, dtype=np.uint64).reshape(-1)
    return dict(_cpp_fastgap.brute_force_gap(idx._core, events))
