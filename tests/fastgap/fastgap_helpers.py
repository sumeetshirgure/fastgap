"""Shared helpers for the fastgap tests: DEM generators and an exact DEM sampler."""

from __future__ import annotations

import os

import numpy as np
import stim

# CLAUDE.md §9 asks for 10^5 shots per configuration; FASTGAP_TEST_SHOTS lowers it for quick runs.
SHOTS = int(os.environ.get("FASTGAP_TEST_SHOTS", "100000"))


def surface_dem(task: str, d: int, p: float, rounds=None) -> stim.DetectorErrorModel:
    circuit = surface_circuit(task, d, p, rounds)
    return circuit.detector_error_model(decompose_errors=True)


def surface_circuit(task: str, d: int, p: float, rounds=None) -> stim.Circuit:
    return stim.Circuit.generated(
        f"surface_code:{task}",
        distance=d,
        rounds=rounds or d,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p,
    )


def phenomenological_dem(d: int, p: float, rounds=None) -> stim.DetectorErrorModel:
    """Inner-patch style phenomenological noise: data errors between rounds plus measurement
    flips, no circuit-level hook errors."""
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=d,
        rounds=rounds or d,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
    )
    return circuit.detector_error_model(decompose_errors=True)


def yoked_dem(d: int, p: float) -> stim.DetectorErrorModel:
    """Two surface-code patches whose logicals are yoked into one observable (L0 = XOR of both)."""
    a = surface_dem("rotated_memory_x", d, p)
    b = surface_dem("rotated_memory_x", d, p)
    out = stim.DetectorErrorModel()
    out += a.flattened()
    out.append(stim.DemInstruction("shift_detectors", [], [a.num_detectors]))
    out += b.flattened()
    return out.flattened()


def random_grid_dem(rng: np.random.Generator, w: int, h: int, diag: float = 0.4,
                    p_lo: float = 0.02, p_hi: float = 0.25, gauge_density: float = 0.5) -> stim.DetectorErrorModel:
    """A random matchable grid DEM with left (observable) and right boundaries, diagonals (odd
    cycles, hence blossoms), duplicated mechanisms, and observable labels scrambled by a random
    gauge so the observable also sits on interior edges."""
    n = w * h
    s = (rng.random(size=n) < gauge_density).astype(np.int64)
    lines = []

    def emit(dets, obs):
        o = obs
        for x in dets:
            o ^= int(s[x])
        for _ in range(2 if rng.random() < 0.15 else 1):
            p = rng.uniform(p_lo, p_hi)
            lines.append(f"error({p}) " + " ".join(f"D{x}" for x in dets) + (" L0" if o else ""))

    for y in range(h):
        for x in range(w):
            i = y * w + x
            if x + 1 < w:
                emit([i, i + 1], 0)
            if y + 1 < h:
                emit([i, i + w], 0)
            if x + 1 < w and y + 1 < h and rng.random() < diag:
                emit([i, i + w + 1], 0)
        emit([y * w], 1)
        emit([y * w + w - 1], 0)
    for i in range(n):
        lines.append(f"detector({i % w}, {i // w}) D{i}")
    return stim.DetectorErrorModel("\n".join(lines))


def toric_dem(n: int, p: float = 0.05) -> stim.DetectorErrorModel:
    """Code-capacity toric code (one sector): an n x n periodic grid with no boundary. The
    observable is carried by the vertical edges that wrap from the last row to the first, so
    every vertical loop around the torus flips it once."""
    lines = []
    for y in range(n):
        for x in range(n):
            i = y * n + x
            lines.append(f"error({p}) D{i} D{y * n + (x + 1) % n}")
            lines.append(f"error({p}) D{i} D{((y + 1) % n) * n + x}" + (" L0" if y == n - 1 else ""))
    return stim.DetectorErrorModel("\n".join(lines))


def with_obs_detector(dem: stim.DetectorErrorModel) -> stim.DetectorErrorModel:
    """Gidney-style gap DEM: the observable L0 becomes an extra last detector (targets kept on
    every error that flipped the observable)."""
    flat = dem.flattened()
    obs_det = flat.num_detectors
    out = stim.DetectorErrorModel()
    for inst in flat:
        if inst.type != "error":
            if inst.type == "detector":
                out.append(inst)
            continue
        targets = []
        for t in inst.targets_copy():
            if t.is_logical_observable_id():
                targets.append(stim.target_relative_detector_id(obs_det))
            else:
                targets.append(t)
        out.append("error", inst.args_copy(), targets)
    out.append("detector", [], [stim.target_relative_detector_id(obs_det)])
    return out


def sample(dem: stim.DetectorErrorModel, shots: int, seed: int):
    dets, obs, _ = dem.compile_sampler(seed=seed).sample(shots, bit_packed=True)
    return np.ascontiguousarray(dets), obs
