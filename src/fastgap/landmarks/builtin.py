"""Built-in landmark strategies (CLAUDE.md §7.2)."""

from __future__ import annotations

import itertools
from typing import Any, Dict, Iterable, List, Optional, Sequence

import numpy as np

from fastgap.landmarks.base import (
    DIST_INF,
    SIDE_L,
    SIDE_NONE,
    SIDE_R,
    InteriorGraph,
    LandmarkSource,
    dedupe,
)


def _coord_array(coords: Dict[int, List[float]], num_detectors: int):
    """(ids, array) of detectors that have coordinates, padded to a common dimension."""
    ids = [d for d in range(num_detectors) if d in coords and len(coords[d]) > 0]
    if not ids:
        return np.zeros(0, dtype=np.int64), np.zeros((0, 0))
    dim = max(len(coords[d]) for d in ids)
    arr = np.full((len(ids), dim), np.nan)
    for i, d in enumerate(ids):
        c = coords[d]
        arr[i, : len(c)] = c
    return np.asarray(ids, dtype=np.int64), arr


def _layer_tolerance(values: np.ndarray) -> float:
    """One lattice spacing along an axis: the smallest positive gap between distinct values."""
    distinct = np.unique(values[~np.isnan(values)])
    if len(distinct) < 2:
        return 0.0
    return float(np.min(np.diff(distinct)))


def _spatial_axes(dim: int) -> List[int]:
    # Stim convention: (x, y, t) or (x, y, t, ...): with 3 or more coordinates the third is time.
    if dim >= 3:
        return [0, 1]
    return list(range(dim))


def _time_axis(dim: int) -> Optional[int]:
    return 2 if dim >= 3 else None


def _jaccard(a: frozenset, b: frozenset) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def _boundary_axis(graph: InteriorGraph, faces_by_axis: Dict[int, List[frozenset]]) -> Optional[int]:
    """The spatial axis whose two faces best match the L and R boundary node sets (Jaccard)."""
    sides = [frozenset(graph.boundary_nodes(SIDE_L).tolist()), frozenset(graph.boundary_nodes(SIDE_R).tolist())]
    if not any(sides):
        return None
    best_axis, best_score = None, 0.0
    for a, faces in faces_by_axis.items():
        score = sum(max(_jaccard(f, s) for s in sides) for f in faces)
        if score > best_score:
            best_axis, best_score = a, score
    return best_axis


class BoundaryOnly:
    """Ablation baseline: no landmarks beyond L and R."""

    name = "boundary_only"

    def params(self):
        return {}

    def select(self, graph, coords, meta, k):
        return []


class Faces:
    """Time faces (first/last round) and the spatial faces not already covered by L/R.

    A face is the set of detectors within one lattice spacing of the bounding-box extreme along
    an axis. The spatial axis whose faces best match the L and R boundary node sets is already
    represented by the h_L / h_R potentials and is skipped.
    """

    name = "faces"

    def __init__(self, layers: float = 1.0):
        self.layers = layers

    def params(self):
        return {"layers": self.layers}

    def faces(self, graph: InteriorGraph, coords) -> List[LandmarkSource]:
        if not coords:
            return []
        ids, arr = _coord_array(coords, graph.num_detectors)
        if len(ids) == 0:
            return []
        dim = arr.shape[1]
        out = []
        t_axis = _time_axis(dim)
        if t_axis is not None:
            t = arr[:, t_axis]
            out.append(LandmarkSource(frozenset(ids[t == np.nanmin(t)].tolist()), "face:t=min"))
            out.append(LandmarkSource(frozenset(ids[t == np.nanmax(t)].tolist()), "face:t=max"))
        spatial: Dict[int, List[LandmarkSource]] = {}
        for a in _spatial_axes(dim):
            x = arr[:, a]
            if np.all(np.isnan(x)):
                continue
            tol = _layer_tolerance(x) * self.layers + 1e-9
            spatial[a] = [
                LandmarkSource(frozenset(ids[sel].tolist()), f"face:x{a}={side}")
                for side, sel in (("min", x <= np.nanmin(x) + tol), ("max", x >= np.nanmax(x) - tol))
            ]
        # The faces along the L-R axis are already represented by the h_L / h_R potentials.
        covered = _boundary_axis(graph, {a: [f.nodes for f in fs] for a, fs in spatial.items()})
        for a, fs in spatial.items():
            if a != covered:
                out += fs
        return dedupe(out)

    def select(self, graph, coords, meta, k):
        return self.faces(graph, coords)[:k]


class FacesCorners(Faces):
    """``faces`` plus the spacetime corners (nearest detector to each bounding-box corner)."""

    name = "faces_corners"

    def corners(self, graph: InteriorGraph, coords) -> List[LandmarkSource]:
        if not coords:
            return []
        ids, arr = _coord_array(coords, graph.num_detectors)
        if len(ids) == 0:
            return []
        dims = [a for a in range(arr.shape[1]) if not np.all(np.isnan(arr[:, a]))][:3]
        filled = np.nan_to_num(arr[:, dims], nan=0.0)
        lo, hi = filled.min(axis=0), filled.max(axis=0)
        out = []
        for corner_bits in itertools.product((0, 1), repeat=len(dims)):
            corner = np.where(np.asarray(corner_bits) == 1, hi, lo)
            best = int(ids[int(np.argmin(np.sum((filled - corner) ** 2, axis=1)))])
            label = "corner:" + "".join("+" if b else "-" for b in corner_bits)
            out.append(LandmarkSource(frozenset([best]), label))
        return dedupe(out)

    def select(self, graph, coords, meta, k):
        return dedupe(self.faces(graph, coords) + self.corners(graph, coords))[:k]


class FPS:
    """Farthest-point sampling in the graph metric, starting from L and R.

    Repeatedly adds the detector farthest from all current landmarks. Unreachable detectors count
    as farthest (one landmark per disconnected piece helps there). ``candidates`` restricts the
    detectors that may be chosen; ``initial`` adds extra already-placed landmark sets.
    """

    name = "fps"

    def __init__(self, candidates: Optional[Sequence[int]] = None, initial: Sequence[LandmarkSource] = ()):
        self.candidates = None if candidates is None else np.asarray(sorted(set(candidates)), dtype=np.int64)
        self.initial = list(initial)

    def params(self):
        return {"restricted": self.candidates is not None, "num_candidates": None if self.candidates is None else len(self.candidates)}

    def select(self, graph, coords, meta, k):
        n = graph.num_detectors
        dist = np.minimum(graph.h_l[:n], graph.h_r[:n]).astype(np.int64)
        for s in self.initial:
            dist = np.minimum(dist, graph.dijkstra(s.nodes)[:n])
        allowed = np.ones(n, dtype=bool)
        if self.candidates is not None:
            allowed[:] = False
            allowed[self.candidates[self.candidates < n]] = True
        out = []
        for i in range(k):
            masked = np.where(allowed, dist, -1)
            v = int(np.argmax(masked))
            if masked[v] <= 0:
                break
            out.append(LandmarkSource(frozenset([v]), f"fps:{i}"))
            dist = np.minimum(dist, graph.dijkstra([v])[:n])
            allowed[v] = False
        return out


class GreedyCover:
    """Greedy max-cover over hop pairs sampled from real shots (Goldberg & Harrelson, SODA 2005).

    From a candidate pool (``faces_corners`` plus ``fps(pool_factor * k)``), greedily adds the
    landmark that most increases the summed lower bound ``sum LB(u, v)`` over hop pairs. Hop
    pairs come from ``meta["hop_pairs"]`` (an (m, 2) array of detector ids) or are sampled by
    decoding ``num_shots`` shots of ``meta["dem"]`` with an ``fps`` index and collecting the hops
    of the returned walks.
    """

    name = "greedy_cover"

    def __init__(self, pool_factor: int = 4, num_shots: int = 2000, seed: int = 0):
        self.pool_factor = pool_factor
        self.num_shots = num_shots
        self.seed = seed

    def params(self):
        return {"pool_factor": self.pool_factor, "num_shots": self.num_shots, "seed": self.seed}

    def sample_hop_pairs(self, meta) -> np.ndarray:
        if "hop_pairs" in meta:
            return np.asarray(meta["hop_pairs"], dtype=np.int64).reshape(-1, 2)
        if "dem" not in meta:
            raise ValueError("greedy_cover needs meta['hop_pairs'] or meta['dem'] to sample hop pairs")
        from fastgap.decoder import GapDecoder
        from fastgap.index import GapIndex

        dem = meta["dem"]
        tmp = GapIndex.from_dem(dem, observable=meta.get("observable", 0), strategy="fps",
                                num_landmarks=meta.get("num_landmarks", 16))
        dec = GapDecoder(tmp)
        sampler = dem.compile_sampler(seed=self.seed)
        dets, _, _ = sampler.sample(self.num_shots)
        pairs = []
        for row in dets:
            r = dec.decode(row, keep_walk=True)
            for a, b, _, exact in r.walk:
                if a >= 0 and b >= 0 and not exact:
                    pairs.append((a, b))
        return np.asarray(pairs, dtype=np.int64).reshape(-1, 2)

    def select(self, graph, coords, meta, k):
        pairs = self.sample_hop_pairs(meta)
        pool = dedupe(FacesCorners().select(graph, coords, meta, 10 ** 6) + FPS().select(graph, coords, meta, self.pool_factor * k))
        if len(pairs) == 0:
            return pool[:k]
        u, v = pairs[:, 0], pairs[:, 1]

        def lb_of(h):
            hu, hv = h[u], h[v]
            both_inf = (hu >= DIST_INF // 2) & (hv >= DIST_INF // 2)
            return np.where(both_inf, 0, np.abs(hu - hv))

        current = np.maximum(lb_of(graph.h_l), lb_of(graph.h_r))
        cand_lb = [lb_of(graph.dijkstra(s.nodes)) for s in pool]
        chosen = []
        remaining = list(range(len(pool)))
        for _ in range(min(k, len(pool))):
            gains = [np.sum(np.maximum(current, cand_lb[i]) - current) for i in remaining]
            j = int(np.argmax(gains))
            if gains[j] <= 0:
                break
            i = remaining.pop(j)
            chosen.append(pool[i])
            current = np.maximum(current, cand_lb[i])
        return chosen


class Grafted:
    """Magic-state-cultivation escape DEM (CLAUDE.md §7.3).

    Landmark sets on: the virtual pair ("doublet") nodes (as one set, and optionally per colour),
    the time slice where the code switches from grafted to matchable, the first and last rounds,
    the surface-code faces; then the rest of the budget is filled with ``fps`` restricted to the
    irregular region.

    Inputs (all optional) come from ``meta``: ``virtual_nodes`` (a list of ids, or a dict colour
    -> ids), ``switch_time`` (time coordinate of the switch slice) and ``irregular_region``
    (detector ids). Without them, detectors whose coordinates have a 5th entry are taken as
    virtual pair nodes, grouped by that entry as their colour, and the irregular region is every
    detector inside their coordinate bounding box grown by ``region_margin``.
    """

    name = "grafted"

    def __init__(self, per_colour: bool = True, region_margin: float = 2.0):
        self.per_colour = per_colour
        self.region_margin = region_margin

    def params(self):
        return {"per_colour": self.per_colour, "region_margin": self.region_margin}

    @staticmethod
    def virtual_groups(coords, meta) -> Dict[str, List[int]]:
        v = meta.get("virtual_nodes")
        if isinstance(v, dict):
            return {str(c): sorted(int(x) for x in ids) for c, ids in v.items()}
        if v is not None:
            return {"all": sorted(int(x) for x in v)}
        groups: Dict[str, List[int]] = {}
        for d, c in (coords or {}).items():
            if len(c) >= 5:
                groups.setdefault(f"{c[4]:g}", []).append(int(d))
        return groups

    def select(self, graph, coords, meta, k):
        groups = self.virtual_groups(coords, meta)
        virtual = sorted({x for ids in groups.values() for x in ids})
        out: List[LandmarkSource] = []
        if virtual:
            out.append(LandmarkSource(frozenset(virtual), "virtual:all"))
            if self.per_colour and len(groups) > 1:
                for c, ids in sorted(groups.items()):
                    out.append(LandmarkSource(frozenset(ids), f"virtual:{c}"))
        ids, arr = _coord_array(coords or {}, graph.num_detectors)
        t_axis = _time_axis(arr.shape[1]) if len(ids) else None
        if t_axis is not None and "switch_time" in meta:
            sel = np.isclose(arr[:, t_axis], float(meta["switch_time"]))
            out.append(LandmarkSource(frozenset(ids[sel].tolist()), "slice:switch"))
        out += Faces().faces(graph, coords)
        out = dedupe(out)[:k]

        region = meta.get("irregular_region")
        if region is None and virtual and len(ids):
            pos = {int(d): arr[i] for i, d in enumerate(ids)}
            vpos = np.asarray([np.nan_to_num(pos[v][:3]) for v in virtual if v in pos])
            if len(vpos):
                lo = vpos.min(axis=0) - self.region_margin
                hi = vpos.max(axis=0) + self.region_margin
                filled = np.nan_to_num(arr[:, : vpos.shape[1]])
                inside = np.all((filled >= lo) & (filled <= hi), axis=1)
                region = ids[inside].tolist()
        if len(out) < k:
            out += FPS(candidates=region, initial=out).select(graph, coords, meta, k - len(out))
        return dedupe(out)[:k]


class Phenomenological(FacesCorners):
    """Inner patches for qLDPC concatenation: a uniform L1 lattice, where ``faces_corners`` is
    usually already near-exact (confirm with the tightness report)."""

    name = "phenomenological"


def looks_like_cultivation(coords, meta) -> bool:
    if meta.get("virtual_nodes") is not None:
        return True
    return any(len(c) >= 5 for c in (coords or {}).values())


def looks_like_regular_lattice(graph: InteriorGraph, coords) -> bool:
    if not coords:
        return False
    have = sum(1 for d in range(graph.num_detectors) if d in coords and len(coords[d]) >= 2)
    return have == graph.num_detectors and graph.num_detectors > 0
