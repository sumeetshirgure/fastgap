"""Landmark-placement framework (CLAUDE.md §7.1).

A landmark is a set of detectors that all sit at distance 0. Its potential
``h(x) = dist(landmark, x)`` is 1-Lipschitz, so ``|h(u) - h(v)| <= dist(u, v)`` for any choice
of landmarks: placement only affects how tight ``gap_lb`` is, never whether it is valid.

The two boundary sides L and R are always added by the framework, outside any strategy.
"""

from __future__ import annotations

from typing import Any, Dict, List, NamedTuple, Optional, Protocol, Sequence, runtime_checkable

import numpy as np

from fastgap import _cpp_fastgap

DIST_INF = int(_cpp_fastgap.DIST_INF)

SIDE_NONE = 0
SIDE_L = 1
SIDE_R = 2


class LandmarkSource(NamedTuple):
    """A landmark: detector ids at distance 0, plus a label recorded in the index config."""

    nodes: frozenset
    label: str


class InteriorGraph:
    """Read-only view of fastgap's detector graph in internal integer units.

    This is exactly the graph the decoder floods (same edges, same weights). Distances never
    pass through the boundary. ``boundary_side[x]`` is ``SIDE_L``, ``SIDE_R`` or ``SIDE_NONE``.
    """

    def __init__(self, core):
        self._core = core
        arrays = core.graph_arrays()
        self.num_nodes: int = int(core.num_nodes)
        self.num_detectors: int = int(core.num_detectors)
        self.offsets: np.ndarray = arrays["offsets"]
        self.neighbors: np.ndarray = arrays["neighbors"]
        self.weights: np.ndarray = arrays["weights"]
        self.boundary_weight: np.ndarray = arrays["boundary_weight"]
        self.boundary_side: np.ndarray = arrays["boundary_side"]
        self.h_l: np.ndarray = core.h_l
        self.h_r: np.ndarray = core.h_r
        self.normalising_constant: float = core.normalising_constant
        self.median_edge_weight: int = int(core.median_edge_weight)

    def dijkstra(self, sources, cap: Optional[int] = None) -> np.ndarray:
        """Multi-source Dijkstra (all sources at distance 0). Unreachable nodes get DIST_INF."""
        init = [(int(s), 0) for s in sources]
        return self._core.dijkstra(init, -1 if cap is None else int(cap))

    def boundary_nodes(self, side: int) -> np.ndarray:
        return np.flatnonzero(self.boundary_side == side)

    def is_reachable(self, dist: np.ndarray) -> np.ndarray:
        return dist < DIST_INF // 2


@runtime_checkable
class LandmarkStrategy(Protocol):
    """A pluggable landmark-placement strategy.

    ``select`` returns at most ``k`` landmark sources. ``meta`` carries optional extra
    information (the stim DEM under ``"dem"``, sampled hop pairs under ``"hop_pairs"``, ...).
    ``params()`` returns a JSON-serialisable dict that is stored (and hashed) in the index.
    """

    name: str

    def select(
        self,
        graph: InteriorGraph,
        coords: Optional[Dict[int, List[float]]],
        meta: Dict[str, Any],
        k: int,
    ) -> List[LandmarkSource]:
        ...


def strategy_params(strategy) -> Dict[str, Any]:
    params = getattr(strategy, "params", None)
    if callable(params):
        return dict(params())
    return {}


def dedupe(sources: Sequence[LandmarkSource]) -> List[LandmarkSource]:
    """Drops empty and repeated landmark sets, keeping the first occurrence."""
    seen = set()
    out = []
    for s in sources:
        if not s.nodes or s.nodes in seen:
            continue
        seen.add(s.nodes)
        out.append(s)
    return out
