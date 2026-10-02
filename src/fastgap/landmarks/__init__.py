"""Landmark strategies (CLAUDE.md §7).

Future experiments add strategies by implementing ``LandmarkStrategy`` and registering them here
with ``register``. The core never needs to change: any set of exact potentials gives a valid
lower bound.
"""

from __future__ import annotations

from typing import Callable, Dict, Union

from fastgap.landmarks.base import (
    DIST_INF,
    SIDE_L,
    SIDE_NONE,
    SIDE_R,
    InteriorGraph,
    LandmarkSource,
    LandmarkStrategy,
    dedupe,
    strategy_params,
)
from fastgap.landmarks.builtin import (
    FPS,
    BoundaryOnly,
    Faces,
    FacesCorners,
    Grafted,
    GreedyCover,
    Phenomenological,
    looks_like_cultivation,
    looks_like_regular_lattice,
)

STRATEGIES: Dict[str, Callable[[], LandmarkStrategy]] = {}


def register(name: str, factory: Callable[[], LandmarkStrategy]) -> None:
    """Registers a strategy factory under ``name`` (usable as ``strategy=name``)."""
    STRATEGIES[name] = factory


for _cls in (BoundaryOnly, Faces, FacesCorners, FPS, GreedyCover, Grafted, Phenomenological):
    register(_cls.name, _cls)


def resolve(strategy: Union[str, LandmarkStrategy], graph: InteriorGraph, coords, meta) -> LandmarkStrategy:
    """Turns a strategy name (or ``"auto"``) into a strategy object (CLAUDE.md §7.2)."""
    if not isinstance(strategy, str):
        return strategy
    if strategy == "auto":
        if looks_like_cultivation(coords, meta):
            return Grafted()
        if looks_like_regular_lattice(graph, coords):
            return FacesCorners()
        return FPS()
    if strategy not in STRATEGIES:
        raise ValueError(f"unknown landmark strategy {strategy!r}; known: {sorted(STRATEGIES)} or 'auto'")
    return STRATEGIES[strategy]()


__all__ = [
    "DIST_INF",
    "SIDE_L",
    "SIDE_NONE",
    "SIDE_R",
    "InteriorGraph",
    "LandmarkSource",
    "LandmarkStrategy",
    "STRATEGIES",
    "register",
    "resolve",
    "dedupe",
    "strategy_params",
    "BoundaryOnly",
    "Faces",
    "FacesCorners",
    "FPS",
    "GreedyCover",
    "Grafted",
    "Phenomenological",
]
