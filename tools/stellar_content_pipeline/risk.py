"""Risk model + route planning (A*/Dijkstra) for the content pipeline.

The in-engine code likely has its own navigation systems. This pipeline focuses on
*design-time / offline* risk estimation for:
- narrative impact (GalNet bulletins)
- reward scaling (missions)
- route feasibility (jump range)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple
import heapq
import math

from .schema import UniverseState, System, Event, Vec3


@dataclass(frozen=True)
class RiskWeights:
    # Base security (law/order). Lower security -> higher risk.
    w_security: float = 0.55

    # Pirate pressure (explicit input). Higher -> higher risk.
    w_pirates: float = 0.30

    # Environmental anomaly (explicit input).
    w_anomaly: float = 0.20

    # Reputation penalty (player-hostile jurisdiction).
    w_rep: float = 0.15

    # Active events: generic multiplier.
    w_event: float = 0.35


@dataclass(frozen=True)
class RiskBreakdown:
    total: float
    security_term: float
    pirates_term: float
    anomaly_term: float
    rep_term: float
    event_term: float

    def to_dict(self) -> Dict[str, float]:
        return {
            "total": self.total,
            "security": self.security_term,
            "pirates": self.pirates_term,
            "anomaly": self.anomaly_term,
            "rep": self.rep_term,
            "events": self.event_term,
        }


def _clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def _rep_to_risk(rep: float) -> float:
    """Map [-100,100] reputation to a risk term in [0,1]."""
    rep = max(-100.0, min(100.0, rep))
    # Positive rep reduces risk; negative rep increases.
    # This curve keeps small negative rep meaningful, but saturates.
    if rep >= 0.0:
        return 0.0
    return _clamp01((-rep) / 75.0)


def _event_risk_contrib(ev: Event) -> float:
    """Convert an Event to an additive risk contribution [0,1]."""
    k = ev.kind.lower().strip()
    sev = _clamp01(ev.severity)
    if "pirate" in k or "raider" in k:
        return _clamp01(0.6 * sev + 0.15)
    if "flare" in k or "storm" in k or "radiat" in k:
        return _clamp01(0.5 * sev + 0.10)
    if "outbreak" in k or "famine" in k:
        return _clamp01(0.35 * sev + 0.05)
    if "skirmish" in k or "war" in k:
        return _clamp01(0.55 * sev + 0.10)
    if "discovery" in k or "anomaly" in k:
        return _clamp01(0.40 * sev + 0.10)
    # Unknown events are still relevant, but less so.
    return _clamp01(0.25 * sev)


class RiskModel:
    def __init__(self, weights: RiskWeights = RiskWeights()) -> None:
        self.weights = weights

    def system_risk(self, state: UniverseState, system_id: int) -> RiskBreakdown:
        sys = state.system_by_id(system_id)
        if sys is None:
            raise KeyError(f"Unknown system_id={system_id}")

        # Security term: security=1 -> 0 risk, security=0 -> 1 risk.
        security_term = _clamp01(1.0 - sys.security)
        # Make low security more punishing than high security is rewarding.
        security_term = security_term ** 1.35

        pirates_term = _clamp01(sys.pirate_activity)
        anomaly_term = _clamp01(sys.anomaly_level)

        rep_term = 0.0
        if sys.controlling_faction_id != 0:
            f = state.faction_by_id(sys.controlling_faction_id)
            if f is not None:
                rep_term = _rep_to_risk(f.player_rep)

        # Sum active events for the system. Use diminishing returns.
        active = [e for e in state.active_events() if e.system_id == sys.id]
        ev_sum = 0.0
        for e in active:
            ev_sum += _event_risk_contrib(e)
        event_term = _clamp01(1.0 - math.exp(-ev_sum))  # saturating curve

        w = self.weights
        total = (
            w.w_security * security_term
            + w.w_pirates * pirates_term
            + w.w_anomaly * anomaly_term
            + w.w_rep * rep_term
            + w.w_event * event_term
        )
        total = _clamp01(total)
        return RiskBreakdown(
            total=total,
            security_term=security_term,
            pirates_term=pirates_term,
            anomaly_term=anomaly_term,
            rep_term=rep_term,
            event_term=event_term,
        )


@dataclass(frozen=True)
class RoutePlan:
    systems: Tuple[int, ...]
    total_distance_ly: float
    jumps: int
    risk_avg: float
    risk_max: float
    eta_days: float

    def to_dict(self) -> Dict[str, object]:
        return {
            "systems": list(self.systems),
            "distance_ly": self.total_distance_ly,
            "jumps": self.jumps,
            "risk_avg": self.risk_avg,
            "risk_max": self.risk_max,
            "eta_days": self.eta_days,
        }


def _grid_key(p: Vec3, cell: float) -> Tuple[int, int, int]:
    return (int(math.floor(p.x / cell)), int(math.floor(p.y / cell)), int(math.floor(p.z / cell)))


class RoutePlanner:
    """Plans routes between systems subject to a maximum jump range.

    Implementation details:
    - Build a neighbor graph using a 3D spatial hash grid to avoid O(N^2) adjacency.
    - Use A* with Euclidean distance heuristic.
    """

    def __init__(self, *, max_jump_ly: Optional[float] = None) -> None:
        self.max_jump_ly = max_jump_ly
        self._neighbors: Dict[int, List[Tuple[int, float]]] = {}
        self._built_for_key: Optional[Tuple[int, float]] = None  # (count, max_jump)

    def _ensure_graph(self, state: UniverseState, max_jump: float) -> None:
        key = (len(state.systems), float(max_jump))
        if self._built_for_key == key:
            return

        systems = state.systems
        cell = max_jump
        buckets: Dict[Tuple[int, int, int], List[System]] = {}
        for s in systems:
            k = _grid_key(s.pos_ly, cell)
            buckets.setdefault(k, []).append(s)

        neigh: Dict[int, List[Tuple[int, float]]] = {s.id: [] for s in systems}
        for s in systems:
            k0 = _grid_key(s.pos_ly, cell)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        k = (k0[0] + dx, k0[1] + dy, k0[2] + dz)
                        for t in buckets.get(k, []):
                            if t.id == s.id:
                                continue
                            d = s.pos_ly.dist(t.pos_ly)
                            if d <= max_jump + 1e-9:
                                neigh[s.id].append((t.id, d))

        # Sort neighbors by distance for deterministic tie-breaking.
        for sid, lst in neigh.items():
            lst.sort(key=lambda x: (x[1], x[0]))

        self._neighbors = neigh
        self._built_for_key = key

    def plan(self, state: UniverseState, start_id: int, goal_id: int, *, max_jump_ly: Optional[float] = None) -> Optional[Tuple[int, ...]]:
        max_jump = float(max_jump_ly if max_jump_ly is not None else (self.max_jump_ly if self.max_jump_ly is not None else state.max_jump_ly))
        if max_jump <= 0.0:
            raise ValueError("max_jump_ly must be > 0")

        self._ensure_graph(state, max_jump)

        start = state.system_by_id(start_id)
        goal = state.system_by_id(goal_id)
        if start is None or goal is None:
            return None
        if start_id == goal_id:
            return (start_id,)

        # A* setup
        open_heap: List[Tuple[float, int]] = []
        heapq.heappush(open_heap, (0.0, start_id))
        came_from: Dict[int, int] = {}
        g_score: Dict[int, float] = {start_id: 0.0}

        goal_pos = goal.pos_ly

        def h(nid: int) -> float:
            n = state.system_by_id(nid)
            if n is None:
                return 0.0
            return n.pos_ly.dist(goal_pos)

        visited = set()

        while open_heap:
            _, current = heapq.heappop(open_heap)
            if current in visited:
                continue
            visited.add(current)

            if current == goal_id:
                # reconstruct
                path = [current]
                while current in came_from:
                    current = came_from[current]
                    path.append(current)
                path.reverse()
                return tuple(path)

            for nb, dist in self._neighbors.get(current, []):
                tentative = g_score[current] + dist
                if tentative < g_score.get(nb, float("inf")) - 1e-12:
                    came_from[nb] = current
                    g_score[nb] = tentative
                    f = tentative + h(nb)
                    heapq.heappush(open_heap, (f, nb))

        return None

    def build_route_plan(
        self,
        state: UniverseState,
        start_id: int,
        goal_id: int,
        *,
        risk_model: Optional[RiskModel] = None,
        max_jump_ly: Optional[float] = None,
        # Travel model (rough; narrative only)
        days_per_jump: float = 0.22,
        days_per_ly: float = 0.015,
    ) -> Optional[RoutePlan]:
        route = self.plan(state, start_id, goal_id, max_jump_ly=max_jump_ly)
        if route is None:
            return None

        # Distances & jumps
        total = 0.0
        for a, b in zip(route, route[1:]):
            sa = state.system_by_id(a)
            sb = state.system_by_id(b)
            if sa is None or sb is None:
                continue
            total += sa.pos_ly.dist(sb.pos_ly)
        jumps = max(0, len(route) - 1)

        rm = risk_model or RiskModel()
        risks = [rm.system_risk(state, sid).total for sid in route]
        risk_avg = sum(risks) / max(1, len(risks))
        risk_max = max(risks) if risks else 0.0

        eta = jumps * max(0.0, days_per_jump) + total * max(0.0, days_per_ly)

        return RoutePlan(
            systems=route,
            total_distance_ly=total,
            jumps=jumps,
            risk_avg=float(risk_avg),
            risk_max=float(risk_max),
            eta_days=float(eta),
        )
