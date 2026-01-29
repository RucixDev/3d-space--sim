"""Mission synthesis.

Given:
- UniverseState
- Bulletins (optional)

Generate mission offers with:
- route (systems list)
- ETA estimate
- risk estimate
- reward scaling

This is designed as a *content* layer. It is not a simulation authority.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple
import hashlib
import random

from .schema import UniverseState, System, Station
from .risk import RiskModel, RoutePlanner, RoutePlan
from .galnet import Bulletin


@dataclass(frozen=True)
class Mission:
    id: str
    kind: str
    title: str
    description: str
    origin_system_id: int
    origin_station_id: int
    dest_system_id: int
    dest_station_id: int
    issuer_faction_id: int
    reward_credits: int
    deadline_day: float
    route: Tuple[int, ...]
    risk_avg: float
    risk_max: float
    eta_days: float
    tags: Tuple[str, ...]

    def to_dict(self) -> Dict[str, object]:
        return {
            "id": self.id,
            "kind": self.kind,
            "title": self.title,
            "description": self.description,
            "origin_system_id": self.origin_system_id,
            "origin_station_id": self.origin_station_id,
            "dest_system_id": self.dest_system_id,
            "dest_station_id": self.dest_station_id,
            "issuer_faction_id": self.issuer_faction_id,
            "reward_credits": self.reward_credits,
            "deadline_day": self.deadline_day,
            "route": list(self.route),
            "risk_avg": self.risk_avg,
            "risk_max": self.risk_max,
            "eta_days": self.eta_days,
            "tags": list(self.tags),
        }


def _stable_id(*parts: str) -> str:
    h = hashlib.sha1()
    for p in parts:
        h.update(p.encode("utf-8"))
        h.update(b"\0")
    return h.hexdigest()[:16]


def _pick(rng: random.Random, seq: Sequence[str]) -> str:
    return seq[rng.randrange(0, len(seq))]


def _best_station(sys: System) -> Optional[Station]:
    # Prefer shipyards / outposts for mission origins.
    if not sys.stations:
        return None
    pri = {"shipyard": 0, "outpost": 1, "city": 2, "refinery": 3}
    return sorted(sys.stations, key=lambda s: (pri.get(s.kind.lower(), 9), s.id))[0]


class MissionSynth:
    def __init__(
        self,
        *,
        seed: Optional[int] = None,
        risk_model: Optional[RiskModel] = None,
        route_planner: Optional[RoutePlanner] = None,
    ) -> None:
        self.seed = seed
        self.risk_model = risk_model or RiskModel()
        self.route_planner = route_planner or RoutePlanner()

    def generate(self, state: UniverseState, bulletins: Sequence[Bulletin] = (), *, max_items: int = 18) -> List[Mission]:
        rng = random.Random(self.seed if self.seed is not None else (hash(int(state.day * 1000.0)) & 0xFFFFFFFF))

        missions: List[Mission] = []

        # 1) Bulletin-driven missions (high relevance).
        for b in bulletins:
            m = self._mission_from_bulletin(state, b, rng)
            if m is not None:
                missions.append(m)
            if len(missions) >= max_items:
                break

        # 2) Add opportunistic missions to fill up.
        if len(missions) < max_items:
            missions.extend(self._opportunistic_missions(state, rng, count=max_items - len(missions)))

        # Deterministic sorting by (reward desc, risk desc, id).
        missions.sort(key=lambda m: (-m.reward_credits, -m.risk_avg, m.id))
        # Deduplicate by id (seed changes can still collide rarely).
        out: List[Mission] = []
        seen = set()
        for m in missions:
            if m.id in seen:
                continue
            out.append(m)
            seen.add(m.id)
            if len(out) >= max_items:
                break
        return out

    # ---------- internal ----------

    def _mission_from_bulletin(self, state: UniverseState, b: Bulletin, rng: random.Random) -> Optional[Mission]:
        sys = state.system_by_id(b.system_id)
        if not sys:
            return None

        # choose destination away from origin
        dest = self._pick_far_system(state, sys.id, rng)
        if dest is None:
            return None

        kind = self._kind_from_tags(b.tags)
        return self._build_mission(state, kind, origin=sys, dest=dest, issuer_faction_id=b.faction_id, rng=rng, context_tags=b.tags)

    def _opportunistic_missions(self, state: UniverseState, rng: random.Random, *, count: int) -> List[Mission]:
        out: List[Mission] = []
        systems = list(state.systems)
        if len(systems) < 2:
            return out

        for _ in range(max(0, count)):
            origin = systems[rng.randrange(0, len(systems))]
            dest = self._pick_far_system(state, origin.id, rng)
            if dest is None:
                continue

            # Bias kind based on origin economy / risk.
            r = self.risk_model.system_risk(state, origin.id).total
            if r > 0.72 and rng.random() < 0.5:
                kind = _pick(rng, ["escort", "bounty", "rescue"])
            else:
                kind = _pick(rng, ["courier", "courier", "trade", "scan"])
            out.append(self._build_mission(state, kind, origin=origin, dest=dest, issuer_faction_id=origin.controlling_faction_id, rng=rng, context_tags=()))
        return out

    def _kind_from_tags(self, tags: Tuple[str, ...]) -> str:
        t = set(x.lower() for x in tags)
        if "pirates" in t or "security" in t:
            return "bounty" if "event" in t else "escort"
        if "humanitarian" in t:
            return "relief"
        if "science" in t:
            return "scan"
        if "advisory" in t:
            return "escort"
        return "courier"

    def _pick_far_system(self, state: UniverseState, from_id: int, rng: random.Random) -> Optional[System]:
        origin = state.system_by_id(from_id)
        if origin is None:
            return None
        # prefer systems at > median distance, but still reachable.
        candidates: List[Tuple[float, System]] = []
        for s in state.systems:
            if s.id == origin.id:
                continue
            d = origin.pos_ly.dist(s.pos_ly)
            candidates.append((d, s))
        if not candidates:
            return None
        candidates.sort(key=lambda x: x[0])
        # choose among top 40% farthest
        start = int(len(candidates) * 0.60)
        far = candidates[start:]
        # random among farthest, but ensure route exists
        for _ in range(16):
            _, s = far[rng.randrange(0, len(far))]
            if self.route_planner.plan(state, origin.id, s.id) is not None:
                return s
        # fallback: any reachable
        for _, s in reversed(candidates):
            if self.route_planner.plan(state, origin.id, s.id) is not None:
                return s
        return None

    def _build_mission(
        self,
        state: UniverseState,
        kind: str,
        *,
        origin: System,
        dest: System,
        issuer_faction_id: int,
        rng: random.Random,
        context_tags: Tuple[str, ...],
    ) -> Mission:
        op = self.route_planner.build_route_plan(state, origin.id, dest.id, risk_model=self.risk_model)
        if op is None:
            # unreachable: degrade to local mission
            dest = origin
            op = self.route_planner.build_route_plan(state, origin.id, origin.id, risk_model=self.risk_model)

        assert op is not None

        o_station = _best_station(origin)
        d_station = _best_station(dest)
        o_sid = o_station.id if o_station else 0
        d_sid = d_station.id if d_station else 0

        # Reward model: distance + risk + kind.
        base = 900 + int(2400 * op.total_distance_ly / max(1.0, state.max_jump_ly))
        base += int(1800 * op.risk_avg)
        if kind in ("bounty", "escort", "rescue"):
            base = int(base * 1.35)
        if kind in ("relief", "scan"):
            base = int(base * 1.15)

        # Mild randomization for variety but deterministic across seed.
        base = int(base * (0.90 + 0.25 * rng.random()))
        base = max(300, base)

        # Deadline model: generous for courier/trade, tighter for rescue/bounty.
        urgency = 1.8
        if kind in ("rescue",):
            urgency = 1.15
        elif kind in ("bounty", "escort"):
            urgency = 1.35
        elif kind in ("scan",):
            urgency = 1.6
        deadline = state.day + op.eta_days * urgency + (0.25 + 0.75 * rng.random())

        f = state.faction_by_id(issuer_faction_id) if issuer_faction_id else None
        issuer_name = f.name if f else "Independent"

        title, desc, tags = self._mission_text(kind, origin, dest, o_station, d_station, issuer_name, rng, context_tags, op)

        mid = _stable_id("mission", kind, origin.name, dest.name, f"{state.day:.3f}", title)
        return Mission(
            id=mid,
            kind=kind,
            title=title,
            description=desc,
            origin_system_id=origin.id,
            origin_station_id=o_sid,
            dest_system_id=dest.id,
            dest_station_id=d_sid,
            issuer_faction_id=issuer_faction_id,
            reward_credits=int(base),
            deadline_day=float(deadline),
            route=tuple(op.systems),
            risk_avg=float(op.risk_avg),
            risk_max=float(op.risk_max),
            eta_days=float(op.eta_days),
            tags=tags,
        )

    def _mission_text(
        self,
        kind: str,
        origin: System,
        dest: System,
        o_station: Optional[Station],
        d_station: Optional[Station],
        issuer_name: str,
        rng: random.Random,
        context_tags: Tuple[str, ...],
        plan: RoutePlan,
    ) -> Tuple[str, str, Tuple[str, ...]]:
        o_loc = o_station.name if o_station else origin.name
        d_loc = d_station.name if d_station else dest.name

        if kind == "courier":
            title = _pick(rng, ["Priority Courier", "Sealed Dispatch", "Contract Delivery"]) + f" — {origin.name} → {dest.name}"
            desc = (
                f"Deliver sealed cargo from {o_loc} in {origin.name} to {d_loc} in {dest.name}. "
                f"Issuer: {issuer_name}. Estimated route: {plan.jumps} jumps / {plan.total_distance_ly:.1f} ly. "
                "Discretion advised."
            )
            tags = ("courier", "trade")
        elif kind == "trade":
            title = _pick(rng, ["Bulk Freight", "Commodity Haul", "Industrial Transfer"]) + f" — {origin.name} → {dest.name}"
            desc = (
                f"Move bulk freight from {o_loc} to {d_loc}. "
                f"Issuer: {issuer_name}. Keep to the declared route to qualify for bonus payout."
            )
            tags = ("trade",)
        elif kind == "escort":
            title = _pick(rng, ["Escort Needed", "Convoy Guard", "Protection Contract"]) + f" — {origin.name} → {dest.name}"
            desc = (
                f"Provide armed escort for a civilian transport departing {o_loc}. "
                f"Destination: {d_loc}. Current risk estimate {plan.risk_avg:.2f} (max {plan.risk_max:.2f}). "
                "Maintain comms discipline and stay within escort envelope."
            )
            tags = ("security", "escort")
        elif kind == "bounty":
            title = _pick(rng, ["Bounty Posted", "Interdiction Target", "Raid Response"]) + f" — {dest.name}"
            desc = (
                f"Track and neutralize hostile actors operating near {dest.name}. "
                f"Staging point: {o_loc}. Issuer: {issuer_name}. "
                f"Threat index: {plan.risk_max:.2f}. Evidence package required for payout."
            )
            tags = ("security", "bounty")
        elif kind == "rescue":
            title = _pick(rng, ["Distress Response", "Rescue Contract", "Emergency Extraction"]) + f" — {dest.name}"
            desc = (
                f"Respond to an active distress signal reported near {dest.name}. "
                f"Launch from {o_loc}. Issuer: {issuer_name}. "
                "Time-sensitive — delay may void salvage rights."
            )
            tags = ("rescue", "urgent")
        elif kind == "relief":
            title = _pick(rng, ["Relief Shipment", "Humanitarian Aid", "Supply Run"]) + f" — {origin.name} → {dest.name}"
            desc = (
                f"Transport relief supplies from {o_loc} to {d_loc}. "
                f"Issuer: {issuer_name}. "
                "Non-lethal contract; defensive armament permitted."
            )
            tags = ("humanitarian", "courier")
        elif kind == "scan":
            title = _pick(rng, ["Survey Mission", "Anomaly Scan", "Science Contract"]) + f" — {dest.name}"
            desc = (
                f"Perform a high-resolution sensor sweep in {dest.name} and return data. "
                f"Departure: {o_loc}. Issuer: {issuer_name}. "
                "Bring extended-range scanners. Avoid contaminating the signal."
            )
            tags = ("science", "scan")
        else:
            title = f"Contract — {origin.name} → {dest.name}"
            desc = f"A general contract issued by {issuer_name}."
            tags = ("general",)

        # Add 'risk' tag if relevant.
        tag_list = list(tags)
        if plan.risk_avg >= 0.7:
            tag_list.append("high-risk")
        elif plan.risk_avg >= 0.45:
            tag_list.append("moderate-risk")

        # Preserve context tags if meaningful.
        for ct in context_tags:
            cl = ct.lower()
            if cl in ("pirates", "conflict", "hazard") and cl not in tag_list:
                tag_list.append(cl)

        return title, desc, tuple(tag_list)
