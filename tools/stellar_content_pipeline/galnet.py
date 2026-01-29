"""GalNet content synthesis.

Generates:
- bulletins.json (structured)
- digest.txt (human-readable)

The output is designed to be simple to ingest by the game UI.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple
import hashlib
import random
import textwrap

from .schema import UniverseState, Event, System
from .risk import RiskModel


@dataclass(frozen=True)
class Bulletin:
    id: str
    created_day: float
    expires_day: float
    system_id: int
    faction_id: int
    severity: float  # 0..1
    subject: str
    body: str
    tags: Tuple[str, ...]

    def to_dict(self) -> Dict[str, object]:
        return {
            "id": self.id,
            "created_day": self.created_day,
            "expires_day": self.expires_day,
            "system_id": self.system_id,
            "faction_id": self.faction_id,
            "severity": self.severity,
            "subject": self.subject,
            "body": self.body,
            "tags": list(self.tags),
        }


def _clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def _stable_id(*parts: str) -> str:
    h = hashlib.sha1()
    for p in parts:
        h.update(p.encode("utf-8"))
        h.update(b"\0")
    return h.hexdigest()[:16]


def _wrap(s: str, width: int = 78) -> str:
    return "\n".join(textwrap.fill(line, width=width) for line in s.splitlines())


def _pick(rng: random.Random, seq: Sequence[str]) -> str:
    return seq[rng.randrange(0, len(seq))]


def _station_econ_tag(econ: str) -> str:
    e = econ.lower().strip()
    if e in ("industrial", "mining", "refinery"):
        return "industry"
    if e in ("agricultural", "food"):
        return "agri"
    if e in ("military", "security"):
        return "mil"
    if e in ("research", "science"):
        return "science"
    return "general"


class GalNetSynth:
    def __init__(self, *, seed: Optional[int] = None, risk_model: Optional[RiskModel] = None) -> None:
        self.seed = seed
        self.risk_model = risk_model or RiskModel()

    def generate_bulletins(self, state: UniverseState, *, max_items: int = 14) -> List[Bulletin]:
        rng = random.Random(self.seed if self.seed is not None else int(state.day * 1000.0) & 0xFFFFFFFF)

        # 1) Event-driven bulletins (high priority).
        bulletins: List[Bulletin] = []
        used_sys = set()
        for ev in state.active_events():
            sys = state.system_by_id(ev.system_id)
            if not sys:
                continue
            used_sys.add(sys.id)
            bulletins.append(self._bulletin_from_event(state, ev, sys, rng))

        # 2) Risk advisories for the riskiest systems without events.
        sys_risks: List[Tuple[float, System]] = []
        for s in state.systems:
            if s.id in used_sys:
                continue
            r = self.risk_model.system_risk(state, s.id).total
            sys_risks.append((r, s))
        sys_risks.sort(key=lambda x: (-x[0], x[1].id))

        for r, s in sys_risks[: max(0, min(4, max_items - len(bulletins)))]:
            if r < 0.70:
                break
            bulletins.append(self._risk_advisory(state, s, r, rng))

        # 3) Ambient 'economy' / 'traffic' items to keep the feed alive.
        if len(bulletins) < max_items:
            ambient = self._ambient_bulletins(state, rng, count=max_items - len(bulletins))
            bulletins.extend(ambient)

        # Sort by severity then stable id.
        bulletins.sort(key=lambda b: (-b.severity, b.id))
        return bulletins[:max_items]

    def digest(self, state: UniverseState, bulletins: Sequence[Bulletin], *, max_lines: int = 40) -> str:
        lines: List[str] = []
        lines.append(f"GALNET DAILY DIGEST — Day {state.day:.2f}")
        lines.append("=" * 32)
        lines.append("")

        if not bulletins:
            lines.append("No bulletins available.")
            return "\n".join(lines)

        # Top 1-3 critical items get expanded.
        critical = [b for b in bulletins if b.severity >= 0.78][:3]
        for b in critical:
            sys = state.system_by_id(b.system_id)
            sys_name = sys.name if sys else f"System {b.system_id}"
            lines.append(f"• {b.subject} [{sys_name}]  (sev {b.severity:.2f})")
            lines.append(_wrap(b.body, width=76))
            lines.append("")

        # Remaining items get one-line.
        rest = [b for b in bulletins if b not in critical]
        if rest:
            lines.append("Other headlines:")
            for b in rest[: max(0, max_lines - len(lines) - 2)]:
                sys = state.system_by_id(b.system_id)
                sys_name = sys.name if sys else f"System {b.system_id}"
                lines.append(f"- {b.subject} [{sys_name}]")
        return "\n".join(lines[:max_lines])

    # ---------- internal generators ----------

    def _bulletin_from_event(self, state: UniverseState, ev: Event, sys: System, rng: random.Random) -> Bulletin:
        kind = ev.kind.lower().strip()
        sev = _clamp01(0.55 + 0.45 * ev.severity)

        faction_id = sys.controlling_faction_id
        f = state.faction_by_id(faction_id) if faction_id else None
        faction_name = f.name if f else "Independent"

        if ev.headline.strip():
            subject = ev.headline.strip()
        else:
            subject = self._event_subject(kind, sys.name, rng)

        if ev.details.strip():
            body = ev.details.strip()
        else:
            body = self._event_body(kind, sys.name, faction_name, ev.severity, rng)

        tags = ["event"]
        if "pirate" in kind or "raider" in kind:
            tags += ["security", "pirates"]
        elif "flare" in kind or "storm" in kind:
            tags += ["hazard", "weather"]
        elif "outbreak" in kind or "famine" in kind:
            tags += ["humanitarian"]
        elif "war" in kind or "skirmish" in kind:
            tags += ["conflict"]
        elif "discovery" in kind or "anomaly" in kind:
            tags += ["science"]
        else:
            tags += ["general"]

        bid = _stable_id("ev", ev.id, sys.name, f"{state.day:.3f}")
        return Bulletin(
            id=bid,
            created_day=float(state.day),
            expires_day=float(state.day + max(0.5, ev.duration_days)),
            system_id=sys.id,
            faction_id=faction_id,
            severity=float(sev),
            subject=subject,
            body=_wrap(body, width=76),
            tags=tuple(tags),
        )

    def _risk_advisory(self, state: UniverseState, sys: System, risk: float, rng: random.Random) -> Bulletin:
        risk = _clamp01(risk)
        sev = _clamp01(0.45 + 0.55 * risk)

        f = state.faction_by_id(sys.controlling_faction_id) if sys.controlling_faction_id else None
        issuer = f.name if f else "GalNet Safety Board"

        band = "HIGH" if risk >= 0.85 else "ELEVATED"
        subject = f"Travel Advisory: {band} RISK in {sys.name}"
        reasons: List[str] = []
        if sys.security < 0.35:
            reasons.append("weak local enforcement")
        if sys.pirate_activity > 0.45:
            reasons.append("increased interdiction reports")
        if sys.anomaly_level > 0.40:
            reasons.append("unstable navigational anomalies")
        if not reasons:
            reasons.append("multiple unverified incidents")

        body = (
            f"{issuer} has issued a {band} travel advisory for the {sys.name} system. "
            f"Pilots are urged to review their routes and avoid unnecessary transits. "
            f"Primary contributing factors: {', '.join(reasons)}. "
            f"If transit is required, travel in convoy and keep comms disciplined."
        )

        tags = ("advisory", "security", "navigation")
        bid = _stable_id("adv", sys.name, f"{state.day:.3f}", f"{risk:.3f}")
        return Bulletin(
            id=bid,
            created_day=float(state.day),
            expires_day=float(state.day + 1.25),
            system_id=sys.id,
            faction_id=sys.controlling_faction_id,
            severity=float(sev),
            subject=subject,
            body=_wrap(body, width=76),
            tags=tags,
        )

    def _ambient_bulletins(self, state: UniverseState, rng: random.Random, *, count: int) -> List[Bulletin]:
        items: List[Bulletin] = []
        systems = list(state.systems)
        if not systems:
            return items

        for _ in range(max(0, count)):
            sys = systems[rng.randrange(0, len(systems))]
            # choose a station if possible
            if sys.stations:
                st = sys.stations[rng.randrange(0, len(sys.stations))]
                econ_tag = _station_econ_tag(st.economy)
                subj = _pick(
                    rng,
                    [
                        f"Market Watch: {st.name} posts new contracts",
                        f"Shipyard Update: {st.name} opens new bays",
                        f"Traffic Notice: {st.name} adjusts approach vectors",
                        f"Local Spotlight: {st.name} celebrates Founders Day",
                    ],
                )
                body = _pick(
                    rng,
                    [
                        f"{st.name} in the {sys.name} system announced a fresh slate of commercial contracts. "
                        "Independent captains are advised to check docking restrictions before arrival.",
                        f"Engineering crews at {st.name} report improved turnaround times after a maintenance sprint. "
                        "Expect short-term congestion at peak hours.",
                        f"Approach control at {st.name} updated lane allocations to reduce near-miss incidents. "
                        "Pilots should confirm their assigned vector on final approach.",
                        f"Civic groups at {st.name} are hosting a public celebration. "
                        "Local authorities warn that ceremonial launches may temporarily disrupt traffic.",
                    ],
                )
                sev = 0.12 + 0.15 * rng.random()
                tags = ("local", econ_tag)
                bid = _stable_id("amb", st.name, sys.name, f"{state.day:.3f}", str(rng.randrange(0, 1_000_000)))
                items.append(
                    Bulletin(
                        id=bid,
                        created_day=float(state.day),
                        expires_day=float(state.day + 0.9),
                        system_id=sys.id,
                        faction_id=sys.controlling_faction_id,
                        severity=float(sev),
                        subject=subj,
                        body=_wrap(body, width=76),
                        tags=tags,
                    )
                )
            else:
                subj = _pick(
                    rng,
                    [
                        f"Astronomy Note: Clear skies over {sys.name}",
                        f"Traffic Bulletin: Low volume through {sys.name}",
                        f"Survey Update: Unmanned probes in {sys.name}",
                    ],
                )
                body = _pick(
                    rng,
                    [
                        f"Long-range observatories report unusually stable sensor conditions in {sys.name}. "
                        "Explorers may find improved scan clarity for the next cycle.",
                        f"Traffic controllers report light transit through {sys.name} and recommend it as an alternate corridor.",
                        f"Autonomous probes in {sys.name} relayed routine telemetry. No actionable anomalies detected.",
                    ],
                )
                sev = 0.08 + 0.12 * rng.random()
                tags = ("ambient",)
                bid = _stable_id("amb2", sys.name, f"{state.day:.3f}", str(rng.randrange(0, 1_000_000)))
                items.append(
                    Bulletin(
                        id=bid,
                        created_day=float(state.day),
                        expires_day=float(state.day + 0.8),
                        system_id=sys.id,
                        faction_id=sys.controlling_faction_id,
                        severity=float(sev),
                        subject=subj,
                        body=_wrap(body, width=76),
                        tags=tags,
                    )
                )

        return items

    def _event_subject(self, kind: str, sys_name: str, rng: random.Random) -> str:
        kind = kind.lower().strip()
        if "pirate" in kind:
            return _pick(rng, [f"Pirate Activity Surges in {sys_name}", f"Convoy Attacked Near {sys_name}", f"Interdictions Reported in {sys_name}" ])
        if "flare" in kind or "storm" in kind:
            return _pick(rng, [f"Solar Flare Warning for {sys_name}", f"Radiation Storm Disrupts Sensors in {sys_name}", f"Navigation Hazards Spike in {sys_name}" ])
        if "outbreak" in kind:
            return _pick(rng, [f"Outbreak Containment Measures in {sys_name}", f"Medical Advisory Issued for {sys_name}" ])
        if "famine" in kind:
            return _pick(rng, [f"Relief Effort Mobilizes for {sys_name}", f"Supply Shortages Reported in {sys_name}" ])
        if "war" in kind or "skirmish" in kind:
            return _pick(rng, [f"Clashes Escalate in {sys_name}", f"Tensions Rise Along {sys_name} Lanes" ])
        if "discovery" in kind or "anomaly" in kind:
            return _pick(rng, [f"New Discovery Reported in {sys_name}", f"Researchers Debate Signals from {sys_name}" ])
        return _pick(rng, [f"Situation Developing in {sys_name}", f"Advisory Posted for {sys_name}" ])

    def _event_body(self, kind: str, sys_name: str, faction_name: str, severity: float, rng: random.Random) -> str:
        sev = _clamp01(severity)
        kind = kind.lower().strip()

        if "pirate" in kind:
            return (
                f"Multiple independent captains reported coordinated pirate activity in the {sys_name} system. "
                f"Local patrols under {faction_name} authority have increased random inspections and urged civilians "
                f"to travel in escorted formations. Incident severity is assessed as {sev:.0%}."
            )
        if "flare" in kind or "storm" in kind:
            return (
                f"Astrophysics monitors flagged a radiation event affecting the {sys_name} system. "
                "Expect intermittent sensor noise and degraded jump accuracy. "
                f"Operators recommend delaying long-range transits until conditions stabilize (impact {sev:.0%})."
            )
        if "outbreak" in kind:
            return (
                f"Medical authorities in {sys_name} announced containment measures following a reported outbreak. "
                "Docking protocols may include additional screening and reduced station capacity. "
                f"Public health impact is estimated at {sev:.0%}."
            )
        if "famine" in kind:
            return (
                f"Supply channels in {sys_name} show stress indicators and rising prices. "
                "Relief shipments are being organized, but captains should anticipate tighter trade controls. "
                f"Economic disruption is estimated at {sev:.0%}."
            )
        if "war" in kind or "skirmish" in kind:
            return (
                f"Armed clashes were reported along established corridors in {sys_name}. "
                "Civilians are advised to avoid conflict zones and monitor emergency broadcasts. "
                f"Operational risk is estimated at {sev:.0%}."
            )
        if "discovery" in kind or "anomaly" in kind:
            return (
                f"Survey teams in {sys_name} reported anomalous readings that merit further investigation. "
                "Research outfits are offering rewards for confirmed data. "
                f"Scientific significance is estimated at {sev:.0%}."
            )
        return (
            f"Authorities in {sys_name} issued a general notice citing ongoing developments. "
            f"{faction_name} officials recommend staying alert and avoiding unnecessary travel until more information is available."
        )
