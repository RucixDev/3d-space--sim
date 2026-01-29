"""Simple timeline simulator for generating *emergent* content arcs.

This module is intentionally lightweight; it's not trying to be the authoritative
universe simulation. Instead it:
- evolves a few scalar fields (pirate_activity, anomaly_level) stochastically
- spawns lightweight Events when thresholds are exceeded
- produces a daily GalNet + mission pack

This is extremely useful during development when the in-engine simulation is
still being built: it creates believable, consistent 'feeds' for UI testing.
"""

from __future__ import annotations

from dataclasses import replace
from typing import List, Optional, Tuple
import random

from .schema import UniverseState, System, Event, Faction
from .risk import RiskModel
from .galnet import GalNetSynth, Bulletin
from .missions import MissionSynth, Mission


def _clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def _spawn_event_id(prefix: str, day: float, sys_id: int, rng: random.Random) -> str:
    return f"{prefix}_{sys_id}_{int(day*1000)}_{rng.randrange(0, 1_000_000):06d}"


class TimelineSimulator:
    def __init__(self, *, seed: Optional[int] = None) -> None:
        self.seed = seed
        self.risk_model = RiskModel()

    def simulate(
        self,
        initial: UniverseState,
        *,
        days: int = 7,
        step_days: float = 1.0,
        bulletins_per_day: int = 12,
        missions_per_day: int = 14,
    ) -> List[Tuple[UniverseState, List[Bulletin], List[Mission]]]:
        if days <= 0:
            return []

        rng = random.Random(self.seed if self.seed is not None else (hash(int(initial.day * 1000.0)) & 0xFFFFFFFF))

        state = initial
        out: List[Tuple[UniverseState, List[Bulletin], List[Mission]]] = []

        for _ in range(days):
            # Generate content for current day
            galnet = GalNetSynth(seed=rng.randrange(0, 2**31 - 1), risk_model=self.risk_model)
            bulletins = galnet.generate_bulletins(state, max_items=bulletins_per_day)

            ms = MissionSynth(seed=rng.randrange(0, 2**31 - 1), risk_model=self.risk_model)
            missions = ms.generate(state, bulletins, max_items=missions_per_day)

            out.append((state, bulletins, missions))

            # Step simulation forward
            state = self._step(state, rng=rng, dt=step_days)

        return out

    def _step(self, state: UniverseState, *, rng: random.Random, dt: float) -> UniverseState:
        # 1) expire events naturally (they are filtered via is_active(day))
        next_day = state.day + dt

        # 2) evolve systems fields in a bounded, mean-reverting way
        new_systems: List[System] = []
        for s in state.systems:
            rb = self.risk_model.system_risk(state, s.id)
            risk = rb.total

            # Pirate activity drifts upward with risk but is pulled down by security.
            drift = 0.15 * (risk - 0.45) - 0.10 * (s.security - 0.5)
            noise = (rng.random() - 0.5) * 0.10
            pirates = _clamp01(s.pirate_activity + (drift + noise) * dt)

            # Anomaly level drifts slowly, spikes with anomaly-related events.
            drift_a = 0.03 * (risk - 0.5)
            noise_a = (rng.random() - 0.5) * 0.04
            anomaly = _clamp01(s.anomaly_level + (drift_a + noise_a) * dt)

            new_systems.append(
                System(
                    id=s.id,
                    name=s.name,
                    pos_ly=s.pos_ly,
                    security=s.security,
                    controlling_faction_id=s.controlling_faction_id,
                    pirate_activity=pirates,
                    anomaly_level=anomaly,
                    stations=s.stations,
                )
            )

        # 3) spawn emergent events
        new_events: List[Event] = list(state.events)
        sys_by_id = {s.id: s for s in new_systems}

        def has_active(kind_sub: str, sys_id: int) -> bool:
            ks = kind_sub.lower()
            for e in new_events:
                if e.system_id == sys_id and ks in e.kind.lower() and e.is_active(day=next_day):
                    return True
            return False

        for s in new_systems:
            # Pirate spikes
            if s.pirate_activity > 0.75 and not has_active("pirate", s.id) and rng.random() < 0.35:
                new_events.append(
                    Event(
                        id=_spawn_event_id("ev_pirates", next_day, s.id, rng),
                        kind="pirate_spike",
                        system_id=s.id,
                        severity=_clamp01(0.65 + 0.35 * rng.random()),
                        start_day=next_day,
                        duration_days=1.5 + 1.5 * rng.random(),
                    )
                )

            # Radiation storms / flares
            if s.anomaly_level > 0.70 and not has_active("flare", s.id) and rng.random() < 0.22:
                new_events.append(
                    Event(
                        id=_spawn_event_id("ev_flare", next_day, s.id, rng),
                        kind="solar_flare",
                        system_id=s.id,
                        severity=_clamp01(0.55 + 0.40 * rng.random()),
                        start_day=next_day,
                        duration_days=0.8 + 0.8 * rng.random(),
                    )
                )

            # Discoveries are rarer but flavorful.
            if s.anomaly_level > 0.55 and not has_active("discovery", s.id) and rng.random() < 0.10:
                new_events.append(
                    Event(
                        id=_spawn_event_id("ev_disc", next_day, s.id, rng),
                        kind="anomaly_discovery",
                        system_id=s.id,
                        severity=_clamp01(0.40 + 0.40 * rng.random()),
                        start_day=next_day,
                        duration_days=2.0 + 3.0 * rng.random(),
                    )
                )

        # 4) produce next UniverseState
        next_state = UniverseState(
            day=next_day,
            factions=state.factions,
            systems=tuple(new_systems),
            events=tuple(new_events),
            max_jump_ly=state.max_jump_ly,
        )
        # Validate again (cheap)
        next_state._validate()
        return next_state

    # ---------- convenience API (used by other tooling) ----------

    def step_state(self, state: UniverseState, *, dt: float = 1.0, rng: Optional[random.Random] = None) -> UniverseState:
        """Advance a UniverseState by *dt* days.

        This is a thin wrapper around the internal :meth:`_step` used by
        :meth:`simulate`, exposed for other tools (e.g., interactive dev shells)
        that need incremental evolution.
        """

        if dt <= 0:
            return state

        if rng is None:
            rng = random.Random(self.seed if self.seed is not None else (hash(int(state.day * 1000.0)) & 0xFFFFFFFF))

        return self._step(state, rng=rng, dt=dt)

    def advance_state(
        self,
        state: UniverseState,
        *,
        days: float,
        step: float = 1.0,
        rng: Optional[random.Random] = None,
    ) -> UniverseState:
        """Advance *state* forward by *days* (can be fractional).

        For deltas larger than *step*, this method iterates in fixed increments
        so the emergent-event logic has multiple chances to trigger.
        """

        if days <= 0:
            return state
        if step <= 0:
            raise ValueError("step must be > 0")

        if rng is None:
            rng = random.Random(self.seed if self.seed is not None else (hash(int(state.day * 1000.0)) & 0xFFFFFFFF))

        remaining = float(days)
        st = state
        while remaining > 1e-9:
            dt = step if remaining > step else remaining
            st = self._step(st, rng=rng, dt=dt)
            remaining -= dt
        return st
