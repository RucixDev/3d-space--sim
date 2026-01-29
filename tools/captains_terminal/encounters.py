"""Risk-driven travel encounters.

This is *not* meant to be a full combat system.
It is a lightweight way to turn route-risk numbers into gameplay texture:
- delays
- minor hull damage
- credit loss / salvage gain

The goal is to make "travel" feel like a decision, not a loading screen.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List
import random

from ..stellar_content_pipeline.risk import RiskBreakdown


@dataclass(frozen=True)
class EncounterResult:
    day_delta: float = 0.0
    credits_delta: int = 0
    hull_delta: int = 0
    log_lines: List[str] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        # dataclasses can't set default mutable values safely; fix at runtime.
        if self.log_lines is None:
            object.__setattr__(self, "log_lines", [])


def _chance_from_risk(r: float) -> float:
    # base 2% even in safe space; scales up aggressively.
    return max(0.02, min(0.55, 0.02 + (r ** 1.6) * 0.55))


def maybe_encounter(rng: random.Random, *, risk: RiskBreakdown) -> EncounterResult:
    """Sample a travel encounter given risk breakdown."""

    if rng.random() > _chance_from_risk(risk.total):
        return EncounterResult()

    # Choose an encounter "theme".
    # Weighted by which terms dominate.
    w_pirates = 0.15 + 1.25 * risk.pirates_term + 0.60 * risk.event_term
    w_anom = 0.12 + 1.10 * risk.anomaly_term + 0.35 * risk.event_term
    w_security = 0.10 + 0.90 * risk.security_term

    tot = w_pirates + w_anom + w_security
    x = rng.random() * tot

    if x < w_pirates:
        return _pirate_encounter(rng, risk)
    if x < w_pirates + w_anom:
        return _anomaly_encounter(rng, risk)
    return _security_encounter(rng, risk)


def _pirate_encounter(rng: random.Random, risk: RiskBreakdown) -> EncounterResult:
    sev = (0.35 + 0.75 * rng.random()) * max(0.15, risk.total)
    sev = min(1.0, sev)

    # Outcome: either you pay them off, or you take damage escaping.
    if rng.random() < 0.55:
        loss = int(80 + 520 * sev + 220 * rng.random())
        delay = 0.05 + 0.22 * sev
        return EncounterResult(
            day_delta=delay,
            credits_delta=-loss,
            hull_delta=0,
            log_lines=[
                "Pirate shakedown in the void.",
                f"You transferred {loss} cr to avoid escalation. (+{delay:.2f} days delay)",
            ],
        )

    dmg = int(2 + 18 * sev + 6 * rng.random())
    delay = 0.06 + 0.28 * sev
    salvage = int(25 + 140 * rng.random()) if rng.random() < 0.15 else 0
    lines = [
        "Ambush: raiders on intercept.",
        f"You burned hard and took {dmg}% hull damage escaping. (+{delay:.2f} days delay)",
    ]
    if salvage > 0:
        lines.append(f"Recovered drifting salvage: +{salvage} cr")

    return EncounterResult(day_delta=delay, credits_delta=salvage, hull_delta=-dmg, log_lines=lines)


def _anomaly_encounter(rng: random.Random, risk: RiskBreakdown) -> EncounterResult:
    sev = (0.25 + 0.85 * rng.random()) * max(0.10, risk.total)
    sev = min(1.0, sev)

    # Mostly time loss; sometimes minor damage.
    delay = 0.08 + 0.42 * sev

    if rng.random() < 0.25:
        dmg = int(1 + 10 * sev + 3 * rng.random())
        return EncounterResult(
            day_delta=delay,
            credits_delta=0,
            hull_delta=-dmg,
            log_lines=[
                "Spatial anomaly: gravimetric shear event.",
                f"Course correction cost time and scraped the hull (-{dmg}% hull, +{delay:.2f} days).",
            ],
        )

    return EncounterResult(
        day_delta=delay,
        credits_delta=0,
        hull_delta=0,
        log_lines=[
            "Sensor bloom: unexpected radiation pocket.",
            f"You throttled down and recalibrated (+{delay:.2f} days).",
        ],
    )


def _security_encounter(rng: random.Random, risk: RiskBreakdown) -> EncounterResult:
    sev = (0.20 + 0.90 * rng.random()) * max(0.10, risk.total)
    sev = min(1.0, sev)

    # Inspections are delays, with a small chance of a fine.
    delay = 0.04 + 0.30 * sev
    if rng.random() < 0.18:
        fine = int(40 + 380 * sev + 80 * rng.random())
        return EncounterResult(
            day_delta=delay,
            credits_delta=-fine,
            hull_delta=0,
            log_lines=[
                "System patrol: random inspection.",
                f"Minor paperwork issue -> fine {fine} cr (+{delay:.2f} days).",
            ],
        )

    return EncounterResult(
        day_delta=delay,
        credits_delta=0,
        hull_delta=0,
        log_lines=[
            "System patrol: transponder verification.",
            f"All clear, but you lost time in holding pattern (+{delay:.2f} days).",
        ],
    )
