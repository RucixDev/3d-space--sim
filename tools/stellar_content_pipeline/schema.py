"""StellarForge content pipeline: schema + validation utilities.

This package is intentionally dependency-free (stdlib only) so it can run anywhere.

The goal: generate *game-ready* narrative content (GalNet bulletins + mission offers)
from a minimal 'universe state' JSON snapshot.

Design philosophy:
- Deterministic output when seed + input are identical.
- Strict, explicit validation with helpful error messages.
- No reliance on the in-engine types; output is JSON that can be ingested by the game.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional, Tuple, Iterable
import json
import math


class SchemaError(ValueError):
    """Raised when the input JSON violates the expected schema."""

    def __init__(self, message: str, *, path: str = "$") -> None:
        super().__init__(f"{message} (at {path})")
        self.path = path


def _type_name(x: Any) -> str:
    return type(x).__name__


def _req(obj: Dict[str, Any], key: str, *, path: str) -> Any:
    if key not in obj:
        raise SchemaError(f"Missing required field '{key}'", path=f"{path}.{key}")
    return obj[key]


def _opt(obj: Dict[str, Any], key: str) -> Any:
    return obj.get(key)


def _as_int(v: Any, *, path: str) -> int:
    if isinstance(v, bool) or not isinstance(v, int):
        raise SchemaError(f"Expected int, got {_type_name(v)}", path=path)
    return v


def _as_float(v: Any, *, path: str) -> float:
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise SchemaError(f"Expected number, got {_type_name(v)}", path=path)
    return float(v)


def _as_str(v: Any, *, path: str) -> str:
    if not isinstance(v, str):
        raise SchemaError(f"Expected string, got {_type_name(v)}", path=path)
    return v


def _as_bool(v: Any, *, path: str) -> bool:
    if not isinstance(v, bool):
        raise SchemaError(f"Expected bool, got {_type_name(v)}", path=path)
    return v


def _as_list(v: Any, *, path: str) -> List[Any]:
    if not isinstance(v, list):
        raise SchemaError(f"Expected list, got {_type_name(v)}", path=path)
    return v


def _as_dict(v: Any, *, path: str) -> Dict[str, Any]:
    if not isinstance(v, dict):
        raise SchemaError(f"Expected object, got {_type_name(v)}", path=path)
    return v


@dataclass(frozen=True)
class Vec3:
    x: float
    y: float
    z: float

    def dist(self, other: "Vec3") -> float:
        dx = other.x - self.x
        dy = other.y - self.y
        dz = other.z - self.z
        return math.sqrt(dx * dx + dy * dy + dz * dz)

    @staticmethod
    def from_json(obj: Any, *, path: str) -> "Vec3":
        o = _as_dict(obj, path=path)
        x = _as_float(_req(o, "x", path=path), path=f"{path}.x")
        y = _as_float(_req(o, "y", path=path), path=f"{path}.y")
        z = _as_float(_req(o, "z", path=path), path=f"{path}.z")
        return Vec3(x=x, y=y, z=z)


@dataclass(frozen=True)
class Faction:
    id: int
    name: str

    # Optional 'baseline' reputation band for the player (-100..100). If omitted, 0.
    player_rep: float = 0.0

    @staticmethod
    def from_json(obj: Any, *, path: str) -> "Faction":
        o = _as_dict(obj, path=path)
        fid = _as_int(_req(o, "id", path=path), path=f"{path}.id")
        name = _as_str(_req(o, "name", path=path), path=f"{path}.name")
        player_rep = _as_float(_opt(o, "player_rep") or 0.0, path=f"{path}.player_rep")
        return Faction(id=fid, name=name, player_rep=float(player_rep))


@dataclass(frozen=True)
class Station:
    id: int
    name: str
    kind: str = "outpost"  # e.g., outpost, city, refinery, shipyard
    economy: str = "mixed"  # e.g., industrial, agricultural, military, research

    @staticmethod
    def from_json(obj: Any, *, path: str) -> "Station":
        o = _as_dict(obj, path=path)
        sid = _as_int(_req(o, "id", path=path), path=f"{path}.id")
        name = _as_str(_req(o, "name", path=path), path=f"{path}.name")
        kind = _as_str(_opt(o, "kind") or "outpost", path=f"{path}.kind")
        economy = _as_str(_opt(o, "economy") or "mixed", path=f"{path}.economy")
        return Station(id=sid, name=name, kind=kind, economy=economy)


@dataclass(frozen=True)
class System:
    id: int
    name: str
    pos_ly: Vec3

    # 0..1 where 1 is maximum security / order.
    security: float = 0.5

    # controlling faction id (0 means unclaimed / independent)
    controlling_faction_id: int = 0

    # extra 'pirate pressure' 0..1 (adds risk)
    pirate_activity: float = 0.0

    # environmental hazards 0..1 (adds risk)
    anomaly_level: float = 0.0

    stations: Tuple[Station, ...] = field(default_factory=tuple)

    @staticmethod
    def from_json(obj: Any, *, path: str) -> "System":
        o = _as_dict(obj, path=path)
        sid = _as_int(_req(o, "id", path=path), path=f"{path}.id")
        name = _as_str(_req(o, "name", path=path), path=f"{path}.name")
        pos = Vec3.from_json(_req(o, "pos_ly", path=path), path=f"{path}.pos_ly")
        security = _as_float(_opt(o, "security") if _opt(o, "security") is not None else 0.5, path=f"{path}.security")
        if not (0.0 <= security <= 1.0):
            raise SchemaError("security must be within [0,1]", path=f"{path}.security")
        controlling = _as_int(_opt(o, "controlling_faction_id") or 0, path=f"{path}.controlling_faction_id")
        pirate = _as_float(_opt(o, "pirate_activity") if _opt(o, "pirate_activity") is not None else 0.0, path=f"{path}.pirate_activity")
        anomaly = _as_float(_opt(o, "anomaly_level") if _opt(o, "anomaly_level") is not None else 0.0, path=f"{path}.anomaly_level")
        pirate = max(0.0, min(1.0, pirate))
        anomaly = max(0.0, min(1.0, anomaly))

        st_list = []
        for i, st in enumerate(_as_list(_opt(o, "stations") or [], path=f"{path}.stations")):
            st_list.append(Station.from_json(st, path=f"{path}.stations[{i}]"))

        return System(
            id=sid,
            name=name,
            pos_ly=pos,
            security=security,
            controlling_faction_id=controlling,
            pirate_activity=pirate,
            anomaly_level=anomaly,
            stations=tuple(st_list),
        )


@dataclass(frozen=True)
class Event:
    """A lightweight, game-agnostic world event.

    `kind` examples:
      - pirate_spike
      - solar_flare
      - faction_skirmish
      - famine
      - outbreak
      - discovery
      - shipyard_strike
    """

    id: str
    kind: str
    system_id: int
    severity: float  # 0..1
    headline: str = ""
    details: str = ""
    start_day: float = 0.0
    duration_days: float = 1.0

    @staticmethod
    def from_json(obj: Any, *, path: str) -> "Event":
        o = _as_dict(obj, path=path)
        eid = _as_str(_req(o, "id", path=path), path=f"{path}.id")
        kind = _as_str(_req(o, "kind", path=path), path=f"{path}.kind")
        system_id = _as_int(_req(o, "system_id", path=path), path=f"{path}.system_id")
        severity = _as_float(_opt(o, "severity") if _opt(o, "severity") is not None else 0.5, path=f"{path}.severity")
        severity = max(0.0, min(1.0, severity))
        headline = _as_str(_opt(o, "headline") or "", path=f"{path}.headline")
        details = _as_str(_opt(o, "details") or "", path=f"{path}.details")
        start_day = _as_float(_opt(o, "start_day") if _opt(o, "start_day") is not None else 0.0, path=f"{path}.start_day")
        duration_days = _as_float(_opt(o, "duration_days") if _opt(o, "duration_days") is not None else 1.0, path=f"{path}.duration_days")
        duration_days = max(0.01, duration_days)
        return Event(
            id=eid,
            kind=kind,
            system_id=system_id,
            severity=severity,
            headline=headline,
            details=details,
            start_day=start_day,
            duration_days=duration_days,
        )

    def is_active(self, *, day: float) -> bool:
        return self.start_day <= day < (self.start_day + self.duration_days)


@dataclass(frozen=True)
class UniverseState:
    """Snapshot input for content generation."""

    day: float
    factions: Tuple[Faction, ...]
    systems: Tuple[System, ...]
    events: Tuple[Event, ...] = field(default_factory=tuple)

    # Optional max jump range (ly). Used by route planning.
    max_jump_ly: float = 8.0

    @staticmethod
    def from_json(obj: Any, *, path: str = "$") -> "UniverseState":
        o = _as_dict(obj, path=path)
        day = _as_float(_req(o, "day", path=path), path=f"{path}.day")

        factions_list = []
        for i, f in enumerate(_as_list(_req(o, "factions", path=path), path=f"{path}.factions")):
            factions_list.append(Faction.from_json(f, path=f"{path}.factions[{i}]"))

        systems_list = []
        for i, s in enumerate(_as_list(_req(o, "systems", path=path), path=f"{path}.systems")):
            systems_list.append(System.from_json(s, path=f"{path}.systems[{i}]"))

        events_list = []
        for i, e in enumerate(_as_list(_opt(o, "events") or [], path=f"{path}.events")):
            events_list.append(Event.from_json(e, path=f"{path}.events[{i}]"))

        max_jump = _as_float(_opt(o, "max_jump_ly") if _opt(o, "max_jump_ly") is not None else 8.0, path=f"{path}.max_jump_ly")
        max_jump = max(0.25, max_jump)

        st = UniverseState(day=day, factions=tuple(factions_list), systems=tuple(systems_list), events=tuple(events_list), max_jump_ly=max_jump)
        st._validate()
        return st

    def _validate(self) -> None:
        # IDs unique
        f_ids = [f.id for f in self.factions]
        if len(set(f_ids)) != len(f_ids):
            raise SchemaError("Faction ids must be unique", path="$.factions")
        s_ids = [s.id for s in self.systems]
        if len(set(s_ids)) != len(s_ids):
            raise SchemaError("System ids must be unique", path="$.systems")

        # system controlling faction ids exist (or 0)
        fset = set(f_ids)
        for s in self.systems:
            if s.controlling_faction_id != 0 and s.controlling_faction_id not in fset:
                raise SchemaError(
                    f"System '{s.name}' references unknown controlling_faction_id={s.controlling_faction_id}",
                    path=f"$.systems[id={s.id}].controlling_faction_id",
                )

        # event system ids exist
        sset = set(s_ids)
        for e in self.events:
            if e.system_id not in sset:
                raise SchemaError(
                    f"Event '{e.id}' references unknown system_id={e.system_id}",
                    path=f"$.events[id={e.id}].system_id",
                )

    def faction_by_id(self, fid: int) -> Optional[Faction]:
        for f in self.factions:
            if f.id == fid:
                return f
        return None

    def system_by_id(self, sid: int) -> Optional[System]:
        for s in self.systems:
            if s.id == sid:
                return s
        return None

    def active_events(self, *, day: Optional[float] = None) -> Tuple[Event, ...]:
        d = self.day if day is None else day
        return tuple(e for e in self.events if e.is_active(day=d))


def load_universe_state_json(path: str) -> UniverseState:
    """Load and validate an input universe state JSON file."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    return UniverseState.from_json(data)


def universe_state_to_dict(state: UniverseState) -> Dict[str, Any]:
    """Convert a UniverseState into a JSON-serializable dict.

    The resulting dict conforms to the same schema expected by
    :func:`UniverseState.from_json`.

    Notes:
    - We intentionally use :func:`dataclasses.asdict` to keep this dependency-free.
    - Tuples are converted to lists (JSON has no tuple type).
    """

    # dataclasses.asdict() keeps tuples as tuples, which is inconvenient for JSON.
    # Convert any tuples we see into lists.
    def _jsonify(x: Any) -> Any:
        if isinstance(x, tuple):
            return [_jsonify(v) for v in x]
        if isinstance(x, list):
            return [_jsonify(v) for v in x]
        if isinstance(x, dict):
            return {k: _jsonify(v) for k, v in x.items()}
        return x

    return _jsonify(asdict(state))


def _atomic_write_text(path: str, text: str) -> None:
    """Atomically replace a text file (best-effort cross-platform)."""

    import os

    tmp_path = f"{path}.tmp"
    with open(tmp_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    os.replace(tmp_path, path)


def save_universe_state_json(state: UniverseState, path: str, *, indent: int = 2) -> None:
    """Write a UniverseState back to disk as JSON.

    This is primarily used by dev tooling (e.g., timeline prototyping, CLI apps).
    """

    data = universe_state_to_dict(state)
    text = json.dumps(data, ensure_ascii=False, indent=indent) + "\n"
    _atomic_write_text(path, text)
