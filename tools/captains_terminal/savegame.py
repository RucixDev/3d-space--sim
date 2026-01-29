"""Savegame + state helpers for Captain's Terminal.

The terminal is intentionally *stdlib only* and stores a single JSON file:
- Player state (credits, hull, location, accepted missions)
- UniverseState snapshot (day, systems, events, faction rep)

Keeping the UniverseState inside the save makes the tool self-contained:
no external database is required, and the file can be attached to bug reports.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Dict, List, Optional
import json
import os
import random

from ..stellar_content_pipeline.schema import UniverseState, load_universe_state_json, universe_state_to_dict


SAVE_VERSION = 1


@dataclass
class PlayerState:
    name: str = "Commander"
    credits: int = 5000
    hull: int = 100
    current_system_id: int = 0

    # Each record:
    # {
    #   "mission": <mission dict>,
    #   "accepted_day": float,
    # }
    active_missions: List[Dict[str, Any]] = field(default_factory=list)

    completed_mission_ids: List[str] = field(default_factory=list)
    failed_mission_ids: List[str] = field(default_factory=list)

    def clamp(self) -> None:
        self.credits = int(max(0, self.credits))
        self.hull = int(max(0, min(100, self.hull)))


@dataclass
class SaveGame:
    version: int
    rng_seed: int
    universe: UniverseState
    player: PlayerState

    # A monotonic counter used to derive deterministic action seeds.
    # This ensures that random outcomes don't reset each time you launch the terminal.
    rng_step: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "version": int(self.version),
            "rng_seed": int(self.rng_seed),
            "rng_step": int(self.rng_step),
            "universe_state": universe_state_to_dict(self.universe),
            "player": {
                "name": self.player.name,
                "credits": int(self.player.credits),
                "hull": int(self.player.hull),
                "current_system_id": int(self.player.current_system_id),
                "active_missions": list(self.player.active_missions),
                "completed_mission_ids": list(self.player.completed_mission_ids),
                "failed_mission_ids": list(self.player.failed_mission_ids),
            },
        }


def _atomic_write_text(path: str, text: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = f"{path}.tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    os.replace(tmp, path)


def save_game_to_path(save: SaveGame, path: str) -> None:
    """Write the savegame JSON to disk (atomic replace)."""

    text = json.dumps(save.to_dict(), ensure_ascii=False, indent=2) + "\n"
    _atomic_write_text(path, text)


def load_game_from_path(path: str) -> SaveGame:
    """Load a savegame JSON from disk (strict enough for dev tooling)."""

    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    version = int(data.get("version", 0))
    if version != SAVE_VERSION:
        raise ValueError(f"Unsupported savegame version={version}; expected {SAVE_VERSION}")

    rng_seed = int(data.get("rng_seed", 0))
    rng_step = int(data.get("rng_step") or 0)

    u_obj = data.get("universe_state")
    if not isinstance(u_obj, dict):
        raise ValueError("savegame missing universe_state object")
    universe = UniverseState.from_json(u_obj)

    p_obj = data.get("player")
    if not isinstance(p_obj, dict):
        raise ValueError("savegame missing player object")

    player = PlayerState(
        name=str(p_obj.get("name") or "Commander"),
        credits=int(p_obj.get("credits") or 0),
        hull=int(p_obj.get("hull") if p_obj.get("hull") is not None else 100),
        current_system_id=int(p_obj.get("current_system_id") or 0),
        active_missions=list(p_obj.get("active_missions") or []),
        completed_mission_ids=list(p_obj.get("completed_mission_ids") or []),
        failed_mission_ids=list(p_obj.get("failed_mission_ids") or []),
    )
    player.clamp()

    return SaveGame(version=version, rng_seed=rng_seed, universe=universe, player=player, rng_step=rng_step)


def create_new_game(
    *,
    universe_path: str,
    rng_seed: Optional[int] = None,
    player_name: str = "Commander",
    starting_system_id: Optional[int] = None,
) -> SaveGame:
    """Create a new SaveGame using a universe_state.json template."""

    universe = load_universe_state_json(universe_path)

    # Seed: accept user-provided, otherwise derive from OS randomness.
    seed = int(rng_seed) if rng_seed is not None else random.SystemRandom().randrange(1, 2**31 - 1)
    rng = random.Random(seed)

    # Pick a starting system with at least one station if possible.
    sys_id: int
    if starting_system_id is not None:
        sys_id = int(starting_system_id)
        if universe.system_by_id(sys_id) is None:
            raise ValueError(f"Unknown starting_system_id={sys_id}")
    else:
        candidates = [s.id for s in universe.systems if len(s.stations) > 0]
        if not candidates:
            candidates = [s.id for s in universe.systems]
        if not candidates:
            raise ValueError("Universe has no systems")
        sys_id = candidates[rng.randrange(0, len(candidates))]

    player = PlayerState(name=player_name, credits=5000, hull=100, current_system_id=sys_id)

    return SaveGame(version=SAVE_VERSION, rng_seed=seed, universe=universe, player=player, rng_step=0)


def load_or_create_game(
    *,
    save_path: str,
    universe_path: str,
    rng_seed: Optional[int] = None,
    player_name: str = "Commander",
    starting_system_id: Optional[int] = None,
    reset_universe: bool = False,
) -> SaveGame:
    """Load a savegame if present, otherwise create a new one."""

    if os.path.exists(save_path) and not reset_universe:
        return load_game_from_path(save_path)

    return create_new_game(
        universe_path=universe_path,
        rng_seed=rng_seed,
        player_name=player_name,
        starting_system_id=starting_system_id,
    )


def clamp_rep(rep: float) -> float:
    return max(-100.0, min(100.0, float(rep)))


def adjust_faction_rep(state: UniverseState, faction_id: int, delta: float) -> UniverseState:
    """Return a new UniverseState with a single faction's rep adjusted."""

    if faction_id == 0:
        return state

    new_factions = []
    for f in state.factions:
        if f.id == faction_id:
            new_factions.append(replace(f, player_rep=clamp_rep(f.player_rep + float(delta))))
        else:
            new_factions.append(f)

    return UniverseState(
        day=state.day,
        factions=tuple(new_factions),
        systems=state.systems,
        events=state.events,
        max_jump_ly=state.max_jump_ly,
    )


def next_action_rng(save: SaveGame) -> random.Random:
    """Return a deterministic RNG for the next *player action*.

    We don't persist Python's internal RNG state directly, since it's not JSON.
    Instead we keep a monotonic counter and derive new seeds from it.
    """

    # Xorshift-ish mixing (enough for gamey randomness, deterministic across runs).
    seed = (save.rng_seed ^ (save.rng_step * 0x9E3779B1)) & 0xFFFFFFFF
    save.rng_step += 1
    return random.Random(seed)
