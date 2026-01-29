"""Captain's Terminal (text-mode prototype gameplay).

This app turns the offline content pipeline into an *interactive loop*:
- read UniverseState
- generate GalNet + mission board for the day
- let the player accept missions
- simulate travel with risk-driven encounters
- advance the UniverseState timeline and update faction reputation

This is aimed at:
- populating underdeveloped gameplay surfaces (missions, comms, news)
- giving a playable 'vertical slice' while the in-engine versions mature

Run:
    python -m tools.captains_terminal
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Dict, Iterable, List, Optional, Tuple
import argparse
import hashlib
import math

from ..stellar_content_pipeline.galnet import GalNetSynth, Bulletin
from ..stellar_content_pipeline.missions import MissionSynth, Mission
from ..stellar_content_pipeline.risk import RiskModel, RoutePlanner
from ..stellar_content_pipeline.schema import UniverseState, System
from ..stellar_content_pipeline.simulator import TimelineSimulator

from .encounters import maybe_encounter
from .savegame import (
    SaveGame,
    load_or_create_game,
    save_game_to_path,
    adjust_faction_rep,
    next_action_rng,
)
from .ui import header, rule, fmt_credits, prompt_choice, prompt_int, wrap_lines


DAYS_PER_JUMP = 0.22
DAYS_PER_LY = 0.015


def _stable_seed(base_seed: int, day: float, tag: str) -> int:
    """Derive a deterministic per-day seed from (base_seed, day, tag)."""

    msg = f"{base_seed}:{day:.3f}:{tag}".encode("utf-8")
    digest = hashlib.sha1(msg).digest()
    # 32-bit seed; avoid 0.
    seed = int.from_bytes(digest[:4], "big") & 0x7FFFFFFF
    return seed if seed != 0 else 1


def _sys_name(state: UniverseState, sid: int) -> str:
    s = state.system_by_id(sid)
    return s.name if s else f"System#{sid}"


def _format_day(day: float) -> str:
    # Keep it simple: day 0.00 is an epoch. Show both integer day and fraction.
    d0 = int(math.floor(day))
    frac = day - d0
    return f"Day {d0} (+{frac:.2f})"


def _print_wrapped(text: str) -> None:
    for ln in wrap_lines(text):
        print(ln)


def _show_status(save: SaveGame) -> None:
    st = save.universe
    pl = save.player

    print(header("CAPTAIN STATUS"))
    print(f"Time: {_format_day(st.day)}")
    print(f"Location: {_sys_name(st, pl.current_system_id)} (id={pl.current_system_id})")
    print(f"Credits: {fmt_credits(pl.credits)}")
    print(f"Hull: {pl.hull}%")
    print(f"Active missions: {len(pl.active_missions)}")

    # Show top 3 faction reps.
    reps = sorted(((f.player_rep, f.name) for f in st.factions), key=lambda x: x[0], reverse=True)
    if reps:
        print("Reputation:")
        for rep, name in reps[:3]:
            sign = "+" if rep >= 0 else ""
            print(f"  {name}: {sign}{rep:.1f}")

    # Active events in current system.
    active_here = [e for e in st.active_events() if e.system_id == pl.current_system_id]
    if active_here:
        print("Alerts:")
        for e in active_here[:3]:
            print(f"  - {e.kind} (sev {e.severity:.2f})")

    print(rule("-"))


def _generate_daily_content(save: SaveGame, *, max_bulletins: int = 10, max_missions: int = 12) -> Tuple[List[Bulletin], List[Mission]]:
    state = save.universe

    rm = RiskModel()

    g_seed = _stable_seed(save.rng_seed, state.day, "galnet")
    m_seed = _stable_seed(save.rng_seed, state.day, "missions")

    galnet = GalNetSynth(seed=g_seed, risk_model=rm)
    bulletins = galnet.generate_bulletins(state, max_items=max_bulletins)

    ms = MissionSynth(seed=m_seed, risk_model=rm)
    missions = ms.generate(state, bulletins, max_items=max_missions)

    return bulletins, missions


def _show_galnet(save: SaveGame, *, max_items: int = 10) -> None:
    bulletins, _ = _generate_daily_content(save, max_bulletins=max_items, max_missions=0)

    print(header("GALNET"))
    if not bulletins:
        print("No bulletins today.")
        print(rule("-"))
        return

    for i, b in enumerate(bulletins, start=1):
        print(f"[{i:02d}] {b.subject}")
        body = (b.body or "").strip()
        if body:
            lines = wrap_lines(body)
            for ln in lines[:3]:
                print(ln)
            if len(lines) > 3:
                print("...")
        print(f"System: {_sys_name(save.universe, b.system_id)}  |  Severity: {b.severity:.2f}")
        print(f"Tags: {', '.join(b.tags) if b.tags else '(none)'}")
        print(rule("-"))


def _mission_short(save: SaveGame, m: Mission) -> str:
    origin = _sys_name(save.universe, m.origin_system_id)
    dest = _sys_name(save.universe, m.dest_system_id)
    return f"{m.kind.upper():7s} {fmt_credits(m.reward_credits):>10s}  ETA {m.eta_days:>4.1f}d  Risk {m.risk_avg:.2f}  {origin} -> {dest}"


def _show_board(save: SaveGame, *, max_items: int = 12) -> None:
    st = save.universe
    pl = save.player

    bulletins, missions = _generate_daily_content(save, max_bulletins=10, max_missions=max_items)

    print(header("MISSION BOARD"))
    print(f"Board date: {_format_day(st.day)}")
    print(f"Docked near: {_sys_name(st, pl.current_system_id)}")
    print(rule("-"))

    if not missions:
        print("No mission offers at this time.")
        print(rule("-"))
        return

    # Filter out missions already accepted/completed/failed (by id).
    accepted_ids = {rec.get("mission", {}).get("id") for rec in pl.active_missions}
    accepted_ids.update(pl.completed_mission_ids)
    accepted_ids.update(pl.failed_mission_ids)

    visible: List[Mission] = [m for m in missions if m.id not in accepted_ids]
    if not visible:
        print("No *new* mission offers today (you've already interacted with them).")
        print(rule("-"))
        return

    for i, m in enumerate(visible, start=1):
        print(f"[{i:02d}] {_mission_short(save, m)}")
        print(f"     Title: {m.title}")
        # Keep one wrapped line of description.
        desc = (m.description or "").strip()
        if desc:
            line = wrap_lines(desc)[0]
            print(f"     {line}")
        print(f"     Deadline: {m.deadline_day:.2f}  Tags: {', '.join(m.tags)}")

    print(rule("-"))
    choice = prompt_int("Accept which mission number? (blank to cancel) > ", min_v=1, max_v=len(visible), allow_blank=True)
    if choice is None:
        return

    m = visible[int(choice) - 1]
    pl.active_missions.append({"mission": m.to_dict(), "accepted_day": float(st.day)})
    print(f"Accepted: {m.title}")


def _resolve_missions(save: SaveGame) -> None:
    """Resolve mission completion/failure based on day + location."""

    st = save.universe
    pl = save.player
    now = float(st.day)
    here = int(pl.current_system_id)

    new_active: List[Dict[str, Any]] = []

    completed: List[Dict[str, Any]] = []
    failed: List[Dict[str, Any]] = []

    for rec in pl.active_missions:
        m = rec.get("mission") or {}
        mid = str(m.get("id") or "")
        if not mid:
            continue

        deadline = float(m.get("deadline_day") or 0.0)
        dest_sys = int(m.get("dest_system_id") or 0)
        issuer = int(m.get("issuer_faction_id") or 0)
        reward = int(m.get("reward_credits") or 0)
        risk_avg = float(m.get("risk_avg") or 0.0)

        if now > deadline + 1e-9:
            failed.append(rec)
            pl.failed_mission_ids.append(mid)
            # Reputation hit.
            save.universe = adjust_faction_rep(save.universe, issuer, delta=-(2.0 + 3.5 * risk_avg))
            continue

        if dest_sys != 0 and here == dest_sys:
            completed.append(rec)
            pl.completed_mission_ids.append(mid)
            pl.credits += reward
            # Reputation gain.
            save.universe = adjust_faction_rep(save.universe, issuer, delta=(2.5 + 4.0 * risk_avg))
            continue

        new_active.append(rec)

    pl.active_missions = new_active
    pl.clamp()

    if completed or failed:
        print(header("MISSION UPDATE"))
        for rec in completed:
            m = rec.get("mission") or {}
            print(f"COMPLETED: {m.get('title') or m.get('id')}  (+{fmt_credits(int(m.get('reward_credits') or 0))})")
        for rec in failed:
            m = rec.get("mission") or {}
            print(f"FAILED:    {m.get('title') or m.get('id')}  (deadline missed)")
        print(rule("-"))


def _travel(save: SaveGame) -> None:
    st = save.universe
    pl = save.player

    systems = sorted(st.systems, key=lambda s: s.id)
    if not systems:
        print("Universe contains no systems.")
        return

    print(header("NAVIGATION"))
    print(f"Current: {_sys_name(st, pl.current_system_id)} (id={pl.current_system_id})")
    print("Destinations:")
    for s in systems[:40]:
        # Keep list from exploding in huge universes.
        sec = int(round(s.security * 100))
        print(f"  {s.id:>3d}: {s.name:<18s}  sec {sec:>3d}%  pirates {s.pirate_activity:.2f}  anom {s.anomaly_level:.2f}")
    if len(systems) > 40:
        print(f"  ... ({len(systems) - 40} more systems)")

    dest_id = prompt_int("Enter destination system id (blank to cancel) > ", min_v=0, max_v=max(s.id for s in systems), allow_blank=True)
    if dest_id is None:
        return

    dest_id = int(dest_id)
    if dest_id == pl.current_system_id:
        print("Already there.")
        return

    if st.system_by_id(dest_id) is None:
        print("Unknown destination.")
        return

    rm = RiskModel()
    rp = RoutePlanner()
    plan = rp.build_route_plan(st, pl.current_system_id, dest_id, risk_model=rm, days_per_jump=DAYS_PER_JUMP, days_per_ly=DAYS_PER_LY)
    if plan is None:
        print("No route found within jump range.")
        return

    print(rule("-"))
    print(f"Route: {' -> '.join(_sys_name(st, sid) for sid in plan.systems)}")
    print(f"Jumps: {plan.jumps}, Distance: {plan.total_distance_ly:.1f} ly")
    print(f"ETA:   {plan.eta_days:.2f} days, Risk avg {plan.risk_avg:.2f} (max {plan.risk_max:.2f})")
    print(rule("-"))

    confirm = prompt_choice("Commit to travel? [y/N] > ", allow_blank=True).lower()
    if confirm not in ("y", "yes"):
        return

    rng = next_action_rng(save)
    sim = TimelineSimulator(seed=None)

    # Simulate segment by segment.
    base_days = 0.0
    extra_days = 0.0
    credit_delta = 0
    hull_delta = 0

    # Keep a small log.
    travel_log: List[str] = []

    for a, b in zip(plan.systems, plan.systems[1:]):
        sa = st.system_by_id(a)
        sb = st.system_by_id(b)
        if sa is None or sb is None:
            continue
        dist = sa.pos_ly.dist(sb.pos_ly)
        seg_days = DAYS_PER_JUMP + dist * DAYS_PER_LY
        base_days += seg_days

        # Encounter check based on arrival system's risk.
        rb = rm.system_risk(st, sb.id)
        enc = maybe_encounter(rng, risk=rb)
        extra_days += enc.day_delta
        credit_delta += int(enc.credits_delta)
        hull_delta += int(enc.hull_delta)
        travel_log.extend(enc.log_lines)

    total_days = base_days + extra_days

    # Apply deltas.
    pl.credits += credit_delta
    pl.hull += hull_delta
    pl.clamp()

    # Advance universe timeline.
    save.universe = sim.advance_state(save.universe, days=total_days, step=1.0, rng=rng)

    # Arrive.
    pl.current_system_id = dest_id

    print(header("ARRIVAL"))
    print(f"Arrived in {_sys_name(save.universe, dest_id)}")
    print(f"Travel time: {total_days:.2f} days (base {base_days:.2f} + delays {extra_days:.2f})")
    if credit_delta != 0:
        print(f"Credits change: {fmt_credits(credit_delta)}")
    if hull_delta != 0:
        print(f"Hull change: {hull_delta}%")
    print(f"Now: {_format_day(save.universe.day)}")
    print(rule("-"))

    if travel_log:
        print("In-flight log:")
        for ln in travel_log[:10]:
            _print_wrapped(ln)
        if len(travel_log) > 10:
            print(f"... ({len(travel_log) - 10} more lines)")
        print(rule("-"))

    _resolve_missions(save)


def _repair(save: SaveGame) -> None:
    st = save.universe
    pl = save.player

    if pl.hull >= 100:
        print("Hull already at 100%.")
        return

    sys = st.system_by_id(pl.current_system_id)
    sec = sys.security if sys else 0.5

    # Safer systems -> better repair economy.
    cost_per = int(round(22 - 12 * sec))
    cost_per = max(8, min(25, cost_per))

    missing = 100 - pl.hull
    total = missing * cost_per

    print(header("REPAIR BAY"))
    print(f"Hull: {pl.hull}% -> 100% (missing {missing}%)")
    print(f"Rate: {fmt_credits(cost_per)} per %")
    print(f"Full repair cost: {fmt_credits(total)}")
    print(f"Available credits: {fmt_credits(pl.credits)}")

    if pl.credits <= 0:
        print("You can't afford repairs.")
        return

    confirm = prompt_choice("Proceed? [y/N] > ", allow_blank=True).lower()
    if confirm not in ("y", "yes"):
        return

    affordable = min(missing, pl.credits // cost_per)
    if affordable <= 0:
        print("Insufficient credits.")
        return

    pl.hull += affordable
    pl.credits -= affordable * cost_per
    pl.clamp()
    print(f"Repaired {affordable}% hull.")


def _wait(save: SaveGame) -> None:
    st = save.universe

    print(header("WAIT"))
    days = prompt_int("Wait how many days? (1..30, blank cancels) > ", min_v=1, max_v=30, allow_blank=True)
    if days is None:
        return

    rng = next_action_rng(save)
    sim = TimelineSimulator(seed=None)
    save.universe = sim.advance_state(save.universe, days=float(days), step=1.0, rng=rng)

    print(f"Advanced time by {days} days. Now {_format_day(save.universe.day)}")
    _resolve_missions(save)


def _help() -> None:
    print(header("COMMANDS"))
    print("status   - show ship + world status")
    print("galnet   - view today's bulletins")
    print("board    - view mission board / accept")
    print("missions - list active missions")
    print("travel   - plot route and travel")
    print("repair   - repair hull (costs credits)")
    print("wait     - advance time")
    print("save     - save to disk")
    print("quit     - save and exit")
    print(rule("-"))


def _show_missions(save: SaveGame) -> None:
    pl = save.player
    st = save.universe

    print(header("ACTIVE MISSIONS"))
    if not pl.active_missions:
        print("No active missions.")
        print(rule("-"))
        return

    now = float(st.day)
    for i, rec in enumerate(pl.active_missions, start=1):
        m = rec.get("mission") or {}
        title = str(m.get("title") or m.get("id") or "(unknown)")
        dest = int(m.get("dest_system_id") or 0)
        deadline = float(m.get("deadline_day") or 0.0)
        reward = int(m.get("reward_credits") or 0)
        remaining = deadline - now

        print(f"[{i:02d}] {title}")
        print(f"     Dest: {_sys_name(st, dest)}  Deadline: {deadline:.2f}  (in {remaining:.2f} days)")
        print(f"     Reward: {fmt_credits(reward)}")

    print(rule("-"))


def _parse_args(argv: Optional[List[str]]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Captain's Terminal - interactive dev gameplay loop")
    p.add_argument("--universe", default="tools/examples/universe_state.json", help="Universe template (used if save doesn't exist)")
    p.add_argument("--save", default="out/captains_terminal/savegame.json", help="Savegame path")
    p.add_argument("--seed", type=int, default=None, help="Base RNG seed (new saves only)")
    p.add_argument("--player", default="Commander", help="Player name (new saves only)")
    p.add_argument("--start-system-id", type=int, default=None, help="Starting system id (new saves only)")
    p.add_argument("--reset-universe", action="store_true", help="Ignore existing save; start a fresh game")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = _parse_args(argv)

    save = load_or_create_game(
        save_path=args.save,
        universe_path=args.universe,
        rng_seed=args.seed,
        player_name=args.player,
        starting_system_id=args.start_system_id,
        reset_universe=args.reset_universe,
    )

    print(header("WELCOME TO CAPTAIN'S TERMINAL"))
    print("A text-mode slice of the 3D space sim: news, missions, travel, reputation.")
    print(f"Save: {args.save}")
    print(rule("-"))
    _help()

    try:
        while True:
            cmd = prompt_choice("cmd> ", allow_blank=True).lower()
            if cmd in ("", "help", "h", "?"):
                _help()
                continue
            if cmd in ("status", "s"):
                _show_status(save)
                continue
            if cmd in ("galnet", "g", "news"):
                _show_galnet(save)
                continue
            if cmd in ("board", "b", "missions-board"):
                _show_board(save)
                continue
            if cmd in ("missions", "m"):
                _show_missions(save)
                continue
            if cmd in ("travel", "t"):
                _travel(save)
                continue
            if cmd in ("repair", "r"):
                _repair(save)
                continue
            if cmd in ("wait", "w"):
                _wait(save)
                continue
            if cmd in ("save", "sv"):
                save_game_to_path(save, args.save)
                print("Saved.")
                continue
            if cmd in ("quit", "q", "exit"):
                save_game_to_path(save, args.save)
                print("Saved. Goodbye.")
                break

            print("Unknown command. Type 'help'.")

    except KeyboardInterrupt:
        print("\nInterrupted.")
        try:
            save_game_to_path(save, args.save)
            print("Saved.")
        except Exception:
            pass

    return 0
