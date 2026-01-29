# StellarForge Content Pipeline (Offline)

This directory contains an **offline, dependency-free** (stdlib only) content generator for:
- **GalNet bulletins + daily digest**
- **Mission offers** (courier, escort, bounty, rescue, relief, scan)
- **Navigation risk table** (CSV)
- Optional **multi-day timeline simulation** (emergent arcs)

It is meant to help push the game forward while in-engine versions of these systems are still evolving.

## Quick start

From the repo root:

```bash
python -m tools.stellar_content_pipeline \
  --input tools/examples/universe_state.json \
  --output out/content \
  --seed 123
```

Outputs will be written to `out/content/`:
- `galnet_bulletins.json`
- `galnet_digest.txt`
- `missions.json`
- `nav_risk.csv`

## Timeline mode (recommended for UI testing)

Generate a 7-day rolling feed:

```bash
python -m tools.stellar_content_pipeline \
  --input tools/examples/universe_state.json \
  --output out/content \
  --seed 123 \
  --simulate-days 7
```

This writes per-day packs into:

- `out/content/timeline/day_0000/`
- `out/content/timeline/day_0001/`
- ...
- `out/content/timeline/timeline_index.json`

## Input format

See `tools/examples/universe_state.json`.

At minimum, the JSON must contain:
- `day` (number)
- `factions` (list of `{id, name, player_rep?}`)
- `systems` (list of `{id, name, pos_ly, security?, controlling_faction_id?, pirate_activity?, anomaly_level?, stations?}`)
- Optional `events` (list of `{id, kind, system_id, severity?, headline?, details?, start_day?, duration_days?}`)
- Optional `max_jump_ly` (number)

The generator is **strict** and will error with clear messages if input is malformed.

## Design notes

- Deterministic outputs when `--seed` is supplied.
- Route planning uses **A\*** over a neighbor graph built with a **3D spatial hash grid** (fast enough for hundreds of systems).
- Rewards scale with **distance + route risk**, with different mission kinds adding multipliers.
- Timeline simulator evolves a few scalar fields and spawns lightweight events for believable arcs.

## Why this exists

The project already has UI surfaces for comms, GalNet, and mission control. This tool provides a
ready-to-ingest content stream so those screens can be populated even while the in-engine simulation
is still under construction.

## Captain's Terminal (interactive prototype gameplay)

If you want a *playable loop* (news → missions → travel → rewards) while the 3D UI is still evolving:

```bash
python -m tools.captains_terminal \
  --universe tools/examples/universe_state.json \
  --save out/captains_terminal/savegame.json
```

Inside the terminal:
- `galnet` to read today's bulletins
- `board` to accept missions
- `travel` to move between systems (with risk-driven encounters)
- `missions` to track deadlines and payouts

This tool writes a single JSON save file containing both the **player state** and a snapshot of the **UniverseState**.
