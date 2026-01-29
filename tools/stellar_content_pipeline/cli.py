"""CLI entrypoint for the Stellar content pipeline.

Usage (from repo root):
  python -m tools.stellar_content_pipeline --input tools/examples/universe_state.json --output out/content --seed 123

Outputs:
  - galnet_bulletins.json
  - galnet_digest.txt
  - missions.json
  - nav_risk.csv

Optional timeline mode:
  python -m tools.stellar_content_pipeline --input ... --output out/content --seed 123 --simulate-days 7

This writes per-day packs into:
  out/content/timeline/day_0000/
  out/content/timeline/day_0001/
  ...
and a timeline_index.json.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from typing import Optional

from .schema import load_universe_state_json, UniverseState
from .risk import RiskModel
from .galnet import GalNetSynth
from .missions import MissionSynth
from .simulator import TimelineSimulator


def _write_json(path: str, data: object) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def _write_text(path: str, text: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def _write_nav_risk_csv(path: str, state: UniverseState, risk_model: RiskModel) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["system_id", "system_name", "risk_total", "risk_security", "risk_pirates", "risk_anomaly", "risk_rep", "risk_events"])
        for s in sorted(state.systems, key=lambda x: x.id):
            rb = risk_model.system_risk(state, s.id)
            w.writerow([s.id, s.name, f"{rb.total:.4f}", f"{rb.security_term:.4f}", f"{rb.pirates_term:.4f}", f"{rb.anomaly_term:.4f}", f"{rb.rep_term:.4f}", f"{rb.event_term:.4f}"])


def _generate_pack(state: UniverseState, *, seed: Optional[int], out_dir: str, max_bulletins: int, max_missions: int) -> None:
    risk_model = RiskModel()
    galnet = GalNetSynth(seed=seed, risk_model=risk_model)
    bulletins = galnet.generate_bulletins(state, max_items=int(max_bulletins))
    digest = galnet.digest(state, bulletins)

    mission_synth = MissionSynth(seed=seed, risk_model=risk_model)
    missions = mission_synth.generate(state, bulletins, max_items=int(max_missions))

    _write_json(os.path.join(out_dir, "galnet_bulletins.json"), {"day": state.day, "items": [b.to_dict() for b in bulletins]})
    _write_text(os.path.join(out_dir, "galnet_digest.txt"), digest + "\n")
    _write_json(os.path.join(out_dir, "missions.json"), {"day": state.day, "items": [m.to_dict() for m in missions]})
    _write_nav_risk_csv(os.path.join(out_dir, "nav_risk.csv"), state, risk_model)


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="stellar-content-pipeline", description="Generate GalNet + missions from a universe state JSON snapshot")
    ap.add_argument("--input", required=True, help="Path to universe_state.json")
    ap.add_argument("--output", required=True, help="Output directory")
    ap.add_argument("--seed", type=int, default=None, help="Deterministic seed (optional)")
    ap.add_argument("--max-bulletins", type=int, default=14, help="Max bulletins")
    ap.add_argument("--max-missions", type=int, default=18, help="Max missions")

    ap.add_argument("--simulate-days", type=int, default=0, help="If >0, generate a multi-day timeline (days)")
    ap.add_argument("--step-days", type=float, default=1.0, help="Timeline step size (days)")
    ap.add_argument("--timeline-subdir", default="timeline", help="Subdirectory under output for timeline packs")

    args = ap.parse_args(argv)

    state = load_universe_state_json(args.input)

    out_dir = args.output
    _generate_pack(state, seed=args.seed, out_dir=out_dir, max_bulletins=args.max_bulletins, max_missions=args.max_missions)

    # Timeline mode
    sim_days = int(args.simulate_days)
    if sim_days > 0:
        sim = TimelineSimulator(seed=args.seed)
        timeline = sim.simulate(
            state,
            days=sim_days,
            step_days=float(args.step_days),
            bulletins_per_day=int(args.max_bulletins),
            missions_per_day=int(args.max_missions),
        )
        tl_root = os.path.join(out_dir, args.timeline_subdir)
        index = {"start_day": state.day, "days": sim_days, "step_days": float(args.step_days), "packs": []}

        for i, (st, bulletins, missions) in enumerate(timeline):
            day_dir = os.path.join(tl_root, f"day_{i:04d}")
            # Reuse generated objects to avoid re-rolling RNG inside synth.
            _write_json(os.path.join(day_dir, "galnet_bulletins.json"), {"day": st.day, "items": [b.to_dict() for b in bulletins]})
            digest = GalNetSynth(seed=args.seed).digest(st, bulletins)
            _write_text(os.path.join(day_dir, "galnet_digest.txt"), digest + "\n")
            _write_json(os.path.join(day_dir, "missions.json"), {"day": st.day, "items": [m.to_dict() for m in missions]})
            _write_nav_risk_csv(os.path.join(day_dir, "nav_risk.csv"), st, RiskModel())

            index["packs"].append({"i": i, "day": st.day, "dir": os.path.relpath(day_dir, out_dir).replace("\\", "/")})

        _write_json(os.path.join(tl_root, "timeline_index.json"), index)

    print(f"Wrote content pack to '{out_dir}'" + (f" (timeline: {sim_days} days)" if sim_days > 0 else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
