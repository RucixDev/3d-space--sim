- Added a **player Missile Warning Receiver (MWR) + Defensive Aids** (HUD + optional automation).
- New **HUD Settings → Combat symbology** options (persisted in `hud_settings.txt`): MWR enable, HUD indicator, evasion arrow, toast warnings, auto countermeasures + TTI threshold, prefer heat sinks.
- Combat HUD now shows an inbound-missile banner with **HEAT/RADAR**, **recommended countermeasure**, and **time-to-impact**, plus an optional jink arrow.
- Optional auto countermeasures share the same inventory + cooldown as manual deployment (prefers CHAFF vs radar seekers; FLARES/HEAT SINK vs heat seekers).
- Orbit Analyzer: added a **Forecast** panel that runs a short-horizon multi-body trajectory prediction and reports closest approaches, predicted impact, and dominant-gravity transitions (`apps/stellar_game/OrbitAnalyzerWindow.*`, `sim::detectDominantBodyTransitions`).

Verify:
- In combat, provoke incoming missiles at the player: confirm the banner appears; toggle arrow/auto-CM in HUD Settings and confirm behavior changes immediately.
- Open Orbit Analyzer → Forecast, click Recompute: confirm closest approach table populates and dominant-body transitions appear on close flybys.
- Tests: run `ctest -R test_hud_settings|test_trajectory_events`.
