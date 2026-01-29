## Round 143 - In-Flight Field Repairs (Jury-Rig Hull Patching)

This round adds a **new survival / exploration mechanic**: you can now perform **in-flight hull repairs** in normal space by consuming basic materials from your cargo hold.

- **Field Repairs** consume **Metals + Machinery** to patch hull in real-time, with a **soft cap** (default: **75% hull**) so full repairs still require station services.
- While active, **shields are forced offline** (power rerouted + heat venting), making the action risky to do under pressure.
- Repairs **pause automatically if you’re drifting too fast** (you need to slow down to resume).
- A compact **HUD overlay** shows progress, current hull/cap, and remaining materials, with a one-click cancel.
- Implemented as a new deterministic core sim module (**quote + step**) with unit tests for easy balancing.

### 🚀 What shipped

- **Core sim: `sim::FieldRepair`**
  - Deterministic quote + step functions (cap-aware, supply-aware).
  - Tunable parameters (rate, recipe, cap, heat-per-hull).
  - Unit tests covering cap clamping, supply limits, and step behavior.

- **stellar_game integration**
  - Action Wheel: **Field Ops** page now includes **Field Repairs** toggle with contextual detail.
  - Command Palette includes **Field Repairs: Toggle**.
  - Shields suppressed during repairs (like silent running).
  - Weapon fire / boost immediately **cancels repairs** for safety.
  - HUD overlay for repair progress + pause status.

### Files changed/added

- `include/stellar/sim/FieldRepair.h`
- `src/sim/FieldRepair.cpp`
- `tests/test_field_repair.cpp`
- `apps/stellar_game/main.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`

---

## Round 142 - Mesh Lab Undo/Redo History + MSVC Build Fixes

This round focuses on **tooling iteration** and **build stability**.

- The **Procedural Mesh Lab** now has a full **Undo/Redo history** for the SDF graph + procedural texture graph (and key mesher/preview knobs), with debounced commits and keyboard shortcuts.
- Fixed several recent MSVC compile breaks:
  - Missing `stellar/core/Hash.h` include (`core::hashCombine` / `core::fnv1a64` not found).
  - Raw OpenGL calls in app code now correctly use loaded `render::gl::*` procs (`glUseProgram` / `glBindVertexArray` / `glBindFramebuffer` not found).
  - Two missing braces in `apps/stellar_game/main.cpp` (unclosed **InputEvents** profiling scope and **Contacts** window scope).

### 🚀 What shipped

- **Procedural Mesh Lab history**
  - Snapshot-based history (graphs + knobs), capped to a configurable max.
  - Debounced *coalesce* window so rapid edits collapse into a single undo step.
  - Undo/redo buttons plus shortcuts: **Ctrl+Z** undo, **Ctrl+Y / Ctrl+Shift+Z** redo.
  - Undo first cancels any pending uncommitted edits back to the last stable baseline.

- **Build fixes**
  - `ProceduralMeshLabWindow.cpp` now includes the missing headers and uses `render::gl::*`.
  - `main.cpp` brace structure fixed so the main loop compiles again.

### Files changed/added

- `apps/stellar_game/ProceduralMeshLabWindow.h`
- `apps/stellar_game/ProceduralMeshLabWindow.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

---

## Round 141 - System Events Now Spawn Real Signals + GalNet Headless Build Fix

This round makes **system events physically manifest in space**: pirate raids, civil unrest, and research breakthroughs now produce extra **distress calls**, **wreckage / derelict salvage signals**, and (for trade booms/busts) **convoy traffic density** changes — driven by the same deterministic SystemEvent layer that powers GalNet bulletins.

It also fixes a const-correctness build error in headless builds (MSVC C2662) caused by GalNet calling `Universe::getSystem()` on a `const Universe`.

### 🚀 What shipped

- **Event-reactive signal generation (core sim)**
  - `generateSystemSignals()` now optionally accepts a `SystemConditionsSnapshot`.
  - When provided, it:
    - scales distress call density for certain events (**Pirate Raid**, **Civil Unrest**)
    - injects short-lived "aftermath" **derelict salvage** sites during chaotic events
    - modulates traffic convoy density for **Trade Boom** / **Trade Bust**
    - reduces distress noise during **Security Crackdowns**
  - Encounter planners are fed the **effective** security knobs (baseline + dynamics + event) so the content feels consistent with the system’s current state.

- **Gameplay integration**
  - The game now passes the current `SystemConditionsSnapshot` into the signal generator when seeding a system, so events immediately show up as extra signal contacts.

- **Headless build fix**
  - Added a `const` overload for `Universe::getSystem(...)` that forwards to the cached non-const implementation (logically const).
  - Fixes MSVC error **C2662** when compiling `sim::GalNet` bulletins.

- **Tests**
  - New `test_signal_event_reactivity` validates Pirate Raid boosts distress density and adds wreckage salvage sites when conditions are provided.

### Files changed/added

- `include/stellar/sim/Universe.h`
- `include/stellar/sim/Signals.h`
- `src/sim/Signals.cpp`
- `apps/stellar_game/main.cpp`
- `tests/test_signal_event_reactivity.cpp`
- `PATCH_NOTES.md`

---

## Round 140 - Actionable Pirate Ultimatums in Comms

This round makes pirate extortion messages interactive: you can now respond directly from the **Comms** inbox.

### ☠️ What shipped

- **Pirate "ULTIMATUM" messages are now actionable**
  - Selecting the ultimatium shows a **Live pirate channel** panel.
  - Clear progress bar: **delivered / required** tribute value and **time left**.

- **One-click response options**
  - **Auto-jettison tribute** uses the existing cargo jettison planner to meet the demanded value with minimal overpay.
  - **Refuse** immediately resolves the negotiation and triggers hostility (same behavior as the Threat HUD).

- **Safer cargo dumping**
  - Auto mode **avoids mission-reserved cargo by default**.
  - Optional toggle: **allow mission cargo** if you’re desperate.
  - Shows a quick **witness warning** when dumping is likely to be observed (may add bounty), matching the Threat HUD heuristic.

### Files changed/added

- `apps/stellar_game/CommsWindow.h`
- `apps/stellar_game/CommsWindow.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

---

## Round 139 - Silent Running (Stealth Mode) + Thermal/Sensor Tradeoffs

This round adds a *high-risk / high-reward* stealth mechanic inspired by space-sim classics: **Silent Running**.

### 🥷 What shipped

- **New gameplay toggle: Silent Running**
  - Default keybind: **Ctrl+X** (rebindable in Controls).
  - Also available via **Action Wheel → Navigation → Silent Running**.
  - While engaged:
    - **Shields are forced offline** (and will not regenerate).
    - Ship **cooling is heavily reduced**, and **baseline heat rises**.

- **Core sim upgrade: ThermalSystem now supports silent running**
  - Adds `ThermalInputs::silentRunning`.
  - Adds `ThermalParams::silentCoolMult` + `ThermalParams::heatPerSilentSec`.
  - Unit test coverage added.

- **Stealth has real systemic consequences**
  - **Radar sensor power** is reduced while silent running (you see less unless you ping).
  - **Cargo scan acquisition** is influenced by a ship “loudness” heuristic:
    - hull class, **speed**, **heat**, and **sensor ping** emissions
    - silent running reduces scan lock likelihood, but rising heat can counteract it over time
  - **Cargo scan duration** is longer while silent running (harder to hold a stable lock).
  - **Missile lock time** now scales with target “signature” and your own emissions:
    - ping helps locks; silent running hurts locks.

- **HUD feedback**
  - Radar shows a clear **SILENT** indicator when the mode is active.

### Files changed/added

- `include/stellar/sim/ThermalSystem.h`
- `src/sim/ThermalSystem.cpp`
- `tests/test_thermal.cpp`
- `apps/stellar_game/ControlsConfig.h`
- `apps/stellar_game/ControlsWindow.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`


## Round 138 - Station Security / Traffic Control (No-Fire Zone + Speeding + Trespass)

This round fleshes out an underdeveloped slice of "station life": approaching the mail-slot now has real consequences.

### 🚨 What shipped

- **New headless sim module: `sim::StationSecurity`**
  - Deterministic, unit-testable station enforcement logic that emits at most one event per tick (priority: weapons → trespass → speeding).
  - Implements:
    - **No-Fire Zone**: weapon discharge inside the station’s no-fire bubble.
    - **Traffic speed envelope**: overspeed is detected via a soft distance-based speed limit with tolerance + persistence timers.
    - **Docking-slot trespass**: entering the mail-slot tunnel without valid clearance.
  - Uses **cooldowns + accumulators** to avoid per-frame spam.
  - Uses a **graduated response** model (strikes): **warning → fine → bounty**.
  - Fine/bounty magnitudes are modulated by the station’s **law profile** (strictness/corruption/fine base).

- **Gameplay integration (stellar_game)**
  - Added station traffic enforcement to the main loop:
    - Violations generate **diegetic Comms** warnings from "Port Authority".
    - Fines are posted to the **Law Ledger** and apply rep penalties.
    - Repeated/serious violations escalate to **bounties**, causing station defenses to engage as intended.

- **Tests**
  - Added `tests/test_station_security.cpp` validating:
    - weapon discharge escalation (warning → fine → bounty)
    - speeding persistence threshold + fine escalation
    - docking-slot trespass detection + clearance suppression

## Round 137

### 🕵️ Black Market: Full Contraband Trading (Buy + Sell) + Stings

- **Implemented true black-market buying** for contraband (previously only black-market selling existed, and buying was intentionally blocked).
- **Unified black-market behavior across the game** by routing illegal trade through `sim::BlackMarket` logic instead of ad-hoc price multipliers.
- **Black market stings are now real events**: deals can be stings, triggering contraband enforcement (confiscation + fines) based on the station’s law profile.
- Added a new **Station Menu → Black Market** page with:
  - live fence availability + risk indicators
  - one-click “Sell all contraband” (excluding mission-reserved cargo)
  - a contraband hold summary with estimated fence payouts
- Updated the **Market UI** so illegal commodities can be **bought and sold via the black market** directly from the market table when a fence is available.
- Updated **trade loop automation** to support **illegal outbound manifests** (black market buys) and **illegal inbound sales**, including stings.

### 🧯 Refactor: Contraband Enforcement Side Effects

- Extracted enforcement side effects (rep, fines/ledger, police heat/alert, smuggle mission failure attribution) into a reusable helper so multiple systems (police scans, stings) share identical consequences.

## Round 136 - ShipScan: scan intel reports + identification memory

This round turns the in-flight **Scan** action into a meaningful, reusable mechanic (not just a mission checkbox).

### Highlights

- **New `sim::ShipScan` module (headless)**
  - Generates a deterministic-but-imperfect **scan report** from a snapshot of target state.
  - Reports include: **threat rating**, **cargo value estimate (with error band)**, and **jammer hints**.
  - Scan quality is driven by **sensor track strength** and is degraded by **jamming**.

- **Gameplay integration**
  - Completing a scan of a ship contact now posts a **Comms** message from **"Ship Scanner"** with a full report.
  - Scan intel is cached (in-session) and shown in:
    - **Radar tooltips**
    - **System Map → Selection details**

- **Identification memory**
  - A completed scan grants a temporary **ID memory window** so a scanned contact stays named for a while even if it later drops below the passive sensor identification threshold.

- **Tests**
  - Added `tests/test_ship_scan.cpp` validating scan quality, threat ordering, cargo error-band behavior, and jammer detection gating.

### Files changed/added

- `include/stellar/sim/ShipScan.h`
- `src/sim/ShipScan.cpp`
- `apps/stellar_game/main.cpp`
- `tests/test_ship_scan.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`

---

## Round 135 - StationServices: economy-backed repair/refuel/restock

This round turns the "quick docked services" (repair, refuel, countermeasures, ordnance) into **real, economy-backed
station transactions** instead of fixed credit-only button presses.

### Highlights

- **New `sim::StationServices` module (headless)**
  - Quotes and applies **hull repair**, **refuel**, **countermeasure restock**, and **missile rearm**.
  - Services consume actual station commodities via `econ::takeInventory(...)`.
  - Prices are based on `econ::quote(...)` (ask price) and then multiplied by the **effective station fee**.
  - Supports **partial fills** when limited by **station stock** or **player credits**.
  - Includes small combinatorial planners for "restock all" and "rearm all" so limited resources are spent sensibly.

- **New shared `sim::CountermeasureLoadout` helpers**
  - Centralizes the flare/chaff/heat-sink caps and provides clamp/restock helpers.
  - Removes duplicated cap logic from the game UI.

- **UI integration**
  - Both **Market Details → Services** and **Station Menu → Services** now use the headless module.
  - Buttons show `[stock]` / `[credits]` hints when purchases are limited.

- **Tests**
  - Added `tests/test_station_services.cpp` covering repair/refuel/restock/rearm constraints.

### Files changed/added

- `include/stellar/sim/CountermeasureLoadout.h`
- `include/stellar/sim/StationServices.h`
- `src/sim/StationServices.cpp`
- `apps/stellar_game/main.cpp`
- `tests/test_station_services.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`

---

## Round 134 - ElectronicWarfare: radar jamming + ghost returns

This round builds on the sensor/radar model by adding a light-weight **Electronic Warfare (EW)** layer.
Some NPCs can field **low-grade radar jammers** that suppress your effective sensor power and inject
plausible "ghost" noise returns on the HUD radar.

### Highlights

- **New `sim::ElectronicWarfare` module**
  - `computeJammingSnr(...)` + `jamming01FromSnr(...)`: inverse-square jammer field mapped to [0..1] via `snr/(1+snr)`.
  - `applyJammingToSensorPower(...)`: reduces effective sensor power (active ping partially pierces jamming).
  - `generateGhostBlips(...)`: deterministic, drifted noise blips driven by `(seed, time, range, jamming)`.
  - Pure/headless and unit-tested.

- **Radar HUD integrates jamming**
  - Aggregates jammer sources from nearby contacts (occluded jammers are heavily attenuated).
  - Displays a "JAM xx%" indicator when jammed.
  - Optional ghost returns drawn as non-selectable signal blips (alpha scales with noise strength).

- **Pirate packs can include a jammer ship**
  - Leaders are more likely to carry a jammer; occasional wing ships do as well.

- **New tests**
  - `tests/test_electronic_warfare.cpp` covers curve sanity + determinism + bounds.

## Round 133 - Station services UI robustness

- Fixed docked station services inventory call signature.
- Fixed countermeasure cap mismatch between new game defaults and HUD init.

## Round 132 - SensorModel: non-omniscient radar, occlusion + active ping

This round targets a very visible but under-developed system in many prototypes: **perfectly omniscient radar**.
We add a small, deterministic **sensor signal model** (range falloff + occlusion + smoothing) and wire it into the
HUD radar so contacts can be *ghosts*, become *identified* with proximity, and pop a tactical **active ping** sweep.

### Highlights

- **New `sim::SensorModel` module**
  - `computeSensorStrength01(...)`: inverse-square strength curve mapped to [0..1] via `snr/(1+snr)`
  - `updateSensorTrack(...)`: exponential smoothing + identification hysteresis (prevents flicker)
  - Pure/headless (no renderer dependencies), so it is easy to unit-test and reuse.

- **Radar HUD is no longer omniscient for ships**
  - Contacts are filtered through the sensor track.
  - Unidentified contacts show as **generic signal blips** until identified.
  - Icon alpha scales with sensor strength (low confidence feels "faint" instead of binary on/off).
  - Planets/stations can **occlude** line-of-sight, attenuating detection.

- **Implemented the existing `sensorPing` control action (default: Ctrl+O)**
  - Emits a short active ping that boosts sensor power for a moment.
  - Draws a sweep ring on the radar and posts an IntegrationHub event.

- **Tests**
  - Added `tests/test_sensor_model.cpp` to validate the strength curve, occlusion attenuation, smoothing half-life, and identification hysteresis.

---

## Round 130 - TrafficLanes: corridor-bundled dual-carriageway lanes + geometry API

This round targets a very visible but under-developed system: **in-system traffic lanes**.
Previously, each convoy picked a lane arc based on its unique id, which meant traffic looked like a set of
independent, slightly-random curves rather than a coherent set of **shared corridors**.

We introduce **corridor-bundled lanes**: lane geometry is now derived from the *station pair* so multiple convoys
between the same endpoints follow the same corridor. We also add an optional **dual-carriageway** behavior so
opposite-direction traffic uses mirrored arcs (reducing head-on overlap and making the lane map read better).

### Highlights

- **New `sim::TrafficLaneGeometry` + `computeTrafficLaneGeometry(...)`**
  - Computes and exposes a stable **laneKey**, a chord-aligned frame (`dir/side/up`), and arc parameters.
  - Enables efficient sampling/evaluation when tools/UI need many points.

- **Corridor-bundled lane mode (default-on)**
  - New `TrafficLaneParams::bundleByStationPair` groups convoys into shared corridors using
    `(systemId, min(from,to), max(from,to))`.

- **Dual-carriageway lanes (default-on)**
  - New `TrafficLaneParams::dualCarriageway` mirrors the corridor for reverse direction, producing distinct
    inbound/outbound lanes.

- **More stable RNG consumption**
  - Lane direction sampling now uses a fixed-cost spherical method (no rejection sampling), making it easier
    to extend the lane model without accidentally perturbing unrelated random draws.

- **Tests**
  - Added new coverage validating corridor bundling, dual-carriageway sign flip, mid-flight arc sanity,
    and legacy (unbundled) behavior.

### Files changed/added

- `include/stellar/sim/TrafficLanes.h`
- `src/sim/TrafficLanes.cpp`
- `src/sim/TrafficConvoyLayer.cpp`
- `tests/test_traffic_lane_geometry.cpp`
- `PATCH_NOTES.md`

---

## Round 129 - FormationPlanner: min-cost convoy escort slot assignment (no crossovers)

This round focuses on an under-developed, very visible behavior: **convoy police escorts flying formation**.
Previously, escorts were assigned slots by an id-ordering heuristic, which could cause ships to **swap sides** and
"cross" through the formation (and occasionally bump) whenever the ordering disagreed with their current positions.

We add a small, headless planner that solves a **min-cost assignment** (Hungarian / Kuhn–Munkres) between escorts
and candidate formation slots. The result is a more readable, stable wing with less crossing and better spacing.

### Highlights

- **New headless `sim::FormationPlanner` (Hungarian assignment)**
  - Builds candidate slot targets using the existing `Formation.h` helpers.
  - Assigns members to slots by minimizing total travel distance, using **integer-quantized costs** for determinism.
  - Optional **sticky slot hints** (penalty-based) to avoid near-tie thrash.

- **Game integration (police convoy escorts)**
  - `stellar_game` now precomputes an optimal slot index for each escort in a follow-wing each frame.
  - Adds a small per-escort **slot memory** so wings keep their structure unless there's a clear benefit to swapping.
  - Formation radius now scales gently with wing size for additional clearance.

- **Tests**
  - New coverage for side-swap avoidance, determinism, and stickiness tie-breaking.

### Files changed/added

- `include/stellar/sim/FormationPlanner.h`
- `src/sim/FormationPlanner.cpp`
- `tests/test_formation_planner.cpp`
- `apps/stellar_game/main.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`

---

## Round 128 - MissileDefense: deterministic evasion vectors (smarter jinks)

This round focuses on an under-developed combat layer: **missile defense beyond detection**.
We add a small, headless helper that computes a **deterministic evasion direction** based on
closest-approach geometry, and wire it into the game so NPCs jink more intelligently.

### Highlights

- **New headless `sim::planMissileEvasion(...)`**
  - Computes a unit vector that moves the target **away from the predicted closest-approach point**
    under a constant relative-velocity estimate.
  - Optional projection into the **plane perpendicular to the line-of-sight** to produce a more
    lateral jink (tends to increase LOS rate against PN-like guidance).
  - Deterministic, seeded tie-breaking for near head-on cases.

- **Game integration**
  - `stellar_game` now uses `planMissileEvasion()` when refreshing an NPC’s evasion direction,
    replacing the ad-hoc `cross(missileVelDir, toTargetDir)` heuristic.

- **Tests**
  - Added coverage for offset-pass behavior, LOS-lateral constraint, and determinism.

### Files changed

- `include/stellar/sim/MissileDefense.h`
- `src/sim/MissileDefense.cpp`
- `tests/test_missile_defense.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

---

## Round 123 - Dynamic resolution render scale (FPS boost)

This round introduces a **render-scale / dynamic resolution** path that can noticeably improve FPS by rendering the 3D scene + post-processing at a lower internal resolution and upscaling in the final compositor.

### Highlights

- **PostFX: render scale + DRS controls**
  - Manual **Render scale** (0.25..1.0)
  - Optional **Dynamic resolution** that nudges the render scale toward a target FPS using a smoothed CPU frame-time signal.
  - Rate-limited scale changes to avoid resize thrash.

- **Bloom now respects render scale**
  - Bloom’s half-res chain is sized relative to the **scene** resolution (not the window), avoiding accidental upsampling work.

- **Retro compositor fixed for scaled scene textures**
  - Retro mode used `texelFetch` under the assumption the scene texture matched window resolution.
  - It now maps output pixels → scene texels correctly when renderScale < 1.

- **Optional upscale sharpening**
  - Cheap unsharp-mask style filter (default off) to offset the blur from upscaling.

### Files changed/added

- `include/stellar/render/PostFX.h`
- `src/render/PostFX.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

---

## Round 121 - FPS Boost: Fast Nebula Updates + Streaming Point Buffers

This round focuses on **real, measurable frame-rate improvements** in the most per-frame expensive procedural background layers (nebula + point sprites).

### Highlights

- **NebulaField: cached turbulence noise (time-sliced refresh)**
  - Previously, the nebula evaluated **fBm Perlin noise for every puff every frame**.
  - Now, each puff stores a cached noise value, and we **refresh only a small batch per frame** using exponential smoothing.
  - Result: similar drifting/breakup look with dramatically lower CPU cost at typical puff counts.

- **PointRenderer: VBO capacity caching + MapBufferRange streaming uploads**
  - Avoids reallocating the point VBO each draw call.
  - Grows the buffer only when needed, then streams updates through `glMapBufferRange` (fallback to `glBufferSubData`).
  - Benefits all point-heavy layers: **starfield**, **nebula**, and **CPU particles**.

### Files changed

- `include/stellar/render/Nebula.h`
- `src/render/Nebula.cpp`
- `include/stellar/render/PointRenderer.h`
- `src/render/PointRenderer.cpp`
- `PATCH_NOTES.md`

---

## Round 119 - Integration Hub Automation Persistence (SaveGame v32)

This round upgrades the **Integration Hub** from a session-only devtool into a persistent **player workflow / macro** system by saving and restoring automation rules through **SaveGame** (quicksave/quickload).

### Highlights

- **SaveGame v32: Integration Hub automation serialization**
  - New keys in the save format:
    - `hubAutomationsEnabled` (bool)
    - `hubRules` (count)
    - `hub_rule` and `hub_action` records (strings base64-encoded like comms)
  - Strings (rule names, tag text, message templates) use the existing base64 token strategy so the save file stays whitespace-safe.

- **Stable enum IDs for persisted automation payloads**
  - `GameEventKind` and `GameActionKind` now have **explicit numeric values** and a "do not reorder" comment.
  - This prevents future enum edits from silently breaking saved automation rules.

- **Quicksave/quickload now preserves Integration Hub rules (v32+)**
  - Quicksave captures the current automation rule-set.
  - Quickload clears runtime queues (events/actions/scheduled actions) and restores automations for saves with `version >= 32`.

- **New test: SaveGame Integration Hub roundtrip**
  - Adds `tests/test_savegame_hub_automations.cpp` to validate that hub rules survive a save/load cycle.

### Files changed

- `apps/stellar_game/GameSignals.h`
- `apps/stellar_game/IntegrationHubWindow.h`
- `apps/stellar_game/main.cpp`
- `include/stellar/sim/SaveGame.h`
- `src/sim/SaveGame.cpp`
- `tests/test_savegame_hub_automations.cpp`
- `PATCH_NOTES.md`

---

## Round 118 - NavAssist Follow Mode (Escort-Friendly Autopilot)

This round adds a **Follow** mode to the Navigation Assist computer and wires it into the **Traffic Escort** HUD so escort contracts can auto-engage a stable formation-follow behavior.

### Highlights

- **NavAssistComputer: new `Follow` mode**
  - Maintains a configurable distance behind a target ship using a smooth pursuit/offset intercept.
  - Uses the existing flight controller pipeline (heading/velocity) to keep behavior consistent with other assist modes.

- **Traffic Escort HUD integration**
  - Adds a one-click **Target + Follow** action for convoy ships.

### Files changed

- `include/stellar/sim/NavAssistComputer.h`
- `src/sim/NavAssistComputer.cpp`
- `apps/stellar_game/main.cpp`

---

## Round 117 - Diagnostics Capture Bundle (Cross-System Debugging Glue)

This round adds a **one-click diagnostics bundle** that cross-integrates the game’s existing systems: the **Integration Hub**, **Flight Recorder**, **Runtime Validation Repro Pack**, and the **Screenshot pipeline**.

The goal is simple: when something looks wrong, you can capture *everything needed* to reproduce or share the issue in one folder — with stable filenames and a manifest.

### Highlights

- **New GameAction: `CaptureBundle`**
  - Creates a unique folder under an output root (default `captures/`) using a sanitized label + timestamp.
  - Writes a `manifest.json` describing the capture.
  - Optionally captures:
    - World screenshot (`world.png`)
    - UI screenshot (`ui.png`)
    - Flight Recorder trace (`flight_trace.json`)
    - Integration Hub trace (`integration_trace.json`)
    - Runtime Validation repro pack snapshot (`repro_pack.json`)
  - Optionally **stops** the Flight Recorder before exports and **copies the bundle path** to clipboard.

- **Integration Hub: new automation action + quick capture button**
  - `CaptureBundle` is available as an automation action kind, with checkboxes for bundle contents and behavior.
  - Event inspector adds a **Bundle** button next to the screenshot buttons for fast “capture everything” moments.
  - Starter rule updated: **Validation Watchdog → Capture bundle** (disabled by default).

- **Runtime Validation: export helper + capture button**
  - Adds `exportReproPackJsonToPath(...)` so other systems can request a repro pack at a specific file path.
  - The Runtime Validation window now includes a **Capture bundle** button next to **Dump now**.

### Files changed

- `apps/stellar_game/GameSignals.h`
- `apps/stellar_game/IntegrationHubWindow.cpp`
- `apps/stellar_game/RuntimeValidationWindow.h`
- `apps/stellar_game/RuntimeValidationWindow.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

---

## Round 116 - Time Trial Ghost Racing (Cross-System Gameplay + Telemetry)

This round cross-integrates **Time Trials**, the existing **Flight Recorder sample format**, and the in-world **ghost rendering pass** to add a proper **best-run “ghost”** you can race against.

### Highlights

- **Best-run ghost replay (per course)**
  - While a run is active, the game records a lightweight telemetry strip of the player ship (position/orientation/velocity).
  - If you set a **new personal best**, that run becomes the stored **best ghost** for that course.

- **In-world rendering (ship + trail)**
  - The best-run ghost is rendered as a translucent ship (teal) using the same instance path as the Flight Recorder ghost.
  - A decimated trail is also drawn in-world for easy line-of-flight comparison.

- **Objective HUD “ghost split”**
  - The Objective HUD now shows a simple **ahead/behind** split versus the ghost by comparing distance to the **next objective**:
    - next gate during the run, or
    - the anchor station during docking finishes.

- **Time Trials UI controls**
  - New **Ghost (best run)** section:
    - enable/disable, ship/trail toggles
    - record rate (Hz) + lead/lag
    - reset best time + ghost for the current course

### Files changed

- `apps/stellar_game/TimeTrialWindow.h`
- `apps/stellar_game/TimeTrialWindow.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

---

## Round 115 - Delayed Automations + Nav Auto-run Phase Events (Cross-System Orchestration)

This round focuses on the *glue code* that makes disparate systems feel like one coherent game:

- the **Integration Hub automation engine** can now **sequence actions over time** (not just instantly), and
- the existing **Nav Auto-run** state machine now emits structured **Integration Hub events** for its key phases.

Together, this enables truly hands-free workflows (travel → supercruise → docking) that are **traceable**, **automatable**, and **recordable** in the Flight Recorder.

### Highlights

- **New concept: scheduled actions (delayed execution)**
  - The Integration Hub now maintains a small **scheduled queue** for actions whose `tRealSec` is in the future.
  - Each frame, due actions are moved into the pending queue and executed normally.
  - Actions tab now shows a **Scheduled** section (with countdown), plus clear/copy utilities.

- **Automation actions now support `delaySec`**
  - Every automation action template can optionally set a **Delay (sec)**.
  - This allows reliable sequencing like:
    - `StopFlightRecorder` → *(+0.15s)* `ExportFlightRecorderTrace`
    - `RequestScreenshot` → *(+0.25s)* `ExportIntegrationTrace`

- **Integration trace export upgraded**
  - Trace JSON schema bumped to **version 4**.
  - Adds optional `scheduledActions` export.
  - Adds `delaySec` field to exported automation actions.

- **Nav Auto-run emits phase events (cross-integrates navigation → automation → camera/comms)**
  - The Nav Auto-run loop now pushes `GameEventKind::NavAssist` events for key transitions:
    - `AutoRunSupercruise`
    - `AutoRunDockingComputer`
    - `AutoRunReplot`
    - `AutoRunStopRange`
    - `AutoRunStopFuel`
    - `AutoRunComplete`
  - These events appear in the Integration Hub timeline and can drive automations.

- **New starter automations (disabled by default)**
  - `AutoRunSupercruise` → `SetCameraRigPreset(Travel)`
  - `AutoRunDockingComputer` → `SetCameraRigPreset(Docking)`
  - `AutoRunStop*` → `TransmitComms` (overlay)

### Files touched
- `apps/stellar_game/IntegrationHubWindow.h`
- `apps/stellar_game/IntegrationHubWindow.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

---

## Round 114 - Mission Deep Links + Objective Quick Actions (Cross-System Integration)

This round cross-integrates **Missions**, **Comms**, the **Objective HUD**, and the **Integration Hub** action pipeline so gameplay intent ("track this mission", "route there", "auto-run", "capture the run") can flow through *one* observability-friendly command path.

### Highlights

- **New Integration Hub action: `SetTrackedMissionId`**
  - Allows UI surfaces and automations to set the currently tracked mission (Objective HUD / Ship-Status mission summary).
  - Payload:
    - `u64a`: mission id (0 clears)

- **Comms mission deep-links**
  - Mission briefings (Comms channel: Mission) now surface actionable buttons when `CommsMessage::relatedId` is present:
    - **Track mission**
    - **Plot mission** (Sync nav to mission’s next stop)
    - **Auto-run mission** (arms the existing nav auto-run pipeline)

- **Objective HUD: one-click mission operations**
  - Added a new **Quick actions** strip:
    - **Plot mission**
    - **Auto-run mission** (also switches camera preset to *Travel*)
    - **Capture + Auto-run** (clears Integration Hub log + starts Flight Recorder + auto-runs)
  - This creates a tight end-to-end loop: *Objective HUD → Integration Hub → Nav Auto-run / Camera Rig / Flight Recorder*.

- **Mission & logistics UI now routes through the Integration Hub**
  - Mission list "Plot route" and "Track" go through the same command pipeline (better traceability + automation).
  - Cargo sourcing helpers now provide **Plot source** + **Auto-run source** via the Integration Hub (cross-integrates economy → navigation → camera).

- **New starter automation rule (disabled by default)**
  - `MissionAccepted` → `SetTrackedMissionId` (tracks newly accepted missions automatically, through the same pipeline).

### Files touched
- `apps/stellar_game/GameSignals.h`
- `apps/stellar_game/IntegrationHubWindow.cpp`
- `apps/stellar_game/CommsWindow.h`
- `apps/stellar_game/CommsWindow.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

---

## Round 112 - Actionable Comms + GoToStation Travel Action (Cross-System Gameplay Integration)

This round focuses on **cross integrating existing systems into a tighter gameplay loop** by making
in-universe **Comms** messages *actionable*.

### Highlights

- **New Integration Hub action: `GoToStation`**
  - A high-level travel command that cross-integrates:
    - route planning (system hops)
    - station targeting (on arrival, or immediately if already in-system)
    - optional **Nav Auto-run** (hands off the trip to the existing auto-run + docking computer flow)
  - Payload:
    - `u64b`: system id
    - `u64a`: station id
    - `b`: arm auto-run

- **Comms → Navigation → Docking (one-click)**
  - If a comms message includes a location, the UI now shows:
    - **Target** (if already in the system) / **Plot route** (if not)
    - **Go & Dock** (arms auto-run so you travel and dock hands-free)
  - This makes mission briefings and transmissions a first-class navigation entry point.

- **Integration Hub automation upgraded: 2× u64 fields**
  - Automation action templates can now drive both `u64a` **and** `u64b` via independent sources
    (constant, event u64a, event u64b).
  - Enables automations for multi-id actions like `GoToStation` without hacks.
  - **Integration trace JSON version bumped to 3** to reflect the new schema.

### Files touched
- `apps/stellar_game/GameSignals.h`
- `apps/stellar_game/IntegrationHubWindow.h`
- `apps/stellar_game/IntegrationHubWindow.cpp`
- `apps/stellar_game/CommsWindow.h`
- `apps/stellar_game/CommsWindow.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

---

## Round 111 - Runtime Validation → Integration Hub → Screenshot Capture (Cross-System Debug Workflow)

This round strengthens the **debug / QA workflow** by cross-integrating three existing systems:

- **Runtime Validation** (NaN/Inf watchdog + smoke checks)
- **Integration Hub** (event/action timeline + automations)
- **Photo Mode / Screenshot pipeline** (world/UI capture)

### Highlights

- **New Integration Hub action: `RequestScreenshot`**
  - Bridges automations and tools into the existing screenshot capture path.
  - Supports UI vs world capture, timestamped filenames, optional pause-for-capture, and copy-path-to-clipboard.

- **Integration Hub: event detail quick actions**
  - Timeline entries now have one-click buttons to request a screenshot (UI/world) with a safe, event-derived basename.

- **Runtime Validation emits through the Integration Hub event stream**
  - Watchdog failures and repro-pack results now flow through `hubPushEvent`, so:
    - automations can react to validation events
    - Flight Recorder marker capture can correlate telemetry with validation spikes

- **Runtime Validation: optional screenshot-on-failure**
  - Configurable “On watchdog failure → request screenshot(s)” controls:
    - UI and/or world capture
    - pause-for-capture
    - timestamp/copy-path
    - out dir + base name

- **Starter automation rule (disabled by default)**
  - `Validation:Watchdog` → RequestScreenshot (UI) + export traces (Integration/Flight)

### Files touched
- `apps/stellar_game/GameSignals.h`
- `apps/stellar_game/IntegrationHubWindow.h`
- `apps/stellar_game/IntegrationHubWindow.cpp`
- `apps/stellar_game/RuntimeValidationWindow.h`
- `apps/stellar_game/RuntimeValidationWindow.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

---

## Round 110 - Integration Hub → Flight Recorder Marker Track (Cross-System Correlation)

This round cross-integrates the **Integration Hub** event stream with the **Flight Recorder**, enabling
rapid correlation between *telemetry* (position, velocity, orientation) and *discrete gameplay/devtools*
events (docking clearance, time trials, validation alerts, mission sync, etc.).

### Highlights
- **Event mirroring hook**
  - `IntegrationHubWindowState` now exposes an optional `onEvent` sink called before an event is moved into
    the hub log. This keeps the Integration Hub dependency-light while still allowing external systems to
    subscribe.
- **Flight Recorder markers**
  - Added a marker ring buffer that captures Integration Hub events during recording (filterable / size limited).
  - New UI section lists markers and lets you click to scrub the replay playhead directly to that moment.
- **Trace export upgrade**
  - Flight recorder trace export can now include markers as **instant events** on a dedicated "Markers" track
    in the Chrome/Perfetto trace JSON (toggleable via "Trace: include markers").

### Files touched
- `apps/stellar_game/IntegrationHubWindow.h`
- `apps/stellar_game/FlightRecorderWindow.h`
- `apps/stellar_game/FlightRecorderWindow.cpp`
- `apps/stellar_game/main.cpp`

---

## Round 109 - Mission Nav Sync (Cross-System Integration)

This round cross-integrates the **Mission Tracker**, the **route planner**, and the existing **Nav Auto-run + Docking Computer** pipeline.

### Highlights
- **New Integration Hub action: `SyncNavToMission`**
  - Given a `MissionId`, it plots the route to the mission’s *next stop* (via leg for passenger/multi-delivery, otherwise final destination).
  - Arms `pendingArrivalTargetStationId` (auto-target on jump completion) and immediately targets the station if you’re already in-system.
  - Optional “arm auto-run” flag to hand the trip off to the existing auto-run flow.
- **Unified codepath** for mission route plotting
  - Mission Board “auto plot on accept” and Mission Tracker “Plot Route” now go through the same action pipeline (less duplication, better observability in Integration Hub).
- **Mission Tracker one-click convenience**
  - Added **Auto-run** button next to Plot Route to arm auto-run for the tracked mission.
- **Starter automation rule (disabled by default)**
  - `MissionAccepted -> SyncNavToMission` for users who want the game to auto-sync nav when accepting missions.

### Files touched
- `apps/stellar_game/GameSignals.h`
- `apps/stellar_game/IntegrationHubWindow.cpp`
- `apps/stellar_game/main.cpp`

---

## Round 107 - Cross Integration: Integration Hub ↔ Comms

This round focuses on **cross integrating existing code and systems** by connecting the Integration Hub action/event layer to the in-universe **Comms** system.

You can now treat Comms as a first-class automation target:

- New `TransmitComms` **GameAction** (subject/body via `subject|body`) that routes actions into the Comms log and optionally pops the incoming overlay.
- New *starter* automation rules (disabled by default) that demonstrate end-to-end integration:
  - **Docking Clearance → Comms** (ATC-style clearance response)
  - **Mission Complete → Comms** (wraps completion events as a transmission)
  - **TimeTrial Finish → Comms** (race control summary)

### Files changed/added

- `apps/stellar_game/GameSignals.h`
- `apps/stellar_game/IntegrationHubWindow.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

---

## Round 98 - Reliability: CI + Sanitizers + Repro Packs

This round continues the push to **make sure existing code is working as intended** by adding:

- A modern **GitHub Actions CI** pipeline (CMake + CTest) that builds headless on Ubuntu/Windows, plus a Clang **ASan+UBSan** job.
- Optional CMake sanitizer toggles (`STELLAR_ENABLE_ASAN`, `STELLAR_ENABLE_UBSAN`) so local dev builds can match CI.
- A **Repro Pack** exporter in the Runtime Validation window that writes a small JSON bundle (ship state + optional Flight Recorder telemetry) when the watchdog trips.

### What’s new

- **CI: CMake + CTest workflows**
  - Headless build + tests on **Ubuntu + Windows**
  - Clang **AddressSanitizer + UndefinedBehaviorSanitizer** build that runs the full test suite
  - Optional Windows render build to catch UI/render integration issues at compile-time

- **Runtime Validation: Repro Pack JSON**
  - Toggle “Dump on watchdog failure” and specify an output filename
  - Optional “Unique name per hit” and “Include Flight Recorder telemetry”
  - Manual “Dump now” button for bug reports

### Files changed/added

- `.github/workflows/c-cpp.yml`
- `CMakeLists.txt`
- `apps/stellar_game/RuntimeValidationWindow.h`
- `apps/stellar_game/RuntimeValidationWindow.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

---

## Round 97 - Reliability: Runtime Validation + Regression Tests

This round focuses on **making sure existing code is working as intended** by adding:

- An in-game **Runtime Validation** dashboard (watchdog + smoke checks) to catch common sim corruption early.
- New **unit tests** for the most recent under-tested systems (Aerodynamics + Time Trials) so regressions get caught immediately.

### What’s new

- **New: Runtime Validation (UI → Runtime Validation)**
  - Continuous **NaN/Inf watchdog** for ship state (pos/vel/orient/angVel) with optional auto-pause.
  - Optional **magnitude watchdog** to detect runaway values without modifying the simulation.
  - **Deterministic smoke checks** you can run and copy as a report:
    - TimeTrial gate pass logic and generator determinism/invariants
    - Aerodynamics sign/invariant sanity checks

- **Tests: regression coverage**
  - Added `stellar_tests` coverage for:
    - `sim::computeAerodynamics(...)`
    - `sim::timeTrialGatePassed(...)` and deterministic course generation

### Files changed/added

- `apps/stellar_game/RuntimeValidationWindow.h` *(new)*
- `apps/stellar_game/RuntimeValidationWindow.cpp` *(new)*
- `tests/test_aerodynamics.cpp` *(new)*
- `tests/test_time_trial.cpp` *(new)*
- `apps/stellar_game/main.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`

---

## Round 96 - 3D Gameplay: Time Trials (Gate Courses)

This round focuses on **moment-to-moment 3D gameplay**: a lightweight, deterministic **gate course** system you can run as a time trial.
Courses are generated from a stable seed (system + station + user seed), then rendered as HUD gate markers and tracked as an objective.

### What’s new

- **New: Time Trials (Windows → Time Trials)**
  - Generate a deterministic station slalom course:
    - gate count, gate radius, course radius, height, jitter, loops
    - stable seed option (system + station) plus a user seed for shareable variants
  - Run phases: **Ready → Running → Finished**, with best time tracking per-course.
  - Quality-of-life controls:
    - Copy shareable **Course Code** to clipboard
    - Reset best time for current course

- **HUD integration**
  - 3D gate markers rendered in the forward HUD (with offscreen clamping).
  - Objective HUD switches to a compact time-trial status panel while a course is active.

### Files changed/added

- `include/stellar/sim/TimeTrial.h` *(new)*
- `src/sim/TimeTrial.cpp` *(new)*
- `apps/stellar_game/TimeTrialWindow.h` *(new)*
- `apps/stellar_game/TimeTrialWindow.cpp` *(new)*
- `apps/stellar_game/main.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`

---

## Round 95 - 3D Camera Rig: Orbit Controls + Horizon Lock + Body Avoidance

This round focuses on the **3D camera**—a notably under-developed gameplay layer compared to the sim/proc systems.
It replaces the hardcoded chase camera with a small but powerful **camera rig** that supports
**orbit controls**, **gravity-aware horizon stabilization**, and **large-body collision avoidance**.

### What’s new

- **New: 3D Camera Rig (Tools → 3D Camera Rig)**
  - **Orbit camera** with `Alt + LMB` orbit, `Alt + MMB` pan, `Alt + wheel` zoom.
  - **Horizon lock** driven by a reference “up” vector:
    - near planets: align to **gravity up**
    - in deep space: fall back to **orbital angular momentum** normal when available
  - **Springy smoothing** for position/orientation/FOV using half-life parameters.

- **Safety: Avoid dominant gravity body**
  - Camera is kept outside the rendered radius of the currently dominant gravity body
    (star/planet) with a small padding, preventing “inside the planet” clips.

- **Better depth precision**
  - Dynamic **near plane** derived from camera distance to the ship improves z precision
    (reduces flicker) without requiring engine-level reverse-Z changes.

- **Cinematic camera compatibility**
  - When the Cinematic Camera tool is enabled, the rig can optionally consume its
    offset/FOV samples (still benefits from smoothing + body avoidance).

### Files changed/added

- `apps/stellar_game/CameraRigWindow.h` *(new)*
- `apps/stellar_game/CameraRigWindow.cpp` *(new)*
- `apps/stellar_game/main.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 94 - 3D Flight Model: Aerodynamic Lift + Stability (Atmospheric)

This round deepens the **flight model** once you leave vacuum: ships now experience **aerodynamic lift**,
**induced drag**, and optional **stability / control-surface moments** when flying through an atmosphere.
The result is a controllable “bank-to-turn” feel—while keeping the sim deterministic and headless-friendly.

### What’s new

- **New `sim::Aerodynamics` module (CPU / deterministic)**
  - Computes lift using the classic relation `L = q · S · CL` and adds induced + stall drag.
  - Provides simple stability:
    - rotational damping proportional to dynamic pressure
    - optional “weathervane” alignment to the velocity vector
    - optional control-surface authority (maps pilot torque input to aero angular accel)

- **Integrated into ship physics**
  - `Ship::stepWithExternalForces(...)` now supports external *angular* acceleration (in addition to linear),
    enabling aerodynamic moments (and future effects like wind gusts).

- **UI: Physics (experimental)**
  - Toggle **Aerodynamic lift + stability** and tune wing area, lift slope, stall angle, induced drag, etc.
  - Optional **Affect NPC ships** toggle.

- **UI: Ship → Status → Gravity / Orbit**
  - Atmosphere section now shows AoA/sideslip, CL, stall %, lift g-load, and extra drag when aero is enabled.

### Files changed/added

- `include/stellar/sim/Aerodynamics.h` *(new)*
- `src/sim/Aerodynamics.cpp` *(new)*
- `include/stellar/sim/Ship.h`
- `src/sim/Ship.cpp`
- `apps/stellar_game/main.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 93 - ImGui Build Plumbing + Build Info Window

This round targets a rough edge in the **render+UI build plumbing** (a notably under-developed part of the repo compared to the sim/proc layers):
fixing the common `imgui.h` include failure when compiling the core library with `STELLAR_ENABLE_IMGUI=ON`.

It also adds a small but high-leverage developer tool: an in-game **Build Info** window with a one-click copy-to-clipboard summary for bug reports.

### What’s new

- **Fix: `imgui.h` include failure in core library builds**
  - When `STELLAR_ENABLE_IMGUI=ON`, the core `stellar` target now links the ImGui target so the ImGui include paths are available while compiling core sources that have optional ImGui-backed utilities (e.g. `ui/TextFx`).
  - Added a configure-time sanity guard: if ImGui is enabled but no ImGui target is available, we disable ImGui with a clear warning (instead of failing later with a missing header error).

- **New: Build Info window (UI → Build Info)**
  - Shows **platform, architecture, compiler, build config**, and key **feature flags**.
  - Shows dependency versions (**SDL** + **Dear ImGui**) and, optionally, **OpenGL runtime strings** (vendor/renderer/version).
  - Includes a **Copy summary to clipboard** button for easy sharing in logs/issues.

### Files changed/added

- `CMakeLists.txt`
- `apps/stellar_game/BuildInfoWindow.h` *(new)*
- `apps/stellar_game/BuildInfoWindow.cpp` *(new)*
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

## Round 92 - Resonant Asteroid Belts: Kirkwood Gaps + Trojan Swarms

This round expands **procedural generation** into one of the most under-developed “world building” layers in the sim:
**minor bodies**. Systems now get deterministic **asteroid belts / debris discs** that express *dynamical structure*:
resonance gaps, resonant ridges, and L4/L5 trojan swarms.

### What’s new

- **New `proc::AsteroidBeltGenerator`**
  - Generates a small **belt plan** per system:
    - **Main belt** between inner terrestrials and the first gas giant (when possible)
    - **Outer belt** beyond the outermost planet (Kuiper-ish)
    - **Trojan swarms** around the dominant planet (two lobes at ±60°)
    - Fallback **debris disk** for sparse layouts
  - Adds curated **mean-motion resonance features**:
    - **Interior “Kirkwood-style” gaps** for the main belt (2:1, 3:1, 5:2, 7:3, ...)
    - **Exterior resonant ridges** for the outer belt (3:2, 2:1, 5:2, ...)
  - Everything is deterministic from *(universeSeed + system.stub.seed + belt id)*.

- **Density field + deterministic sampling**
  - A belt is defined by a **continuous density function** `asteroidBeltDensity01(a, theta)`.
  - Point clouds are generated with a **Mitchell-style best-candidate sampler** using a low-discrepancy (Halton) candidate stream.
  - Output reads “blue-noise-ish”, while still respecting gaps/ridges.

- **Procedural System Lab UI: Minor Bodies**
  - New panel listing generated belts, showing:
    - span + thickness
    - controlling planet (if any)
    - resonance list (gap vs ridge)
  - Live **radial density plot** (mean over azimuth)
  - 2D **belt-plane scatter** with optional resonance ring overlay

- **Tests**
  - Added a CPU-only test to validate determinism + verify resonance dips behave as expected.

### Files changed/added

- `include/stellar/proc/AsteroidBeltGenerator.h` *(new)*
- `src/proc/AsteroidBeltGenerator.cpp` *(new)*
- `apps/stellar_game/ProceduralSystemLabWindow.h`
- `apps/stellar_game/ProceduralSystemLabWindow.cpp`
- `tests/test_asteroid_belts.cpp` *(new)*
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 91 - Nebula 2.0: Coherent Filaments + Importance-Sampled Puffs

This round focuses on **procedural sky generation** by upgrading the cheap point-sprite nebula layer into something that reads like
actual gas clouds: **filaments, cavities, and large coherent masses**—still deterministic, still light enough for real-time.

### What’s new

- **Deterministic importance sampling for nebula puffs**
  - Generates a **low-discrepancy (Halton) candidate pool** on the sphere.
  - Evaluates a coherent **3D density field** (domain-warped ridged fBm + large blobs + cavity carving).
  - Keeps the highest-density candidates, producing natural-looking structure without unbounded rejection loops.

- **Density-aware appearance**
  - Puff **alpha and size** scale with the sampled density so filaments read stronger and voids stay wispy.
  - Coherent HSV palette fields keep neighboring puffs in related color bands.

- **Smoother, more coherent animation**
  - Update pass uses animated **3D fBm turbulence** for radial jitter + alpha breakup (instead of 2D-only drift).
  - Dense regions wobble slightly less, helping the silhouette feel stable.

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - Added a CPU-only test verifying determinism (same seed -> same hashed point cloud) and broad sanity bounds.

### Files changed/added

- `include/stellar/render/Nebula.h`
- `src/render/Nebula.cpp`
- `tests/test_nebula_field.cpp` *(new)*
- `PATCH_NOTES.md`

## Round 90 - Resonant Moon Systems: Mutual-Hill Spacing + Irregular Moons

This round deepens **procedural system generation** by making moon systems feel less like random scatter and more like *miniature dynamical architectures*.

### What’s new

- **Moon generator upgraded to v2** (does **not** perturb other system procedural streams)
  - Places moons using a **best-candidate** search in log-orbit space.
  - Enforces a stability-friendly spacing heuristic using **mutual Hill radius** separation.
  - Uses a more physics-inspired inner cutoff based on an estimated **Roche limit** + conservative radius buffer.

- **Resonant chains (gentle snapping pass)**
  - If adjacent moons are already near a simple rational period ratio (e.g., 2:1, 3:2, 4:3), the outer moon is nudged toward that resonance.
  - Produces more readable “Galilean-style” systems without hard-locking everything into perfect commensurability.

- **Irregular moons for gas giants**
  - Optional small, **high-inclination / higher-eccentricity** moons placed in the outer Hill-stable band.
  - Adds variety and “real solar system” flavor (including frequent retrograde captures) without any new simulation complexity.

- **Procedural System Lab UI**
  - Moon table gains:
    - **P ratio** (period ratio vs previous moon around the same parent), annotated when close to a resonance.
    - **dRH** (separation in **mutual Hill radii**) vs the previous moon.

### Files changed/added

- `src/proc/SystemGenerator.cpp`
- `apps/stellar_game/ProceduralSystemLabWindow.cpp`
- `PATCH_NOTES.md`

## Round 89 - Starfield 2.0: Spherical Blue-Noise Relaxation + Kelvin Color Temperature

- **Star placement** got a major quality upgrade.
  - Added 3 distribution modes:
    - `UniformRandom` (legacy)
    - `Fibonacci` (fast quasi-uniform lattice)
    - `RelaxedFibonacci` (default): runs a deterministic, repulsion-based relaxation pass using a coarse 3D grid to approximate **spherical blue-noise**.
  - Added `Starfield::Settings` controls for band shaping, random rotation, tangent jitter, and relaxation iterations/strength.
- **Star appearance** is more physically plausible:
  - Stars sample an approximate **blackbody temperature** (Kelvin) from a brightness-biased spectral-type mix (M..O).
  - Temperature is converted to approximate sRGB; a tiny dust-reddening term avoids an overly-clean palette.
- **Tests**: added `test_starfield` to sanity-check unit-length directions and basic nearest-neighbor uniformity (CPU-only; runs in headless builds).

### Files changed/added

- `include/stellar/render/Starfield.h`
- `src/render/Starfield.cpp`
- `tests/test_starfield.cpp` *(new)*
- `PATCH_NOTES.md`

## Round 88 - FrameGraph GraphViz DOT Export + Profiler UI

- Added `render::FrameGraph::toDot(...)`: export the compiled FrameGraph as a GraphViz DOT file for quick visual debugging in tools like GraphViz/xdot.
  - Optional clustering by transient **physical texture id** to make aliasing behavior visible.
  - Optional lifetime (firstUse..lastUse) and size annotations.
  - Optional invisible schedule edges to encourage a stable left-to-right layout.
- The **Profiler** window GPU FrameGraph panel gains a new **Export graph (GraphViz DOT)** section:
  - Copy DOT to clipboard
  - Write a `.dot` file to disk
  - Toggle what metadata is included (backbuffer/external/lifetimes/etc.)
- Extended `test_frame_graph` to validate DOT export emits the expected nodes/edges.

### Files changed/added

- `include/stellar/render/FrameGraph.h`
- `src/render/FrameGraph.cpp`
- `apps/stellar_game/ProfilerWindow.h`
- `apps/stellar_game/ProfilerWindow.cpp`
- `tests/test_frame_graph.cpp`
- `PATCH_NOTES.md`

## Round 87 - Cinematic Camera: Centripetal Catmull-Rom + Arc-Length Constant Speed

- Added `math::Spline` utilities: uniform + centripetal Catmull-Rom evaluation and a tiny arc-length table (s->u inversion) for constant-speed sampling.
- **Cinematic Camera** window now supports:
  - **Centripetal** mode (reduces overshoot/self-intersections around sharp corners)
  - **Constant speed (arc-length)** mode with an adjustable table resolution
  - lazy per-segment arc-length caching so per-frame sampling stays cheap
- Added `test_spline` to validate endpoint invariants and ensure arc-length reparameterization reduces speed variance.

### Files changed/added

- `include/stellar/math/Spline.h` *(new)*
- `apps/stellar_game/CinematicCameraWindow.h`
- `apps/stellar_game/CinematicCameraWindow.cpp`
- `tests/test_spline.cpp` *(new)*
- `PATCH_NOTES.md`

## Round 86 - BeatTracker: Spectral-Flux Onset Detection + BPM UI

- Added `dsp::BeatTracker`: a tiny, dependency-free beat/onset tracker built on spectral flux with an adaptive threshold + peak picking.
- The **Audio Analyzer** window gains a new **Beat / Onset** section:
  - live onset + threshold plots for tuning sensitivity
  - beat pulse (shader-friendly), BPM estimation, and confidence indicator
  - UI controls for tau/sensitivity/min interval
- Added `test_beat_tracker` with a deterministic synthetic spectrum stream to validate beat triggering + tempo convergence.

### Files changed/added

- `include/stellar/dsp/BeatTracker.h` *(new)*
- `src/dsp/BeatTracker.cpp` *(new)*
- `apps/stellar_game/AudioAnalyzerWindow.h`
- `apps/stellar_game/AudioAnalyzerWindow.cpp`
- `tests/test_beat_tracker.cpp` *(new)*
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 85 - glTF Export Power-Ups: Tangents + MSFT_lod Packed LOD Chains

- `render::exportMeshToGltf` can now optionally export per-vertex tangents (`TANGENT` / VEC4) for correct normal mapping in standard glTF viewers.
- Added `render::exportMeshLodsToGltf` which packs a LOD chain into a single `.gltf`/`.bin` using the `MSFT_lod` extension (instead of emitting one file per LOD).
- Procedural Mesh Lab export UI gains toggles for **Include tangents** and **Pack LODs into one glTF (MSFT_lod)**.
- Added `test_gltf_export` to validate tangent packing/layout and basic `MSFT_lod` JSON emission.

## Round 84 - Fence-based Async GPU Readback + Hitch-free Screenshots

- Upgraded `render::AsyncTextureReadback` to optionally use OpenGL sync objects (`glFenceSync`/`glClientWaitSync`) for non-blocking readiness checks, avoiding GPU stalls when mapping PBOs (fallback to frame-delay heuristics when unavailable).
- Extended readbacks to support `ReadbackSource::Framebuffer` using `GL_PIXEL_PACK_BUFFER` + `glReadPixels` (backbuffer/FBO capture, not just textures).
- Screenshot capture is now queued and written out asynchronously during the readback poll step (reduces frame hitching); clipboard copy and Photo Mode last-path update happen on completion.
- Added sync-object entry points to the OpenGL loader (`glFenceSync`, `glClientWaitSync`, `glDeleteSync`).

## Round 83 - Multi-thread Profiler + Multi-track Chrome Trace Export

- Made `core::Profiler` thread-safe so `STELLAR_PROFILE_SCOPE` can be used inside worker threads without data races.
- `ProfilerEvent` now stores a hashed `threadId` and `ProfilerFrame` stores `mainThreadId` (thread that began the frame) for better diagnostics.
- Chrome trace exporter now supports `splitThreads` (default ON): span events are emitted on separate tids per thread and named via `thread_name` metadata (Perfetto / chrome://tracing).
- Profiler UI adds a **Split threads** toggle in the Export Trace section.
- Added `test_profiler_mt` to validate multi-thread capture.

## Round 82 - Mission Board Batch Route Preview (single-source Dijkstra)

- Added `sim::NavRouteBatch` (single-source shortest paths) that computes best routes from one origin to all systems in a nearby set using the same implicit jump-range graph as A*.
- Mission Board route preview now uses `NavRouteBatch` to compute per-offer (via/to) jumps/distance/fuel with far fewer solves (at most a handful of batch passes vs N x A* per offer).
- Added `test_nav_route_batch` to validate the batch solver against existing A* expectations for hop/distance/fuel cost modes.

## Round 81 - FluidSim2D Multigrid Projection + Diagnostics

- Added a multigrid V-cycle accelerator for the pressure projection step in `proc::FluidSim2D` (still Stable Fluids-style boundaries). This dramatically improves convergence at higher grid sizes without cranking iteration counts.
- Added `FluidSim2DStats` (pressure residual + divergence metrics) and surfaced it in the Procedural Fluid Lab UI.
- Added projection controls to the lab (enable/disable multigrid + smoothing/coarse solve knobs).

# Patch Notes

## Round 80 — JobSystem TaskGroup + help-while-wait parallelFor (lower overhead, safer nesting)

The `core::JobSystem` is used all over the simulation/proc stack (hyperlane routing, scanners, route planning),
and it already provided a clean `submit()` + `parallelFor()` interface.

However, `parallelFor()` previously used `submit()` + `std::future` internally, which creates extra allocations
and bookkeeping overhead. More importantly, when `parallelFor()` is invoked from inside a worker job
(nested parallelism), **blocking waits** can become fragile unless the waiting thread can still help execute
queued work.

This round improves existing code by adding a lightweight **TaskGroup** + **help-while-wait** primitive and
rewriting `parallelFor()` to use it.

### What changed

- `JobSystem` gains:
  - `enqueue()` — fire-and-forget scheduling with no future
  - `tryRunOne()` — run a queued task inline (used for help-while-wait)
  - `TaskGroup` — submit N tasks and wait without futures
  - `isWorkerThread()` — tiny TLS-based helper for future oversubscription heuristics
- `parallelFor()` now:
  - keeps the user callable by reference (old semantics) while still supporting move-only functors
  - schedules helper loops through `TaskGroup` (avoids future allocations)
  - waits using help-while-wait semantics (safer nested usage)
- `test_jobs` expanded to cover `TaskGroup` and nested `parallelFor()` inside a job.

### Files changed/added

- `include/stellar/core/JobSystem.h`
- `tests/test_jobs.cpp`
- `PATCH_NOTES.md`

## Round 79 — FuzzySearch upgrade (DP matcher, CamelCase acronyms, better highlighting)

The UI uses fuzzy search in a lot of places (command palette, controls, logs, market dashboards, Comms inbox).
The previous implementation used a greedy subsequence scan per token, which was fast but often produced
**non-intuitive matches** (especially for CamelCase acronyms) because it didn't search for the *best* alignment.

This round improves existing code by replacing the matcher with a lightweight, deterministic **dynamic-programming**
scorer/backtracker (in the "fzy" family of fuzzy finders):

- Still **ASCII-only** (as intended for these UI labels), but now prefers:
  - **CamelCase anchors** (FlightRecorderWindow → `frw`),
  - **word-boundary anchors** (separators),
  - **tight consecutive runs** (better highlight + less surprising result ordering).
- Returns **best-path match positions** for more faithful highlighting.

### Files changed/added

- `src/ui/FuzzySearch.cpp`
- `tests/test_fuzzy_search.cpp`
- `PATCH_NOTES.md`

## Round 78 — FluidSim2D BFECC Dye Advection (MacCormack correction) + UI Controls

Semi-Lagrangian advection is extremely stable, but it's also **very diffusive**: dye/smoke loses crispness quickly.
This round improves existing code by upgrading the FluidSim2D dye transport to an optional
**BFECC / MacCormack-style corrected advection** pass with a monotonic clamp to avoid ringing/overshoot.

### What changed

- `FluidSim2DParams` gains:
  - `dyeAdvectionCorrection` (0..1 blend between plain semi-Lagrangian and full BFECC)
  - `dyeAdvectionClamp` to clamp corrected values to the source stencil min/max
- `FluidSim2D` adds `advectBFECC(...)` and uses it for dye channels during `step()`
  - forward advect writes directly to the destination buffer
  - backward advect uses a single scratch buffer to estimate and compensate truncation error
- Procedural Fluid Lab now exposes **Advection Quality** controls (BFECC amount + clamp toggle).
- `test_fluid_sim2d` updated to exercise the BFECC path and assert dye remains bounded.

### Files changed/added

- `include/stellar/proc/FluidSim2D.h`
- `src/proc/FluidSim2D.cpp`
- `apps/stellar_game/ProceduralFluidLabWindow.cpp`
- `tests/test_fluid_sim2d.cpp`
- `PATCH_NOTES.md`



## Round 77 — Save/Load Comms Inbox (Base64 tokenization + id-stable import)

This round improves existing code by making the new diegetic Comms system **persist across quicksave/quickload**,
without destabilizing the line-oriented save format.

### What changed

- Added `stellar::core::base64Encode/base64Decode` (RFC 4648, whitespace-tolerant decode) used for safe string tokens in SaveGame.
- SaveGame now stores `comms` with `comms_msg` lines (id, time, channel, link ids, unread/pinned, and base64 payloads for from/subject/body).
  - Uses a `~` sentinel for empty strings to keep parsing token-based.
- Game quicksave/load now exports/imports the comms log:
  - Uses new `CommsLog::replace(...)` to restore messages and recompute the internal id allocator.
  - Clears overlay queue/selection on load to avoid dangling references.
- Added deterministic tests:
  - `test_base64` round-trip + spot encodings
  - `test_savegame_comms` ensures comms messages survive save/load including newlines + markup

### Files changed/added

- `include/stellar/core/Base64.h` *(new)*
- `src/core/Base64.cpp` *(new)*
- `include/stellar/sim/SaveGame.h`
- `src/sim/SaveGame.cpp`
- `include/stellar/sim/Comms.h`
- `src/sim/Comms.cpp`
- `apps/stellar_game/main.cpp`
- `tests/test_base64.cpp` *(new)*
- `tests/test_savegame_comms.cpp` *(new)*
- `CMakeLists.txt`
- `PATCH_NOTES.md`


## Round 62 — Smuggling Scanner (black market route suggestions w/ risk + availability modeling)

This round adds a missing gameplay/tooling primitive: a **headless smuggling route scanner**.
The core sim already had deterministic contraband rules + black market pricing, but there was no reusable module that could answer:
"What illegal goods should I haul, where, and how risky is it?"

### What changed

- Added a new sim module `stellar::sim::SmugglingScanner`:
  - Buys legally at the origin station (official ask + origin fee).
  - Sells at the destination station via the black market (black-market bid, includes fence cut).
  - Estimates sting probability, fine size, and computes clean/stung/expected profit.
- Added three black market availability modes so UIs can choose the behavior they want:
  - **TodayOnly**: only return opportunities where the fence is available today.
  - **Expected**: scale score by `access01` (useful for "I can wait for a good day").
  - **Ignore**: treat as always available (useful for tooling/tests).
- Added a mean/variance option for scoring:
  - `ExpectedProfit - riskLambda * stdDev` (riskLambda is configurable).
- Added a deterministic unit test that searches for a valid legal→illegal commodity pair, forces a strongly profitable setup via inventory edits, then validates determinism.

### Files changed/added

- `include/stellar/sim/SmugglingScanner.h`
- `src/sim/SmugglingScanner.cpp`
- `tests/test_smuggling_scanner.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`


## Round 61 — Deterministic preview LOD sampler (jittered grid + unbiased stub reservoir)

The **Procedural Galaxy Lab** preview was still using a brute-force sector scan with an early-break at `maxStubs`.
That approach becomes *very* expensive at large view radii and also biases the preview (it fills from the first sectors in scan order).
This round adds a **deterministic, adaptive LOD sampler** that stays fast when zoomed out while keeping the preview stable under panning.

### What changed

- Replaced the brute-force `min..max` sector loops in `rebuildPreview()` with an **adaptive jittered-grid sampler**:
  - Picks an XY **stride** based on the view bounds + expected systems/sector, targeting ~O(`maxStubs`) sector generations.
  - Anchors sampling blocks to the **world grid** so results are stable as you pan.
  - Uses a deterministic hash (`SplitMix64`) of `(blockCoord, z, stride, pass, seed)` to choose a representative sector per block.
  - Skips blocks outside the view circle using a fast **AABB→circle min-distance** test.
  - Always samples a small neighborhood around the view center at full resolution, then refines with up to 2 additional passes (halving stride) if the sample is sparse.
- Added an **unbiased stub reservoir**:
  - Maintains a max-heap of the best `maxStubs` stubs keyed by a deterministic per-stub hash (avoids scan-order bias).
  - Materializes the final set sorted by key for stable ordering (helpful for deterministic coloring/selection).
- Added preview diagnostics:
  - `previewSectorsGenerated`, `previewCandidateStubs`, and `previewSectorStrideXY` displayed in the UI.

### Files changed/added

- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `PATCH_NOTES.md`


## Round 60 — Heatmap dirty-flag caching + optional animation (missing required glue)

The GPU galaxy heatmap is an **analytic** background pass, so rendering it every frame wastes GPU time
when the view/params are static. This round adds a **stable render-key hash** so the heatmap is only
re-rendered when an input that affects it changes (a classic “dirty flag” approach), with an optional
animation mode when you *do* want it to update continuously.

### What changed

- Added a stable **heatmap render key** (`heatmapRenderKey`) that hashes:
  - view (center/zoom/Z slice), heatmap controls (mode/exposure/contours/LOD), and galaxy parameters.
  - render target dimensions.
- `renderGalaxyHeatmap(...)` now **skips the offscreen render pass** when the key and dimensions match the
  last render (and animation is off), reusing the previous texture.
- Added UI controls:
  - **Animate**: forces per-frame re-render using `iTime`.
  - **Force Re-render** button: invalidates the cache and redraws next frame.
- Stored cache state in `ProceduralGalaxyLabWindowState` (`heatmapLastHash/LastW/LastH`).

### Files changed/added

- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `PATCH_NOTES.md`

## Round 59 — Unclamped internal ShaderToy uniforms + robust GL state restore (missing required glue)

This round plugs two subtle but important gaps in the **ShaderToy-based procedural rendering** stack:

1) **Engine-driven uniforms** (galaxy radii, world-space centers, etc.) were going through the *same* parameter pathway as UI-edited shader params, which meant they could be **accidentally clamped** to a tiny range (often `[0..1]`).
2) The GPU heatmap pass was doing minimal OpenGL state restoration, which could subtly **perturb later rendering** (especially texture bindings / clear color) on some drivers.

### What changed

- Added **`ShaderToyParamSet::setRawValue(...)`**:
  - Writes a parameter by name **without clamping** to `minValue/maxValue`.
  - Still normalizes unused components to zero based on the parameter type (scalar/vec2/vec3/color).
  - Intended for internal, engine-synced uniforms that naturally exceed UI ranges.

- Refactored parameter normalization:
  - Split “sanitize value shape” from “clamp to range”, so both clamped and unclamped updates share the same component rules.

- Procedural Galaxy Lab heatmap glue:
  - Switched the internal `g*` heatmap uniforms to use **`setRawValue`**.
  - Removed the old “`max==min` disables clamping” hack from the heatmap param schema.

- Added a small **RAII OpenGL state guard** around the heatmap offscreen pass:
  - Restores framebuffer, viewport, program, VAO, enable flags, **clear color**, active texture, and the **2D bindings for texture units 0..3**.
  - Uses the project’s SDL-loaded `gl::ActiveTexture` entry point (avoids calling `::glActiveTexture` directly, which is not exported on Windows’ legacy OpenGL import libs).

### Files changed/added

- `include/stellar/render/ShaderToyParams.h`
- `src/render/ShaderToyParams.cpp`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `PATCH_NOTES.md`

## Round 58 — Fix broken test harness integration + add CMake “test-lint” guardrails (missing required code)

This round targets an under-developed but *critical* area: **project correctness via tests**.

The bespoke test harness (in `tests/test_harness.h`) expects each `tests/test_*.cpp` translation unit to define a
single entry point function:

- `int test_xxx();`

The runner in `tests/main.cpp` provides the only `main()`.

Two tests in the repository (`test_galaxy_morphology.cpp`, `test_system_moons.cpp`) still defined their own `main()`,
which caused:

- **linker failures** (multiple definition of `main`), and
- an **incomplete registry** (generated registry referenced missing `test_*` symbols).

### What changed

- **Converted** `tests/test_galaxy_morphology.cpp` and `tests/test_system_moons.cpp` to the required convention:
  - `int test_galaxy_morphology()`
  - `int test_system_moons()`
  - Both now use the shared `CHECK(...)` macro and return `failures ? 1 : 0`.

- Added a lightweight **configure-time “test-lint” sanity check** to `tests/CMakeLists.txt`:
  - Fails CMake configuration if any `test_*.cpp` contains `int main(`.
  - Fails CMake configuration if any `test_*.cpp` is missing its required `int test_xxx(` entry point.

This prevents the same class of breakage from creeping back in later.

### Files changed/added

- `tests/CMakeLists.txt`
- `tests/test_galaxy_morphology.cpp`
- `tests/test_system_moons.cpp`
- `PATCH_NOTES.md`

## Round 57 — Heatmap now respects Center Z + Warp/Flare + hard disc thickness (missing shader glue)

This round is a **completeness / correctness** pass for the GPU galaxy heatmap:
previous rounds added shader-driven LOD and many morphology knobs, but the heatmap still behaved like a
*2D-only* approximation.

The CPU generator is truly 3D: it has a warped midplane, flared thickness, and the preview slice is centered
at **Center Z** with a **Z Half-Range** filter. The heatmap now follows those same semantics.

### What changed

- **Center Z is now wired into the heatmap**:
  - Added `gViewCenterZLy` uniform and use it to evaluate the Z-slice `[centerZ - zHalf, centerZ + zHalf]`.

- **Warp + flare are now reflected in the heatmap density**:
  - Added uniforms for warp/flare parameters (`gWarp*`, `gFlare*`).
  - Implemented shader helpers `thicknessHalf()` and `warpZLy()` mirroring the CPU morphology helpers.

- **Correct vertical slice integration with hard disc clipping**:
  - The CPU generator hard-clips the disc to `|z - warpZ| <= halfThickness` and applies an exponential vertical falloff inside it.
  - The heatmap now computes a **closed-form integral** of `exp(-|z-warpZ|/halfT)` over the intersection of:
    - the preview slice, and
    - the clipped disc slab.
  - This prevents the heatmap from over-estimating density when `Z Half-Range` exceeds the actual (flared) thickness.

- **Blob fields sampled at the correct Z**:
  - Cluster/void blob fields are now evaluated at `z = Center Z` (instead of implicitly `z = 0`), so moving the slice in Z remains coherent.

### Files changed/added

- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `PATCH_NOTES.md`

## Round 56 — Fix missing heatmap uniform wiring (no-clamp params) + ImGui legend ID collisions

This round is a **stability / completeness** pass: it generates the missing glue that was required for the
GPU galaxy heatmap (Rounds 54–55) to behave correctly.

### What changed

- **Fixed heatmap parameter clamping bug**:
  - The heatmap builds its uniform schema programmatically using `ShaderToyParamSet`.
  - `ShaderToyParamSet::setValue()` clamps values to `[minValue, maxValue]` when `max > min`.
  - The default `maxValue` is `1`, which unintentionally clamped *world-scale* uniforms (radii, cell sizes, etc.)
    to `1.0`, breaking the heatmap.
  - Heatmap params now explicitly disable clamping by setting `max==min` (see comment in code).

- **Fixed ImGui ID collisions** in the legend:
  - Several `ImGui::ColorButton()` calls reused the same `"##"` identifiers inside loops.
  - Added `ImGui::PushID/PopID` so every swatch has a unique ID.

### Files changed/added

- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `PATCH_NOTES.md`

## Round 55 — Derivative-driven dynamic LOD for procedural galaxy heatmap (shader AA + blob skipping)

This round targets the requested under-developed area: **procedural generation + custom rendering**, with a focus on
**dynamic LOD implemented inside shaders**.

The GPU galaxy heatmap introduced in Round 54 looked great, but at large zoom levels high-frequency FBM detail and thin
contour lines could shimmer/alias, and the expensive cluster/void blob scans were sometimes wasted on sub-pixel features.

### What changed

- Added **derivative-driven procedural LOD** helpers to the ShaderToy GLSL micro-library:
  - `featurePx()` estimates how large a feature is in screen pixels using `dFdx/dFdy`.
  - `fbm2_lod(...)` / `fbm3_lod(...)` automatically **band-limit** FBM octaves as they become sub-pixel.
  - Optional **energy preservation** keeps contrast more constant if desired.

- Upgraded the **Procedural Galaxy Lab heatmap shader**:
  - Spiral and density noise now use **`fbm2_lod`** when enabled.
  - Cluster/void blob fields can now **skip** the neighborhood scan when the blob radius is sub-pixel (big perf win when zoomed out).
  - Contour lines are now **`fwidth()` anti-aliased** so they don’t sparkle during pan/zoom.

- Added new **Heatmap LOD UI controls**:
  - Toggle Dynamic LOD.
  - Pixel thresholds (Lo/Hi) defining when details fade out.
  - “Energy Preserve” knob.
  - Toggle “Skip tiny blob scans”.

### Files changed/added

- `src/render/ShaderToy.cpp`
- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `PATCH_NOTES.md`

## Round 54 — GPU galaxy heatmap preview (analytic density field + contours)

This round targets the requested under-developed area: **procedural generation + custom rendering**.
The Procedural Galaxy Lab previously rendered only *discrete stubs* (points) via ImGui draw lists.
That made it hard to reason about the *continuous* density field you are tuning (spiral arms, noise, clusters/voids, bar/ring, etc.).

This patch adds a fast **GPU heatmap** rendered to an offscreen texture using the existing **ShaderToy** pipeline,
then composited *behind* the point preview.

### What changed

- Added a **GPU heatmap background** to the **Procedural Galaxy Lab**:
  - Offscreen **`render::RenderTarget2D`** sized to the preview canvas.
  - A dedicated **ShaderToy snippet** that evaluates an analytic density model and outputs a false-color heatmap.
  - Optional **contour lines** to make gradients/filaments easier to see at a glance.

- Added new UI controls:
  - Enable/disable heatmap.
  - Heatmap **mode**: Density / Spiral Arms / Clusters / Voids / Morphology.
  - Resolution (Full/Half/Quarter) to trade quality for performance.
  - Exposure + gamma + contour count/width.

- The heatmap uses a **GPU-side streaming-safe blob field** for cluster/void visualization.
  It is deterministic and visually correlated with the CPU generator, but not guaranteed to match it *exactly*.

### Files changed/added

- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `PATCH_NOTES.md`

## Round 53 — Procedural “cosmic voids” (streaming-safe negative-density bubbles)

This round targets the under-developed area requested: **procedural generation**.
The galaxy generator can already add *positive* structure (spiral arms, clusters, bar/ring, etc.).
This patch adds a new complementary tool: **large-scale void bubbles** that carve *negative space* and create
distinctive cavities/filaments without precomputing or storing global state.

### What changed

- Added **`proc::GalaxyVoids`**: a deterministic, streaming-safe coarse-cell field that produces smooth bubble
  influences (`void01` in 0..1) with jittered centers/radii.

- Integrated voids into **`proc::GalaxyGenerator`**:
  - New parameters in `proc::GalaxyParams`: `voidStrength`, `voidCellSizeLy`, `voidChancePerCell`, `voidRadiusLy`,
    `voidRadiusJitter`, `voidStrengthJitter`, `voidFalloffPower`.
  - Density modulation (when enabled):
    - `mean *= clamp01(1 - voidStrength * void01)`
  - Works in **both** generation modes:
    - legacy per-sector Poisson sampling
    - streaming-safe minimum-separation (blue-noise-ish) placement
  - Default settings keep `voidStrength=0`, so the legacy distribution and deterministic regression signatures remain unchanged.

- Upgraded the **Procedural Galaxy Lab**:
  - New “Cosmic Voids” parameter block.
  - New **Color: void** visualization mode + legend.
  - Hover tooltips show void influence and the strongest void bubble’s ID/radius.

### Tests

- Added `test_galaxy_voids` to validate determinism, value ranges, and “disabled” behavior.

### Files changed/added

- `include/stellar/proc/GalaxyVoids.h` *(new)*
- `src/proc/GalaxyVoids.cpp` *(new)*
- `include/stellar/proc/GalaxyGenerator.h`
- `src/proc/GalaxyGenerator.cpp`
- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `tests/test_galaxy_voids.cpp` *(new)*
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 52 — Deterministic convoy escort formations (headless Formation helpers)

This round focuses on an under-developed gameplay cue: **NPC group motion**.
Convoy escorts used to chase the exact same follow point, causing “stacking” and jittery convergences.
Now escorts keep a readable, deterministic **wing formation** around their assigned convoy.

### What changed

- Added **`sim::Formation`** helpers (headless, deterministic):
  - Formation patterns: **Trail**, **Wedge**, **Diamond**, **Ring**.
  - A small `Basis` builder (`makeBasisFromForward`) to generate a safe orthonormal frame.
  - Seeded per-slot jitter (optional) to avoid unnaturally perfect symmetry without introducing drift.

- Upgraded police **convoy escort follow** behaviour:
  - Escorts are assigned stable **slot indices** per `followId`.
  - Each escort now chases its own **formation anchor** around the followed contact.
  - The result is clearer convoy silhouettes and less "all NPCs collapse onto one point" behaviour.

### Tests

- Added `test_formation` to validate determinism, basic spacing sanity, and basis construction.

### Files changed/added

- `include/stellar/sim/Formation.h` *(new)*
- `apps/stellar_game/main.cpp`
- `tests/test_formation.cpp` *(new)*
- `PATCH_NOTES.md`

## Round 51 — Ship HUD nav markers + build fixes (MSVC)

This round targets the most under-developed area introduced recently: **the procedural Ship HUD**, along with a set of
blocking **MSVC build issues** that snuck in during the last patch.

### Added

- **Attitude instrument nav markers**:
  - **Retrograde marker** (opposite your current velocity) so you can quickly orient for braking / rendezvous.
  - **Nearest-contact bearing marker** (ship-local yaw/pitch to the closest contact) for fast situational awareness without a full target-lock system.

### Fixes

- Fixed MSVC error where `ShaderToyParamSet::applyToShader()` required a non-const `ShaderProgram&` even though uniforms are set through const methods.
- Removed accidental **duplicate definitions** in `ShaderToyGraph` (`userHeader_` member + `setUserHeader()` implementation).
- Added missing `RenderTarget2D::isInited()` accessor used by procedural render systems.
- Split the large `GpuSurfaceCache` GLSL fragment source into multiple literals to avoid MSVC’s per-literal size limit (`C2026: string too big`).
- Fixed Ship HUD telemetry compilation errors by using `Ship::positionKm()` / `Ship::velocityKmS()` getters consistently, and by referencing the correct fuel variables.

### Files changed/added

- `include/stellar/render/ShaderToyParams.h`
- `src/render/ShaderToyParams.cpp`
- `include/stellar/render/ShaderToyGraph.h`
- `src/render/ShaderToyGraph.cpp`
- `include/stellar/render/RenderTarget.h`
- `src/render/GpuSurfaceCache.cpp`
- `apps/stellar_game/ShipHudOverlay.h`
- `apps/stellar_game/ShipHudOverlay.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`



## Round 50 — Procedural cratered planets (seam-free impact synthesis + ejecta rays)

This round focuses on an under-developed procedural generation area: **planet surface detail**.
The existing noise-based albedo/height fields looked good at a distance, but lacked the kind of “readable” structure
(craters, rims, ejecta) that sells scale when you fly close.

### What changed

- Added a **multi-scale crater field** for **Rocky / Desert / Ice** surfaces:
  - Seam-free (evaluated on the **unit sphere**).
  - Uses a lightweight **Worley-style nearest-feature lookup** over a jittered integer grid.
  - Shapes each impact into a **bowl + raised rim + ejecta ring**.
  - Large impacts optionally include **directional ejecta rays** (streaks) in the outer ring.

- The crater field affects both:
  - **Albedo** (darkened bowls, brighter rims/ejecta)
  - **Height** (feeds directly into the procedural **normal map** pipeline)

- Implemented in **both** paths so visuals remain consistent:
  - CPU surface synthesis (`render::ProceduralPlanet`)
  - GPU surface cache synthesis (`render::GpuSurfaceCache` fragment shader)

### Files changed/added

- `src/render/ProceduralPlanet.cpp`
- `src/render/GpuSurfaceCache.cpp`
- `PATCH_NOTES.md`



## Round 49 — Shader parameters pipeline (DSL → uniforms → UI → serialization)

**Big focus:** making procedural shaders *controllable* (and shareable) by introducing a ShaderToy-style parameter pipeline.

### Added
- **Comment-driven parameter DSL** parsed from shader code across all passes:
  - `// @group <name>` ... `// @endgroup`
  - `// @param <type> <name> ["Label"] ...`
  - Types: `float`, `int`, `bool`, `vec2`, `vec3`, `color`
- **Auto-injected uniform block**: the graph now injects the parsed `uniform` declarations into every pass wrapper *before* the `#line 1` reset, keeping compile error line numbers aligned to user code.
- **Procedural Shader Lab → Parameters panel**:
  - live sliders / toggles / vector drags / color picker
  - filter box, per-group layout, “Reset to defaults”, and an “Advanced” mode (shows GLSL type + uniform name + copy uniform block).
- **.stoy serialization**: optional `PARAMS_BEGIN`/`PARAMS_END` block storing UI overrides, loaded back automatically.

### Tweaks
- Reaction Diffusion preset now exposes key values (F/K/du/dv/paint radius/edge gain/etc.) via the new `// @param` system.

### Fixes
- Fixed a broken `uiWindows.add(...)` block in `apps/stellar_game/main.cpp` that could prevent compilation.
- Added `HudWidgetId::Ship` to `include/stellar/ui/HudLayout.h` to match the existing Ship HUD layout logic.

## Round 48 — Procedural Shader Lab v2: multi-pass ShaderToy graph (Buffers A–D), channel routing, feedback & reaction-diffusion preset

This round pushes the most under-developed area: **procedural generation shaders**.
Instead of a single fullscreen snippet, the lab now supports a **ShaderToy-like multi-pass workflow** (Buffer A–D + Image)
so you can build feedback-driven simulations (reaction diffusion, trails, ping-pong filters, etc.) directly in-engine.

### What's new

- **Multi-pass ShaderToy graph runner (`ShaderToyGraph`)**
  - Renders **Buffer A → Buffer B → Buffer C → Buffer D → Image** every frame.
  - Each pass can sample up to **4 inputs** (`iChannel0..3`) routed from any buffer.
  - Feedback “just works”: when a buffer samples itself, it reads the **previous frame** (classic ping-pong).

- **ShaderToy wrapper upgraded**
  - Added common uniforms: `iTimeDelta`, `iFrame`, `iPass`.
  - Added channels: `iChannel0..3` + `iChannelResolution[4]`.
  - Null channels automatically bind a **1×1 black** texture so sampling is always defined.

- **Procedural Shader Lab UI overhaul**
  - Edit any pass (Image / Buffer A–D) with per-pass:
    - enable/disable (buffers)
    - resolution scaling (1× / 1⁄2× / 1⁄4×)
    - channel routing drop-downs for `iChannel0..3`
  - Preview any buffer output or the final Image output.
  - New sim controls: pause, single-step, fixed time-step, and buffer reset.

- **New multi-pass preset: Reaction Diffusion (Gray–Scott)**
  - Buffer A performs the simulation with self-feedback.
  - Image visualizes the evolving chemical fields.
  - Mouse painting injects perturbations into the sim.

- **Graph save/load format (`.stoy`)**
  - Save & load the full graph (all pass codes + routing + scales) as a single text file.
  - Live reload works for `.stoy` files too.

### Files touched

- `include/stellar/render/ShaderToy.h`
- `src/render/ShaderToy.cpp`
- `include/stellar/render/ShaderToyGraph.h` *(new)*
- `src/render/ShaderToyGraph.cpp` *(new)*
- `apps/stellar_game/ProceduralShaderLabWindow.h`
- `apps/stellar_game/ProceduralShaderLabWindow.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`

---

## Round 44 — Procedural Sky v2: seeded Milky Way glow, star clustering, blackbody colors + dithered nebula compositing

This round shifts focus to an under-developed part of *rendering*: **procedural background fidelity**.
The Procedural Sky shader is now closer to “space-sim-grade” — with a seeded Milky Way band, coherent star clustering,
more believable star colors, and less raymarch banding.

### What's new

- **Seeded Milky Way layer (optional)**
  - Adds a seam-free galactic band with **clumps** and **dust lanes**.
  - Galactic orientation is deterministic per sky seed (so it varies between systems while staying stable).

- **Star clustering + diffraction spikes**
  - A coherent cluster field biases star probability and brightness into **sparkly clumps** (mostly along the galactic plane).
  - Rare bright stars can render subtle **diffraction spikes**.

- **Blackbody-ish star colors**
  - Stars sample a temperature distribution (cooler stars are more common) and approximate blackbody RGB for richer variety.

- **Nebula raymarch improvements**
  - Adds per-pixel jitter to reduce banding at low step counts.
  - Uses cheap front-to-back compositing for more “volumetric” depth and to prevent over-bright fog.

- **VFX settings persistence + UI**
  - Added Milky Way + cluster parameters to VFX settings (now **v3**) and exposed them in the VFX UI.

### Files touched

- `include/stellar/render/ProceduralSky.h`
- `src/render/ProceduralSky.cpp`
- `include/stellar/ui/VfxSettings.h`
- `src/ui/VfxSettings.cpp`
- `apps/stellar_game/main.cpp`
- `tests/test_vfx_settings.cpp`
- `PATCH_NOTES.md`

---

## Round 43 — Procedural Ship HUD v4: segment-display readouts (14‑seg vector font) + HUD settings fidelity

This round continues to invest in the **procedural Ship HUD** — specifically the most under-developed part of the cockpit feel: **how values are presented**.

Instead of rendering the important readouts (speed, shield %, FSD state, etc.) as normal UI text, the HUD can now render them in a **procedural 14‑segment-style display** (vector strokes), with optional glitch‑driven segment dropouts.

### What's new

- **Segment-display value readouts (optional)**
  - New `Ship HUD → Segment readouts` toggle.
  - Renders value strings using a stylized **14‑segment / starburst** display approximation (digits + A‑Z + punctuation).
  - Includes a subtle **soft glow** pass to keep it readable against dark backgrounds.
  - When Ship HUD **Glitch FX** is enabled, the seg-display can deterministically **drop segments** per time-slice to read as “signal instability” (not random flicker).

- **HUD settings + command palette integration**
  - Added `shipHudSegmentText` to `HudSettings` (**v5**) and persisted it to `hud_settings.txt`.
  - Added a HUD settings UI checkbox and a command palette action: **Toggle Ship HUD segment readouts**.

- **Fix: HUD settings dirty tracking missed Ship HUD fields**
  - `hudSettingsEquivalent(...)` now correctly compares all Ship HUD fields (detail level, glitch/decor/glyphs/microtext, seed/nonce, and segment text).
  - This makes “Unsaved changes” detection and Save/Discard behavior consistent when tweaking the Ship HUD.

### Files touched

- `apps/stellar_game/ShipHudOverlay.h/.cpp`
- `apps/stellar_game/main.cpp`
- `include/stellar/ui/HudSettings.h` / `src/ui/HudSettings.cpp`
- `PATCH_NOTES.md`

---

## Round 42 — Procedural Ship HUD v3: squarified treemap layout + Attitude (navball-lite) instrument

This round keeps focus on the **procedural Ship HUD** and pushes it into a more pilotable, “ship-grade” cockpit surface:

- The **layout generator** has been upgraded from a greedy BSP split to a **squarified treemap** packer, producing
  more readable, squarer panels (especially important for circular instruments).
- A new **Attitude** instrument adds a “navball-lite” display: a horizon + pitch ladder + roll, with **prograde**
  (velocity) and **gravity** direction markers when available.

### What's new

- **Squarified treemap panel packing**
  - Replaced the previous recursive split layout with a deterministic squarified treemap algorithm.
  - Produces significantly better aspect ratios for panels without sacrificing determinism.
  - Adds subtle procedural variety by randomizing row placement and fill direction while staying seed-stable.

- **New instrument: `ShipHudInstrument::Attitude`**
  - Adds an attitude indicator inspired by real-world artificial horizons and spaceflight FDAIs:
    - Horizon line + pitch ladder (±10/20/30°)
    - Roll handling via horizon rotation
    - Prograde marker (ship velocity direction in ship-local angles)
    - Gravity marker (down vector in ship-local angles when gravity is present)
  - The Attitude panel is automatically included at **detail level ≥ 2**.

- **Smarter reference frame selection**
  - Near meaningful gravity: the indicator uses **gravity up** (−g) so the horizon is “local.”
  - In deep space: it falls back to an **orbit-plane up** reference (pos × vel) to avoid a misleading “star gravity horizon.”

### Files touched

- `include/stellar/ui/ProceduralShipHud.h` / `src/ui/ProceduralShipHud.cpp`
- `apps/stellar_game/ShipHudOverlay.h/.cpp`
- `apps/stellar_game/main.cpp`
- `tests/test_ship_hud_plan.cpp`

---

## Round 41 — Procedural Ship HUD v2: themed skins, glyph atlas, circuit decor, panel dropouts + overlay refactor

This round revisits an under-developed area: the **procedural Ship HUD**.
The previous version proved the concept (procedural layout + basic gauges), but it lacked the
**rich procedural “skin” layer** that makes a HUD feel like a manufactured instrument panel.

### What's new

- **Dedicated Ship HUD renderer module**
  - Added `apps/stellar_game/ShipHudOverlay.h/.cpp`.
  - `main.cpp` now computes telemetry once and hands rendering to `ui::drawShipHudOverlay(...)`.
  - The overlay contains a small ring-buffer helper (`ShipHudHistory`) so sparklines stay consistent.

- **Seed-derived HUD theme variant (`ProceduralShipHudPlan::themeVariant`)**
  - Each generated plan now includes a stable, deterministic theme selector.
  - Renderers can use this to vary line weight, corner rounding, tick density, and decoration style.

- **Procedural panel decor (opt-in)**
  - Deterministic circuit traces / hatch / radar / barcode motifs per panel.
  - Designed to be subtle (alpha-controlled) and scale with detail level.

- **Procedural glyph icons + microtext (opt-in)**
  - Vector glyph per instrument (speed, shield, heat, fuel, pips, fsd, throttle, target, cargo, g-force).
  - Small deterministic microtext “calibration/data” blocks for authenticity.

- **Improved gauge rendering**
  - Added tick marks (major/minor) to arcs/dials.
  - Added more consistent value formatting and “danger” behavior.

- **Panel-level glitch dropouts (opt-in)**
  - If glitch FX is enabled, individual panels can briefly lose signal and show procedural noise.
  - Dropouts are deterministic per time-slice to read as “signal loss” rather than random flicker.

- **HUD settings expanded + persisted (HudSettings v4)**
  - New toggles: decor, glyphs, microtext, panel dropouts.
  - New scalar: decor alpha.
  - All settings are saved/loaded in `HudSettings`.

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - Added `test_ship_hud_plan.cpp`: determinism + packing/overlap sanity checks for ship HUD plans.

## Round 39 — Procedural rings v2: ringlets, micro-divisions, spiral density waves, spokes + SystemLab preview

This round targets a previously **underdeveloped visual proc-gen module**: **planetary rings**.
The old generator produced nice broad bands, but it lacked the *fine structure* that makes rings feel alive.

### What's new

- **Richer, physically-inspired ring structure** in `generateRingTexture(...)`
  - Multi-scale radial banding (coarse + fine) with periodic warp (seam-free at u=0/1).
  - Dozens of **ringlets** (narrow high-density stripes), including occasional **arc-like segments** (azimuth-localized).
  - A mix of **major divisions** and many **micro-gaps** (thin carved lanes).
  - Stylized **spiral density waves** (integer azimuthal modes → seamless wrapping).
  - Anisotropic **wake-like streak noise** plus rare **spokes** (radial dark features).
  - Warm-dusty vs cold-icy palettes, with subtle radial tinting and ringlet highlights.

- **Procedural System Lab: Rings inspector + CPU preview**
  - New **Rings** panel shows per-planet predicted ring presence/variant (mirrors the logic in `main.cpp`).
  - Click a planet to preview its ring texture (CPU raster via ImDrawList, no GL texture needed).
  - Includes a radial mean-density plot and a few tweak knobs (chance/opacity, preview resolution, variant).

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - New `test_rings.cpp`: determinism, seam continuity, and basic alpha/color sanity checks.
  - `tests/CMakeLists.txt`: excludes `test_rings.cpp` when `STELLAR_ENABLE_RENDER` is disabled.

### Files changed

- `src/render/ProceduralRings.cpp`
- `apps/stellar_game/ProceduralSystemLabWindow.h`
- `apps/stellar_game/ProceduralSystemLabWindow.cpp`
- `tests/test_rings.cpp`
- `tests/CMakeLists.txt`

## Round 38 — Procedural resource fields v3: structural sub-features + density-aware blue-noise + heatmap

This round revisits the resource-field generator and pushes it into *believable belt structure* territory:
**gaps, clumps, streaks, and density-driven placement/yields** — while staying fully deterministic.

### What's new

- **Deterministic structural features per field** (derived from field id)
  - Torus belts: angular **gaps** (resonant “lanes”), multiple **hotspot clumps**, optional low-amplitude **spokes**.
  - Sheet fields: filament **streaks/streams**.
  - Cluster pockets: internal **hotspots/subclusters**.
  - Exposed via a new `ResourceFieldFeature` list on `ResourceFieldPlan`.

- **Density field + adaptive blue-noise placement**
  - New `sim::resourceFieldDensity01(...)` returns a stable density value in `[0..1]` for any point in a field.
  - The Mitchell-style best-candidate selection now uses **density-aware clearance**:
    - tighter packing in hotspots
    - looser packing in voids/gaps
  - `ResourceAsteroid` stores `density01`, and `baseUnits` scales with density (richer cores, poorer edges).

- **Procedural System Lab: stronger inspection tools**
  - Fixed star-class display (`Star::cls`).
  - Resource field table now shows a compact **structure summary**.
  - Selected field view lists features + shows mean asteroid density.
  - Scatter preview adds an optional **density heatmap** overlay + optional density-colored points.

- **API polish**
  - Added `resourceFieldLayoutName(...)` and `resourceFieldFeatureKindName(...)`.
  - Added `filterFeaturesForField(...)` helper.

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - `test_resource_field.cpp` now validates:
    - feature determinism + presence per field
    - asteroid `density01` range + consistency with `resourceFieldDensity01()`

### Files changed

- `include/stellar/sim/ResourceField.h`
- `src/sim/ResourceField.cpp`
- `apps/stellar_game/ProceduralSystemLabWindow.h`
- `apps/stellar_game/ProceduralSystemLabWindow.cpp`
- `tests/test_resource_field.cpp`
- `PATCH_NOTES.md`

## Round 37 — Procedural resource fields v2: orbital-plane belts + QMC blue-noise placement + inspector

This round targets an under-developed procedural loop: **mining sites / resource fields**.
The project already had deterministic field generation, but the distributions could look clumpy and the
field orientation was essentially arbitrary. This update makes fields feel more *astrophysical* and adds
a dedicated inspector so you can iterate on the generator quickly.

### What's new

- **Low-discrepancy sampling utility (`core::LowDiscrepancy`)**
  - Header-only Halton helpers used for procedural candidate sampling.

- **Resource field placement upgraded to deterministic “best-candidate” blue-noise**
  - Each asteroid picks the best of a fixed Halton candidate set (Mitchell-style) to maximize clearance.
  - Dramatically reduces visible clumping while keeping determinism and avoiding unbounded rejection loops.

- **Optional orbital-plane alignment for belts/sheets**
  - `sim::generateResourceFields(...)` now accepts a *preferred plane normal*.
  - `sim::generateSystemSignals(...)` computes the anchor station's orbit normal and passes it through,
    so belts/sheets align with the local “ecliptic”.

- **Procedural System Lab upgraded with a Signals/Resource Fields inspector**
  - Live generation controls (resource field count + optional signal categories).
  - Resource field table (kind/layout/radii/arc/richness/commodities/id).
  - Selected-field details: yield mix + a 2D asteroid scatter preview in field-local coordinates.
  - Clipboard export for quick debugging.

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - `test_resource_field.cpp` now also validates preferred-plane alignment.

### Files added/changed

- `include/stellar/core/LowDiscrepancy.h` *(new)*
- `include/stellar/sim/ResourceField.h`
- `src/sim/ResourceField.cpp`
- `src/sim/Signals.cpp`
- `apps/stellar_game/ProceduralSystemLabWindow.h`
- `apps/stellar_game/ProceduralSystemLabWindow.cpp`
- `apps/stellar_game/main.cpp`
- `tests/test_resource_field.cpp`

## Round 36 — Planet surface proc-gen v2: tectonic plates, craters, biomes, storms

This round targets the most under-developed part of the pipeline: **planetary surface generation**.
The old surfaces were intentionally “good enough” noise; this update makes planets read as *worlds*:
large-scale coherent structure (plates/continents), mid-scale features (mountain belts/cracks), and
micro detail (craters/speckle) — while staying **deterministic** and cache-friendly.

### What’s new

- **Tectonic plate field (spherical Voronoi)** for rocky / desert / ocean / ice
  - Deterministically generated plate sites on the unit sphere.
  - Plate boundaries produce mountain belts / ridges.
  - Domain warping breaks up hard Voronoi edges into organic coastlines and belts.

- **Crater field stamping** (rocky / desert / ice)
  - Deterministic crater catalogs per body (direction, angular radius, depth, rim).
  - Craters affect both **height** (for normals) and **albedo** (interior darkening + rim brightening).

- **Ocean worlds get real biomes**
  - Continents emerge from the plate field + macro noise.
  - A lightweight lat/elevation/moisture model assigns **desert / savanna / grass / forest / jungle / tundra / snow**.
  - Beaches and shallow-water tinting make coastlines pop.

- **Gas giants + clouds: storms & vortices**
  - Banding is now domain-warped (less “stripey”).
  - Deterministic vortices create large storm spots + swirling structure.

- **Major CPU-side performance win for normal maps**
  - Normal map generation now precomputes a height map once and then does finite differences.
  - This is ~4× cheaper for heavy surfaces (plates + crater loops), and reduces hitches when new bodies stream in.

- **GPU surface synthesis updated**
  - The `GpuSurfaceCache` fragment shader now mirrors the new plate/crater/biome/storm logic so the default
    GPU path gets the same visual upgrade.

### Files changed/added

- `src/render/ProceduralPlanet.cpp`
- `src/render/GpuSurfaceCache.cpp`

## Round 35 — Procedural moons + System Lab

This round fills a major “system-level proc-gen” gap by generating **planet-centric moons**
with reasonable orbital stability constraints (Hill sphere), and gives you a new in-game
lab window to inspect the generated architecture.

### What's new

- **Procedural moons in `proc::SystemGenerator`**
  - Each planet can spawn 0–N moons based on type/mass and available stable orbit space.
  - Moon orbits are constrained to a fraction of the host planet’s **Hill radius** and are
    spaced outward with log-ish separation.
  - **Determinism-preserving design:** moons use a *separate RNG stream per planet*
    seeded from `SystemStub.seed`, so existing downstream draws (stations, names, etc.)
    are not perturbed.

- **New sim primitive: `sim::Moon` + `StarSystem::moons`**
  - Lightweight procedural-only data: `id`, `name`, `type`, `massEarth`, `radiusEarth`,
    planet-centric `OrbitElements`, and `parentPlanetIndex`.

- **New helpers in `sim::Units`**
  - `moonPosKm(parentPlanet, moon, tDays)` and `moonVelKmS(...)` compose parent + moon
    states for rendering/physics.

- **World rendering now draws moons**
  - Moons get the same planet surface pipeline (procedural surface + normals + optional clouds +
    atmosphere rim/volumetric) with conservative “air retention” heuristics to avoid tiny moons
    looking like mini-Earths.
  - For now moons are **ambience-only** (not targetable/scan targets) to keep the gameplay/UI
    surface area stable while the new body class settles.

- **New window: Procedural System Lab**
  - Browse the **current system** or type any system id.
  - Inspect star parameters, planet table, and a dedicated moons table with Hill-radius fractions.

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - Added `test_system_moons.cpp` to validate determinism and basic orbital sanity.

### Files added/changed

- `include/stellar/sim/Celestial.h`
- `include/stellar/sim/System.h`
- `include/stellar/sim/Units.h`
- `src/proc/SystemGenerator.cpp`
- `apps/stellar_game/main.cpp`
- `apps/stellar_game/ProceduralSystemLabWindow.h` *(new)*
- `apps/stellar_game/ProceduralSystemLabWindow.cpp` *(new)*
- `CMakeLists.txt`
- `tests/test_system_moons.cpp` *(new)*

## Round 34 — Galaxy morphology field (bar + ring + warp + flare)

This round adds a **streaming-safe macro morphology field** that layers cleanly with existing inhomogeneity controls
(spiral arms, density noise, sparse clusters) and the min-separation (blue-noise-ish) placement option.

### What's new

- **New module: `proc::GalaxyMorphology`**
  - Deterministic sampling at any world position:
    - `bar01` (elliptical bar mask)
    - `ring01` (gaussian-ish inner ring mask)
    - `warpZLy` midplane offset (sinusoidal in azimuth with radial ramp)
    - local `thicknessHalfLy` with optional radial **flare**
    - combined `densityMul` (bar + ring)

- **`proc::GalaxyGenerator` integrates morphology**
  - Bar + ring now participate in the **inhomogeneous thinning path** and the **min-separation** placement path:
    - `mean *= (1 + barStrength*bar01) * (1 + ringStrength*ring01)`
  - Warp + flare adjust **disc bounds** and the base vertical falloff:
    - bounds: `abs(z - warpZ) <= halfThickness(r)`
    - vertical falloff is evaluated relative to the warped midplane
  - All new knobs are **disabled by default** to preserve deterministic regression signatures.

- **Procedural Galaxy Lab upgrades**
  - New “Galaxy Morphology” panel (bar, ring, warp, flare) with live regeneration.
  - Optional **Morph guides** overlay (bar ellipse + ring circle).
  - New color mode: **Color: morph** (visualizes bar+ring influence).
  - Tooltips show per-system morphology diagnostics (`bar01`, `ring01`, `warpZLy`, thickness, multiplier).

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - Added `test_galaxy_morphology.cpp` to validate determinism and sane ranges for masks/warp/flare.

### Files added/changed

- `include/stellar/proc/GalaxyMorphology.h` *(new)*
- `src/proc/GalaxyMorphology.cpp` *(new)*
- `include/stellar/proc/GalaxyGenerator.h`
- `src/proc/GalaxyGenerator.cpp`
- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `tests/test_galaxy_morphology.cpp` *(new)*
- `CMakeLists.txt`
- `PATCH_NOTES.md`


## Round 31 — Galaxy hazards / space weather (procedural generation)

This round adds a **galaxy-scale hazard field** (nebulae + storms + rift classification) and hooks it into hyperlane
generation so the lane network can naturally form safer corridors and avoid turbulent space.

### What's new

- **New module: `proc::GalaxyHazards`**
  - Deterministic 3D hazard sampling at any position:
    - `nebula01` (sensor occlusion / “cover”)
    - `storm01` (navigation disruption)
    - combined `hazard01` + `GalaxyHazardKind` classification
  - Optional **time drift** (linear domain drift in a seed-derived direction) for evolving “space weather”.
  - Optional **region-biased art-direction** so hazards feel consistent with `proc::GalaxyRegions`.

- **Hyperlanes can now be hazard-aware**
  - New `HyperlaneParams::hazards` section (opt-in):
    - biases MST topology cost using integrated hazard along candidate edges
    - increases lane `risk01`
    - reduces lane `bandwidth01`

- **Procedural Trade Systems Lab upgrades**
  - Added hazard controls under “Lane generation parameters”.
  - Selection panel shows hazard kind + scalar values when hazards are enabled.
  - Top trade corridors list now includes per-edge hazard estimates.

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - Added `test_galaxy_hazards.cpp` to validate determinism, ranges, time drift, and segment sampling.

### Files added/changed

- `include/stellar/proc/GalaxyHazards.h` *(new)*
- `src/proc/GalaxyHazards.cpp` *(new)*
- `include/stellar/proc/Hyperlanes.h`
- `src/proc/Hyperlanes.cpp`
- `apps/stellar_game/ProceduralTradeSystemsLabWindow.cpp`
- `tests/test_galaxy_hazards.cpp` *(new)*


## Round 30 — Procedural interstellar traffic field (trade flows over hyperlanes)

This round adds a **macro interstellar traffic estimator** that turns the procedural economy + hyperlane graph into a
usable “where is the galaxy busy?” signal.

The output is a deterministic **traffic value per system** and a **corridor flow value per hyperlane edge**, computed
by routing pairwise (or sampled) trade potentials along shortest paths.

### What's new

- **New module: `proc::TradeFlow`**
  - Computes a headless `TradeFlowNetwork`:
    - node `traffic` / `traffic01`
    - edge `flow` / `flow01`
  - Uses a simple gravity-style trade potential based on economic mass + commodity complementarity.
  - Routes flow along the hyperlane network using the existing `proc::HyperlaneRouter`.
  - Automatically switches to a **sampling mode** for larger neighborhoods (configurable sources/partners).

- **`proc::HyperlaneRouter` gains path reconstruction**
  - New `buildPathTo(targetId, outPath)` API to extract the chosen path after `compute()`.
  - This is used by the flow model to distribute traffic along edges.

- **Procedural Trade Systems Lab upgrades**
  - New “Interstellar traffic field” controls (sampling + gravity exponent).
  - Selection panel now shows:
    - the selected system’s traffic value
    - top traffic hubs
    - top trade corridors (with bandwidth/risk/distance info)

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - New unit test validates that a toy flow routes through the cheapest hyperlane path and produces expected
    node/edge normalization.
  - Fixes `test_trade_routes_hyperlanes.cpp` to match the registry convention (`int test_trade_routes_hyperlanes()`).

### Files added/changed

- `include/stellar/proc/TradeFlow.h` *(new)*
- `src/proc/TradeFlow.cpp` *(new)*
- `include/stellar/proc/HyperlaneRouter.h` *(path extraction)*
- `src/proc/HyperlaneRouter.cpp` *(path extraction)*
- `apps/stellar_game/ProceduralTradeSystemsLabWindow.h/.cpp` *(traffic UI + caching)*
- `tests/test_trade_flow.cpp` *(new)*
- `tests/test_trade_routes_hyperlanes.cpp` *(fix name + harness usage)*
- `CMakeLists.txt` *(register new proc source)*

## Round 29 — Hyperlane-aware procedural trade routing (graph + econ)

This round **connects the hyperlane layer to the procedural trade system** so macro routes can be constrained by connectivity and scored by an **effective travel cost along lanes** (optionally accounting for bandwidth and risk).

### What's new

- **New module: `proc::HyperlaneRouter`**
  - Deterministic single-source shortest-path (Dijkstra-style) routing over a `proc::HyperlaneNetwork`.
  - Produces per-destination metrics:
    - travel cost (used for route scoring)
    - physical lane distance (sum of edge lengths)
    - compounded risk
    - bottleneck bandwidth (min-edge)
    - hop count

- **`proc::TradeRoutes` hyperlane-aware overloads**
  - New overloads of `suggestTradeRoutes(...)` accept a `HyperlaneNetwork` + `HyperlaneTravelParams`.
  - Routes that are **unreachable** in the lane graph are dropped.
  - `TradeRoute` now also carries direct distance + lane metrics so UI/debug can show “cost vs lane vs direct”.

- **`proc::TradeIntel` hyperlane-aware overload**
  - New `buildTradeIntel(...)` overload that computes route economics on hyperlane routes.
  - Route economics records mirror the new lane metrics.

- **Procedural Trade Systems Lab**
  - Added a toggle to compute macro trade routes **via hyperlanes**.
  - Added lane generator controls (neighbor radius/K, min degree, extra edges, region cell size).
  - Added lane travel-cost controls (risk weight, bandwidth bias, min bandwidth factor).
  - Route bullets now show cost/lane/direct/hops/bottleneck when hyperlane routing is enabled.

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - Added a unit test that validates hyperlane connectivity can change the best route (unreachable systems are excluded).

### Files added/changed

- `include/stellar/proc/HyperlaneRouter.h` *(new)*
- `src/proc/HyperlaneRouter.cpp` *(new)*
- `include/stellar/proc/TradeRoutes.h` *(extend TradeRoute + new overloads)*
- `src/proc/TradeRoutes.cpp` *(hyperlane routing integration)*
- `include/stellar/proc/TradeIntel.h` *(hyperlane-aware buildTradeIntel + extra metrics)*
- `src/proc/TradeIntel.cpp` *(hyperlane-aware buildTradeIntel implementation)*
- `apps/stellar_game/ProceduralTradeSystemsLabWindow.h/.cpp` *(hyperlane toggle + lane params + richer route display)*
- `tests/test_trade_routes_hyperlanes.cpp` *(new)*
- `CMakeLists.txt` *(register new proc source)*

## Round 28 - Procedural Hyperlane Network
- Added **proc::Hyperlanes**: a deterministic, sparse **hyperlane/trade corridor** graph built from a local set of system stubs (k-NN candidates -> MST backbone -> min-degree augmentation -> hash-based extra edges).
- Wired hyperlane visualization into **Procedural Galaxy Lab** (toggle + tunable parameters, bandwidth/risk color coding, tooltip degree stats).
- Added **test_hyperlanes** for connectivity/determinism regression coverage.


## Round 27 — Procedural galaxy regions (Worley/Voronoi biomes)

This round adds a **galaxy-scale region layer** based on a deterministic 3D cellular (Worley/Voronoi) partition.
It’s a foundational procedural system you can later plug into: economy, security, mission generation, and rendering.

### What's new

- **New module: `proc::GalaxyRegions`**
  - Constant-time sampling (fixed 3×3×3 neighbor checks) that returns:
    - stable region id + seed
    - stylized region kind (Core / Inner Disc / Outer Rim / Nebula / Cluster / Rift)
    - deterministic region name
    - an **edge factor** (`edge01`) using Worley f2–f1 logic (useful for boundaries)

- **Procedural Galaxy Lab upgraded**
  - New **Color by region** toggle + legend.
  - Hover tooltip shows region name/kind/edge.
  - Region cell-size slider for quickly exploring scale.

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - New unit test for determinism + basic spatial coherence.

### Files added/changed

- `include/stellar/proc/GalaxyRegions.h` *(new)*
- `src/proc/GalaxyRegions.cpp` *(new)*
- `apps/stellar_game/ProceduralGalaxyLabWindow.h/.cpp` *(region visualization + tooltip)*
- `tests/test_galaxy_regions.cpp` *(new)*
- `CMakeLists.txt` *(register new proc source)*

## Round 26 — Black market fences for contraband

This round tackles an underdeveloped loop called out in the README: **crime + contraband**.
Contraband used to be mostly a punishment (scan → confiscation → fine). Now it can be an *economy*.

### What's new

- **New module: `sim::BlackMarket` (headless)**
  - Deterministic, per-station black market profile:
    - whether a fence is available today
    - access/risk metrics
    - pricing multipliers (bid/ask) + fence cut
    - sting chance (with optional heat modulation)
  - Utility to apply those multipliers to an official `econ::MarketQuote`.
  - A first gameplay transaction primitive: `sellToBlackMarket(...)`
    - If stung, it reuses the existing **PoliceScan** contraband enforcement math to confiscate + fine.

- **`stellar_sandbox` upgraded**
  - New `--blackMarket` mode prints the fence profile + sample illegal-commodity quotes.
  - `--law` output now correctly prints the **station-local** illegal list (not just faction-wide).

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - New unit test validates determinism, bounds, and basic sell vs sting behavior.

### Files added/changed

- `include/stellar/sim/BlackMarket.h` *(new)*
- `src/sim/BlackMarket.cpp` *(new)*
- `apps/stellar_sandbox/main.cpp` *(black market tooling)*
- `tests/test_black_market.cpp` *(new)*
- `CMakeLists.txt` *(register new sim source)*

## Round 25 — Procedural macro market pricing + trade intel

This round turns the existing **TradeProfile + TradeRoutes** macro signal into *actionable economic intel*.
Everything here remains purely procedural — it does **not** query live station inventories/prices — but it now
produces plausible per-system price surfaces and rough profit estimates for suggested routes.

### What's new

- **New module: `proc::TradePrices`**
  - Deterministic macro **price multipliers** per commodity derived from a system's TradeProfile.
  - Export-heavy goods become cheaper; import-heavy goods become more expensive.
  - Hub/population dampen extremes to keep major markets more stable.

- **New module: `proc::TradeIntel`**
  - Enriches `proc::TradeRoutes` outputs with estimated buy/sell prices and profit numbers.
  - Provides simple **round-trip (2-leg) loop** detection: origin → other → origin.

- **Procedural Trade Systems lab window upgraded again**
  - Adds optional macro economics display (spread, fees, cargo capacity knobs).
  - Shows estimated per-route profit and the best round-trip loops found in the current route set.

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - Adds a focused unit test covering multiplier behavior, determinism, and basic loop profitability.

### Files added/changed

- `include/stellar/proc/TradePrices.h` *(new)*
- `src/proc/TradePrices.cpp` *(new)*
- `include/stellar/proc/TradeIntel.h` *(new)*
- `src/proc/TradeIntel.cpp` *(new)*
- `apps/stellar_game/ProceduralTradeSystemsLabWindow.h/.cpp` *(macro econ + loop display)*
- `tests/test_trade_prices.cpp` *(new)*
- `CMakeLists.txt` *(register new proc sources)*

## Round 24 — Procedural macro trade routes (profile-driven lane signal)

This round extends the "procedural trade profiles" work into a *galaxy-level* trade-lane signal you can query cheaply
without touching live station inventories/prices.

### What's new

- **New module: `proc::TradeRoutes`**
  - Generates deterministic **export** and **import** route suggestions from a system to nearby candidate systems.
  - Uses `proc::TradeProfile` export/import scores, hub/population/wealth, and distance falloff to estimate a
    "trade gravity" score.
  - Returns a compact list of routes (direction + best-fit commodity + score + distance + risk proxy).

- **Procedural Trade Systems lab window upgraded**
  - Click any system row to inspect its top exports/imports and **macro trade routes**.
  - Adds basic tuning knobs (max routes, distance exponent, optional max distance).

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - Adds a focused unit test validating basic ordering/selection logic and max-route truncation.

### Files added/changed

- `include/stellar/proc/TradeRoutes.h` *(new)*
- `src/proc/TradeRoutes.cpp` *(new)*
- `apps/stellar_game/ProceduralTradeSystemsLabWindow.h/.cpp` *(selection + macro route display)*
- `tests/test_trade_routes.cpp` *(new)*
- `CMakeLists.txt` *(register new proc source)*

## Round 23 — TradeProfile-driven station/economy generation

This patch makes the procedural TradeProfile *actually shape* generated systems and station economies.

### What changed

- System generation now derives station archetype weights from the system's TradeProfile.
  Mining, agriculture, industrial, and tech biases are reflected in the station mix.

- Station economy models are tuned per TradeProfile:
  - export-favored commodities get higher production and larger stock/capacity
  - import-favored commodities get higher consumption and lower stock

- Station fees and local volatility incorporate TradeProfile hub/lawlessness.

- Signature/versioning fixes:
  - Universe signature now includes stationCount properly.
  - test_signature updated for new signature inputs.

### Files changed

- include/stellar/proc/TradeProfile.h
- src/proc/TradeProfile.cpp
- include/stellar/proc/TradeEconomy.h
- src/proc/TradeEconomy.cpp
- src/proc/SystemGenerator.cpp
- src/proc/GalaxyGenerator.cpp
- include/stellar/sim/Signature.h
- src/sim/Signature.cpp
- tests/test_signature.cpp
- apps/stellar_game/ProceduralTradeSystemsLabWindow.cpp


## Round 32 — Streaming-safe min-separation star placement (blue-noise-ish)

This round focuses on the *most under-developed* part of the procedural pipeline: **controlling spatial system spacing**.

The legacy in-sector Poisson sampling can produce very tight clumps (or near-collisions) that look noisy and can make navigation/visualization harder.

### What changed

- Added an optional `GalaxyParams::minSystemSeparationLy` parameter.
  - `0` keeps the legacy behavior (so existing deterministic regression signatures remain stable).
  - `> 0` enables a deterministic, streaming-safe "cell + priority" placement scheme that enforces an approximate minimum distance between systems.

- The new placement method is coherent across *sector boundaries* without requiring any multi-sector state.

- Updated Procedural Galaxy Lab UI so you can tweak the parameter interactively.

- Added `test_galaxy_minsep` regression test to ensure:
  - min-separation is respected
  - generation is deterministic

### Files changed/added

- `include/stellar/proc/GalaxyGenerator.h`
- `src/proc/GalaxyGenerator.cpp`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `tests/test_galaxy_minsep.cpp` *(new)*


## Round 33 - Streaming Procedural Star Clusters

- Added **`proc::GalaxyClusters`**: a streaming-safe, deterministic "cluster influence" field (coarse 3D cell grid with jittered centers + smooth radial falloff).
- Integrated optional cluster density into **`proc::GalaxyGenerator`** (works with both the min-separation sampler and the inhomogeneous thinning path).
  - No behavior change when `clusterStrength == 0` (legacy distribution path remains unchanged).
- Expanded **Procedural Galaxy Lab** with:
  - Cluster generation controls (strength, cell size, chance/cell, radius, jitter, falloff)
  - "Color: cluster" mode + legend + tooltip cluster details
- Added `tests/test_galaxy_clusters.cpp` deterministic regression test.


## Round 45 — Procedural night-side city lights (emissive planet pass)

This round focuses on an under-developed rendering cue: **making the night side of planets feel alive**.

### What changed

- Added a new procedural surface kind: **`SurfaceKind::CityLights`**.
  - CPU path: `ProceduralPlanet::generateSurfaceTexture()` now supports it.
  - GPU path: `GpuSurfaceCache` shader now supports it.
  - The lights distribution is *sparse*, land-biased, coastline-biased, and includes a mix of clusters + filamentary “corridors”.

- Extended **`MeshRenderer`** with an **optional emissive texture**:
  - Emissive is additive on top of lighting.
  - A configurable terminator fade reveals emissive mostly on the night side using `smoothstep(start,end,-dot(N,L))`.

- Integrated city lights into the main world render path:
  - Deterministic per-planet probability (seeded), scaled by surface type.
  - Added UI controls in **World visuals**:
    - enable/disable
    - intensity
    - chance
    - terminator fade start/end

- Added city lights thumbnails to the **Surface generator preview** panel.

### Files changed

- include/stellar/render/ProceduralPlanet.h
- src/render/ProceduralPlanet.cpp
- src/render/GpuSurfaceCache.cpp
- include/stellar/render/MeshRenderer.h
- src/render/MeshRenderer.cpp
- apps/stellar_game/main.cpp
- src/ui/VfxSettings.cpp


## Round 46 — Procedural star corona (animated emissive shell)

This round tackles an underdeveloped rendering cue: **making the system star feel alive**.

### What changed

- Added a new renderer: **`render::StarCoronaRenderer`**.
  - Draws an **additive emissive shell** (usually a slightly scaled UV sphere).
  - Uses **seam-free 3D FBM noise** sampled on the unit sphere to break up the limb glow.
  - Adds animated **streamer rays** and **prominence lobes** aligned to a deterministic per-star axis.
  - Tuned for HDR/bloom friendliness (cheap but cinematic).

- Integrated into the world forward pass:
  - Drawn right after the star surface.
  - Controlled via the **World Visuals** window (enable, intensity, shell scale, rim power, noise, animation speed, streamer + prominence shaping, re-roll).
  - Treated as optional: if the shader fails to compile, the game logs a warning and continues.

### Files changed/added

- `include/stellar/render/StarCoronaRenderer.h` *(new)*
- `src/render/StarCoronaRenderer.cpp` *(new)*
- `CMakeLists.txt`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`


## Round 47 — Procedural Shader Lab (live-edit GLSL sketches)

This round focuses on an under-developed part of procedural generation: **iterating on shaders fast**.
Instead of baking everything into C++ strings or rebuilding for each tweak, you now get a **ShaderToy-style**
workflow directly inside the game.

### What changed

- Added **`render::ShaderToy`** (a tiny fullscreen shader runner):
  - Wraps a user snippet into a valid **GLSL 330** fragment shader.
  - Injects **common procedural uniforms** (`iResolution`, `iTime`, `iMouse`, `iSeed`) plus an optional **camera** (`iCam*`).
  - Injects a small, reusable GLSL “micro-library”:
    - integer-hash based value noise (2D/3D)
    - `fbm2/fbm3`, `worley2`, `warp2`
    - `palette()` + a couple SDF helpers (`sdSphere/sdBox`)
    - `rayDirFromUv()` + a simple tonemap
  - Uses `#line 1` before the snippet so driver error logs point at your **snippet line numbers**.

- Added a new in-game window: **Procedural Shader Lab**.
  - Live snippet editor + **compile** button + error log.
  - Built-in presets:
    - Template
    - Domain Warped Nebula
    - Voronoi Circuits
    - Raymarch Tunnel
  - Renders to a preview `RenderTarget2D` and displays it in ImGui.
  - Supports **mouse interaction** in the preview (fills `iMouse`) and a simple **orbit camera** for 3D sketches.
  - File workflow:
    - load/save shader snippet to a path
    - optional “live reload” (polls `last_write_time`)

### Files changed/added

- `include/stellar/render/ShaderToy.h` *(new)*
- `src/render/ShaderToy.cpp` *(new)*
- `apps/stellar_game/ProceduralShaderLabWindow.h` *(new)*
- `apps/stellar_game/ProceduralShaderLabWindow.cpp` *(new)*
- `apps/stellar_game/main.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 62 — Zero-allocation ShaderToy param updates + Heatmap param handles

This round focuses on a subtle but important under-developed part of the procedural pipeline: **parameter plumbing**.
The galaxy heatmap (and ShaderToy-driven tools in general) update dozens of uniforms per redraw; doing that through
name-based lookups can quietly create per-frame overhead and makes schema drift harder to notice.

### What changed

- **`render::ShaderToyParamSet` now supports C++20 heterogeneous lookup** (transparent `std::string_view` key support)
  so `findIndex()` no longer needs to allocate temporary `std::string` objects.
- Added **index-based setters**:
  - `setValue(int index, ...)` (clamped)
  - `setRawValue(int index, ...)` (unclamped, but still sanitizes unused lanes)
- **Procedural Galaxy Lab heatmap** caches all internal uniform indices once at init-time and then updates by index.
  - Faster per-frame updates.
  - Fails fast with a clear error if the internal schema ever drifts.
- Added a small **unit test** covering the new by-index setters and the clamped vs. raw behavior.

### Files changed/added

- `include/stellar/render/ShaderToyParams.h`
- `src/render/ShaderToyParams.cpp`
- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `tests/test_shader_toy_params.cpp` *(new)*
- `PATCH_NOTES.md`

## Round 63 — Headless-friendly ShaderToyParams module (fixes tests & unlocks offline tooling)

This round addresses a **missing integration piece**: `ShaderToyParams` is largely CPU-only (parsing, schema, clamping),
but it previously lived behind the OpenGL build flag. That meant **headless builds (no OpenGL)** could not link code that
uses `ShaderToyParamSet` (including the new unit test), even though those parts don’t need a renderer.

### What changed

- **`src/render/ShaderToyParams.cpp` is now compiled in headless builds** (moved into the always-on core source list).
- Added **feature macros** exported from CMake:
  - `STELLAR_ENABLE_RENDER` (0/1)
  - `STELLAR_ENABLE_IMGUI` (0/1)
  This enables clean, compile-time gating for optional rendering code.
- `ShaderToyParamSet::applyToShader()` is now **a no-op in headless builds**, avoiding any accidental link dependency on the
  OpenGL shader backend when it isn’t present.

### Why it matters

- Fixes **link failures** when building tests without OpenGL.
- Unlocks future **offline / CI tooling** that can parse `.stoy` graphs and validate parameter schemas in a headless build.

### Files changed/added

- `CMakeLists.txt`
- `src/render/ShaderToyParams.cpp`
- `PATCH_NOTES.md`

## Round 64 — Streaming-friendly preview sector cache (LRU) + adjustable preview cap

This round focuses on one of the most under-developed parts of procedural exploration: **interactive pan/zoom performance**.
Previously, the Galaxy Lab preview regenerated every visited sector each time you nudged the view. That’s correct but wasteful
because sector generation is deterministic.

### What changed

- Added a tiny, header-only `core::LruCache` utility (hash map + recency list) suitable for deterministic procedural caches.
- **Procedural Galaxy Lab** now caches generated sector stubs in an **LRU sector cache**:
  - Keys: sector coordinates `(x,y,z)`.
  - Invalidates automatically when the generator *context* changes (`seed`, `GalaxyParams`, or `factionCount`).
  - Capacity is user-controllable, with a one-click **Clear Cache** button.
  - Diagnostics are exposed: hits/misses/evictions and an approximate cached-stub memory estimate.
- Added a missing usability control: **Max Preview Systems** is now adjustable from the UI.
- Added a unit test for the new LRU cache.

### Files changed/added

- `include/stellar/core/LruCache.h` *(new)*
- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `tests/test_lru_cache.cpp` *(new)*
- `PATCH_NOTES.md`

## Round 65 — Parallel preview generation + cache-friendly insertion API

This round targets a common pain point when exploring big galaxies: **preview regeneration stalls** when the view covers lots of
sectors (large radius or high Max Preview Systems).

### What changed

- **Procedural Galaxy Lab** preview generation is now **parallelized on the CPU** using the existing `core::JobSystem`:
  - Sampled sectors are partitioned into chunks and processed concurrently.
  - Each worker keeps a local bottom‑K sample; the main thread merges those into the global bottom‑K.
  - The selection remains **deterministic** (stable hash key per system), regardless of worker count.
- The preview **sector cache** now stores **`shared_ptr` to immutable stub vectors**, so cached sectors can be safely shared
  across workers without copying.
- `core::LruCache` gained a missing but extremely useful API for expensive values:
  - `getOrInsert(key, value)` / `getOrInsert(key, value, onEvict)` — insert a precomputed value while preserving LRU behavior.
  - `insertOrAssign(key, value)` — replace an existing entry (touching MRU).
- UI additions under **Preview Performance**:
  - **Parallel Preview** toggle.
  - **Threads (0=auto)** control.
  - A small diagnostic showing how many worker threads were used in the last rebuild.

### Files changed/added

- `include/stellar/core/LruCache.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `tests/test_lru_cache.cpp`
- `PATCH_NOTES.md`

## Round 66 — Hyperlane overlay in Galaxy Lab + dynamic LOD line rendering

This round exposes one of the most “hidden” procedural systems in the repo: **hyperlane generation**. The procedural module
(`proc::generateHyperlaneNetwork`) already existed and was tested, but it wasn’t wired into the main exploration tooling.

### What changed

- **Procedural Galaxy Lab** gained a new **Hyperlanes (Navigation Graph)** overlay:
  - Generates a deterministic, sparse lane network (**MST + kNN + extra edges**) over the current preview slice.
  - Optional **hazard modulation** controls (nebula/storm fields) are now tweakable from the UI.
  - Edges can be colored by **risk** (default) or by **bandwidth**.
- Added **dynamic LOD for lane drawing** to keep the map readable at all zoom levels:
  - Caps rendered edges (manual or auto based on canvas size).
  - Skips edges that become **sub‑pixel** via a minimum screen‑length threshold.
  - Thickness and opacity scale with bandwidth.
- Added **caching + deterministic node subset selection**:
  - Hyperlanes are computed on at most **Max Nodes** (bottom‑K by stable hash key) so generation cost is bounded.
  - The overlay rebuilds only when its **input hash** changes (seed, preview stubset hash, hyperlane params, node cap).
- Added a small usability improvement: when hovering a system, the tooltip shows **hyperlane degree** and average bandwidth/risk,
  and (optionally) highlights incident edges.

### Files changed/added

- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `PATCH_NOTES.md`

## Round 67 - Hyperlane Route Inspector (interactive pathfinding + overlay)

### What changed
- **Procedural Galaxy Lab:** added a hyperlane **Route Inspector** that turns the generated lane graph into an interactive navigation tool:
  - **Shift+Click** a system to set **route start**.
  - **Ctrl+Click** a system to set **route end**.
  - **Right-click** a system to **clear** the current route.
  - Uses the existing **`proc::HyperlaneRouter` (Dijkstra)** to compute lowest-cost paths over the cached hyperlane subset.
  - Adds controls for **riskWeight / bandwidthBias / minBandwidthFactor**, plus overlay **opacity/width**.
  - Draws a **highlighted route polyline** (risk-colored segments, bandwidth-weighted width) on top of the hyperlane LOD pass.
  - Exposes **route metrics** (cost, distance, compound risk, bottleneck bandwidth, hops) + **copy-to-clipboard** route export.
- **Tests:** added Catch2 coverage for `HyperlaneRouter` (path choice + metrics composition).

### Files changed/added
- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `tests/test_hyperlane_router.cpp`
- `PATCH_NOTES.md`

## Round 68 — K-shortest Hyperlane Routing (Yen + A*), multi-route overlay, and tests

### What changed
- **Core / proc:** added `proc::plotKHyperlaneRoutesAStarCost(...)` — a **loopless K-shortest path enumerator** for hyperlane graphs:
  - Uses **A*** point-to-point search with a straight-line heuristic (admissible because hyperlane cost is always ≥ Euclidean distance).
  - Wraps that solver with **Yen’s algorithm** to enumerate **K alternative routes** (ordered by total travel cost, deterministic tie-breaks).
  - Returns full `HyperlanePathMetrics` for each alternative (cost, distance, compound risk, bottleneck bandwidth, hops).

- **Procedural Galaxy Lab:** upgraded the Hyperlane Route Inspector to support **multiple alternative routes**:
  - New **K Routes** slider.
  - A selectable table of route candidates (hops / cost / dist / risk / BW).
  - Optional **Draw alternatives** mode: renders all K routes with adjustable opacity/width, while the selected route stays bold on top.

- **Tests:** added a deterministic unit test that validates ordering + metrics for a small hand-built hyperlane graph.

### Files changed/added
- `include/stellar/proc/HyperlaneKRoutes.h` *(new)*
- `src/proc/HyperlaneKRoutes.cpp` *(new)*
- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `tests/test_hyperlane_k_routes.cpp` *(new)*
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 69 — Hyperlane Traffic Heatmap (Betweenness Centrality) + Chokepoint Overlay

This round turns hyperlanes into something you can *reason* about: a lightweight **graph analytics layer** that estimates which
lanes are likely to be “busy corridors” (and which systems are critical connectors), then exposes that in the Galaxy Lab.

### What changed

- **Core / proc:** added `proc::estimateHyperlaneBetweennessCentrality(...)`:
  - Computes **edge + node betweenness centrality** over the hyperlane graph using a Brandes-style accumulation step.
  - Uses **weighted Dijkstra** per source (weights match `HyperlaneRouter` travel cost), so “traffic” reflects **fastest routes**.
  - Supports **deterministic approximation via source sampling** (`sampleSources`) to keep it interactive on large previews.

- **Procedural Galaxy Lab:** Hyperlane overlay gained a new **Color Mode**:
  - **Risk** (existing)
  - **Bandwidth** (existing)
  - **Traffic (betweenness)** *(new)*
    - Adjustable **Traffic Samples (0=exact)**.
    - Optional **Highlight chokepoints** (top-N highest betweenness edges are redrawn thicker on top).
    - Extra **Traffic Width Boost** and **Traffic Alpha Boost** controls.
    - Displays a small cache stat (compute time + max centrality).

- **Route overlay integration:** when Traffic mode is active, the route polylines also use the same traffic coloring
  (still bandwidth-whitened) so you can see whether the chosen route uses major corridors or quieter backways.

- **Tests:** added a deterministic unit test for betweenness on a simple 4-node chain.

### Files changed/added

- `include/stellar/proc/HyperlaneCentrality.h` *(new)*
- `src/proc/HyperlaneCentrality.cpp` *(new)*
- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `tests/test_hyperlane_centrality.cpp` *(new)*
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 70 — Hyperlane Structural Criticality (Bridges + Articulation Points) + Criticality Overlay

This round adds a second, complementary graph analytics layer for hyperlanes: **structural fragility**.
Where betweenness centrality answers “what lanes are *busy*?”, this answers “what lanes/systems are *single points of failure*?”.

### What changed

- **Core / proc:** added `proc::analyzeHyperlaneVulnerability(...)`:
  - Finds **bridges** (cut-edges) and **articulation points** (cut-vertices) on the undirected hyperlane graph.
  - Computes per-bridge **cut size** (how many nodes end up on the smaller side if the bridge fails) and per-articulation **impact**
    (how many nodes fall outside the largest remaining component if the system disappears).
  - Works per connected component so “impact” stays meaningful even if the lane graph is not fully connected.

- **Procedural Galaxy Lab:** Hyperlane overlay gained a new **Color Mode**:
  - **Criticality (bridges/articulation)** *(new)*
    - Highlights fragile **bridge lanes** (thicker + hotter color for higher cut impact).
    - Draws a ring overlay around the top **articulation systems** (radius/color scale with impact).
    - Adds UI controls for **Top Bridges / Top Articulation**, plus width/alpha boosts.

- **Route overlay integration:** when Criticality mode is active, route segments are colored by the same criticality metric
  (so you can immediately see if the chosen path relies on bridge lanes).

- **Tests / build hygiene:**
  - Converted the Round 69 centrality test from Catch2 style to the repo’s **`test_harness`** entrypoint format.
  - Added a deterministic unit test for bridge/articulation detection on simple graphs.

### Files changed/added

- `include/stellar/proc/HyperlaneVulnerability.h` *(new)*
- `src/proc/HyperlaneVulnerability.cpp` *(new)*
- `apps/stellar_game/ProceduralGalaxyLabWindow.h`
- `apps/stellar_game/ProceduralGalaxyLabWindow.cpp`
- `tests/test_hyperlane_centrality.cpp` *(rewritten to harness style)*
- `tests/test_hyperlane_vulnerability.cpp` *(new)*
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 71 — Procedural Fluid Lab (Stable Fluids + Curl Noise Forcing)

This round adds an **interactive procedural fluid simulation** playground aimed at quickly generating
convincing “space smoke” / nebula textures and experimenting with turbulence controls.

### What changed

- **Core / proc:** added `proc::FluidSim2D` *(new)*:
  - 2D incompressible solver in the “Stable Fluids” family (semi-Lagrangian advection, implicit diffusion, pressure projection).
  - **RGB dye** channels for colorful swirls.
  - Optional **vorticity confinement** to re-inject small-scale rolling motion.
  - Optional divergence-free **curl-noise forcing** (stream function based) for controllable procedural turbulence.
  - Safety clamps for extreme brush input (max speed / max dye).

- **Game UI:** added **Procedural Fluid Lab** window *(new)*:
  - Paint dye + push velocity by dragging directly on the preview.
  - Live-tweak viscosity/diffusion/dissipation, vorticity, curl-noise parameters, and safety clamps.
  - **PNG export** of the current simulated texture (useful for baking animated sheets, nebula masks, etc.).

- **Tests:** added a deterministic solver smoke test that checks:
  - no NaNs/Infs during stepping
  - divergence stays reasonably small

### Files changed/added

- `include/stellar/proc/FluidSim2D.h` *(new)*
- `src/proc/FluidSim2D.cpp` *(new)*
- `apps/stellar_game/ProceduralFluidLabWindow.h` *(new)*
- `apps/stellar_game/ProceduralFluidLabWindow.cpp` *(new)*
- `apps/stellar_game/main.cpp`
- `tests/test_fluid_sim2d.cpp` *(new)*
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 72 — Text Animation Modifier (TextFx)

### What’s new
- Added a lightweight **markup-driven text animation/modifier system** (`stellar::ui::textfx`) that can be used anywhere you draw UI text.
- New **Text Animation Lab** window to author/preview effects in real time (Visual → Text Animation Lab).
- **Toast overlay** now supports TextFx markup (existing bracketed gameplay text like `[CONTRABAND]` stays literal; only recognized tags affect rendering).

### Markup examples
- `[wave amp=6 freq=0.25 speed=2]Hello[/wave]`
- `[pulse min=0.2 max=1 speed=2][color #ff4444]WARNING[/color][/pulse]`
- `[grad #ff00aa #00ccff]NEBULA[/grad]`
- `[scramble amount=0.9 rate=28 set=hex]HACKING...[/scramble]`
- `[type cps=28 fade=0.06]Incoming transmission...[/type]`

### Notes
- Unknown tags are preserved literally (important for existing strings that use brackets).
- The core compiler/evaluator is **ImGui-free**; ImGui rendering wrappers compile only when `STELLAR_ENABLE_IMGUI` is enabled.

### Files changed/added

- `include/stellar/ui/TextFx.h` *(new)*
- `src/ui/TextFx.cpp` *(new)*
- `apps/stellar_game/TextAnimationLabWindow.h` *(new)*
- `apps/stellar_game/TextAnimationLabWindow.cpp` *(new)*
- `apps/stellar_game/main.cpp`
- `CMakeLists.txt`
- `tests/test_text_fx.cpp` *(new)*
- `PATCH_NOTES.md`

## Round 73 — Shader Lab External iChannel Textures (Checker/Noise + Live Fluid/Flow)

This round upgrades the **Procedural Shader Lab** so your shader graphs can sample **non-feedback textures**
in iChannel0..3 — similar to how ShaderToy exposes image/noise inputs.

### What changed

**Core / render:**
- `render::ShaderToyGraph` now supports **External0..External3** channel sources.
  The app can bind textures via `ShaderToyGraph::setExternalTexture(slot, tex)` and route them per-pass.

**Game UI (Procedural Shader Lab):**
- Adds an **External Textures (iChannel)** panel with live previews and controls.
- The lab binds a curated default set:
  - **External0:** Checker pattern
  - **External1:** Seeded, tiling RGB noise
  - **External2:** Live CPU **fluid dye** texture (from `proc::FluidSim2D`)
  - **External3:** Live **flow/velocity** texture derived from the same fluid sim
    (velocity encoded into RG, where 0.5 = 0).
- New **Inkflow (External fluid flow)** preset demonstrates advecting a feedback buffer using the external flow field.

### Files changed/added

- `include/stellar/render/ShaderToyGraph.h`
- `src/render/ShaderToyGraph.cpp`
- `apps/stellar_game/ProceduralShaderLabWindow.h`
- `apps/stellar_game/ProceduralShaderLabWindow.cpp`
- `PATCH_NOTES.md`

## Round 74 — Audio Analyzer / Oscilloscope + FFT Spectrum

This round focuses on the **procedural audio subsystem**, adding tooling to *see* what the synth is producing.

### What’s new

**Core (headless-friendly DSP):**
- New `stellar::dsp::AudioAnalyzer` module implementing a **radix-2 FFT** analyzer with a **Hann window**, plus log-banded spectrum output.
- Produces:
  - time-domain waveform window
  - linear spectrum magnitude
  - log-frequency bands with optional smoothing and dB normalization

**Game (UI tooling):**
- New **Audio Analyzer** window (Menu → Audio → Audio analyzer… or Windows → Audio Analyzer)
  - oscilloscope waveform
  - spectrum view (dB or linear)
  - quick SFX preview buttons
  - export the most recent capture to a **mono 16-bit WAV**

**Audio engine upgrade:**
- `AudioEngine` now maintains a lightweight ring-buffer capture of the **mixed output** (mono) for visualizers/tools.

### Files changed/added

- `include/stellar/dsp/AudioAnalyzer.h` *(new)*
- `src/dsp/AudioAnalyzer.cpp` *(new)*
- `apps/stellar_game/AudioAnalyzerWindow.h` *(new)*
- `apps/stellar_game/AudioAnalyzerWindow.cpp` *(new)*
- `apps/stellar_game/AudioEngine.h`
- `apps/stellar_game/AudioEngine.cpp`
- `apps/stellar_game/main.cpp`
- `tests/test_audio_analyzer.cpp` *(new)*
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 75 — Procedural Mission Briefings + Risk-Sorted Mission Board

This round builds out a historically under-developed area: **missions** had solid mechanics, but little contract-style
presentation and no consistent way to reason about **risk** at a glance.

### What changed

- **Sim (headless / deterministic):** added `stellar::sim::MissionBriefing`:
  - `computeMissionRisk(...)` estimates **overall / danger / law** risk (0..1) from:
    - `effectiveSystemSecurityProfile` (security, piracy, traffic, contest)
    - straight-line **distance (ly)** to the mission's next objective
    - jurisdiction **law profile** + **black market profile** for smuggling
  - `generateMissionBriefing(...)` produces a deterministic contract pack:
    - short **contract code**
    - procedural **contact handle**
    - title + synopsis + bullet list (optionally using `ui::textfx` markup)

- **Game (Mission Board UX):**
  - Added offer **sorting**: reward, deadline, risk, distance, reward/ly.
  - Added a per-offer **risk summary line** (tier + component breakdown).
  - Added a **Brief** popup with:
    - TextFx-rendered contract title/synopsis/bullets
    - risk progress bars (overall, danger, law)
    - copy-to-clipboard (strips markup)
    - accept & plot directly from the contract view

- **Tests:** new deterministic coverage for the briefing generator.

### Files changed/added

- `include/stellar/sim/MissionBriefing.h` *(new)*
- `src/sim/MissionBriefing.cpp` *(new)*
- `tests/test_mission_briefing.cpp` *(new)*
- `apps/stellar_game/main.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 76 — Comms Inbox + Incoming Transmission Overlay

This round focuses on a long-missing piece of *diegetic UX*: the game had great systems (pirates, police scans,
missions) but **no in-universe comms history**. When something happened, it was easy to miss the toast, and there
was no way to re-open a contract or recall an ultimatum.

### What’s new

- **Sim:** new lightweight `stellar::sim::CommsLog` + `CommsMessage` module
  - bounded message log (defaults to 256)
  - unread + pinned support
  - small helpers that generate *markup-ready* transmissions for:
    - pirate ultimatums
    - authority bounty demands
    - corrupt-scan bribe offers
    - mission contract briefings

- **Game:** new **Comms / Inbox** window
  - filter + channel selection + unread-only
  - message list + details view (TextFx markup rendering)
  - convenience actions: mark read/unread, pin, plot route to the message’s linked destination

- **HUD:** new **incoming transmission overlay**
  - non-interactive, animated preview using TextFx
  - queues transmissions so you can see what arrived without opening a window

### Files changed/added

- `include/stellar/sim/Comms.h` *(new)*
- `src/sim/Comms.cpp` *(new)*
- `apps/stellar_game/CommsWindow.h` *(new)*
- `apps/stellar_game/CommsWindow.cpp` *(new)*
- `tests/test_comms.cpp` *(new)*
- `apps/stellar_game/main.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 77 — TextFx Performance Pass (LRU Cache + Sweep-Line Span Evaluation)

This round focuses on **improving existing code** in one of the most “every-frame” paths: TextFx markup rendering.
As the Comms inbox + overlay leaned on TextFx more heavily, the previous implementation was doing extra work:

- Recompiling the same markup strings every frame in UI code.
- Evaluating each glyph by scanning *all spans* (`O(glyphs * spans)`), even when only a few spans were active.

### What’s new

- **TextFx Program now stores `glyphCount`**
  - computed during compilation (no need to re-scan UTF-8 in hot paths)

- **New `ui::textfx::ProgramCache` (bounded, LRU-ish)**
  - caches compiled markup Programs
  - supports heterogeneous lookup (`std::string_view`) to avoid allocations on cache hits

- **TextFx Draw optimized with a sweep-line active span set**
  - maintains the set of spans that apply to the current glyph as we draw left-to-right
  - reduces per-frame span scanning overhead significantly for longer marked-up strings

- **Comms UI updated to use TextFx caching throughout**
  - overlay + inbox no longer recompiles / re-strips markup every frame
  - “Copy plain” uses cached plain text instead of recompiling

- **Build / headless support**
  - Moved **Starfield** + **Nebula** generators into the CPU/core build so procedural sky code can compile and be tested
    in environments without OpenGL (better CI coverage).

- **Tests**
  - added coverage for `Program::glyphCount` correctness and basic `ProgramCache` behavior

### Files changed/added

- `include/stellar/ui/TextFx.h`
- `src/ui/TextFx.cpp`
- `apps/stellar_game/CommsWindow.cpp`
- `apps/stellar_game/TextAnimationLabWindow.cpp`
- `tests/test_text_fx.cpp`
- `PATCH_NOTES.md`

## Round 99 — Station Docking Pads (Procedural Layout + Clearance Pad Assignment)

This round expands the *station docking* loop beyond the single “center hangar point” by introducing a
**procedural docking pad layout** per station and wiring it into the docking clearance + docking flow.

### What’s new

- **Procedural docking pad layout (deterministic)**
  - New `sim::DockingPads` module generates a stable set of pad poses in station-local space.
  - Uses a farthest-point sampling pass ("blue-noise-ish") so pads spread nicely without obvious grid artifacts.

- **Clearance now assigns a pad number**
  - `DockingClearanceState` stores `assignedPad` (1-based).
  - When clearance is granted, the station returns an assigned pad number (stable for that clearance).

- **Docked ship snaps to assigned pad**
  - While docked, the ship is attached to the station at the assigned pad pose (instead of a fixed center point).
  - The pad is recovered automatically on load by selecting the nearest pad to the saved docked ship position.

- **Undock orientation fix**
  - Undocking now orients the ship to face away from the station along the slot axis.

### Files changed/added

- `include/stellar/sim/DockingPads.h` *(new)*
- `src/sim/DockingPads.cpp` *(new)*
- `include/stellar/sim/DockingClearanceService.h`
- `src/sim/DockingClearanceService.cpp`
- `apps/stellar_game/main.cpp`
- `tests/test_docking_clearance.cpp`
- `tests/test_docking_pads.cpp` *(new)*
- `CMakeLists.txt`
- `PATCH_NOTES.md`

## Round 100 — 3D Distance & FOV Toolkit (Dolly Zoom Lock + Physical Lens Readout)

This round focuses on **distance + FOV**: getting the camera to behave predictably as you change distance,
while also exposing the key projection numbers you usually end up re-deriving by hand.

### What’s new

- **New `math::fov` utilities** (header-only)
  - Vertical↔horizontal FOV conversion for any aspect ratio.
  - Dolly-zoom FOV math (keep framing constant while distance changes).
  - Angular diameter (exact via `asin(r/d)`), view-height-at-distance, and units-per-pixel helpers.
  - Simple physical camera equivalence: focal length ↔ FOV given a sensor size.

- **Camera rig: distance-aware FOV tools**
  - **Telemetry**: viewport aspect, vertical/horizontal FOV, camera distance, view height at depth, and scale (U/px).
  - **Dolly zoom lock**: keeps the subject’s framing stable while you orbit-zoom.
    - Supports auto-capture and manual capture of the reference framing.
  - **Physical lens readout** (sensor + focal length → FOV) with a one-click “apply to base FOV”.
  - **Dominant body framing helper**: compute angular diameter as seen from the camera and set base FOV to frame it.

### Fixes

- Fixed a couple build-breaking issues inside the camera rig module (`kPi` constant + cinematic FOV gating).

### Files changed/added

- `include/stellar/math/Fov.h` *(new)*
- `apps/stellar_game/CameraRigWindow.h`
- `apps/stellar_game/CameraRigWindow.cpp`
- `tests/test_fov_math.cpp` *(new)*
- `PATCH_NOTES.md`

## Round 101 - Procedural generation rendering

### Rendering upgrades

- **GPU surface cache crater shader fixed**: renamed GLSL variable `out` to avoid a reserved keyword collision, re-enabling the GPU surface cache when the shader previously failed to compile.
- **Star corona shader fixed**: removed an `n` redefinition by renaming the scalar noise term, restoring the renderer when it was disabled due to fragment compilation errors.
- **RenderTarget2D: new descriptor**
  - Optional depth attachment (depthless targets for fullscreen procedural baking).
  - Optional mipmap allocation + `generateMips()` helper for post-bake refresh.
- **ProceduralGraphBaker quality pipeline**
  - Optional mip generation for baked textures (better minification / fewer moire artifacts).
  - Optional lightweight dithering to reduce 8-bit banding.
  - Now uses a depthless render target when baking 2D procedural graphs.

### UI

- Procedural Lab: added toggles for mip generation and dithering.
- Procedural Mesh Lab: same controls for baked 2D surface textures; stats now include mip generation time.

### Files changed/added

- `include/stellar/render/RenderTarget.h`
- `src/render/RenderTarget.cpp`
- `include/stellar/render/ProceduralGraph.h`
- `src/render/ProceduralGraph.cpp`
- `src/render/GpuSurfaceCache.cpp`
- `src/render/StarCoronaRenderer.cpp`
- `apps/stellar_game/ProceduralLabWindow.h`
- `apps/stellar_game/ProceduralLabWindow.cpp`
- `apps/stellar_game/ProceduralMeshLabWindow.h`
- `apps/stellar_game/ProceduralMeshLabWindow.cpp`
- `PATCH_NOTES.md`

## Round 102 - Seamless Triplanar + Micro-Detail Normals for Procedural Materials

### Procedural graph baking

- **ProceduralGraphBaker**: added an optional *"Pack height in alpha"* mode that writes the scalar graph output `t` into the baked texture's alpha channel. This keeps the existing RGB albedo workflow intact while enabling downstream uses like height-based micro normals.

### Raymarch preview shading

- **SdfRaymarcher**: upgraded the albedo projection from a single dominant-axis selection to an optional full **triplanar blend** (with sharpenable weights).
- Added an optional **micro-normal** perturbation path that derives a surface-gradient from a height channel (alpha or luminance), then projects that gradient onto the geometric normal's tangent plane for stable detail shading.

### UI

- **Procedural Lab**: new bake-quality toggle: *Pack height in alpha*.
- **Procedural Mesh Lab**:
  - new bake-quality toggle: *Pack height in alpha*.
  - new raymarch controls: *Triplanar blend*, *Tri sharpness*, *Height from alpha*, *Micro normal strength*, and *Micro normal step*.

### Files changed/added

- `include/stellar/render/ProceduralGraph.h`
- `src/render/ProceduralGraph.cpp`
- `include/stellar/render/SdfRaymarcher.h`
- `src/render/SdfRaymarcher.cpp`
- `apps/stellar_game/ProceduralLabWindow.h`
- `apps/stellar_game/ProceduralLabWindow.cpp`
- `apps/stellar_game/ProceduralMeshLabWindow.h`
- `apps/stellar_game/ProceduralMeshLabWindow.cpp`
- `PATCH_NOTES.md`

## Round 104

### Cross-system Integration Hub (Actions + Events)

This round adds a new **Integration Hub** debug window that acts as a lightweight bridge between otherwise loosely-coupled systems:

- **GameActionQueue**: systems (e.g., Time Trial) can request cross-system actions (target station, engage docking computer, set camera mode) without reaching into `main.cpp` state directly.
- **GameEventLog**: systems can emit structured events for debugging and for attaching to repro packs.
- **JSON trace export**: write a compact `stellar_integration_trace` JSON to disk for bug reports.

### Time Trial cross-integration refactor

- Time Trial now emits **actions/events** (via sinks) instead of relying on one-off request booleans.
- The central game loop drains and executes actions in order, recording an execution event.

### Runtime Validation + Repro Packs

- Runtime repro packs can optionally embed a recent slice of Integration Hub events (time window + cap).
- Watchdog hits and repro-pack writes emit Validation events into the Integration Hub log.


## Round 105

### Run Capture: Time Trial ↔ Flight Recorder ↔ Integration Trace

This round deepens **cross-system integration** by wiring together gameplay + devtools so a time-trial run can produce a shareable, reproducible trace bundle.

- **Time Trial → Flight Recorder automation**
  - New toggles in the Time Trial window:
    - Clear Integration Hub log on arm (clean run trace)
    - Auto start Flight Recorder on first gate pass
    - Auto stop Flight Recorder on finish
    - Auto export traces on finish
  - When enabled, Time Trials now auto-drive the Flight Recorder and trace exporters via the Integration Hub action queue.

- **Integration Hub: expanded action vocabulary**
  - Added actions for:
    - Start/Stop Flight Recorder
    - Export Flight Recorder trace
    - Export Integration trace
    - Clear Integration Hub log

- **Main loop: action executor upgrades**
  - Action executor now performs the new trace/recording actions and emits structured events for them.
  - `EngageDockingComputer` now honors the `engage=false` action form (disengage).

- **Fix**
  - Fixed a bad reference in the docking-clearance request event path (was referencing a non-existent `dockingClearanceSvc`).

### Files changed/added

- `apps/stellar_game/GameSignals.h`
- `apps/stellar_game/IntegrationHubWindow.cpp`
- `apps/stellar_game/TimeTrialWindow.h`
- `apps/stellar_game/TimeTrialWindow.cpp`
- `apps/stellar_game/FlightRecorderWindow.h`
- `apps/stellar_game/FlightRecorderWindow.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

## Round 106

### Integration Hub Automations: Event → Action Rules (IFTTT-style)

This round turns the Integration Hub into a real *cross-system glue layer* by adding a small **automation/rule engine**:

- New **Automation** tab in the Integration Hub window.
- Define rules that match incoming **GameEvents** (by kind + tag match mode) and automatically enqueue **GameActions**.
- Includes safety rails:
  - Per-rule cooldown
  - Per-frame max action budget to avoid accidental action spam / loops
- Supports simple template expansion for action strings (e.g. export paths):
  - `{kind}`, `{tag}`, `{u64a}`, `{u64b}`, `{tRealSec}`, `{tRealMs}`, `{tSimDays}`, `{rule}`
- Ships with a few **starter preset rules** (disabled by default) so you can enable and tweak quickly.

### Missions emit Integration Hub gameplay events

To make automations useful for gameplay, core mission milestones now emit structured events:

- `MissionAccepted`
- `MissionComplete`
- `MissionFailed`
- `MissionLegComplete`

These include IDs in `u64a/u64b` plus a human-readable message.

### Mission Briefing: optional risk hints & reputation cues

`MissionBriefingParams` now includes:

- `includeRiskHints`
- `includeReputationCues`

When enabled, the generated briefing appends short, context-sensitive hints derived from the computed risk model.

### Fixes / cleanup

- Added missing `GameEventKind::Debug` (was referenced by the main loop).
- Moved `timeRealSec` into the main time block so lambdas can consistently timestamp emitted events.

### Files changed/added

- `apps/stellar_game/GameSignals.h`
- `apps/stellar_game/IntegrationHubWindow.h`
- `apps/stellar_game/IntegrationHubWindow.cpp`
- `apps/stellar_game/main.cpp`
- `include/stellar/sim/Mission.h`
- `include/stellar/sim/MissionBriefing.h`
- `src/sim/MissionBriefing.cpp`
- `PATCH_NOTES.md`

## Round 118 - Nav Assist Follow (formation) + Escort HUD integration

Adds a new **Nav Assist** mode designed for the *most annoying part of escort gameplay*: staying comfortably in formation with a moving convoy.

- **Nav Assist: Follow**
  - Chases a *moving follow point behind the target* based on its travel direction.
  - Uses exponential smoothing on the follow direction to avoid jitter when target velocity is noisy.
  - Can optionally keep the ship **facing the target** while translating (default), which makes convoy escorting feel much more natural.
- New default keybind: **Home** (`NavAssistFollow`).
- **Traffic Escort HUD** gets a one-click **Target + Follow** button to immediately lock onto the convoy and engage Follow mode (stays safely inside contract max range).

### Files changed/added

- `apps/stellar_game/ControlsConfig.h`
- `apps/stellar_game/ControlsWindow.cpp`
- `apps/stellar_game/main.cpp`
- `include/stellar/sim/NavAssistComputer.h`
- `src/sim/NavAssistComputer.cpp`
- `PATCH_NOTES.md`


## Round 122 - Star corona perf pass (GPU + CPU)

This round targets a real-world perf hotspot: **procedural star corona shading** (lots of fragment work + per-frame instance uploads).

### What changed

- **Shader early-out**: fragments with negligible rim contribution now return immediately, skipping the expensive 3D FBM noise path.
- **Quality knob**: `noiseOctaves` (1..8) lets you trade detail for speed.
- **Streaming instance uploads**: instance VBO updates now prefer `glMapBufferRange` with invalidation and fall back to an orphan + `glBufferSubData` path.
- **VAO state caching**: instanced attribute pointers are configured once per mesh VAO (instead of every draw) to reduce GL driver overhead.

### New UI controls

In **World visuals → Star corona**:

- `Corona noise octaves` *(perf)*
- `Corona rim early-out` *(perf)*

### Files changed/added

- `include/stellar/render/StarCoronaRenderer.h`
- `src/render/StarCoronaRenderer.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

## Round 124 - Fast HDR buffers (R11G11B10F) to improve FPS

This round targets a common bottleneck on integrated GPUs: **memory bandwidth**.
The HDR scene buffer and bloom targets previously used **RGBA16F** (64 bpp). We
now support a packed HDR format, **R11G11B10F** (32 bpp), which can significantly
reduce bandwidth and improve frame rate.

### What changed

- **New PostFX setting:** HDR buffer format
  - `RGBA16F (quality)`
  - `R11G11B10F (fast)`
- **Automatic fallback:** if the driver doesn't support `R11G11B10F` as an FBO
  color attachment, PostFX automatically falls back to `RGBA16F` at runtime.
- **Bloom transient textures** now match the active HDR format so the whole
  post pipeline benefits from the reduced bandwidth.
- Added a UI control under **PostFX → Performance / Render Scale** and shows the
  **active** runtime format.

### Files changed/added

- `include/stellar/render/Gl.h`
- `include/stellar/render/PostFX.h`
- `src/render/PostFX.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

## Round 125 - Async auto exposure (PBO) + bloom bandwidth cuts

This round targets two **frame-time spikes** that can show up on mid-range and integrated GPUs:

1) **Auto exposure GPU readback stalls**: the previous implementation used a synchronous `glGetTexImage` each frame (even though it reads only a 1×1 texture). On many drivers this forces a CPU↔GPU sync point.
2) **Redundant bloom clears**: bloom passes rendered full-screen, but also cleared the color buffer beforehand—doubling the memory bandwidth for no visual gain.

### What changed

- **Async auto exposure readback** (optional, enabled by default):
  - Uses a tiny **pixel-pack-buffer (PBO) ring** + **fence sync**.
  - Consumes the newest ready sample **without stalling**.
  - Adds a small **~1–2 frame latency** to exposure adaptation (visually fine, vastly smoother).
- **Removed redundant `glClear` calls** in bloom bright-pass and blur passes.

### New UI control

In **PostFX → Tonemap / Output → Auto exposure**:

- `Async readback (no stalls)` *(recommended)*

### Files changed/added

- `include/stellar/render/PostFX.h`
- `src/render/PostFX.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`


## Round 126 - CPU-side render caches + streamed line/mesh buffers

This round focuses on **CPU-side per-frame work** and **driver stalls** that show up as FPS drops on integrated GPUs:

- **Orbit line sampling** (planets + stations) was being re-generated every frame.
- **Asteroid spin parameters** were being re-derived from RNG every frame.
- **LineRenderer / MeshRenderer** were using upload paths that can trigger buffer reallocations and implicit sync.

### What changed

- **Orbit line cache**: planet + station orbit line vertices are now built **once per system** and reused each frame.
- **Asteroid spin cache**: deterministic spin axis/phase/rate is now computed **once per asteroid id** and cached.
- **LineRenderer streaming VBO**: grows on demand, then uses `glMapBufferRange` (with invalidate) for low-stall updates.
- **MeshRenderer streaming instance VBO**: same growth + mapping approach to reduce per-draw driver overhead.

### Files changed/added

- `include/stellar/render/LineRenderer.h`
- `src/render/LineRenderer.cpp`
- `include/stellar/render/MeshRenderer.h`
- `src/render/MeshRenderer.cpp`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`


## Round 127 - Fewer GL driver calls + point streaming ring buffers (FPS)

This round targets two common sources of frame time spikes in OpenGL apps:

1) **GL driver overhead** from repeatedly rebinding per-instance attribute layouts.
2) **Implicit sync/stalls** when streaming dynamic point-sprite data into a buffer
   the GPU is still consuming.

### What changed

- **MeshRenderer instance attribute layout is now bound once per mesh**.
  The instance attribute pointers (locations 3..6) are VAO state. Re-specifying
  them every draw call costs CPU time on some drivers. We now bind them lazily
  once per mesh and reuse the state.

- **PointRenderer now uses a small VAO/VBO ring (triple-buffer) for streaming**.
  Each draw rotates to the next buffer so the CPU never overwrites the buffer
  potentially still in-flight on the GPU. This reduces driver synchronization
  stalls when drawing large starfields/nebula/particle sprites.

### Files changed/added

- `include/stellar/render/PointRenderer.h`
- `src/render/PointRenderer.cpp`
- `include/stellar/render/MeshRenderer.h`
- `src/render/MeshRenderer.cpp`
- `PATCH_NOTES.md`


## Round 128 - Objective HUD + Integration Hub build fixes + safer HUD layout padding

This round fixes a hard build break introduced by stale Objective HUD / Integration Hub glue
code and improves the robustness of HUD overlay positioning.

### What changed

- **Fixed Objective HUD emitting an invalid `GameEvent` initializer**.
  A stale brace-initializer referenced a non-existent `EventAttrs` type, causing a large
  cascade of syntax errors. Objective HUD now emits a normal `GameEvent` via a helper.

- **Added `makeActionSyncNavToMission()` helper + migrated call sites**.
  `SyncNavToMission` actions were being created with brittle aggregate initializers.
  The new helper prevents field-order mistakes and avoids MSVC narrowing diagnostics.

- **HUD layout helpers now support per-widget padding**.
  `hudSetNextWindowPosFromLayout()` / `hudCaptureWindowPosToLayout()` accept optional
  padding in pixels so HUD windows can be anchored inside a safe edge margin.
  This unblocks the Objective HUD fallback tracker’s padded placement.

- **Fixed incorrect `projectToScreenAny()` usage in the offscreen waypoint indicator**.
  Updated the call to match the current signature and removed invalid `ImVec2 +=` usage
  (not supported by ImGui’s vector type).

- **Declared missing runtime state for HUD toggling and mining toast timestamps**.
  Restores compilation and keeps future work on per-commodity mining-toasts straightforward.

### Files changed/added

- `apps/stellar_game/GameSignals.h`
- `apps/stellar_game/main.cpp`
- `PATCH_NOTES.md`

## OpenAI Patch Round 2 - Integration Hub Trace Import + Replay

This patch expands the **Integration Hub** into a proper “trace lab”:
you can now **import** trace JSON into a safe staging buffer, selectively apply it to the current hub,
and **replay** recorded actions with timing controls.

### Fixes

- Fixed a stray `} // namespace` in `IntegrationHubWindow.cpp` that prematurely closed `stellar::game`
  and caused a cascade of MSVC syntax errors.

### New

- Added a small, dependency-free JSON parser: `include/stellar/core/SimpleJson.h` (header-only).
- Added Integration Hub **Import** tab:
  - Load a `stellar_integration_trace` JSON from file or clipboard into a staging buffer.
  - Apply events/actions/rules selectively (with clear warnings before pushing *pending* actions).
- Added **Replay scheduler**:
  - Schedule actions from imported trace (or current action history).
  - Timed replay (preserve original trace deltas), speed multiplier, lead-in delay.
  - Filters: exclude Toasts and file-output actions by default.
  - Optional: temporarily disable automations during replay and automatically restore afterwards.

## OpenAI Patch Round 142 - Scan Intel UX + Test Harness Hardening

This patch focuses on **stability + observability**:

### Fixes

- Fixed headless test build breakages by expanding the bespoke test harness:
  - Added `CHECK_EQ()`.
  - Fixed missing `failures` counter in `test_signal_event_reactivity`.
  - Removed the accidental Catch2 dependency in `test_system_event_economy`.

- Fixed `IntegrationHubWindow.cpp` accidentally calling `atomicWriteFile()` with a `std::string`
  payload (triggered the `WriteFn must be invocable` static assertion).

- Fixed MSVC build errors in the Contacts UI:
  - Updated threat estimation to use the new `sim::ShipScanInput` shape.
  - Removed references to stale fields (`weaponSecondary`, armor/resist/bounty scan inputs).
  - Removed assignment to a `const` local (`scanRangeKm`).

- Fixed missing include for `sim::stationPosKm()` usage in `CopilotWindow.cpp`.

### New

- **Scan intel tooltip expansion**: Contacts now show key scan report fields (quality/threat/health/cargo/EW)
  directly in the hover tooltip.
- **Integration Hub event for scans**: Completing a contact scan now emits a `Gameplay/ContactScan` event
  so the rule engine and exported traces can react to scans.

### Files changed/added

- `apps/stellar_game/main.cpp`
- `apps/stellar_game/IntegrationHubWindow.cpp`
- `apps/stellar_game/CopilotWindow.cpp`
- `tests/test_harness.h`
- `tests/test_signal_event_reactivity.cpp`
- `tests/test_system_event_economy.cpp`
- `PATCH_NOTES.md`

## OpenAI Patch Round 143 - Mission Control: Itinerary Planner + Runner

This patch tackles an **underdeveloped gameplay workflow**: once you've stacked multiple missions,
it's hard to decide *where to go next* and easy to waste time re-plotting routes.

### New

- **Mission itinerary planner (sim)**: a deterministic, heuristic planner that builds a "best next stops"
  itinerary for a selected set of active missions.
  - Plans using either **min-hops** or **min-distance** cost models.
  - Groups by *station objective* (or optionally by *system*) and scores stops via reward/cost with
    deadline urgency and risk penalties.
  - Designed to be fast enough for UI, and safe even when route batch reachability is incomplete
    (falls back to straight-line hop estimates).

- **Mission Control window (game UI)**:
  - Mission selection tools (select all/none/tracked) + optional visibility toggles for completed/failed.
  - Plan summary table with one-click **Route/Go** buttons per stop.
  - Per-mission **briefing renderer** (markup title/synopsis/bullets) with risk + objective metadata.

- **Itinerary runner** (optional automation):
  - Arms from the current plan and then advances stops automatically as mission objectives update.
  - Can auto-plot the next stop and auto-track a representative mission on each leg.
  - Uses the existing Integration Hub action pipeline for go-to + camera rig travel preset.

### Files changed/added

- `include/stellar/sim/MissionPlanner.h` *(new)*
- `src/sim/MissionPlanner.cpp` *(new)*
- `apps/stellar_game/MissionControlWindow.h` *(new)*
- `apps/stellar_game/MissionControlWindow.cpp` *(new)*
- `apps/stellar_game/main.cpp`
- `CMakeLists.txt`
- `PATCH_NOTES.md`
