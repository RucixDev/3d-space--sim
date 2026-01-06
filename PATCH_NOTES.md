## 2026-01-06 (Patch) - Traffic Convoys: Schedule Hydration + Signal Expiry Consistency

This round tightens **robustness** in the traffic/convoy replay pipeline, especially for older or
hand-edited saves where shipment schedule metadata may be missing.

- Added `shipmentScheduleLooksValid(...)` and `hydrateShipmentScheduleFromId(...)` so convoy replay
  can deterministically recover missing `departDay/arriveDay/distKm/speedKmS` from shipment ids.
- `TrafficConvoyLayer` ledger replay now hydrates schedules before evaluation, and skips broken
  station references to avoid spawning "(0,0,0)" convoys.
- Clarified `SignalSite.expireDay` semantics (expired when `timeDays > expireDay`) and aligned
  convoy filtering to match.
- New regression coverage in `test_traffic_convoy_layer`: replay a schedule-less shipment and
  assert deterministic mid-flight visibility.


## 2026-01-06 (Patch) - Game: Ledger-driven Traffic Convoys + TrafficLedger Schedule Fixes

This round focuses on **bug fixing + integration polish** for the traffic stack:

- The game already runs `simulateNpcTradeTraffic(...)` to nudge markets.
- We already had `TrafficLedger` to record the *actual* shipments chosen by that nudger.
- We already had `Signals` + `TrafficConvoyLayer` to expose those shipments as **moving convoy signals**.

### Key fixes & integrations

- **Game integration:** `stellar_game` now maintains a per-system `TrafficLedger` and passes it into:
  - `simulateNpcTradeTraffic(...)` (records shipments while markets are nudged)
  - `generateSystemSignals(..., trafficLedger)` and the convoy refresh path

  Result: visible **Traffic Convoy** signals now reflect the *real* commodity flows driving the station economies.

- **Bug fix:** `TrafficLedger::makeNpcTradeShipment(...)` schedule metadata now uses the same **distance-at-arrival refinement** logic as the traffic lane prototype, so `distKm/speedKmS/arriveDay` are internally consistent even when stations move during the trip.

- **Safety guard:** `TrafficLedger::record(...)` now ignores duplicate shipment ids (protects against accidental double-simulation in callers).

### Tests

- `test_traffic_ledger` now sanity-checks that a shipment's stored `distKm` and `speedKmS` match the computed endpoints/timeline.


## 2026-01-06 (Patch) - Signals: Ledger-backed Traffic Convoy Signals (Integration)

This round continues the **integration + bug-fixing** theme by closing the loop between:

- `sim::simulateNpcTradeTraffic(...)` (ambient inventory nudger)
- `sim::TrafficLedger` (records the chosen shipments)
- `sim::Signals` (in-system signal generation)

### Key change

- `generateSystemSignals(...)` now accepts an **optional** `TrafficLedger*`.
  - When provided, `SignalKind::TrafficConvoy` sites are derived from the ledger shipments via
    `sim::TrafficConvoyLayer` (so visible convoys can reflect the *actual* trade flows that
    moved station inventories).
  - If the ledger has no matching shipments, signals fall back to the existing deterministic
    `sim::TrafficLanes` convoy prototype (so UI/gameplay doesn’t go empty).

### Tests

- Extended `test_signals` to validate that ledger-backed convoy sites map back to recorded shipments
  (ids + commodity + units) and become active when sampled mid-flight.


## 2026-01-06 (Patch) - Game: Moving Traffic Convoy Signals (Integration + Bug Fixes)

This round focuses on **bug fixing and integrating existing systems** by making the new
`sim::Signals` *TrafficConvoy* sites work end-to-end inside `stellar_game`:

- **Enabled** traffic convoy signal generation on system entry (uses the existing `sim::Signals` + `sim::TrafficLanes`).
- **Fixed a missing enum mapping** in the game that previously would have mis-typed convoy sites when enabled.
- **Made convoys actually move in-game:**
  - Periodic refresh to **add/remove** active convoys as time advances.
  - Per-frame state evaluation so convoy positions/velocities stay current.
- **Fixed persistence bug:** convoy ids are deterministic, but convoys should *not* be written into `resolvedSignalIds`.
- **UI polish:** System map + selection panel show convoy commodity and basic manifest details; signal scan grants data for convoy scans too.


## 2026-01-05 (Patch) - Signals: Traffic Convoy Sites + Stable Daily Anchoring

This round improves the **headless signal generator** and tightens determinism:

- **New signal kind:** `SignalKind::TrafficConvoy`
  - Optional payload: `TrafficConvoy` + evaluated `TrafficConvoyState` at the query time.
  - Driven by the existing `sim::TrafficLanes` prototype (deterministic lane traffic).
- **New generation knobs:** `SignalGenParams::includeTrafficConvoys` and `SignalGenParams::trafficLaneParams`.
- **Bug fix:** static signal placement now uses a stable *mid-day anchor time* inside the current integer day,
  preventing small-dt drift when querying signals multiple times within the same day.
- **Sandbox tooling:** `stellar_sandbox --signals --signalTraffic`
  - Prints convoy sites with route/cargo/progress details.
  - JSON output includes a `trafficConvoy` object per site.
- **Tests:** extended `test_signals` coverage for convoy signal generation.

## 2026-01-05 (Patch) - Lambert Planner: Incremental Porkchop Stepper + In-game Auto-search

This round focuses on **interactive navigation tooling** and makes the Lambert porkchop search
usable inside real-time UI without freezing the frame.

Key changes:

- **New incremental API:** `sim::LambertPorkchopStepper`
  - Evaluates the same fixed-size `(departure, TOF)` grid as `searchLambertPorkchop(...)`, but in
    **bounded chunks** via `step(maxCells)`.
  - Maintains **best-K** candidates and optionally captures the full grid (`storeGrid`).
- **Refactor:** `searchLambertPorkchop(...)` now runs the stepper to completion internally
  (keeps behavior deterministic while removing duplicate looping logic).
- **New test:** `test_lambert_planner_stepper`
  - Verifies the stepper reproduces the one-shot helper results and reports correct progress.

- **Game UI:** Transfer planner (Lambert) gained an **Auto-search (porkchop)** section
  - Incremental progress bar + best-K table.
  - One-click **Apply**: arms the maneuver node at the candidate departure time and fills RTN Δv.
  - Optional coarse RK4 miss-distance validation on apply.

## 2026-01-05 (Patch) - Traffic Convoy Layer: Replay NPC Trade Shipments as Moving Convoys

This round closes the loop between **ambient economy traffic** and **visible in-system convoys**.

We already had:
- `sim::simulateNpcTradeTraffic(...)` (background market nudging)
- `sim::TrafficLedger` (records the shipments chosen by that nudger)
- `sim::TrafficLanes` (a deterministic lane + convoy kinematics prototype)

Now we can turn those recorded shipments into *actual moving convoy trajectories* for tooling and future gameplay hooks.

Key changes:

- **New core module:** `stellar/sim/TrafficConvoyLayer.*`
  - `convoyFromShipment(...)` converts a `TrafficShipment` to a `TrafficConvoy` (lossless metadata copy).
  - `generateTrafficConvoysFromLedger(...)` returns `TrafficConvoyView` items with evaluated `pos/vel/progress` at the query time.
  - `sampleTrafficConvoyPathKm(...)` provides a simple polyline sampler for future UI/renderer lane visualization.

- **Sandbox tooling:** `stellar_sandbox --ledgerConvoys`
  - Runs the background traffic sim for the chosen system, records shipments, then prints them as moving convoys.
  - Supports JSON output (`--json`) and backfill control (`--ledgerBackfill <d>`).

- **New test:** `test_traffic_convoy_layer`
  - Verifies determinism and endpoint correctness for ledger-derived convoy playback.

## 2026-01-05 (Patch) - Orbital Navigation: Lambert Porkchop Transfer Planner

This round adds a **navigation / astrodynamics tooling layer** on top of the existing Lambert solver:

- **New `LambertPlanner` module** (`include/stellar/sim/LambertPlanner.h`, `src/sim/LambertPlanner.cpp`)
  - Deterministic "porkchop" search over a grid of **departure time** × **time-of-flight**.
  - Returns **best-K** candidates sorted by a configurable score:
    - `depart`: minimize only departure Δv (useful for flybys)
    - `arrive`: minimize arrival relative speed (soft arrival)
    - `total`: minimize departure Δv + arrival relative speed (rendezvous-friendly)
    - `weighted`: custom blend
  - Optional full-grid capture for tooling/plotting.

- **Sandbox integration:** `stellar_sandbox --lambert`
  - Search transfers **from the chosen station** (`--fromSys`, `--fromStation`) to a target
    station/planet (`--lambertTarget`, `--lambertIdx`).
  - JSON output supports emitting the full porkchop grid with `--json --lambertGrid`.

- **New test:** `test_lambert_planner.cpp`
  - Verifies planner wiring/determinism using a canonical circular-orbit case.

## 2026-01-05 (Patch) - Supercruise: Corridor-Aware Drop Window + Diagnostics

This round upgrades the **supercruise guidance** to be more robust when approaching targets:
it now tracks **lateral (cross-track) motion** and provides a simple **predicted miss distance**
metric, which can be used by UI tooling and helps avoid dropping while sliding past the destination.

Key changes:

- `SupercruiseHud` now reports:
  - `lateralKmS`: relative speed orthogonal to the ship→destination line
  - `missKm`: predicted closest-approach distance under constant relative velocity

- `SupercruiseParams` adds `maxLateralFrac` (default 0.6) so the safe-drop window can require the
  approach to be reasonably "on-axis" (prevents weird sideways drops).

## 2026-01-05 (Patch) - Traffic Ledger (Economic Shipments → Gameplay Hook)

This round adds a small but powerful missing piece for the economy/traffic stack: a **traffic ledger**
that can record the *actual* station-to-station shipments chosen by the existing ambient model
(`sim::simulateNpcTradeTraffic`).

The goal is to bridge the gap between "invisible inventory deltas" and concrete gameplay/UI targets:

- Tooling can print **daily shipment manifests**.
- The game can spawn **visible convoys** that correspond to real commodity flows.
- Tests can verify **determinism** of shipment generation across seeds.

Key changes:

- **New core module:** `stellar/sim/TrafficLedger.*`
  - `TrafficShipment` records `{from,to,commodity,units}` plus deterministic schedule metadata
    (`departDay`, `arriveDay`, `distKm`, `speedKmS`).
  - `TrafficLedger` is a lightweight in-memory log with pruning + queries.

- **Traffic simulation hook:** `sim::simulateNpcTradeTraffic(...)`
  - Added an optional `TrafficLedger*` parameter to record shipments as the simulation runs.
  - Recording is **side-effect free**: schedule metadata is derived from shipment ids so the traffic
    RNG sequence (and inventory results) are not perturbed.

## 2026-01-05 (Patch) - Traffic Lanes (Ambient Convoys Prototype)

This round starts building a **real traffic layer**: deterministic, in-system **trade convoys** that
travel between stations on predictable "lanes".

The goal is to support the kind of space-sim flavor mentioned in the README (convoys you can
eventually **interdict** or **protect**) while keeping the feature **headless + deterministic** so it
can be tested and used by tooling.

Key changes:

- **New core module:** `stellar/sim/TrafficLanes.*`
  - Generates a daily schedule of convoys (seed + system + day stamp).
  - Convoys pick **export-ish** commodities from the origin station and bias toward destinations that
    **consume** them.
  - Provides a lightweight evaluator that returns a convoy's **pos/vel/progress** at a given time.
  - Adds a smooth **lane arc** offset (zero at endpoints) derived from the convoy id — no extra saved
    state required.

- **Sandbox tooling:** `stellar_sandbox`
  - **New flags:** `--convoys`, `--convoyAll`, `--convoyWindow <d>`.
  - JSON output includes a `convoys` section with schedule + evaluated state.

- **Tests:** added `test_traffic_lanes` for determinism and endpoint correctness.

## 2026-01-05 (Patch) - Market Dashboard (Economy Analytics UI)

This round adds a new **Market Dashboard** window to make trade planning easier without
digging through individual station screens.

Key changes:

- **New window:** *Market Dashboard*
  - Found under **Windows → Market Dashboard** (and in the command palette via the window registry).
  - Scan scope:
    - **Current system** (fast, focused)
    - **Nearby systems** within a configurable radius (includes a max-systems cap)

- **Commodity-centric market view:**
  - Shows **effective** bid/ask prices (including station fees), **inventory**, and a simple **trend**
    computed from the station's recent price samples.

- **Quick arbitrage summary:**
  - Displays best **buy** and best **sell** locations in the current scan and estimates profit per unit.
  - Includes a “how many units could you move” estimate based on **hold mass** and **credits**.

- **Integrated routing/targeting:**
  - **Route** button plots a route to the system and sets a **pending arrival target station**.
  - If already in-system, it can **target** the station immediately.

## 2026-01-05 (Patch) - Floating Origin Rendering (Precision @ Huge Distances)

This round tackles a classic space-sim problem: **GPU float precision jitter** when you get far from the
system origin. Even if simulation runs in doubles, most render paths (instance buffers, uniforms) end up
as floats, so tiny local offsets (ships/stations) start to wobble when absolute coordinates get large.

Key changes:

- **Floating origin / camera-relative rendering (game):**
  - A per-frame **render origin** in km (`gRenderOriginKm`) is subtracted from all world positions before
    converting to render units.
  - Default behavior follows the **player ship**, keeping it near `(0,0,0)` in render space.
  - **New CVar:** `r.floating_origin.enabled` (archived)

- **VFX correctness under origin shifts:**
  - `ParticleSystem` now supports `shiftOrigin(deltaU)` so already-spawned particles remain stable when
    the origin moves.
  - Beam visuals are shifted in the same way so laser lines remain consistent.

- **Lighting + atmosphere fixes:**
  - The star is no longer assumed to be at render-space origin when floating origin is enabled.
  - Atmosphere shader now takes an explicit **sun position** uniform (`uSunPos`) and the game sets it
    each frame.

## 2026-01-05 (Patch) - Workspaces 2.1: Per-Workspace Layout Profiles

This round continues improving **UI systems** by making **workspaces** and **layout profiles**
play nicely together.

Key changes:

- **New core helper:** `stellar/ui/UiLayoutProfiles.*`
  - Headless utilities to discover layout profiles in `ui_layouts/`, sanitize filenames,
    generate unique `.ini` paths, and safely copy/rename/delete profile files.
  - Used by the game UI to keep file operations consistent and safer.

- **Workspaces window upgrades:**
  - Each workspace can now be **bound** to a specific layout profile (Dear ImGui `.ini`) from a list.
  - **Assign unique layout file**: gives the workspace its own dedicated ini (cloned from the assigned
    layout if possible, otherwise seeded from current runtime layout).
  - **Save current layout → assigned**: snapshots the current runtime docking/window layout into the
    workspace’s assigned ini (auto-assigns a unique ini if needed).
  - **Duplicate** now auto-creates a dedicated layout file for the copy so layouts don’t conflict.
  - **New from current** can optionally create a dedicated layout file and activates immediately.
  - **Rename** can optionally rename the workspace’s layout file (when safe and not shared).
  - **Delete** can optionally delete the workspace’s layout file (when safe and not shared), and now
    properly activates a replacement workspace when deleting the active one.

## 2026-01-04 (Patch) - Multi-Viewport UI (Detachable Windows)

This round improves **UI systems** with optional **multi-viewport** support (a.k.a. *platform windows*).

When enabled, you can drag ImGui panels outside the main window and they will become real OS windows —
ideal for multi-monitor cockpit setups.

- **UI Settings → Multi-monitor (Viewports):**
  - Enable/disable multi-viewport at runtime.
  - Options:
    - hide detached windows from the taskbar
    - disable auto-merge into the main window
    - borderless (no OS decorations)
- **Persistence:** stored in `ui_settings.txt` alongside theme/scale/docking.
- **Renderer integration:** adds the required `UpdatePlatformWindows()` / `RenderPlatformWindowsDefault()` pass.

## 2026-01-04 (Patch) - UI Developer Tools: CVar Browser + Palette Console Mode

This round improves **UI systems** with a new developer-facing runtime-tuning panel and faster command execution.

- **New window:** **CVars** (Debug)
  - Searchable list of all core CVars (name + help), optional *relevance sort*.
  - Edit bool/int/float/string CVars in-place (read-only CVars are disabled).
  - Quick-set line for `name = value` plus **Save/Load** to a config file.
  - Quality-of-life: reset-to-default + copy-to-clipboard (`name = value`).

- **Command Palette upgrade:** prefix the query with **`>`** to enter *Console mode*
  - `> your.command args` runs the command directly.
  - Shows recent console history as runnable entries.

## 2026-01-04 (Patch) - Atmospheric Drag + Reentry Heating (Experimental)

This round adds an **experimental atmosphere physics** layer: a deterministic, headless
atmosphere model (per-planet) plus optional in-game **drag + heat** integration.

- **New core module:** `stellar/sim/Atmosphere.*`
  - Deterministic planet atmosphere parameters derived from `PlanetType` + rough mass/radius.
  - Exponential density falloff: ρ(h) = ρ0 · exp(-h/H).
  - Drag based on **dynamic pressure**: q = ½ ρ v², with a = q · CdA / m.
  - Optional heating rate proportional to q (kPa) for a simple “reentry risk” signal.
- **Game integration:**
  - Toggle + tuning sliders in **World Visuals → Physics (experimental)**.
  - Live readout in **Ship / Status → Gravity / Orbit** (altitude, density, dynamic pressure, drag, heat).

- **Tests:** new `test_atmosphere` covers basic in-atmosphere drag directionality and vacuum cutoff.


## 2026-01-04 (Patch) - Lambert Transfer Planner (Star-Centric 2-Body)

This round adds an experimental **Lambert solver** (universal variables, 0-rev) and a first-pass
**transfer planner** in the Trajectory / Maneuver Planner UI.

The goal is a “Kerbal-ish” interplanetary helper: pick a planet/station target, choose a time-of-flight,
and automatically compute the **maneuver-node Δv** that reaches the target’s future position (star-centric).
Once computed, the existing **Maneuver Computer** can be used to execute the burn.

- **New core module:** `stellar/sim/LambertSolver.*`
  - Universal-variable Lambert solver using Stumpff functions + bisection on the z-parameter.
  - Supports short/long-way solutions and prograde/retrograde selection via a reference plane normal.
  - Units: km / km/s / km^3/s^2 (consistent with the gravity/orbit tooling added earlier).
- **Game integration (Ship / Status → Trajectory / Maneuver Planner):**
  - New **Transfer planner (Lambert)** section under Maneuver Node:
    - Requires a planet or station target.
    - Choose **TOF (hours)** + long-way + prograde.
    - Computes Δv in RTN (relative to the **star** at burn time) and auto-sets Reference Body to **Star**.
    - Optional **coarse RK4 validation** to estimate miss distance at arrival.

- **Tests:** new `test_lambert_solver` validates canonical circular-orbit 90° (short-way) and 270° (long-way) cases.


## 2026-01-04 (Patch) - Maneuver Computer (Auto-Execute Maneuver Node)

This round adds an experimental **maneuver computer** that can *arm* the current maneuver node from the Trajectory / Maneuver Planner
and automatically execute it as a **continuous burn** (a pragmatic approximation of the planner’s instantaneous Δv).

- **New core module:** `stellar/sim/ManeuverComputer.*`
  - Captures an absolute-time maneuver plan (node time + world-space Δv).
  - Centers the burn by starting ~½ the estimated duration early.
  - Uses a minimal attitude controller to align the ship’s **forward (+Z)** axis to the burn vector.
  - Throttle shaping on the last step to reduce Δv overshoot.
  - Tracks achieved Δv by projecting the ship’s *actual* velocity change onto the burn direction.
- **Game integration (Ship / Status → Trajectory / Maneuver Planner):**
  - New **Maneuver computer** section with:
    - **Arm / Abort / Clear** controls.
    - Live telemetry (time-to-node, alignment error, Δv remaining, burn timing).
    - Tuning sliders (alignment tolerance, face gain, lead time, Δv tolerance, abort timeout).
    - Optional manual override disengage (deadzone).
- **Tests:** new `test_maneuver_computer` covers centered-burn timing, Δv convergence, and miss/abort behavior.


## 2026-01-04 (Patch) - Trajectory Preview + 1-Node Maneuver Planner (RK4)

This round adds a deterministic **trajectory predictor** (RK4 integrator) and a first-pass **maneuver node** workflow.
It’s meant to be a “Kerbal-lite” planning layer for the existing spaceflight model: predict the next N minutes,
plan an instantaneous delta-v, and visualize the post-burn branch.

- **New core module:** `stellar/sim/TrajectoryPredictor.*`
  - Headless RK4 propagation in km / km/s.
  - Supports a single **instantaneous maneuver node** (exact-time event split).
- **Game integration (Ship / Status):**
  - New **Trajectory / Maneuver Planner (experimental)** panel:
    - Toggle world-space trajectory line.
    - Gravity model chooser: ballistic / effective (game settings) / physical (scale=1).
    - Reference body selector: auto dominant / star / specific planet.
    - Maneuver node inputs in **RTN frame** (Radial, Tangential/along-track, Normal).
    - “Compute: Circularize at node” helper to auto-fill delta-v.
- **World visualization:**
  - Cyan pre-burn path, magenta post-burn path.
  - RTN axis cross drawn at the maneuver node.
  - Min-altitude marker (amber / red when impact predicted).
- **Tests:** new `test_trajectory_predictor` covers ballistic motion, exact maneuver event sampling, and circular-orbit stability.

## 2026-01-04 (Patch) - Procedural Planet Rings + Blend-State Hygiene

This round adds a purely-visual **planet ring system** (procedural annulus textures + meshes) and fixes a subtle
render-state leak where additive blending could carry over into opaque passes.

- **New render generator:** `stellar/render/ProceduralRings.*`
  - Deterministic RGBA ring textures with banding, gaps, and azimuthal clumping.
  - LRU cache (`RingTextureCache`) with adjustable resolution (World Visuals).
- **New mesh helper:** `render::Mesh::makeRing(...)`
  - Builds a double-sided annulus in the XZ plane with UVs designed for the ring generator.
- **Game integration:**
  - Rings are generated per-planet (seeded, type-weighted probability) and aligned to orbit planes with optional tilt.
  - Toggle + tuning in **World Visuals → Secondary layers → Planet rings**.
  - Ring texture preview added to the **Surface generator preview** panel.
- **Render fix:** explicit `glDisable(GL_BLEND)` between transparent background layers and opaque geometry to prevent
  blend-state bleed into planet/cube passes.

## 2026-01-04 (Patch) - Experimental Newtonian Gravity + Live Orbit Diagnostics

This round introduces an **opt-in Newtonian gravity layer** for normal-space flight plus a lightweight
**two-body orbit solver** for debugging and future gameplay (safer supercruise drops, close approaches,
"periapsis warning" HUD, etc.). Gravity is disabled by default to preserve the existing arcade feel.

- **New core module:** `stellar/sim/Gravity.*`
  - Deterministic gravity acceleration from the system **star + planets** (km/s^2) with softening + optional clamp.
  - Convenience helpers for **dominant gravity body** selection (useful for orbit UI).
- **New core module:** `stellar/sim/OrbitalMechanics.*`
  - State-vector to orbit scalars: eccentricity, semi-major axis, periapsis/apoapsis, and period (when bound).
- **Ship physics upgrade:** `Ship::stepWithExternalAccel(...)`
  - Allows environmental accelerations (gravity now, more later) to be integrated consistently with sub-stepping.
- **Game integration:**
  - Toggle: **World Visuals → Physics (experimental) → Newtonian gravity**.
  - Applies gravity to **player + NPC contacts** in normal space.
  - Adds **Ship / Status → Gravity / Orbit** readout (dominant body, g-levels, peri/apo, period).
- **Tests:** new `test_gravity` covers gravity direction/magnitude and orbit-solver sanity (circular + elliptical).


## 2026-01-03 (Patch) - Normal-Space Nav Assist (Approach + Match Velocity)

This round adds a **normal-space navigation assist computer** that sits between raw manual thrust and the
full **DockingComputer**. It provides two practical modes aimed at moment-to-moment piloting:

- **Nav Assist: Approach** — closes to a configurable standoff distance and stabilizes there.
- **Nav Assist: Match Velocity** — matches a target's velocity and holds the current separation (great for rendezvous).

Highlights:

- **New core module:** `stellar/sim/NavAssistComputer.*`
  - Deterministic, headless, testable guidance built on the existing `FlightController` (including lead/intercept).
  - Two independent tuning profiles (approach vs. match velocity), plus optional auto-disengage on arrival.
- **Game integration:** new keybinds in `stellar_game`
  - **Nav Assist: Approach** (default: PageUp)
  - **Nav Assist: Match Velocity** (default: PageDown)
  - Safe behavior: disengages on strong manual input and automatically drops when leaving normal-space.
- **Tests:** new `test_nav_assist` validates approach convergence and moving-target velocity matching.


## 2026-01-03 (Patch) - Parallel Trade Run Planner + Sandbox Multi-threaded Scanning

This round upgrades the **multi-leg trade run planner** with a **parallel execution path** backed by the engine
`JobSystem`, making large station scans dramatically faster when `--threads` is enabled in `stellar_sandbox`.

- **New planner API:** `stellar::sim::planTradeRunsParallel(JobSystem&, ...)`
  - Snapshots station data + station economy states on the calling thread (avoids `Universe` cache thread-safety issues).
  - Computes per-stop outgoing leg manifests in parallel (cargo-manifest planner is the hot path).
  - Preserves deterministic ordering / tie-breaks so results match the serial planner.
- **Trade run routing micro-optimization:** if two systems are within `jumpRangeLy`, the planner skips A* and emits the
  direct 1-hop route (same result, less CPU).
- **Sandbox upgrade:** `stellar_sandbox --tradeRun` automatically uses the parallel planner when `--threads` creates a pool.
- **New tests:** `test_trade_runs_parallel` validates serial/parallel equivalence and determinism.


## 2026-01-02 (Patch) - Missile Countermeasures (Flares/Chaff) + Decoy-aware Seeker Guidance

This round adds a first-pass **countermeasure system** (flares/chaff) and upgrades missile guidance
so seeker missiles can be **fooled by decoy targets** exposed as `CombatTargetKind::Decoy`.

- **New sim module:** `stellar/sim/Countermeasures.*`
  - Deterministic flare/chaff burst spawner + lifetime integration.
  - Exposes decoys as `SphereTarget` entries with `decoyHeat` / `decoyRadar` signatures (strength decays with TTL).
- **Combat upgrade:** decoy-aware missile seekers
  - New `MissileSeekerType` + seeker params on `sim::Missile` (`seekerFovCos`, `decoyResistance`).
  - `stepMissiles()` compares locked-target score vs best decoy score and biases steering accordingly.
- **Game integration:** deploy flares in `stellar_game`
  - New keybind action: **Deploy Countermeasure** (default: Backspace).
  - Active countermeasures are appended to the projectile/missile target list so they can attract/hit.
- **Tests:** extended `test_combat` for countermeasure integration and decoy steering.


## 2026-01-01 (Patch) - Multi-Leg Trade Run Planner (Beam Search) + Sandbox --tradeRun

This round adds a **multi-leg trade run planner** that can suggest profitable chains like
`A -> B -> C -> ...` (not necessarily closed loops). It uses the existing cargo-manifest
planner per leg and (optionally) constrains each leg by jump-range reachability.

- **New core module:** `stellar/sim/TradeRunPlanner.*`
  - Beam-search expansion of the best outgoing trade legs at each stop (tunable width + fanout).
  - Supports ranking by **total profit**, **profit per ly**, **profit per hop**, or **profit per route cost**.
  - Optional jump-range reachability using the existing A* route planner (per-leg hop routes stored).
  - Deterministic ordering and tie-breaks for stable CLI + test behavior.
- **Sandbox upgrade:** `stellar_sandbox --tradeRun` prints (or emits JSON) for best trade runs from the chosen station.
  - Knobs: `--runLegs`, `--runLimit`, `--runBeam`, `--runLegCandidates`, `--runMinLegProfit`, `--runMinProfit`, `--runScore`.
  - By default, trade runs **do not** enforce `--jr` reachability unless `--jr` is explicitly provided.
- **New tests:** `test_trade_runs` forces a profitable 2-leg chain (Food A→B, Ore B→C) and checks determinism.

## 2026-01-01 (Patch) - K-Shortest Jump Routes (Yen) + Sandbox --kRoutes

This round upgrades the galaxy jump route planner to support **multiple alternative loopless routes**
(using Yen's K-shortest paths algorithm) rather than only the single best path.

- **New NavRoute API:** `plotKRoutesAStarCost(...)` and `plotKRoutesAStarHops(...)` with a simple `KRoute` result.
  - Returns up to **K** routes ordered by total cost (with deterministic tie-breaks).
  - Uses the same underlying A* cost model as the existing single-route planner.
- **Sandbox upgrade:** `stellar_sandbox --route --kRoutes <n>` prints (or emits JSON) for a ranked list of routes.
  - Works with `--routeCost hops|dist|fuel` so you can explore alternatives for different navigation styles.
- **New tests:** `test_k_routes` builds a small grid graph with 3 distinct loopless routes and checks ordering,
  determinism, and route metric consistency.

## 2026-01-01 (Patch) - Trade Loop Scanner + Sandbox --tradeLoop

This round adds a multi-leg trade loop scanner to find profitable closed routes that start and end at an origin station.
It expands the existing trade manifest planner into **2-leg** (A→B→A) and **3-leg** (A→B→C→A) loop search.

- **New core module:** `stellar/sim/TradeLoopScanner.*`
  - Builds per-leg cargo manifests via `econ::bestManifestForCargo(...)` (net-of-fees) without mutating economies.
  - Deterministic ordering with tie-breaks on system/station ids.
  - Tunable search fanout (`--loopLegCandidates`) and result cap (`--loopLimit`).
- **Sandbox upgrade:** `stellar_sandbox --tradeLoop` prints (or emits JSON) for the best loops from the chosen station.
  - New knobs: `--loopLegs`, `--loopLimit`, `--loopLegCandidates`, `--loopMinLegProfit`, `--loopMinProfit`.
- **New tests:** `test_trade_loops` forces a profitable round-trip (Food outbound, Ore inbound) and checks determinism.

## 2026-01-01 (Patch) - Headless Signals Planner + ResourceField Integration

This round makes the in-system 'signal site' layer first-class in the core sim and sandbox tooling.
The main goal is to have deterministic, test-covered generation of: resource fields, daily derelicts,
distress calls, and mission salvage sites.

- **New core module:** `stellar/sim/Signals.*` generates deterministic in-system signal sites.
  - Supports resource fields (persistent), daily derelicts (1-day TTL), distress calls (configurable TTL),
    and mission salvage sites derived from active `SaveGame` missions.
  - Adds lightweight JSON- and text-friendly `SignalSite` structs and `signalKindName(...)` helper.
- **New core helpers:** `stellar/sim/Units.h` centralizes AU↔km conversion and station/planet orbital helpers.
- **Resource fields now built + tested:** wired `src/sim/ResourceField.cpp` into the core library and
  enabled the existing `test_resource_field`.
- **Sandbox upgrade:** `stellar_sandbox --signals` prints (or emits JSON) for the generated sites, with
  knobs for field count and distress density.
- **New tests:** `test_signals` covers determinism, resolved derelict marking, and mission salvage site surfacing.

## 2026-01-01 (Patch) - Clamp Utilities + UI/String Safety Pass

This round hardens the MSVC build further by removing a few remaining "gotcha" patterns
that can trigger overload-resolution and pointer-arithmetic errors on Windows toolchains.

- **New core utility:** `stellar/core/Clamp.h` with `core::clamp(...)` and `core::clampCast(...)`.
  - Used to clamp nav route mode + commodity filter indices in `stellar_game`.
  - Used by SaveGame parsing for `navRouteMode`.
- **Navigation persistence:** replaced the last `std::clamp` usage in nav-route save/load with `core::clamp*`.
- **Quaternion/vector rotation:** switched one remaining `Quatd * Vec3d` call site to `Quatd::rotate(...)`
  for portability (and clearer intent).
- **Trade UI:** combo preview now uses a local `std::string` buffer (avoids `std::string_view` → `const char*`).
- **Trade toasts:** rebuilt buy toasts with explicit `std::string` assembly to avoid any `const char* + const char*`
  pitfalls on MSVC.

## 2026-01-01 (Patch) - Remaining MSVC Compile Fixes + SaveGame/Nav Sync

This round completes the set of MSVC compile fixes you reported and syncs code with the existing patch notes:

- **SaveGame nav persistence implemented**: `SaveGame` now actually contains the nav persistence fields (`navRoute`, hop, mode, auto-run, fuel-range constraint, pending arrival station) and `SaveGame.cpp` serializes/deserializes them (robust against stale route counts).
- **Tests updated**: `test_savegame` now exercises nav persistence and includes a stale-count robustness case for `navRoute`.
- **Math QoL**: added `Quatd * Vec3d` convenience operator (rotates the vector), matching usage in `stellar_game`.
- **UI string/pointer fixes**: resolved several `const char*` vs `std::string_view` / pointer-concatenation issues in `stellar_game` (commodity codes, trade error toasts, mission text).
- **Supercruise distress**: fixed a missing `jurisdiction` identifier by passing the system faction id directly into distress planning.
- **SDL/GL lifetime**: restored the intended GL-scope closing brace so GL resources are destroyed before `SDL_GL_DeleteContext`.


## 2026-01-01 (Patch) - stellar_game Build Fixes (ImGui 1.91.2 + MSVC)

This round focuses on resolving a large batch of compile errors in `apps/stellar_game`:

- **Event-loop scoping fix:** removed an extra closing brace in `main.cpp` that prematurely ended the SDL keydown
  handler and caused cascading scope/identifier errors (`key`, `dtReal`, `playerDamageContact`, etc.).
- **ImGui compatibility (v1.91.2 + `IMGUI_DISABLE_OBSOLETE_FUNCTIONS`):**
  - Updated renamed symbols:
    - `ImGuiCol_TabActive` → `ImGuiCol_TabSelected`
    - `ImGuiSelectableFlags_AllowItemOverlap` → `ImGuiSelectableFlags_AllowOverlap`
  - Switched the fetched ImGui tag to **`v1.91.2-docking`** so docking-only flags (e.g. `ImGuiWindowFlags_NoDocking`)
    are available.
- **Trade helper fixes:** aligned the Trade Helper UI with current `sim::TradeOpportunity` fields
  (`buyPrice`, `sellPrice`, `unitsPossible`, `unitMassKg`, `netProfitTotal`, etc.) and fixed the commodity filter flag name.
- **Bookmarks init order:** load bookmarks after `bookmarksPath`/`bookmarks` are declared.
- **Math:** added `double * Vec3d` operator overload.


## 2026-01-01 (Patch) - Navigation Persistence + Smarter Auto-Run

This round makes long-haul travel feel better by persisting plotted routes and route planner settings
through save/load, and by making auto-run more resilient.

- **Save/Load**:
  - SaveGame version bumped to **21**.
  - Persisted nav route, route cursor, auto-run toggle, route mode, and fuel-range constraint.
  - Persisted pending arrival station helper (auto-target on arrival).
- **Game**:
  - Quickload restores navigation state (with validation/repair against current system).
  - Galaxy map defaults selection to the **next hop** after load when a route is active.
  - Auto-run now:
    - waits briefly after load (avoids surprise immediate jumps),
    - attempts a one-shot replot if the next hop is out of range (e.g., cargo/range changed),
    - stops and targets the nearest station if fuel is insufficient.
- **Tests**:
  - Expanded `test_savegame` to cover nav persistence and robustness against stale `nav_route` counts.

## 2025-12-31 (Patch) - Distress Encounters 2.0 (Rescue Requests + Deterministic Plans)

This round expands supercruise distress signals into an interactive rescue loop:
scanning a distressed ship can now transfer requested supplies and pay out a small
reward (+ optional rep), while still preserving the chance of pirate ambush.

- **New core module**: `stellar::sim::Distress` (`include/stellar/sim/Distress.h` + `src/sim/Distress.cpp`)
  - Deterministic distress plan generation (victim request, reward, ambush params).
- **Game**:
  - Distress signals now generate a stable plan at spawn time (so scans show consistent info).
  - Distress sites can spawn a **"Distress Vessel"** contact that requests a commodity.
  - Completing a normal contact scan within 15,000 km transfers cargo toward the request.
  - Upon completion: grants credits + optional faction reputation and marks the signal "completed" for UI.
  - Target HUD now shows additional detail for distress signals and distress victims.
- **Tests**: added `test_distress`.

## 2025-12-30 (Patch) - Ballistics Lead Solver + NPC Weapons/Loadouts (Combat Integration)

This round pushes combat toward a more unified, systems-driven model: NPCs now use real weapons (beams/projectiles), loadout-derived stats, and the same power distributor mechanics as the player.

### New: `sim::Ballistics`
- Added `stellar::sim::Ballistics` (`include/stellar/sim/Ballistics.h` + `src/sim/Ballistics.cpp`):
  - Robust intercept-time solver for constant-speed projectiles that inherit shooter velocity.
  - Helper `solveProjectileLead(...)` that returns both **time-to-hit** and the **lead point**.
- Added unit tests: `test_ballistics`.

### Game: HUD lead indicator uses shared ballistics
- Replaced the HUD lead quadratic in `stellar_game` with `sim::solveProjectileLead(...)`.
  - Keeps the HUD consistent with how projectiles are actually simulated.

### Game: NPC combat now uses shared Combat + Loadout + Distributor
- Contacts gained a lightweight “combat loadout”:
  - `ShipHullClass`, module Mks (thrusters/shields/distributor), weapon type, AI skill.
  - Per-contact distributor config/state + pips.
- NPC shield regen is now:
  - Driven by **loadout-derived regen rate**.
  - Scaled by **SYS pips**.
  - Limited by and consumes **SYS capacitor**, like the player.
- Pirate + police firing now routes through `sim::tryFireWeapon(...)`:
  - Beams raycast and apply immediate damage via shared `applyDamage(...)`.
  - Projectile weapons spawn `sim::Projectile` instances (with proper lead aiming via `sim::Ballistics`).
- Projectile hit resolution now supports **NPC → NPC ship** hits (not just NPC → player).

## 2025-12-30 (Patch) - Power distributor (Pips + Capacitors)

This round adds a 3-channel power distributor system (ENG/WEP/SYS) inspired by classic space sims.

- **New core module**: `stellar::sim::PowerDistributor` (pips normalization, capacitor recharge, boost consumption, weapon capacitor cost helper).
- **Gameplay integration** (player ship):
  - Boost now consumes **ENG capacitor** (partial boost supported per frame when energy is low).
  - Weapon firing now requires and consumes **WEP capacitor**.
  - Shield regen is now limited by and consumes **SYS capacitor**, and scales with SYS pips.
- **HUD**: optional capacitor rings around the reticle (toggled by `hudShowDistributorRings`).
- **Ship/Status UI**: new Power distributor section to view capacitors + adjust pips.
- **Save game**: pips + capacitor fractions are persisted (SaveGame version bumped to 17).
- **Tests**: added unit tests for normalization, recharge redistribution, boost consumption, and weapon cost heuristics.


## 2025-12-30 (Patch) - Combat (Core) + Weapon/Projectile Extraction


- **Sim:** added `sim::Combat` (`stellar/sim/Combat.h` + `src/sim/Combat.cpp`) to move combat math out of `main.cpp`:
  - Shield-first `applyDamage(...)`
  - Ray-sphere intersection + nearest-hit raycast with optional aim-cone filter
  - Projectile stepping with segment collision checks + hit event emission
  - A shared `tryFireWeapon(...)` helper that spawns beams/projectiles and returns hit metadata
- **Game:** refactored player firing + projectile updates to use the shared combat module.
  - Kinetic projectiles now also collide with asteroids (sparks only; mining still requires the mining laser).
- **Tests:** added `test_combat` to lock down damage, ray hits, nearest-target selection, and projectile hit emission.

## 2025-12-30 (Patch) - ShipLoadout (Core) + Balance Extraction
- **Sim:** added `sim::ShipLoadout` (`stellar/sim/ShipLoadout.h`) to centralize gameplay tuning tables:
  - Hull definitions (Scout / Hauler / Fighter)
  - Mk module definitions (Thrusters / Shields / Distributors)
  - Weapon definitions (beam/pulse/cannon/rail/mining)
- **Sim:** added a shared derived-stat helper `computeShipDerivedStats(...)` so the ship stat formula lives in one place.
- **Game:** refactored `stellar_game` to use the shared definitions (removes a large block of balance data from `main.cpp`).
- **Tests:** added `test_ship_loadout` to validate monotonic upgrades (Mk3 > Mk1) and lock down the derived-stat formula.

## 2025-12-30 (Patch) - Deterministic Per-Faction Law Profiles (Smuggling)
- **Sim:** added `sim::LawProfile` (`sim/Law.h` + `sim/Law.cpp`): deterministic per-faction tuning for:
  - Cargo-scan *strictness* (scan frequency multiplier)
  - Enforcement *fine schedule* (base + rate)
  - *Corruption* scalar (affects bribe offer chance)
  - Rep penalties for contraband compliance vs. evasion
- **Game:** wired the profile into the smuggling loop:
  - Contraband fine + rep penalty now come from the faction's law profile.
  - Police bribe chance is multiplied by faction corruption.
  - Cargo-scan start rate is multiplied by faction scan strictness.
- **Tests:** added `test_law` to lock down determinism, bounds, and monotonicity.

## 2025-12-30 (Patch) - FlightController (Core) + Autopilot/NPC AI Refactor
- **Sim:** added `sim::FlightController` (new module) to centralize "approach/chase + face direction" logic:
  - Produces a normalized `sim::ShipInput` from a moving target point + velocity.
  - Conservative speed profiling (distance gain + stopping-distance clamp) for smoother arrivals.
  - Optional boost selection (disabled by default) and optional roll alignment to an up vector.
- **Game:** migrated two in-game controllers onto the shared module:
  - Player station **Autopilot** now uses `sim::approachTarget(...)` (better damping near the slot).
  - NPC "chaseTarget" behavior is now powered by the same controller (consistent feel across roles).
- **Tests:** added `test_flight_controller` to validate convergence + input range safety.

## 2025-12-30 (Patch) - Ship Boost Caps + stellar_game Compile Fixes (Convoys)
- **Sim:** upgraded `sim::Ship` to support configurable **boost acceleration caps**:
  - New APIs: `setMaxLinearAccelBoostKmS2(...)` / `setMaxAngularAccelBoostRadS2(...)` plus getters.
  - Defaults preserve legacy behavior (boost = **1.8x linear**, **1.4x angular**) unless explicitly overridden.
- **Game:** fixed several `stellar_game` build breaks in convoy/escort code caused by stale Ship/Contact fields:
  - Updated convoy + escort spawns to use `Ship` setters (no more direct member access).
  - Added `Contact::tradeCargoValueCr` and kept it in sync for normal traders (used for ambush scaling).
- **Tests:** added `test_ship` to lock down boost-cap behavior.

## 2025-12-29 (Patch) - Universe Nearby Query Accelerator + MSVC .c_str() Fix
- **Fix (MSVC):** removed invalid `.c_str()` calls on `econ::CommodityDef::name` in the target info panel (`CommodityDef::name` is `const char*`).
- **Sim/Perf:** rewrote `Universe::queryNearby(...)` to avoid scanning the entire sector bounding box when `maxResults` is small:
  - Best-first sector expansion using a safe AABB distance lower bound.
  - Early-out once no remaining sector can beat the current worst of the kept top-N.
  - Preserves deterministic ordering (distance, then `SystemId`).
- **Tests:** added `test_query_nearby` which cross-checks `Universe::queryNearby(...)` against a brute-force sector scan for multiple cases.

## 2025-12-29 (Patch) - Build Fixes (MSVC)
- **Fix:** corrected an invalid numeric literal in `ProceduralLivery.cpp` (previously `0xDECALull`), replacing it with a stable hashed seed salt (`hashCombine(seed, fnv1a64("decal"))`).
- **Fix:** removed a duplicate `lightPos_` member declaration in `MeshRenderer.h` that caused MSVC `C2086` redefinition errors.
- **Fix:** added missing `stellar/math/Vec2.h` include in `apps/stellar_game/main.cpp` so `math::Vec2d` is defined for System Map 2.0.
- **Fix:** replaced invalid `ImGui::GetMouseDelta(...)` call with `ImGui::GetMouseDragDelta(...)` in the Hangar preview (ImGui 1.91+).
- **Fix:** removed dependency on `IM_PI` (not exposed in this ImGui build) and used `math::kPi` for HUD ring math.
- **Fix:** corrected `std::string_view` UI usage (commodity name/code) by using span-safe ImGui calls (`%.*s` + `TextUnformatted(text, end)`).
- **Fix:** fixed Galaxy list jump-range scope error (`jrMaxLy` → `galaxyPlanMaxLy`).

## 2025-12-29 (Patch) - Hangar 3D Preview + Procedural Ship Livery System
- **UI/Render:** added a new **Hangar / Livery** window (Windows → Hangar / Livery):
  - A real-time **3D ship preview** rendered into an offscreen OpenGL **FBO render target** and shown in the UI via `ImGui::Image`.
  - Mouse controls: **drag to rotate**, **scroll to zoom**, optional **auto-spin**.
- **Render:** introduced `render::RenderTarget2D` (RGBA8 + depth) for general-purpose **render-to-texture previews**.
- **Render:** added `render::ProceduralLivery` – a new procedural **UV texture** generator for ships/stations with multiple patterns:
  - Solid, Stripes, Hazard (diagonal warning stripes), Camo, Hex, and Digital.
  - Optional **decal** (band + dot), plus **wear/dirt** and subtle fabric-like noise.
- **World:** the player ship now renders in a separate cube pass with the generated livery texture (toggleable: “Apply in world”).
- **Render:** `MeshRenderer` now supports a configurable point-light position (`setLightPos`) so world lighting stays “star at origin” while UI previews can use a nicer offset light.
- **UX:** livery settings persist to `livery.txt` with **Save / Load / Reset** plus optional **auto-save on close**.

## 2025-12-28 (Patch) - System Map 2.0 (Orrery) + Time Scrub + Atlas Icons
- **UI:** rebuilt the Scanner → **System Map** into an **infinite-canvas "orrery"**:
  - Cursor-anchored **mouse wheel zoom**, RMB/MMB **pan**, **follow ship**, and a right-click **context menu**.
  - Added a subtle **grid/scale** overlay (concentric AU rings + scale bar) to make distances readable at any zoom level.
- **UI/Render:** System Map now draws **all icons from the HUD `SpriteAtlas`**, keeping the procedural look consistent and avoiding a pile of texture binds.
  - Selected objects get a pulsing **selection ring**.
  - The current nav target also shows a **supercruise drop-distance ring** (station/planet/signal) to help approach planning.
- **UI:** added an **orbit preview time scrub** (**t ± days**) plus optional animation (days/sec) so you can forecast planetary/station motion **without changing sim time**.
- **UI:** added a right-side **Selection/Quick Targets** panel with:
  - Per-target info + actions: **Scan (K)**, **Supercruise (H)**, **Autopilot to station (P)**, **Center map on selection**.
  - Quick target lists for **Planets / Stations / Signals** (atlas icons included).

## 2025-12-28 (Patch) - Planet Atmospheres (Limb Glow) + Procedural Cloud Shells
- **Render:** added `render::AtmosphereRenderer` (additive limb glow) for lightweight **planet atmospheres**.
  - Fresnel-ish rim term + optional **day-side boost** and **forward-scatter** highlight when looking toward the star.
  - Rendered as an instanced, slightly larger UV-sphere shell with additive blending (no gameplay impact).
- **Render:** added a new procedural surface kind: **Clouds** (`render::SurfaceKind::Clouds`) generating a seam-free alpha mask on the unit sphere.
- **Render:** `MeshRenderer` now optionally outputs **alpha from texture** (`setAlphaFromTexture`, `setAlphaMul`) enabling translucent shells (cloud layer) without a new mesh pipeline.
- **World:** planets now render optional **cloud layers** (alpha blended) and optional **atmospheres** (additive) on top of their surface.
  - Per-planet tint/strength is derived from planet type + seed; clouds also get a subtle deterministic **rotation**.
- **UI:** **World Visuals** window gains a new **Secondary layers** section for tuning:
  - Atmosphere: intensity, rim power, shell scale, day-side boost, forward scatter, tint-with-star.
  - Clouds: opacity, shell scale, spin rate (requires procedural surfaces).
- **UI:** Scanner planet tooltips can now show a **Clouds** preview below the surface preview.

## 2025-12-28 (Patch) - HUD Layout Editor (Drag/Drop) + Persistent Overlay Placement
- **UI:** added a **HUD Layout system** (`hud_layout.txt`) that persists normalized overlay placement (pivot-based) plus per-widget scale and enable flags.
- **UI:** added a new **HUD Layout** window (HUD → Layout) with:
  - **Edit mode** (drag HUD panels to reposition)
  - Per-widget **enabled** toggles + **scale** sliders + **anchor/pivot** presets
  - Quick actions: **Save/Load/Reset**
  - Keyboard shortcuts: **Ctrl+H** edit mode, **Ctrl+S** save, **Ctrl+L** load, **Ctrl+R** reset
- **UI:** while in edit mode, the HUD draws a **safe-area guide** and **pivot markers** for HUD widgets.
- **UI:** Radar / Objective / Threat / Jump HUD overlays now use the persisted layout for position + scale.

## 2025-12-28 (Patch) - Hyperspace Jump Tunnel (PostFX) + Jump HUD
- **Render/PostFX:** added a new **hyperspace / jump tunnel** screen-space effect (swirl + rings/sparkles) composited in the PostFX pass.
  - Controlled via new `render::PostFXSettings` parameters: `hyperspace`, `hyperspaceTwist`, `hyperspaceDensity`, `hyperspaceNoise`, `hyperspaceIntensity`.
  - Designed to play nicely with existing bloom/tonemap (the effect is added **pre-tonemap**).
- **Gameplay/UI:** **FSD** now supports a clean **abort during charge**: press **J** again while charging to cancel (fuel is only consumed when the charge completes).
- **UI:** added a lightweight top-center **Jump HUD overlay** (toggleable) showing destination + charge/hyperspace progress + ETA.
- **UI:** Post FX panel gains a new **Hyperspace / Jump FX** section with auto-FSD controls and manual tuning sliders.

## 2025-12-28 (Patch) - Procedural Nebula Background Layer (Point Sprites)
- **Render:** added `render::NebulaField`, a background **nebula/gas layer** built from large **additive point sprites**.
  - Optional **galactic-plane banding** ("band power"), mild **turbulence**, and tunable **parallax** (0..1) for depth.
  - HDR-friendly **intensity** control (lets bloom/PostFX give the nebula a soft glow).
- **Render:** added a new procedural **cloud sprite** generator (`makeCloudSpriteTextureRGBA`) for soft, noisy alpha puffs.
- **UI:** **VFX Lab** gains a new **Nebula** section: density, variant, banding, radii, parallax, intensity, opacity, size range, and turbulence.

## 2025-12-28 (Patch) - Procedural Planet Surfaces + Star-Relative Lighting
- **Render:** added a new procedural **surface texture** generator for spherical bodies (`render::ProceduralPlanet`) plus an LRU-ish **runtime GL texture cache** (`render::SurfaceTextureCache`).
  - Generates 2:1 **equirectangular** albedo textures for **Rocky / Desert / Ocean / Ice / Gas Giant / Star** surfaces.
  - Deterministic from a seed; safe to call from the main thread and cached to avoid re-uploading.
- **Render:** upgraded `MeshRenderer` lighting to use a simple **star-at-origin point light** (per-fragment) so planets/stations get a more believable terminator.
  - Added `MeshRenderer::setUnlit(bool)` to render emissive/unlit bodies (stars, debug).
- **World:** star + planets now optionally render with procedural surface textures (toggleable).
  - Added a high-level HDR knob: **star intensity** to let the star drive bloom in the PostFX pass.
- **UI:** new **World Visuals** window (Windows → World Visuals) to:
  - Toggle procedural surfaces + UI previews.
  - Adjust surface texture resolution.
  - Preview the surface generator output live.
- **UI:** Scanner tooltips for **Star** and **Planets** optionally show a surface preview (matches the world seed).

## 2025-12-28 (Patch) - Combat HUD Suite (Procedural Reticle + Lead + Flight Path Marker)
- **HUD → Combat:** added a new **Combat HUD** package with high-visibility flight/combat symbology:
  - **Procedural reticle** (`HudReticle`) drawn from the HUD atlas (tintable alpha-mask sprite).
  - **Weapon cooldown rings** around the reticle (primary/secondary, tinted by weapon color).
  - **Heat ring** around the reticle (grows with ship heat).
  - **Flight path marker** (`HudVelocity`) that shows your current velocity vector (optionally relative to the **local frame**).
  - **Projectile lead indicator** (`HudLead`) for contact targets using an analytic moving-target intercept; hold **Shift** to show time-to-impact.
- **Render:** expanded procedural sprite generation with new HUD kinds: **HudReticle**, **HudLead**, **HudVelocity**.
- **UX:** firing now tracks the **last fired** weapon so the lead marker matches the weapon you’re actually using.

## 2025-12-28 (Patch) - Galaxy Map 2.0 (Pan/Zoom + Atlas Icons)
- **UI:** revamped the Galaxy window into a **two-column layout**: interactive map + system list on the left, selection + route planner on the right.
- **Map:** upgraded the mini-map to an **infinite-canvas** style viewer:
  - **Mouse wheel** zoom (keeps the world point under the cursor stable).
  - **MMB/RMB drag** panning.
  - **Double-click** a system to center the view.
  - **Right-click** context menu (select, center view, plot route, engage FSD).
- **Map visuals:** systems are now rendered with the **HUD procedural atlas** (star icon tinted by **star class**), with optional **faction glyph** overlays, plus a pulsing **next-hop** highlight when a route is plotted.

## 2025-12-28 (Patch) - HUD Icon Atlas + Offscreen Target Indicator
- **Render:** added a lightweight **procedural icon atlas** (`render::SpriteAtlas`) that packs many `(kind, seed)` sprites into **one GL texture**.
  - Incremental updates via `glTexSubImage2D` (loader support added) so new icons can be uploaded without recreating the whole texture.
- **HUD:** Radar + Tactical overlay + target label icons now draw from the **atlas** (fewer texture binds, fewer per-size sprite-cache entries).
- **HUD:** added an optional **offscreen target indicator** (edge-of-screen arrow + icon + distance). Toggle it from **HUD → Offscreen target indicator**.
- **Sprite Lab:** displays the current HUD atlas (size/capacity) and provides a **Clear HUD atlas** button.

## 2025-12-28 (Patch) - Tactical Overlay + Textured Point Sprites
- **UI/HUD:** added a **Tactical overlay** (**`**) that projects in-world objects onto the screen with **procedural icons**.
  - Hover markers for a tooltip; **MMB** targets the hovered marker (avoids conflicting with weapon fire).
  - Range / marker cap / per-type filters + label toggle (HUD menu).
- **Render:** upgraded `PointRenderer` with a **textured point-sprite** path (samples a small RGBA sprite via `gl_PointCoord`).
- **VFX:** starfield + particle passes can now switch between **circle points** and **textured sprites** (VFX Lab toggles).
- **Art:** added small **procedural radial sprites** for stars + particles (deterministic from the sim seed).

## 2025-12-28 (Patch) - HUD Radar + System Map + Expanded Procedural Icons
- **UI/HUD:** added a **Radar** overlay (**R**) with clickable blips (target stations/planets/contacts/signals/cargo/asteroids).
- **UI/Scanner:** added a new **System Map** tab to the scanner (**F6**) featuring a 2D orbit view with **click-to-target** icons and zoom/filter toggles.
- **Render:** expanded procedural icon generation with new sprite kinds: **Cargo**, **Asteroid**, **Signal**.
- **UI:** scanner lists for signal sources, asteroids, and floating cargo now include procedural icons + hover previews.
- **UI/HUD:** target marker now shows the target’s procedural icon beside the distance label.

## 2025-12-28 (Patch) - HDR PostFX (bloom + tonemap)
- **Render:** added an HDR offscreen render target and a lightweight **PostFX** pipeline (bright-pass + ping-pong Gaussian blur + composite).
- **Render:** added film-ish finishing controls (exposure/gamma, vignette, subtle grain, chromatic aberration) and an optional screen-space **warp streak** effect.
- **UI:** new **Post FX** tuning window (**F12**) to live-adjust bloom + tonemap.
- **Fix:** scoped GL-backed objects so they destruct **before** `SDL_GL_DeleteContext` (prevents GL deletes after context teardown).

## 2025-12-28 (Patch) - Starfield + particle VFX pass
- **Render:** added a deterministic procedural **Starfield** background (point sprites) that follows the camera and twinkles.
- **Render:** added a lightweight CPU **ParticleSystem** (point sprites) for **thrusters**, **impact sparks**, and **explosions**.
- **Render:** `PointRenderer` now supports **alpha** per vertex and an **additive blend mode** for glowy VFX.
- **Game:** VFX hooks added for thruster plume, weapon impacts, contact destruction, and player death.
- **UI:** new **VFX Lab** window (**F11**) for toggling starfield/particles and tuning basic counts/intensity.

## 2025-12-28 (Patch) - Procedural UI sprite generation + icon-rich panels
- **Render:** added a small **procedural sprite generator** that can produce deterministic **RGBA icon sprites** (Commodity / Faction / Mission / Station / Planet / Star) from a `(kind, seed)` pair.
- **Render:** added `Texture2D::createRGBA(...)` for uploading raw RGBA8 pixel buffers to OpenGL textures (used by the sprite generator + future procedural art).
- **UI:** Market, Missions, Contacts, and System Scanner now show **procedural icons** (with hover tooltips that preview larger variants).
- **UI:** new **Sprite Lab** window (**F10**) to preview sprite kinds, tweak seed/size, browse a commodity icon atlas, and inspect system icon samples.

## 2025-12-27 (Patch) - Ship-fit mission boards + route-aware multi-hop tuning
- **Mission Board:** cargo-bearing offers (Delivery / Multi-hop / Smuggle) now size their payload to the player's **current free cargo capacity**, so generated jobs are accept-able immediately without dumping cargo.
- **Passengers:** passenger offers now respect available **seat capacity** (no more 6-person offers on a 2-seat ship).
- **Multi-hop routing:** Multi-delivery via stops are chosen to avoid **extreme detours**, and payouts/deadlines now scale with the **full route distance** (origin→via→destination).
- **Salvage:** salvage offer quantities are capped to ensure the requested goods can **fit in your hold** at completion time.
- **Tests:** expanded `test_missions` with ship-capacity sanity checks for mission-board generation.

## 2025-12-27 (Patch) - Economy-aware delivery missions (profit + legality fixes)
- **Mission Board (core):** Delivery and Multi-delivery cargo is now picked using **origin/destination market prices** (demand-aware), instead of being fully random.
- **Legality:** "legal" delivery missions will no longer ask you to deliver **contraband** into a destination jurisdiction where it's illegal (those are now reserved for **Smuggle** jobs).
- **Rewards:** Delivery payouts now scale with the cargo's expected **acquisition cost/value**, preventing negative-profit jobs for high-value commodities.
- **Smuggling:** Smuggle offers now only appear when the destination faction actually has contraband rules; otherwise the offer downgrades to a normal courier job.
- **Tests:** expanded `test_missions` with lightweight checks for delivery legality/profitability and smuggling legality.

## 2025-12-27 (Patch) - Salvage recovery missions + mission derelict sites
- **Missions:** added a new **Salvage** mission type.
  - Salvage jobs spawn a **mission derelict** signal in-system; recover the requested goods and return to the station for payout.
  - Completion requires that you have **visited the site** (tracked via `m.scanned`) so you can't instantly complete from dock.
- **Signals + Scanner:** mission salvage signals and cargo pods are tagged **[MISSION]** in the scanner UI.
- **Objective + Tracker:** Objective HUD and the Ship/Status mission tracker can now point you at the **mission site** (`Target site`).
- **Core (MissionLogic):** Mission Board can generate Salvage offers and dock completion now supports Salvage.

## 2025-12-27 (Patch) - Pirate extortion / tribute (combat + cargo loop)
- **Pirates:** pirate packs may now try **extortion** before opening fire.
  - A **Threat HUD** appears showing the demanded cargo **value** and **time remaining**.
  - Use **Ship/Status → Cargo management → Jettison** to drop cargo pods; jettisoned cargo counts toward the demand by **base value**.
  - If you satisfy the demand, pirates **disengage and flee** (tribute pods are removed to simulate them scooping it).
  - If you ignore the demand or **attack** the pirates, they immediately **turn hostile**.
- **Contacts UI:** pirates involved in an active demand are tagged **[DEMAND]**.

## 2025-12-27 (Patch) - Cargo jettison + contraband dump (smuggling QoL)
- **Ship / Status:** added a **Cargo management** panel to **jettison cargo** (spawns floating pods behind your ship).
- **Smuggling QoL:** added a one-click **Dump contraband** "panic button" to jettison all illegal goods in the current jurisdiction.
- **Law:** dumping cargo near authorities (station comms range, during a scan, or near police) is treated as a **crime** and will generate a **bounty** + reputation hit.
- **Safety:** mission cargo is reserved by default to prevent accidental mission sabotage (toggle: **Allow dumping mission cargo**).

## 2025-12-27 (Patch) - Objective HUD overlay + Mission Board route previews
- **stellar_game:** added an in-flight **Objective HUD overlay** (top-right) that shows your **tracked mission**, next stop hint, and (if present) a quick **route remaining** summary.
- **Ship/Status:** added a toggle for the Objective HUD overlay.
- **Missions (Active):** added a **Plot route** button to quickly plot to the current mission leg and open the Galaxy map.
- **Mission Board:** optional **route preview** line per offer (A*), showing estimated **jumps / distance / fuel** to the *next stop* under your current route-planner settings.

## 2025-12-27 (Patch) - Mission tracker HUD + auto-plot next leg + save persistence
- **stellar_game:** added a **Mission tracker** section to **Ship / Status** showing the tracked mission with quick actions: **Select**, **Plot Route**, **Target** (if in-system), and **Untrack**.
- **stellar_game:** optional **Auto-plot next leg**: when a *tracked* **Multi-delivery** completes leg 1/2, the nav route is automatically re-plotted to the final destination.
- **Missions UI:** Active missions list now supports **Track/Untrack** and shows a **[TRACKED]** tag.
- **SaveGame:** persisted `trackedMissionId` (version bumped to **13**) so the tracker survives save/load.


## 2025-12-27 (Patch) - Shared mission completion + bounty event helpers + Galaxy system list (QoL)
- **Core (MissionLogic):** added `tryCompleteBountyScan(...)` and `tryCompleteBountyKill(...)` so bounty missions can be completed via shared, testable logic.
- **stellar_game:** mission deadline ticking + docked completion now delegate to `sim::MissionLogic` (so prototype behavior matches tooling/tests and avoids duplicated rules).
- **Tests:** expanded `test_missions` to cover deadline expiry, multi-hop delivery leg progression, and bounty scan/kill completions.
- **Galaxy UI:** added a searchable/sortable **System list** under the map with an **IN/OUT of jump range** indicator for faster selection.

## 2025-12-27 (Patch) - In-game route planner modes + safer auto-run (stellar_game)
- **Galaxy UI:** route plotting now supports **Min jumps**, **Min distance**, or **Min fuel** (A*).
- **Galaxy UI:** added an option to constrain route edges to **current-fuel jump range** (vs max range) for "what can I do right now?" planning.
- **Galaxy UI:** jump range overlay now shows both **max range** and **current-fuel range** rings.
- **QoL:** auto-run route now **stops cleanly** (with a single toast) if the next hop is out of range or fuel is insufficient.

## 2025-12-27 (Patch) - Cost-aware route planning (core + stellar_sandbox)
- **Core:** added `sim::plotRouteAStarCost(...)` to plan routes using a weighted cost model (per jump + per ly).
- **Tooling:** `stellar_sandbox --route` now supports `--routeCost hops|dist|fuel` (plus `--fuelBase` / `--fuelPerLy`) and reports total cost + fuel estimate.


## 2025-12-27 (Patch) - Alert decay + Lay Low station option
- **Bugfix/behavior:** local security `policeHeat` now decays over time (exponential decay), so "Alert" cools down naturally as time passes.
- **New dock service:** **Lay low (3h)** in **Security Services**: pay a fee to advance time by 3 hours and reduce local alert/response.

# Patch Notes (Dec 27, 2025 - Police Alert Fix + Port Security Services)

- **Bugfix:** contraband enforcement, police bribes, and black-market contraband sales now raise **local security alert** (`policeHeat`) instead of ship thermal heat.
- **UI clarity:** status labels now distinguish **Alert** (local security) from **Ship heat** (thermal/overheat).
- **New station option:** **Bribe port authority (reduce alert)** lowers local alert for a fee (does not clear bounties).
- **Black market consequence:** selling illegal goods applies a small reputation penalty and raises alert over time.

---

# Patch Notes (Dec 27, 2025 - New Commodities + Save Compatibility)

## Economy: 2 new commodities (new)
- Added two higher-value commodities:
  - **Weapons** (`ARMS`)
  - **Stimulants** (`STIM`)
- Updated station economy models so the new goods are actually produced/consumed:
  - Agricultural: produces small amounts of Stimulants
  - Industrial: produces Weapons
  - TradeHub/Research/Shipyard: increased demand for these goods

## Smuggling / contraband: more variety (upgrade)
- Contraband rules now include the new goods.
  - Factions may ban **Stimulants** as a vice.
  - Factions may ban **Weapons** as controlled tech.

## SaveGame: commodity-count forward compatibility (important)
- SaveGame version bumped to **12**.
- **Older saves remain loadable** even when commodity counts change:
  - `cargo` lines and station override `inventory` lines now parse the remainder of their line and
    default missing commodity slots to **0**.

---

# Patch Notes (Dec 27, 2025 - Bribes & Contraband Scans)

## New: Bribe & comply window for contraband scans (police)
- If a **police** cargo scan completes and detects contraband, they may offer a short "bribe or comply" window.
  - **C**: Pay bribe to keep your cargo.
  - **I**: Comply: cargo is confiscated and you are fined (unpaid fines become bounty).
  - **Run**: Leave scan range to keep cargo but receive a bounty and reputation hit.
- The window is visible in **Ship / Status** with a countdown timer.
- Police ships will generally **hold fire** during the bribe/compliance window (unless already hostile).

---

## Shipyard: Smuggling compartments (new)
- Added **Smuggling compartments Mk1–Mk3** as a Shipyard upgrade.
  - Reduces the *chance* of cargo scans when you are carrying contraband.
  - Makes scans take slightly longer (more time to break range if you try to run).
- Persisted via SaveGame (`smuggleHoldMk`). Save version bumped to **11** (older saves default to Mk0).

## Smuggling: meaningful failure case (new)
- If a cargo scan **confiscates contraband** that is tied to an active **Smuggle** mission, that mission now **fails immediately**.
  - Applies the usual mission-failure reputation hit with the issuing faction.

## Missions UI QoL
- Active cargo missions now show a quick **Cargo: have / need** line (Delivery / Multi-hop / Smuggle).

---

# Patch Notes (Dec 27, 2025 - Hotfix)

## Smuggle missions: completion fix (SDL prototype)
- Fixed the SDL prototype's dock-completion logic to include **Smuggle** missions.
  - Completing a Smuggle mission now consumes cargo and pays out as expected.
  - Smuggle deliveries still **do not** feed the destination's *legal* market inventory.

## Market QoL: reserve mission cargo
- The Market sell flow now reserves cargo required by active **Delivery**, **Multi-hop**, and **Smuggle** missions.
  - The Cargo column shows reserved units.
  - If all units are reserved, the Sell button is disabled with a tooltip explaining why.

---

# Patch Notes (Dec 27, 2025)

## Missions: Smuggling (new)
- Added **Smuggle** missions (new `MissionType::Smuggle`).
  - Smuggling jobs grant the contraband cargo from the contact on acceptance (no legal market inventory/credit checks).
  - Completion consumes the cargo and pays out normally, but **does not** feed the destination's legal market inventory.

## Simulation: shared contraband rules (refactor)
- Introduced `sim::Contraband` (`stellar/sim/Contraband.h`) to centralize deterministic contraband rules:
  - `illegalCommodityMask`, `isIllegalCommodity`, `illegalCommodityListString`, `hasIllegalCargo`
- The SDL prototype now uses the shared contraband helpers so illegality rules stay consistent across UI + headless modules.

## Tooling + tests
- `stellar_sandbox` now recognizes the **Smuggle** mission type name.
- Added a unit test to cover Smuggle mission acceptance + completion invariants.

---

# Patch Notes (Dec 26, 2025)

## Tooling: stellar_sandbox JSON + mission board (new)
- `stellar_sandbox` gained a small CLI parser and JSON emitter for scripting/balancing:
  - `--json` (stdout) and `--out <path>` (file) outputs nearby systems + optional trade ideas
- New headless mission workflow helpers:
  - `--missions` to print deterministic offers
  - `--acceptOffer <idx>` to simulate acceptance with cargo/credits/seats constraints
  - `--autoComplete` (+ `--advanceDays`) to test completion logic without the renderer
  - `--load/--save` to round-trip SaveGame files from tooling

## Navigation: jump route planner (new)
- Moved the game's A* hop-count jump planner into core: `sim::plotRouteAStarHops` (`stellar/sim/NavRoute.h`).
- `stellar_sandbox` gained `--route` mode to plot jump paths headlessly (supports `--json` output).
- Added unit test `test_nav` to cover route finding, unreachable cases, and determinism.

## Core utilities (new)
- Added `stellar/core/Args.h` (tiny arg parser) and `stellar/core/JsonWriter.h` (minimal JSON writer)


## Encounters: squads + pursuit escalation (new)
- Pirates can now spawn in **packs** (leader + wingmen) instead of always solo.
- Police can spawn with **wingmen**, and their response scales with your **local bounty** and a soft pursuit "heat".
- Pirates may **break off and flee** when badly damaged (adds a bit of morale/variety to fights).
- Contacts UI now shows squad role tags: `[LEAD]` / `[WING]`, plus local `Heat`.

## Passenger missions + cabins (new)
- Added a new mission type: **Passenger** (party size stored in `Mission::units`).
- Added `SaveGame::passengerSeats` with Shipyard upgrades to increase capacity.
- Added seat-capacity validation on mission acceptance + completion support.
- SaveGame version bumped to **10** (backwards-compatible: old saves default to 2 seats).

## Mission system scaffolding (new)
- Added a reusable, **headless** mission module: `sim::MissionLogic`.
  - Deterministic mission-board generation persisted in `SaveGame`.
  - Shared accept logic (cargo checks, inventory/credit handling, station fees).
  - Shared completion/deadline helpers for Courier/Delivery/Multi-hop.
  - Game UI now delegates mission-offer generation + acceptance to the shared module.
- Added lightweight helpers:
  - `econ::cargoMassKg(...)` (shared cargo mass calc)
  - `sim::Reputation` helpers (rep + fee shaping)

## Tooling / CI
- Added `CMakePresets.json` for one-command headless builds.
- Updated GitHub Actions to build via CMake+Ninja and run `ctest`.

## Tests
- Added `test_missions` (mission board determinism + accept/complete path).

This patch adds a low-cost **ambient NPC trade traffic** layer that moves commodities between stations, nudging inventories and prices over time (even while you fly). It also persists these "traffic day stamps" in the save file so markets remain deterministic across reloads.

## Build fixes
- Added missing asteroid mining fields used by gameplay code:
  - `AsteroidNode::yield`
  - `AsteroidNode::chunkAccumulator`
- Added missing signal-site one-shot flag:
  - `SignalSource::fieldSpawned`
- Added missing NPC combat stat fields referenced by encounter logic:
  - `Contact::shieldMax`, `Contact::shieldRegenPerSec`, `Contact::hullMax`
- Fixed commodity loop using the correct enum/constant (`econ::kCommodityCount` instead of `CommodityId::COUNT`).
- Updated `spawnCargoPod` to accept a 5th parameter for scatter velocity (fixes "term does not evaluate to a function taking 5 arguments").
- Fixed string concatenation with `const char*` that caused pointer arithmetic ("Scooped" toast).
- Fixed invalid `.c_str()` calls on `CommodityDef::name` (it is `const char*`).

## Small gameplay/quality tweaks
- Ensured spawned pirates/traders/police initialize `shieldMax/hullMax` consistently.
- Added a simple NPC shield regen loop (uses `Contact::shieldRegenPerSec`).

## Economy + tooling upgrades
- Route planner now filters out *infeasible* trades (no source inventory or no destination capacity).
- Added `bestRoutesForCargo(...)` which considers cargo mass limits and market fees, and returns net trip profit.
- `stellar_sandbox` gained a headless `--trade` scan mode so you can quickly find profitable routes without running the full SDL/OpenGL prototype.

## Tests + CI
- Restored a procedural determinism test (`test_proc`) against the current streaming `Universe`.
- Added unit tests for:
  - route planner feasibility and cargo/fee math (`test_route_planner`)
  - SaveGame round-trip serialization (`test_savegame`)
- GitHub Actions now runs a headless `make distcheck` build (render/ImGui disabled) to keep CI fast and reliable.

## Ambient NPC trade traffic (new)
- Added `sim::simulateNpcTradeTraffic(...)` (core library) to model background station-to-station trade flows.
- Traffic simulation is:
  - deterministic per (universe seed, system id, day stamp)
  - bounded in CPU cost (backfills a limited number of days on large time jumps)
- Added save/load support for traffic stamps (`SaveGame` version bumped to **10**).
- Updated `test_savegame` to cover the new serialization.

## Streaming fix: SystemId encodes sector (new)
- Fixed `GalaxyGenerator::makeSystemId(...)` to **pack sector coords + local index** into the 64-bit id.
  - This restores the intended behavior where `Universe::getSystem(id)` can decode the sector and fetch the correct stub **without requiring a hint**.
- Added a regression check to `test_streaming` to ensure `getSystem(id)` returns the correct stub position/fields on a fresh Universe.
- Improved `Universe::getSystem(...)` fallback stub generation to create a **plausible galaxy-disc position** (instead of collapsing to {0,0,0}).

## UI: Window Registry + Window Manager (new)
- Added `stellar/ui/WindowRegistry` (headless): a small registry that binds window keys to runtime `bool` open flags.
- UI workspaces now capture/apply window visibility via the registry, reducing "forgot to wire the new window" bugs.
- Added a new in-game **Window Manager** panel (Windows → Window Manager...) with:
  - search + group filter
  - bulk open/close
  - reset-to-defaults
- Command Palette window toggles are now generated from the registry (includes shortcut hints where available).
- Added unit test `test_window_registry`.
## UI: In-game CPU Profiler (new)
- Added `stellar/core/Profiler` (headless): a lightweight per-frame CPU scope profiler with nested events.
- Added a new in-game **Profiler** window (Windows → Profiler) featuring:
  - rolling frame-time plot
  - flame graph timeline view
  - per-scope aggregate table with filtering + one-click copy-to-clipboard
- Instrumented the game loop with coarse scopes: `Frame`, `InputEvents`, `Sim`, `Present`.
- Added unit test `test_profiler`.

# Patch Notes (Jan 5, 2026)

## UI: Photo Mode + Screenshot Capture (new)
- Added a new in-game **Photo Mode** window (Windows → Photo Mode) for quick capture workflows.
- Supports 2 capture stages:
  - **World only** (captured before any ImGui UI is drawn)
  - **With UI** (captured after ImGui draws)
- Screenshots are written as **PNG** with:
  - configurable output directory (default: `screenshots/`)
  - timestamped, unique filenames (auto-suffix if a name already exists)
  - optional **copy path to clipboard** for fast sharing
- Added console command: `screenshot [ui|world] [basename] [dir]`.

## UI: In-game Log Viewer (new)
- Added `stellar/ui/LogBuffer`: a thread-safe ring buffer that can be fed via a `core::LogSink`.
- Added a new **Log** window (Windows → Log, or UI → Log) with:
  - per-level filters (Trace/Debug/Info/Warn/Error)
  - fuzzy text filter (optional relevance sorting)
  - copy selected line / copy filtered lines
  - export filtered lines to a text file (auto-creates parent folders)
- Added unit test `test_log_buffer`.

## Build/test fix
- Updated `test_atmosphere` to use the current `Star` field names (`massSol`, `radiusSol`).

# Patch Notes (Jan 5, 2026) — Round 6

## Traffic lanes + convoys: smoother arcs + schedule robustness
- **Lane arcs now ease-in/out**: the sideways lane offset uses `sin^2(pi*phase)` so convoys have ~zero lateral velocity at **depart/arrive** (they no longer “kick” sideways when spawning or docking).
- **Arc magnitude clamping tightened**: `arcMinKm` is treated as a preferred minimum *only when allowed* by the distance-based cap (prevents comically large arcs for short station-to-station hops).
- **Schedule generation made more robust**:
  - uses a stable min/max duration clamp
  - applies speed clamps safely
  - performs a single refinement step using destination position at arrival time
  - handles near-zero station separation without producing zero-duration trips

## Portability fixes
- Replaced the non‑standard `M_PI` usage with C++20 `std::numbers::pi_v<double>` in:
  - `sim/TrafficLanes` lane arc math
  - `sim/Interdiction` escape direction sampling

## Tests
- Added regression checks ensuring **endpoint lateral speed ~ 0** for:
  - `test_traffic_lanes`
  - `test_traffic_convoy_layer`

# Patch Notes (Jan 6, 2026) — Round 11

## Save/load: persist recent NPC trade shipments
- Added `SaveGame::trafficShipments` to persist recent NPC trade shipments in a save-friendly format.
- `stellar_game` now restores per-system `TrafficLedger` state from `SaveGame::trafficShipments` on load.
  - Fixes a mismatch where markets were advanced via `trafficStamps`, but convoy signals would fall back to
    deterministic lanes after a save/load (because no new day had advanced to re-record shipments).
- `simulateNpcTradeTraffic(...)` now prunes the optional `TrafficLedger` even when no new day is simulated,
  keeping memory bounded after loads / time tweaks.

## Tests
- Extended `test_savegame` to cover `trafficShipments` round-tripping and robustness against stale section counts.


# Patch Notes (Jan 6, 2026) — Round 12

## Supercruise: safer drop logic at the boundary
- Fixed a long-standing edge-case where a **manual drop exactly at the drop radius** could be flagged as an
  **emergency drop** due to a strict `<` comparison. Safe-drop readiness now treats `dist == dropRadius` as valid.
- `guideSupercruise(...)` no longer forces a drop in the **degenerate `dist≈0` case** when `interdicted=true`
  (keeps interdiction gating consistent with the normal drop path).

## Traffic convoys: cleaner inactive state
- When convoys are queried with `includeInactive=true`, convoys outside their active schedule window are now treated as
  **stationary at their endpoint** (zero velocity), avoiding confusing non-zero velocities in UI/nav.

## Game: robustness
- Added a defensive **de-duplication** pass for `TrafficConvoy` signals by id during refresh, ensuring older saves/patches
  that accidentally left duplicates cannot cause convoy lists to grow.

## Tests
- Extended `test_supercruise` with boundary + interdiction regressions.


# Patch Notes (Jan 6, 2026) — Round 13

## Save/load: traffic shipment robustness
- **De-dupe** `trafficShipments` by shipment id during save-file parsing (guards against corrupted/hand-edited files and older patches that may have accidentally duplicated entries).
- Traffic shipment parsing now treats schedule metadata (depart/arrive/dist/speed) as **optional** while still requiring the core identifiers to be present.

## Game: TrafficLedger restore correctness
- Restored shipments are ingested via `TrafficLedger::record(...)` (now supports move ingestion) so duplicates can't surface as duplicated convoy signals.

## Signals: deterministic distress TTL
- Distress call `expireDay` is now **anchored to the integer day stamp** (instead of `timeDays + ttl`), preventing repeated queries within a day from “extending” deterministic distress calls.

## Tests
- Extended `test_savegame` with a regression that ensures duplicated `trafficShipments` lines are de-duped and do not break parsing of later keys.
## 2026-01-06 (Patch) - Nav Assist + Traffic Convoy Target QoL

## Game: Nav Assist respects moving convoy signals
- When your **current target is a Traffic Convoy signal**, nav assist now uses the convoy's **current velocity** (instead of treating it as stationary).
  - This fixes "match velocity" feeling wrong for convoy targets.
  - It also improves approach guidance by allowing intercept-course logic to lead the target properly.

## UI: Status panel shows convoy details
- The **Status** window's target summary now displays Traffic Convoy **cargo, route, speed, and ETA** when a convoy signal is targeted.
