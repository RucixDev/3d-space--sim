# Patch Notes

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
