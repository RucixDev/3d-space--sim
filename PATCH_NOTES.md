## 2026-01-01 (Patch) - Build Fix: Missing ResourceField Sources + CMake Source Lists

This patch fixes the CMake configure-time errors:
- `Cannot find source file: src/sim/ResourceField.cpp`
- `Cannot find source file: test_resource_field.cpp`
- `No SOURCES given to target ...`

Changes:
- **Added missing files** for the new mining/resource-field scaffolding:
  - `include/stellar/sim/WorldIds.h`
  - `include/stellar/sim/ResourceField.h`
  - `src/sim/ResourceField.cpp`
  - `tests/test_resource_field.cpp`
- **CMake**: restored/updated `CMakeLists.txt` and `tests/CMakeLists.txt` so the `stellar` library and `stellar_tests` targets always have valid source lists.
- **Core**: ensured `src/sim/Distress.cpp` is compiled into the `stellar` library (required by the game build).
- **Tests**: wired in `test_resource_field` and `test_distress` to the test runner.

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