# Patch Notes

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
