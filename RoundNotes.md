- **Round 2:** Windows build fixes (UI loop/save/menu). `NOMINMAX` + `WIN32_LEAN_AND_MEAN`, Time Trial ghost restore compile fix, and UI updates for `SystemConditionsSnapshot` + `CStrView` conversions.
- **Round 3:** Orbit Analyzer maneuver planner now wires cleanly to the in-game RTN maneuver node (dv along/normal/radial + ref-body choice) and compiles (bindings + function signature aligned).
- `apps/stellar_game/OrbitAnalyzerWindow.cpp`: fixed gravity body kind enum usage and corrected impact/closest-approach readouts (impact altitude computed from distance-to-body).
- `apps/stellar_game/main.cpp`: fixed Orbit Analyzer call wiring, GalNet digest -> Integration Hub event push, and removed undefined supercruise/hyperspace guards.
- `apps/stellar_game/GalNetInboxWindow.cpp`: updated to current `CommsLog` API (items/markRead/markPinned) and resolves system names via `Universe`.

Verify:
- Build the `game-release` preset.
- In-game: open Orbit Analyzer, click a planner action (e.g. Circularize) and confirm maneuver node + trajectory preview update; check Forecast impact/approach tables.
- In-game: open GalNet Inbox, select bulletins, mark read/unread/pin, and request an update; enable watchlist digest and confirm an event appears in Integration Hub.
