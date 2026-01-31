- **Round 2:** Fixed Windows build breaks in the new save/menu + UI loop.
- `CMakeLists.txt`: define `NOMINMAX` (and `WIN32_LEAN_AND_MEAN`) for all Windows targets to prevent `min/max` macro conflicts.
- `apps/stellar_game/main.cpp`: fix Time Trial ghost restore (qualify `game::FlightRecorderSample`) so the save snapshot loader compiles.
- `apps/stellar_game/CommsWindow.cpp`: updated to `SystemConditionsSnapshot.event` fields.
- `apps/stellar_game/LogbookWindow.cpp`, `apps/stellar_game/MarketDashboardWindow.cpp`: explicit `CStrView` → string conversions.

Verify:
- Build the `game-release` preset on Windows.
- In-game: open Comms during an active system event; save/load (slot or quicksave); confirm Time Trial ghost loads.