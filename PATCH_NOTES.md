# Patch Notes (Dec 26, 2025)

This patch focuses on fixing the current MSVC compile errors in `apps/stellar_game/main.cpp` and smoothing out a few early gameplay systems that were already being referenced by the code.

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

