# StellarForge — C++ Procedural Space Simulation Framework (Starter)

This is a **starter C++20 framework** aimed at building a space simulation game with **deterministic procedural generation** (galaxy → star systems → planets).  
The core simulation/proc library is intentionally lightweight and self-contained.

This repo now also includes a **real-time prototype app** using **SDL2 + OpenGL** (optional; controlled via CMake options).

It’s inspired by the *kind* of open-ended space adventure pioneered by games like **Pioneer** (exploration/trading across huge numbers of systems), but **all code here is original** and designed as a clean foundation you can extend.

## What you get in this template

- **CMake** project layout (library + sandbox app)
- Deterministic PRNG (**SplitMix64**) with stable seeding for repeatable worlds
- **Procedural generation** modules:
  - Galaxy generator (disc-like distribution)
  - Star + planet generator (simple astrophysics-inspired heuristics)
  - Name generator (syllable-based)
  - Value noise + fBm utilities for future terrain/fields
- A simple **Keplerian orbit** solver (elliptical orbit position in AU)

New in this build:

- **Streaming universe**: systems are generated **on-demand** (infinite indices) with an **LRU cache**
- **Renderer layer**: a minimal OpenGL point renderer (GL function loading via SDL's `SDL_GL_GetProcAddress`)
- **Player ship**: simple input + physics integration
- **Factions + markets**: deterministic faction generation + per-system market generation
- **Persistence**: lightweight save file (`savegame.txt`) for player/system/time/cargo

Early gameplay loops added in this pass:

- **Supercruise** (high-speed in-system travel) + **7-second rule** approach assist
- **Docking corridors** + **request docking** clearance (traffic can deny)
- **Station services gating** (market/refuel/repairs/missions require docking)
- **FSD hyperspace jump** (fuel + cooldown + arrival near destination station)
- **NPC traffic** (traders affect station inventories) + **simple interdictions**
- **Fuel scoop + heat** loop (scoop near star, too close overheats)
- **Missions**: courier, delivery cargo, bounty scan

## Build

By default the project builds:

- `stellar` (core simulation/proc library)
- `stellar_render` (OpenGL helper library)
- `stellar_sandbox` (CLI)
- `stellar_game` (SDL2/OpenGL prototype)
- tests

### Linux / macOS
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/apps/stellar_sandbox/stellar_sandbox --help
```

### Windows (Visual Studio / Developer PowerShell)
```powershell
cmake -S . -B build
cmake --build build --config Release
.\build\apps\stellar_sandbox\Release\stellar_sandbox.exe --help
```

If you only want the headless library + sandbox/tests (no SDL/OpenGL), configure with:

```bash
cmake -S . -B build -DSTELLAR_ENABLE_RENDER=OFF -DSTELLAR_ENABLE_IMGUI=OFF
```

(Legacy alias still supported: `-DSTELLAR_BUILD_GAME=OFF`)

## Run the sandbox

Sample N systems and print one system:
```bash
./stellar_sandbox --seed 12345 --systems 50 --pick 7
```

Print multiple systems:
```bash
./stellar_sandbox --seed 42 --systems 10 --list
```

Show orbital positions at a given time:
```bash
./stellar_sandbox --seed 42 --systems 5 --pick 0 --time 120
```

## Run the real-time prototype

The `stellar_game` executable is built by default.

### Dependencies

- **SDL2** (development headers)
- **OpenGL** (system OpenGL)

If SDL2 isn't installed, CMake can fetch it automatically:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSTELLAR_FETCH_SDL2=ON
cmake --build build -j
```

If you prefer system packages:

- Debian/Ubuntu: `sudo apt install libsdl2-dev`
- macOS: `brew install sdl2`

Then:

```bash
./build/apps/stellar_game/stellar_game
```

The game will auto-load `savegame.txt` if present.

This prototype now auto-saves periodically. Use `rm savegame.txt` to reset.

### Controls (prototype)

- Flight:
  - W/S: forward/back thrust
  - A/D: strafe left/right
  - Space/Ctrl: up/down
  - Arrow keys: pitch/yaw
  - Q/E: roll
  - Shift: boost, X: brake, F: toggle dampers
- Targeting / travel:
  - T: cycle stations, Y: cycle planets, G: clear target
  - B: toggle supercruise, PageUp/PageDown: throttle, Z: toggle supercruise assist
  - 1..5: time acceleration (Pioneer-style; only when safe)
- Docking:
  - L: request docking
  - P: toggle approach autopilot
  - U: undock (when docked)
- Ship systems:
  - O: toggle fuel scoop
- UI:
  - F1 help, F2 nav, F3 station/services, F4 galaxy map, F5 missions

## Next steps you can add easily

- Proper scene graph + render batching
- UI (ImGui) for market/cargo/ship stats
- Hyperjump / travel rules using galaxy positions
- More detailed economy: production chains, station inventories, price history
- Persistence for visited systems, player reputation, mission log
