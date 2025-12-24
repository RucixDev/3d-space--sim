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
- **Player ship**: simple input + physics integration (AU/day units)
- **Factions + markets**: deterministic faction generation + per-system market generation
- **Persistence**: lightweight save file (`savegame.txt`) for player/system/time/cargo

- **In-system time acceleration**: PageUp/PageDown time compression with safety clamps (hold CTRL to override)
- **Docking corridor + docking state**: stations have physical orbits; dock/undock with **G** (trading requires docking)
- **Supercruise-like mode**: toggle with **J**, throttle with **W/S** (simplified high-speed travel)

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
cmake -S . -B build -DSTELLAR_BUILD_GAME=OFF
```

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

The game will auto-load `savegame.txt` if present; press **F5** to save.

## Next steps you can add easily

- Proper scene graph + render batching
- UI (ImGui) for market/cargo/ship stats
- Hyperjump / travel rules using galaxy positions
- More detailed economy: production chains, station inventories, price history
- Persistence for visited systems, player reputation, mission log
