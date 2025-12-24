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
- **Player ship**: flight model + navigation modes (**normal flight / supercruise / docked**)
- **Targeting / nav**: select a **station or planet** as a navigation target (HUD marker)
- **Supercruise**: high-speed travel mode with simplified physics
- **Supercruise assist**: optional **"7-second" approach throttle hold** + **auto-drop** near target
- **Docking**: station approach corridor + speed limit + alignment requirement
- **Docking clearance**: **request docking** (not always granted; must be within range)
- **Approach autopilot**: station approach autopilot that lines you up with the corridor and respects limits
- **Time compression UX**: PageUp/PageDown time acceleration with safety clamps (CTRL to force)
- **Factions + markets**: deterministic faction generation + per-system market generation
- **Persistence**: lightweight save file (`savegame.txt`) for player/system/time/cargo

Gameplay / UX additions (Elite/Pioneer-style progression):

- **Targeting system**: select a **station or planet** as your nav target (HUD marker)
- **Supercruise**: a high-speed travel mode with simplified physics
  - Optional **Supercruise Assist**: **7-second-rule** style auto-throttle + auto-drop
- **Docking approach corridor**: visual corridor line + corridor constraints
- **Approach autopilot**: holds corridor + respects speed limit (when a station is targeted)
- **Request docking clearance**: docking is not always allowed; you must request permission
- **Time compression hotkeys** with safety clamps (hold CTRL to force)

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

### Controls (prototype)

UI / windows:

- **TAB**: toggle Galaxy window
- **F1**: toggle Flight window
- **F2**: toggle Economy window
- **F3**: toggle Nav window
- **F4**: toggle HUD overlay

Flight:

- **WASD + R/F**: translate (normal flight)
- **Arrow keys**: pitch/yaw
- **Q/E**: roll
- **LShift**: boost
- **X**: brake (normal) / throttle cut (supercruise)
- **Z/C**: dampers on/off
- **V**: camera toggle (chase/cockpit)

Navigation:

- **T**: cycle target (stations + planets)
- **Y**: clear target
- **P**: toggle **Approach Autopilot** (normal) / **Supercruise Assist** (supercruise)
- **J**: toggle supercruise

Docking:

- **L**: request docking clearance (must be within range)
- **G**: dock / undock

Time / sim:

- **PageUp / PageDown**: time compression
- **Home**: 1x time
- **End**: max time
- Hold **CTRL** to force time compression beyond safety clamps
- **Space**: pause
- **F5**: save
- **F9**: load

The game will auto-load `savegame.txt` if present; press **F5** to save.

## Next steps you can add easily

- Proper scene graph + render batching
- UI (ImGui) for market/cargo/ship stats
- Hyperjump / travel rules using galaxy positions
- More detailed economy: production chains, station inventories, price history
- Persistence for visited systems, player reputation, mission log
