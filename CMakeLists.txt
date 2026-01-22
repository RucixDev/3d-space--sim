cmake_minimum_required(VERSION 3.20)

project(stellar_forge LANGUAGES CXX)

include(FetchContent)
include(CTest)

# Legacy option preserved for older build scripts / README snippets.
# Prefer STELLAR_ENABLE_RENDER + STELLAR_ENABLE_IMGUI going forward.
option(STELLAR_BUILD_GAME "(DEPRECATED) Build the SDL/OpenGL game prototype" ON)

option(STELLAR_ENABLE_RENDER "Build OpenGL/SDL rendering module and the game app" ON)
option(STELLAR_ENABLE_IMGUI  "Enable Dear ImGui UI (requires STELLAR_ENABLE_RENDER)" ON)

if (NOT STELLAR_BUILD_GAME)
  set(STELLAR_ENABLE_RENDER OFF CACHE BOOL "" FORCE)
  set(STELLAR_ENABLE_IMGUI OFF CACHE BOOL "" FORCE)
endif()

# ImGui cannot be enabled without the render path.
if (NOT STELLAR_ENABLE_RENDER)
  set(STELLAR_ENABLE_IMGUI OFF CACHE BOOL "" FORCE)
endif()


option(STELLAR_FETCH_SDL2 "Fetch SDL2 if not found via find_package (render builds only)" ON)
option(STELLAR_FETCH_IMGUI "Fetch Dear ImGui (render+imgui builds only)" ON)

# Tooling / correctness options
option(STELLAR_ENABLE_ASAN "Enable AddressSanitizer (Clang/GCC only)" OFF)
option(STELLAR_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer (Clang/GCC only)" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- Dependencies (conditional) ----
set(STELLAR_SDL_TARGET "")
set(IMGUI_TARGET "")

if (STELLAR_ENABLE_RENDER)
  find_package(OpenGL QUIET)
  if (NOT OpenGL_FOUND)
    message(WARNING "OpenGL not found. Disabling STELLAR_ENABLE_RENDER/STELLAR_ENABLE_IMGUI. Configure with -DSTELLAR_ENABLE_RENDER=OFF for headless builds, or install OpenGL development packages.")
    set(STELLAR_ENABLE_RENDER OFF CACHE BOOL "" FORCE)
    set(STELLAR_ENABLE_IMGUI OFF CACHE BOOL "" FORCE)
  else()

  # SDL2
  find_package(SDL2 CONFIG QUIET)
  if (NOT SDL2_FOUND AND STELLAR_FETCH_SDL2)
    message(STATUS "SDL2 not found via find_package. Fetching SDL2...")
    set(SDL_TEST OFF CACHE BOOL "" FORCE)
    set(SDL_TESTS OFF CACHE BOOL "" FORCE)
    set(SDL_SHARED OFF CACHE BOOL "" FORCE)
    set(SDL_STATIC ON CACHE BOOL "" FORCE)

    FetchContent_Declare(
      sdl2
      GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
      GIT_TAG release-2.30.9
    )
    FetchContent_MakeAvailable(sdl2)
  endif()

  # Determine SDL target name
  if (TARGET SDL2::SDL2)
    set(STELLAR_SDL_TARGET SDL2::SDL2)
  elseif (TARGET SDL2-static)
    set(STELLAR_SDL_TARGET SDL2-static)
  elseif (TARGET SDL2)
    set(STELLAR_SDL_TARGET SDL2)
  endif()

  if (STELLAR_SDL_TARGET STREQUAL "")
    message(FATAL_ERROR "SDL2 target not found. Install SDL2 or enable STELLAR_FETCH_SDL2.")
  endif()

  # Dear ImGui (optional)
  if (STELLAR_ENABLE_IMGUI AND STELLAR_FETCH_IMGUI)
    FetchContent_Declare(
      imgui
      GIT_REPOSITORY https://github.com/ocornut/imgui.git
      GIT_TAG v1.91.2-docking
    )
    FetchContent_MakeAvailable(imgui)

    add_library(imgui_lib STATIC
      ${imgui_SOURCE_DIR}/imgui.cpp
      ${imgui_SOURCE_DIR}/imgui_demo.cpp
      ${imgui_SOURCE_DIR}/imgui_draw.cpp
      ${imgui_SOURCE_DIR}/imgui_tables.cpp
      ${imgui_SOURCE_DIR}/imgui_widgets.cpp
      ${imgui_SOURCE_DIR}/backends/imgui_impl_sdl2.cpp
      ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
    )
    target_include_directories(imgui_lib PUBLIC
      ${imgui_SOURCE_DIR}
      ${imgui_SOURCE_DIR}/backends
    )
    target_link_libraries(imgui_lib PUBLIC ${STELLAR_SDL_TARGET} OpenGL::GL)
    target_compile_definitions(imgui_lib PUBLIC IMGUI_DISABLE_OBSOLETE_FUNCTIONS)
    set(IMGUI_TARGET imgui_lib)
  endif()

  # Sanity: if ImGui is enabled but no target was created, try to bind to a user-provided target
  # (common names) or disable ImGui to avoid confusing missing-header build errors.
  if (STELLAR_ENABLE_IMGUI AND IMGUI_TARGET STREQUAL "")
    if (TARGET imgui::imgui)
      set(IMGUI_TARGET imgui::imgui)
    elseif (TARGET imgui)
      set(IMGUI_TARGET imgui)
    endif()
  endif()

  if (STELLAR_ENABLE_IMGUI AND IMGUI_TARGET STREQUAL "")
    message(WARNING "STELLAR_ENABLE_IMGUI=ON but Dear ImGui was not provided. Disabling STELLAR_ENABLE_IMGUI. Set STELLAR_FETCH_IMGUI=ON (recommended) or provide an ImGui target before configuring.")
    set(STELLAR_ENABLE_IMGUI OFF CACHE BOOL "" FORCE)
  endif()
  endif() # OpenGL_FOUND
endif()

# ---- Library ----
set(STELLAR_CORE_SOURCES
  src/core/Assert.cpp
  src/core/Hash.cpp
  src/core/Log.cpp
  src/core/CVar.cpp
  src/core/Profiler.cpp
  src/core/ChromeTrace.cpp
  src/core/Base64.cpp

  src/ui/Bookmarks.cpp
  src/ui/FuzzySearch.cpp
  src/ui/LogBuffer.cpp
  src/ui/HudLayout.cpp
  src/ui/HudSettings.cpp
  src/ui/VfxSettings.cpp
  src/ui/ProceduralUi.cpp
  src/ui/TextFx.cpp
  src/ui/ProceduralShipHud.cpp
  src/ui/AudioSettings.cpp
  src/ui/Livery.cpp
  src/ui/UiLayoutProfiles.cpp
  src/ui/UiSettings.cpp
  src/ui/UiWorkspaces.cpp
  src/ui/WindowRegistry.cpp

  src/proc/Noise.cpp
  src/proc/FluidSim2D.cpp
  src/dsp/AudioAnalyzer.cpp
  src/dsp/BeatTracker.cpp
  src/proc/GalaxyHazards.cpp
  src/proc/GalaxyRegions.cpp
  src/proc/GalaxyClusters.cpp
  src/proc/GalaxyVoids.cpp
  src/proc/Hyperlanes.cpp
  src/proc/HyperlaneRouter.cpp
  src/proc/HyperlaneKRoutes.cpp
  src/proc/HyperlaneCentrality.cpp
  src/proc/HyperlaneVulnerability.cpp
  # Procedural geometry generators (CPU-side; safe for headless builds).
  src/render/ProceduralAsteroid.cpp
  src/render/SdfMesher.cpp
  src/render/SdfDualContouring.cpp
  src/render/MeshSimplify.cpp
  src/render/GltfExport.cpp
  src/render/SdfGraph.cpp
  src/render/ProceduralArtifact.cpp
  # Procedural render presets (CPU-side; safe for headless builds).
  src/render/ProceduralPostFX.cpp
  src/render/ProceduralCompositor.cpp
  src/render/AutoExposure.cpp
  src/render/ShaderToyParams.cpp
  # Procedural sky layers (CPU-side; safe for headless builds).
  src/render/Starfield.cpp
  src/render/Nebula.cpp
  src/proc/NameGenerator.cpp
  src/proc/GalaxyGenerator.cpp
  src/proc/GalaxyMorphology.cpp
  src/proc/SystemGenerator.cpp
  src/proc/AsteroidBeltGenerator.cpp
  src/proc/TradeProfile.cpp
  src/proc/TradeEconomy.cpp
  src/proc/TradeRoutes.cpp
  src/proc/TradePrices.cpp
  src/proc/TradeIntel.cpp
  src/proc/TradeFlow.cpp
  src/proc/MiePhase.cpp
  src/proc/PhaseMultipleScattering.cpp

  src/sim/Faction.cpp
  src/sim/FactionProfile.cpp
  src/sim/Orbit.cpp
  src/sim/Gravity.cpp
  src/sim/Atmosphere.cpp
  src/sim/Aerodynamics.cpp
  src/sim/TimeTrial.cpp
  src/sim/OrbitalMechanics.cpp
  src/sim/TrajectoryPredictor.cpp
  src/sim/TrajectoryAnalysis.cpp
  src/sim/ManeuverComputer.cpp
  src/sim/ManeuverProgramComputer.cpp
  src/sim/LambertSolver.cpp
  src/sim/LambertPlanner.cpp
  src/sim/Ship.cpp
  src/sim/Ballistics.cpp
  src/sim/Combat.cpp
  src/sim/Countermeasures.cpp
  src/sim/MissileDefense.cpp
  src/sim/FlightController.cpp
  src/sim/NavAssistComputer.cpp
  src/sim/SupercruiseComputer.cpp
  src/sim/Interdiction.cpp
  src/sim/EncounterDirector.cpp
  src/sim/FormationPlanner.cpp
  src/sim/SquadTactics.cpp
  src/sim/SecurityModel.cpp
  src/sim/SystemSecurityDynamics.cpp
  src/sim/SystemEvents.cpp
  src/sim/SystemEventEconomy.cpp
  src/sim/SystemConditions.cpp
  src/sim/DockingPads.cpp
  src/sim/Docking.cpp
  src/sim/DockingComputer.cpp
  src/sim/DockingClearanceService.cpp
  src/sim/ThermalSystem.cpp
  src/sim/FuelScoopSystem.cpp
  src/sim/PowerDistributor.cpp
  src/sim/Law.cpp
  src/sim/PoliceScan.cpp
  src/sim/LawLedger.cpp
  src/sim/BlackMarket.cpp
  src/sim/Extortion.cpp
  src/sim/NpcCombatAI.cpp
  src/sim/CargoJettisonPlanner.cpp
  src/sim/Distress.cpp
  src/sim/Derelict.cpp
  src/sim/ResourceField.cpp
  src/sim/Mining.cpp
  src/sim/TrafficEscort.cpp
  src/sim/Signals.cpp
  src/sim/Industry.cpp
  src/sim/IndustryService.cpp
  src/sim/IndustryScanner.cpp
  src/sim/Warehouse.cpp
  src/sim/ShipyardService.cpp
  src/sim/Insurance.cpp
  src/sim/Logbook.cpp
  src/sim/SaveGame.cpp
  src/sim/Traffic.cpp
  src/sim/TrafficTransit.cpp
  src/sim/TrafficLedger.cpp
  src/sim/TrafficLanes.cpp
  src/sim/TrafficConvoyLayer.cpp
  src/sim/TrafficConvoyEncounter.cpp
  src/sim/SmugglingScanner.cpp
  src/sim/TradeScanner.cpp
  src/sim/TradeLoopScanner.cpp
  src/sim/SalvageSweep.cpp
  src/sim/TradeRunPlanner.cpp
  src/sim/MarketAnalysis.cpp
  src/sim/Universe.cpp
  src/sim/NavRoute.cpp
  src/sim/NavRouteBatch.cpp
  src/sim/MissionAssist.cpp
  src/sim/MissionBriefing.cpp
  src/sim/Comms.cpp
  src/sim/MissionLogic.cpp

  src/econ/Commodity.cpp
  src/econ/Economy.cpp
  src/econ/Market.cpp
  src/econ/RoutePlanner.cpp
)

set(STELLAR_RENDER_SOURCES
  src/render/Gl.cpp
  src/render/Shader.cpp
  src/render/ShaderToy.cpp
  src/render/ShaderToyGraph.cpp
  src/render/GaussianSurfels.cpp
  src/render/Texture.cpp
  src/render/FrameGraph.cpp
  src/render/TextureReadback.cpp
  src/render/GpuSurfaceCache.cpp
  src/render/GpuParticleField.cpp
  src/render/ProceduralSprite.cpp
  src/render/ProceduralGraph.cpp
  src/render/SdfRaymarcher.cpp
  src/render/ProceduralPlanet.cpp
  src/render/ProceduralRings.cpp
  src/render/ProceduralLivery.cpp
  src/render/AtmosphereRenderer.cpp
  src/render/VolumetricAtmosphereRenderer.cpp
  src/render/ProceduralSky.cpp
  src/render/StarCoronaRenderer.cpp
  src/render/StarfieldGpuRenderer.cpp
  src/render/ParticleSystem.cpp
  src/render/PostFX.cpp
  src/render/RenderTarget.cpp
  src/render/Mesh.cpp
  src/render/Camera.cpp
  src/render/MeshRenderer.cpp
  src/render/PointRenderer.cpp
  src/render/LineRenderer.cpp
)

add_library(stellar STATIC
  ${STELLAR_CORE_SOURCES}
  $<$<BOOL:${STELLAR_ENABLE_RENDER}>:${STELLAR_RENDER_SOURCES}>
)

add_library(stellar::stellar ALIAS stellar)

target_include_directories(stellar PUBLIC
  ${PROJECT_SOURCE_DIR}/include
)

target_compile_features(stellar PUBLIC cxx_std_20)

# Expose build options to C/C++ for feature-gated code paths.
target_compile_definitions(stellar PUBLIC
  STELLAR_ENABLE_RENDER=$<BOOL:${STELLAR_ENABLE_RENDER}>
  STELLAR_ENABLE_IMGUI=$<BOOL:${STELLAR_ENABLE_IMGUI}>
)

# ---- Tooling: sanitizers (optional) ----
# These flags help catch memory errors and undefined behavior early when running tests.
add_library(stellar_sanitizers INTERFACE)

if ((STELLAR_ENABLE_ASAN OR STELLAR_ENABLE_UBSAN) AND (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU"))
  set(_san_list "")
  if (STELLAR_ENABLE_ASAN)
    list(APPEND _san_list address)
  endif()
  if (STELLAR_ENABLE_UBSAN)
    list(APPEND _san_list undefined)
  endif()
  list(JOIN _san_list "," _san_csv)

  target_compile_options(stellar_sanitizers INTERFACE
    -fsanitize=${_san_csv}
    -fno-omit-frame-pointer
    -fno-sanitize-recover=all
  )
  target_link_options(stellar_sanitizers INTERFACE
    -fsanitize=${_san_csv}
  )
elseif (STELLAR_ENABLE_ASAN OR STELLAR_ENABLE_UBSAN)
  message(WARNING "Sanitizers requested but compiler '${CMAKE_CXX_COMPILER_ID}' does not support -fsanitize flags here. Ignoring STELLAR_ENABLE_ASAN/UBSAN.")
endif()

# Propagate sanitizer flags from the core library to dependents (tests/apps).
target_link_libraries(stellar PUBLIC stellar_sanitizers)


if (STELLAR_ENABLE_RENDER)
  target_link_libraries(stellar PUBLIC
    OpenGL::GL
    ${STELLAR_SDL_TARGET}
  )

  if (STELLAR_ENABLE_IMGUI AND NOT IMGUI_TARGET STREQUAL "")
    # The core library contains small ImGui-backed utilities (e.g. TextFx rendering), so
    # it must see ImGui include paths when STELLAR_ENABLE_IMGUI is on.
    target_link_libraries(stellar PUBLIC ${IMGUI_TARGET})
  endif()
  if (WIN32)
    target_link_libraries(stellar PUBLIC gdi32)
  endif()
endif()

# ---- Apps ----
add_executable(stellar_sandbox apps/stellar_sandbox/main.cpp)
target_link_libraries(stellar_sandbox PRIVATE stellar::stellar)
set_target_properties(stellar_sandbox PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/apps
)

if (STELLAR_ENABLE_RENDER AND STELLAR_ENABLE_IMGUI)
  add_executable(stellar_game
    apps/stellar_game/main.cpp
    apps/stellar_game/ShipHudOverlay.cpp
    apps/stellar_game/ActionWheel.cpp
    apps/stellar_game/AudioEngine.cpp
    apps/stellar_game/AudioAnalyzerWindow.cpp
    apps/stellar_game/CommandPalette.cpp
    apps/stellar_game/ControlsWindow.cpp
    apps/stellar_game/CommsWindow.cpp
    apps/stellar_game/SensorContacts.cpp
    apps/stellar_game/MarketDashboardWindow.cpp
  apps/stellar_game/LogbookWindow.cpp
    apps/stellar_game/TradePlannerWindow.cpp
    apps/stellar_game/ConsoleWindow.cpp
    apps/stellar_game/LogWindow.cpp
    apps/stellar_game/CVarWindow.cpp
    apps/stellar_game/ProfilerWindow.cpp
    apps/stellar_game/BuildInfoWindow.cpp
    apps/stellar_game/RuntimeValidationWindow.cpp
  apps/stellar_game/IntegrationHubWindow.cpp
    apps/stellar_game/FlightRecorderWindow.cpp
    apps/stellar_game/CinematicCameraWindow.cpp
    apps/stellar_game/CameraRigWindow.cpp
    apps/stellar_game/TimeTrialWindow.cpp
    apps/stellar_game/OrbitAnalyzerWindow.cpp
    apps/stellar_game/TrafficLanesWindow.cpp
    apps/stellar_game/SystemConditionsWindow.cpp
    apps/stellar_game/PorkchopPlot.cpp
    apps/stellar_game/ProceduralLabWindow.cpp
    apps/stellar_game/ProceduralFluidLabWindow.cpp
    apps/stellar_game/TextAnimationLabWindow.cpp
    apps/stellar_game/ProceduralMeshLabWindow.cpp
    apps/stellar_game/ProceduralShaderLabWindow.cpp
    apps/stellar_game/ProceduralGalaxyLabWindow.cpp
    apps/stellar_game/ProceduralSystemLabWindow.cpp
    apps/stellar_game/ProceduralTradeSystemsLabWindow.cpp
    apps/stellar_game/SpectralMieLabWindow.cpp
    apps/stellar_game/GaussianSurfelReconstructionLabWindow.cpp
    apps/stellar_game/Screenshot.cpp
    apps/stellar_game/PhotoModeWindow.cpp
  )
  target_link_libraries(stellar_game PRIVATE stellar::stellar ${IMGUI_TARGET})
  set_target_properties(stellar_game PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/apps
  )
endif()

# ---- Tests ----
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()
