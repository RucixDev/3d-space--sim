cmake_minimum_required(VERSION 3.20)

project(StellarForge
  VERSION 0.1.0
  DESCRIPTION "Procedural space simulation starter framework"
  LANGUAGES CXX
)

option(STELLAR_BUILD_SANDBOX "Build the sandbox executable" ON)
option(STELLAR_BUILD_TESTS "Build tests" ON)
option(STELLAR_BUILD_GAME "Build the real-time SDL2/OpenGL prototype" ON)
option(STELLAR_FETCH_SDL2 "Fetch SDL2 via FetchContent if not found" ON)

# Library
add_library(stellar)
add_library(stellar::stellar ALIAS stellar)

target_compile_features(stellar PUBLIC cxx_std_20)

target_include_directories(stellar
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# Warnings / compile options
if (MSVC)
  target_compile_options(stellar PRIVATE /W4 /permissive-)
else()
  target_compile_options(stellar PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
endif()

# Sources
target_sources(stellar PRIVATE
  src/core/Hash.cpp
  src/core/Log.cpp

  src/math/Vec3.cpp

  src/proc/GalaxyGenerator.cpp
  src/proc/NameGenerator.cpp
  src/proc/Noise.cpp
  src/proc/SystemGenerator.cpp

  src/sim/Orbit.cpp
  src/sim/Faction.cpp
  src/sim/Market.cpp
  src/sim/SaveGame.cpp
  src/sim/Ship.cpp
  src/sim/Universe.cpp
)

# Nice folder structure for IDEs
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES
  src/core/Hash.cpp
  src/core/Log.cpp
  src/math/Vec3.cpp
  src/proc/GalaxyGenerator.cpp
  src/proc/NameGenerator.cpp
  src/proc/Noise.cpp
  src/proc/SystemGenerator.cpp
  src/sim/Orbit.cpp
  src/sim/Faction.cpp
  src/sim/Market.cpp
  src/sim/SaveGame.cpp
  src/sim/Ship.cpp
  src/sim/Universe.cpp
)

# Rendering helper library (no windowing; uses GL function pointers loaded by the app).
add_library(stellar_render)
add_library(stellar::render ALIAS stellar_render)

target_compile_features(stellar_render PUBLIC cxx_std_20)
target_include_directories(stellar_render
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

if (MSVC)
  target_compile_options(stellar_render PRIVATE /W4 /permissive-)
else()
  target_compile_options(stellar_render PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
endif()

target_sources(stellar_render PRIVATE
  src/render/Gl.cpp
  src/render/Shader.cpp
  src/render/PointRenderer.cpp
)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES
  src/render/Gl.cpp
  src/render/Shader.cpp
  src/render/PointRenderer.cpp
)

if (STELLAR_BUILD_SANDBOX)
  add_subdirectory(apps/stellar_sandbox)
endif()

if (STELLAR_BUILD_GAME)
  # SDL2: prefer a system / vcpkg package, else FetchContent if enabled.
  find_package(SDL2 CONFIG QUIET)
  if (NOT SDL2_FOUND)
    find_package(SDL2 QUIET)
  endif()

  if (NOT SDL2_FOUND AND STELLAR_FETCH_SDL2)
    include(FetchContent)

    # Build SDL without its test suite.
    set(SDL_TEST OFF CACHE BOOL "" FORCE)
    set(SDL_TESTS OFF CACHE BOOL "" FORCE)
    set(SDL_SHARED ON CACHE BOOL "" FORCE)
    set(SDL_STATIC OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
      SDL2
      GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
      GIT_TAG release-2.32.10
    )
    FetchContent_MakeAvailable(SDL2)
  endif()

  if (NOT SDL2_FOUND AND NOT TARGET SDL2::SDL2)
    message(FATAL_ERROR "SDL2 not found. Install SDL2 development libraries or configure with -DSTELLAR_FETCH_SDL2=ON.")
  endif()

  find_package(OpenGL REQUIRED)
  add_subdirectory(apps/stellar_game)
endif()

if (STELLAR_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
