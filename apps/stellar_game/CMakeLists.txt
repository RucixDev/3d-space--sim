add_executable(stellar_game
  main.cpp
  ActionWheel.cpp
  AudioEngine.cpp
  CinematicCameraWindow.cpp
)

target_compile_features(stellar_game PRIVATE cxx_std_20)

# Avoid SDL's main() macro indirection. This keeps the entry point predictable across platforms.
target_compile_definitions(stellar_game PRIVATE SDL_MAIN_HANDLED)

target_link_libraries(stellar_game
  PRIVATE
    stellar::stellar
    stellar::render
)

# SDL2 target name differs depending on how it's provided (system package vs FetchContent).
if (TARGET SDL2::SDL2)
  target_link_libraries(stellar_game PRIVATE SDL2::SDL2)
elseif (TARGET SDL2)
  target_link_libraries(stellar_game PRIVATE SDL2)
endif()

target_link_libraries(stellar_game PRIVATE OpenGL::GL)

if (MSVC)
  target_compile_options(stellar_game PRIVATE /W4 /permissive-)
else()
  target_compile_options(stellar_game PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
endif()
