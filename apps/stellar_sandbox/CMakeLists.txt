add_executable(stellar_sandbox
  main.cpp
)

target_compile_features(stellar_sandbox PRIVATE cxx_std_20)
target_link_libraries(stellar_sandbox PRIVATE stellar::stellar)

if (MSVC)
  target_compile_options(stellar_sandbox PRIVATE /W4)
else()
  target_compile_options(stellar_sandbox PRIVATE -Wall -Wextra -Wpedantic)
endif()
