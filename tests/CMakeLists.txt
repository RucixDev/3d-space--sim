# Test runner
#
# This directory intentionally uses a tiny bespoke harness (see test_harness.h)
# instead of pulling a large third-party test framework.
#
# To reduce maintenance burden, the list of test translation units is globbed and
# a small registry translation unit is generated at CMake configure-time.

file(GLOB STELLAR_TEST_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp"
)

# Keep ordering stable across platforms/generators.
list(SORT STELLAR_TEST_SOURCES)

# Render-only tests.
# (These depend on the optional OpenGL rendering backend being compiled into the library.)
if (NOT STELLAR_ENABLE_RENDER)
  list(FILTER STELLAR_TEST_SOURCES EXCLUDE REGEX "test_frame_graph\\.cpp$")
  list(FILTER STELLAR_TEST_SOURCES EXCLUDE REGEX "test_surface_normalmap\\.cpp$")
endif()

# ---- Generate a registry TU that maps file base-names -> test function symbols ----
# Convention: each test file test_xxx.cpp defines: int test_xxx();
set(STELLAR_TEST_GEN_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(MAKE_DIRECTORY "${STELLAR_TEST_GEN_DIR}")

set(STELLAR_TEST_REGISTRY_CPP "${STELLAR_TEST_GEN_DIR}/generated_test_registry.cpp")

set(_registry "#include \"test_registry.h\"\n\n#include <cstddef>\n\nnamespace stellar::test {\n\n")

# Forward declarations.
foreach(src ${STELLAR_TEST_SOURCES})
  get_filename_component(name "${src}" NAME_WE)
  string(APPEND _registry "int ${name}();\n")
endforeach()

string(APPEND _registry "\nstatic const Case kCases[] = {\n")
foreach(src ${STELLAR_TEST_SOURCES})
  get_filename_component(name "${src}" NAME_WE)
  string(APPEND _registry "  {\"${name}\", &${name}},\n")
endforeach()
string(APPEND _registry "};\n\nconst Case* cases(std::size_t* outCount) {\n  if (outCount) *outCount = sizeof(kCases) / sizeof(kCases[0]);\n  return kCases;\n}\n\n} // namespace stellar::test\n")

file(WRITE "${STELLAR_TEST_REGISTRY_CPP}" "${_registry}")

add_executable(stellar_tests
  main.cpp
  ${STELLAR_TEST_SOURCES}
  ${STELLAR_TEST_REGISTRY_CPP}
)

target_include_directories(stellar_tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(stellar_tests PRIVATE cxx_std_20)

target_link_libraries(stellar_tests PRIVATE stellar::stellar)

add_test(NAME stellar_tests COMMAND stellar_tests)

set_target_properties(stellar_tests PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
)
