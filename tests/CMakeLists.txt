add_executable(stellar_tests
  main.cpp
  test_orbit.cpp
  test_streaming.cpp
  test_query_nearby.cpp
  test_universe_cache.cpp
  test_proc.cpp
  test_economy.cpp
  test_traffic.cpp
  test_route_planner.cpp
  test_manifest_planner.cpp
  test_savegame.cpp
  test_missions.cpp
  test_law.cpp
  test_police_scan.cpp
  test_faction_profile.cpp
  test_security_model.cpp
  test_distress.cpp
  test_resource_field.cpp
  test_signals.cpp
  test_docking.cpp
  test_docking_clearance.cpp
  test_thermal.cpp
  test_nav.cpp
  test_k_routes.cpp
  test_args.cpp
  test_signature.cpp
  test_industry.cpp
  test_industry_service.cpp
  test_industry_scanner.cpp
  test_warehouse.cpp
  test_trade_scanner.cpp
  test_trade_loops.cpp
  test_trade_loops_parallel.cpp
  test_trade_runs.cpp
  test_trade_runs_parallel.cpp
  test_ship.cpp
  test_ship_loadout.cpp
  test_shipyard_service.cpp
  test_insurance.cpp
  test_power_distributor.cpp
  test_flight_controller.cpp
  test_nav_assist.cpp
  test_intercept_course.cpp
  test_interdiction.cpp
  test_combat.cpp
  test_ballistics.cpp
  test_supercruise.cpp
  test_encounter_director.cpp
  test_bookmarks.cpp
  test_fuzzy_search.cpp
  test_log_sinks.cpp
  test_cvars.cpp
  test_jobs.cpp
  test_ui_settings.cpp
  test_hud_settings.cpp
  test_ui_workspaces.cpp
)

target_compile_features(stellar_tests PRIVATE cxx_std_20)
target_link_libraries(stellar_tests PRIVATE stellar::stellar)

add_test(NAME stellar_tests COMMAND stellar_tests)

set_target_properties(stellar_tests PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
)
