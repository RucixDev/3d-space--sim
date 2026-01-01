#include <iostream>

int test_orbit();
int test_streaming();
int test_query_nearby();
int test_universe_cache();
int test_proc();
int test_economy();
int test_traffic();
int test_route_planner();
int test_manifest_planner();
int test_savegame();
int test_missions();
int test_law();
int test_police_scan();
int test_docking();
int test_docking_clearance();
int test_thermal();
int test_nav();
int test_args();
int test_signature();
int test_industry();
int test_industry_service();
int test_industry_scanner();
int test_warehouse();
int test_trade_scanner();
int test_ship();
int test_ship_loadout();
int test_power_distributor();
int test_flight_controller();
int test_intercept_course();
int test_combat();
int test_ballistics();
int test_supercruise();
int test_encounter_director();
int test_bookmarks();
int test_fuzzy_search();
int test_ui_settings();
int test_hud_settings();
int test_distress();
int test_resource_field();
int test_log_sinks();

int main() {
  int fails = 0;

  fails += test_orbit();
  fails += test_streaming();
  fails += test_query_nearby();
  fails += test_universe_cache();
  fails += test_proc();
  fails += test_economy();
  fails += test_traffic();
  fails += test_route_planner();
  fails += test_manifest_planner();
  fails += test_savegame();
  fails += test_missions();
  fails += test_law();
  fails += test_police_scan();
  fails += test_docking();
  fails += test_docking_clearance();
  fails += test_thermal();
  fails += test_nav();
  fails += test_args();
  fails += test_signature();
  fails += test_industry();
  fails += test_industry_service();
  fails += test_industry_scanner();
  fails += test_warehouse();
  fails += test_trade_scanner();
  fails += test_ship();
  fails += test_ship_loadout();
  fails += test_power_distributor();
  fails += test_flight_controller();
  fails += test_intercept_course();
  fails += test_combat();
  fails += test_ballistics();
  fails += test_supercruise();
  fails += test_encounter_director();
  fails += test_bookmarks();
  fails += test_fuzzy_search();
  fails += test_ui_settings();
  fails += test_hud_settings();
  fails += test_distress();
  fails += test_resource_field();
  fails += test_log_sinks();

  if (fails == 0) {
    std::cout << "[stellar_tests] ALL PASS\n";
    return 0;
  }

  std::cerr << "[stellar_tests] FAILS=" << fails << "\n";
  return 1;
}
