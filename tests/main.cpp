#include <iostream>

int test_orbit();
int test_streaming();
int test_query_nearby();
int test_universe_cache();
int test_proc();
int test_economy();
int test_traffic();
int test_route_planner();
int test_savegame();
int test_missions();
int test_nav();
int test_args();
int test_signature();

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
  fails += test_savegame();
  fails += test_missions();
  fails += test_nav();
  fails += test_args();
  fails += test_signature();

  if (fails == 0) {
    std::cout << "[stellar_tests] ALL PASS\n";
    return 0;
  }

  std::cerr << "[stellar_tests] FAILS=" << fails << "\n";
  return 1;
}
