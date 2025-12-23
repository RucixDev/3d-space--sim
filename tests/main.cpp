#include <iostream>

int test_orbit();
int test_streaming();
int test_economy();

int main() {
  int fails = 0;

  fails += test_orbit();
  fails += test_streaming();
  fails += test_economy();

  if (fails == 0) {
    std::cout << "[stellar_tests] ALL PASS\n";
    return 0;
  }

  std::cerr << "[stellar_tests] FAILS=" << fails << "\n";
  return 1;
}
