#include "stellar/ui/FuzzySearch.h"

#include "test_harness.h"

#include <string>

int test_fuzzy_search() {
  int failures = 0;

  {
    const auto r = stellar::ui::fuzzyMatch("abc", "abc");
    CHECK(r.score > 0);
    CHECK(r.positions.size() == 3);
    CHECK(r.positions[0] == 0);
    CHECK(r.positions[1] == 1);
    CHECK(r.positions[2] == 2);
  }

  {
    // Subsequence across separators.
    const auto r = stellar::ui::fuzzyMatch("abc", "a_b_c");
    CHECK(r.score > 0);
    CHECK(r.positions.size() == 3);
    CHECK(r.positions[0] == 0);
    CHECK(r.positions[1] == 2);
    CHECK(r.positions[2] == 4);
  }

  {
    // Multi-token match should enforce token order.
    const auto r = stellar::ui::fuzzyMatch("dock clear", "Docking Clearance Request");
    CHECK(r.score > 0);
    CHECK(!r.positions.empty());
    // Should hit a word boundary at least once.
    bool hasBoundaryHit = false;
    for (int p : r.positions) {
      if (p == 0) {
        hasBoundaryHit = true;
        break;
      }
    }
    CHECK(hasBoundaryHit);
  }

  {
    const auto r = stellar::ui::fuzzyMatch("zzz", "Docking Clearance");
    CHECK(r.score < 0);
    CHECK(r.positions.empty());
  }

  {
    // Whitespace-only query is treated as empty.
    const auto r = stellar::ui::fuzzyMatch("   \t\n", "Anything");
    CHECK(r.score == 0);
    CHECK(r.positions.empty());
  }

  return failures;
}
