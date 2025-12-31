#pragma once

#include <string_view>
#include <vector>

namespace stellar::ui {

// Lightweight ASCII-only fuzzy matcher for UI search/filter.
//
// Notes:
//  - "positions" are byte offsets into the input text (safe for ASCII UI strings).
//  - If query is empty/whitespace, score is 0 and positions is empty.
//  - If no match is found, score is -1.
struct FuzzyMatchResult {
  int score{0};
  std::vector<int> positions;
};

// Compute a fuzzy match score and the matched character positions.
FuzzyMatchResult fuzzyMatch(std::string_view query, std::string_view text);

// Convenience: score-only.
inline int fuzzyMatchScore(std::string_view query, std::string_view text) {
  return fuzzyMatch(query, text).score;
}

} // namespace stellar::ui
