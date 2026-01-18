#include "stellar/ui/FuzzySearch.h"

#include <algorithm>
#include <cctype>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace stellar::ui {
namespace {

static bool isSpace(unsigned char c) { return std::isspace(c) != 0; }

static bool isWordBoundary(char c) {
  // Common separators in command palettes, file paths, UI labels, etc.
  return c == ' ' || c == '_' || c == '-' || c == '/' || c == '\\' || c == '.' || c == ':' || c == '(' || c == ')'
         || c == '[' || c == ']' || c == '{' || c == '}' || c == '#' || c == '|' || c == ',' || c == ';' || c == '=';
}

static bool isAsciiLower(char c) { return c >= 'a' && c <= 'z'; }
static bool isAsciiUpper(char c) { return c >= 'A' && c <= 'Z'; }
static bool isAsciiDigit(char c) { return c >= '0' && c <= '9'; }

static bool isCamelBoundary(char prev, char cur) {
  // lower->Upper is a natural boundary in CamelCase.
  if (isAsciiLower(prev) && isAsciiUpper(cur)) return true;
  // Digit->Alpha is also a useful boundary in many identifiers.
  if (isAsciiDigit(prev) && (isAsciiLower(cur) || isAsciiUpper(cur))) return true;
  return false;
}

static char lowerAscii(char c) { return (char)std::tolower((unsigned char)c); }

static std::string_view trimAscii(std::string_view s) {
  std::size_t b = 0;
  while (b < s.size() && isSpace((unsigned char)s[b])) ++b;
  std::size_t e = s.size();
  while (e > b && isSpace((unsigned char)s[e - 1])) --e;
  return s.substr(b, e - b);
}

static void pushToken(std::vector<std::string_view>& out, std::string_view s, std::size_t b, std::size_t e) {
  if (e > b) out.push_back(s.substr(b, e - b));
}

static constexpr int kNegInf = std::numeric_limits<int>::min() / 4;

static inline int addScore(int a, int b) {
  if (a <= kNegInf / 2) return kNegInf;
  return a + b;
}

struct TokenMatch {
  int score{kNegInf};
  int first{-1};
  int last{-1};
  std::vector<int> positions;
};

static TokenMatch matchToken(std::string_view pat, std::string_view text, std::size_t start) {
  TokenMatch out{};
  pat = trimAscii(pat);
  if (pat.empty()) {
    out.score = 0;
    return out;
  }

  if (start >= text.size()) return out;

  const std::size_t n = text.size() - start;
  const std::size_t m = pat.size();
  if (m > n) return out;

  // Scoring constants (integer scores for determinism).
  // The shape is inspired by common fuzzy finders: strong preference for
  // consecutive matches and boundaries, mild penalty for gaps.
  constexpr int kMatch = 12;
  constexpr int kGap = -1;
  constexpr int kConsecutive = 14;
  constexpr int kBoundary = 14;
  constexpr int kCamel = 14;
  constexpr int kUpperText = 6;
  constexpr int kExactCase = 2;

  std::vector<char> lowerText(n);
  std::vector<int> bonus(n, 0);
  for (std::size_t j = 0; j < n; ++j) {
    const char c = text[start + j];
    lowerText[j] = lowerAscii(c);

    const bool atStart = (start + j) == 0;
    const char prev = atStart ? '\0' : text[start + j - 1];

    int b = 0;
    if (atStart || isWordBoundary(prev)) b += kBoundary;
    if (!atStart && isCamelBoundary(prev, c)) b += kCamel;
    bonus[j] = b;
  }

  // DP tables for last row only (O(m*n) time, O(n) working memory).
  // We keep a full prev-index table for backtracking the best path.
  std::vector<int> prevIdx(m * n, -1);
  auto prevAt = [&](std::size_t i, std::size_t j) -> int& { return prevIdx[i * n + j]; };

  std::vector<int> dPrev(n, kNegInf);
  std::vector<int> mPrev(n, kNegInf);
  std::vector<int> mPrevIdx(n, -1);

  std::vector<int> dCur(n, kNegInf);
  std::vector<int> mCur(n, kNegInf);
  std::vector<int> mCurIdx(n, -1);

  // ---- Row 0 (first pattern char) ----
  const char p0 = lowerAscii(pat[0]);
  for (std::size_t j = 0; j < n; ++j) {
    if (lowerText[j] == p0) {
      int caseBonus = 0;
      if (isAsciiUpper(text[start + j])) caseBonus += kUpperText;
      if (pat[0] == text[start + j]) caseBonus += kExactCase;
      // Leading gap penalty: prefer earlier matches.
      dCur[j] = addScore(kMatch + bonus[j] + caseBonus, (int)j * kGap);
      prevAt(0, j) = -1;
    }

    const int skip = (j > 0) ? addScore(mCur[j - 1], kGap) : kNegInf;

    // Tie-break: prefer "matching now" if equal, so later matches win ties,
    // which tends to be better for future characters.
    if (dCur[j] >= skip) {
      mCur[j] = dCur[j];
      mCurIdx[j] = (dCur[j] > kNegInf / 2) ? (int)j : -1;
    } else {
      mCur[j] = skip;
      mCurIdx[j] = (j > 0) ? mCurIdx[j - 1] : -1;
    }
  }

  dPrev.swap(dCur);
  mPrev.swap(mCur);
  mPrevIdx.swap(mCurIdx);
  std::fill(dCur.begin(), dCur.end(), kNegInf);
  std::fill(mCur.begin(), mCur.end(), kNegInf);
  std::fill(mCurIdx.begin(), mCurIdx.end(), -1);

  // ---- Rows 1..m-1 ----
  for (std::size_t i = 1; i < m; ++i) {
    const char pc = lowerAscii(pat[i]);

    for (std::size_t j = 0; j < n; ++j) {
      if (lowerText[j] == pc) {
        int caseBonus = 0;
        if (isAsciiUpper(text[start + j])) caseBonus += kUpperText;
        if (pat[i] == text[start + j]) caseBonus += kExactCase;

        int best = kNegInf;
        int bestPrev = -1;

        // Non-consecutive transition from the best previous match (with gap penalty).
        if (j > 0 && mPrev[j - 1] > kNegInf / 2) {
          const int cand = addScore(mPrev[j - 1], kGap);
          best = cand;
          bestPrev = mPrevIdx[j - 1];
        }

        // Consecutive transition (prefer on ties).
        if (j > 0 && dPrev[j - 1] > kNegInf / 2) {
          const int cand = addScore(dPrev[j - 1], kConsecutive);
          if (cand >= best) {
            best = cand;
            bestPrev = (int)(j - 1);
          }
        }

        if (best > kNegInf / 2) {
          dCur[j] = addScore(best, kMatch + bonus[j] + caseBonus);
          prevAt(i, j) = bestPrev;
        }
      }

      const int skip = (j > 0) ? addScore(mCur[j - 1], kGap) : kNegInf;
      if (dCur[j] >= skip) {
        mCur[j] = dCur[j];
        mCurIdx[j] = (dCur[j] > kNegInf / 2) ? (int)j : -1;
      } else {
        mCur[j] = skip;
        mCurIdx[j] = (j > 0) ? mCurIdx[j - 1] : -1;
      }
    }

    dPrev.swap(dCur);
    mPrev.swap(mCur);
    mPrevIdx.swap(mCurIdx);
    std::fill(dCur.begin(), dCur.end(), kNegInf);
    std::fill(mCur.begin(), mCur.end(), kNegInf);
    std::fill(mCurIdx.begin(), mCurIdx.end(), -1);
  }

  // Pick the best end position for the last character.
  int bestScore = kNegInf;
  int bestJ = -1;
  for (std::size_t j = 0; j < n; ++j) {
    const int s = dPrev[j];
    if (s > bestScore) {
      bestScore = s;
      bestJ = (int)j;
    } else if (s == bestScore && s > kNegInf / 2 && bestJ >= 0 && (int)j < bestJ) {
      // Tie-break: prefer earlier end for a more "prefixy" highlight.
      bestJ = (int)j;
    }
  }

  if (bestJ < 0 || bestScore <= kNegInf / 2) return out;

  // Backtrack.
  out.positions.resize(m);
  int j = bestJ;
  for (int i = (int)m - 1; i >= 0; --i) {
    out.positions[(std::size_t)i] = (int)(start + (std::size_t)j);
    if (i > 0) {
      const int pj = prevAt((std::size_t)i, (std::size_t)j);
      if (pj < 0) return TokenMatch{}; // should not happen; treat as no match.
      j = pj;
    }
  }

  out.score = bestScore;
  out.first = out.positions.front();
  out.last = out.positions.back();
  return out;
}

} // namespace

FuzzyMatchResult fuzzyMatch(std::string_view query, std::string_view text) {
  FuzzyMatchResult out{};
  query = trimAscii(query);
  if (query.empty()) {
    out.score = 0;
    return out;
  }

  // Tokenize query on whitespace.
  std::vector<std::string_view> tokens;
  tokens.reserve(4);
  {
    std::size_t b = 0;
    bool inTok = false;
    for (std::size_t i = 0; i < query.size(); ++i) {
      if (isSpace((unsigned char)query[i])) {
        if (inTok) {
          pushToken(tokens, query, b, i);
          inTok = false;
        }
      } else if (!inTok) {
        b = i;
        inTok = true;
      }
    }
    if (inTok) pushToken(tokens, query, b, query.size());
  }

  if (tokens.empty()) {
    out.score = 0;
    return out;
  }

  int totalScore = 0;
  std::size_t searchStart = 0;
  out.positions.clear();
  out.positions.reserve(query.size());

  for (std::string_view tok : tokens) {
    if (tok.empty()) continue;

    const TokenMatch m = matchToken(tok, text, searchStart);
    if (m.score <= kNegInf / 2) {
      out.score = -1;
      out.positions.clear();
      return out;
    }

    totalScore += m.score;
    out.positions.insert(out.positions.end(), m.positions.begin(), m.positions.end());
    if (m.last >= 0) searchStart = (std::size_t)m.last + 1;
  }

  // Ensure positions are sorted and unique.
  std::sort(out.positions.begin(), out.positions.end());
  out.positions.erase(std::unique(out.positions.begin(), out.positions.end()), out.positions.end());

  out.score = totalScore;
  return out;
}

} // namespace stellar::ui
