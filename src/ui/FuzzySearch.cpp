#include "stellar/ui/FuzzySearch.h"

#include <algorithm>
#include <cctype>
#include <string>

namespace stellar::ui {

static bool isSpace(unsigned char c) { return std::isspace(c) != 0; }

static bool isWordBoundary(char c) {
  return c == ' ' || c == '_' || c == '-' || c == '/' || c == '\\' || c == '.' || c == ':' || c == '(' || c == ')'
         || c == '[' || c == ']' || c == '{' || c == '}' || c == '#';
}

static char lowerAscii(char c) {
  return (char)std::tolower((unsigned char)c);
}

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

  // ASCII-lowercase view of text is computed on the fly to avoid allocations.
  const auto getLower = [&](std::size_t i) -> char { return lowerAscii(text[i]); };

  int totalScore = 0;
  std::size_t searchStart = 0;
  out.positions.clear();
  out.positions.reserve(query.size());

  for (std::string_view tok : tokens) {
    if (tok.empty()) continue;

    int score = 0;
    int first = -1;
    int last = -1;
    int consecutive = 0;
    std::size_t ti = searchStart;

    for (std::size_t qi = 0; qi < tok.size(); ++qi) {
      const char qc = lowerAscii(tok[qi]);

      // Find next matching char.
      std::size_t pos = std::string::npos;
      for (std::size_t j = ti; j < text.size(); ++j) {
        if (getLower(j) == qc) {
          pos = j;
          break;
        }
      }
      if (pos == std::string::npos) {
        out.score = -1;
        out.positions.clear();
        return out;
      }

      out.positions.push_back((int)pos);

      if (first < 0) first = (int)pos;

      if (last >= 0 && (int)pos == last + 1) {
        ++consecutive;
        score += 18;
      } else {
        consecutive = 0;
        score += 6;
      }

      // Word boundary bonus (start of string or after separator).
      if (pos == 0 || isWordBoundary(text[pos - 1])) score += 12;

      last = (int)pos;
      ti = pos + 1;
    }

    if (first >= 0 && last >= 0) {
      score -= (last - first); // spread penalty
      if (first == 0) score += 30; // prefix bonus
    }

    searchStart = (last >= 0) ? (std::size_t)(last + 1) : searchStart;
    totalScore += score;
  }

  // Ensure positions are sorted and unique (can duplicate when tokens overlap in weird input).
  std::sort(out.positions.begin(), out.positions.end());
  out.positions.erase(std::unique(out.positions.begin(), out.positions.end()), out.positions.end());

  out.score = totalScore;
  return out;
}

} // namespace stellar::ui
