#pragma once

#include <string>
#include <string_view>

namespace stellar::game {

// Parameters used to generate a screenshot file path.
struct ScreenshotRequest {
  // If true, the intent is to capture after UI rendering.
  // (The request callback decides which stage to capture at.)
  bool includeUi{false};

  // Output directory (relative or absolute).
  std::string outDir{"screenshots"};

  // Base file name without extension.
  std::string baseName{"stellarforge"};

  // If true, append a timestamp to the file name.
  bool timestamp{true};
};

// Keep only filename-safe ASCII characters: [A-Za-z0-9_-].
// Whitespace is converted to '_'.
// If the result is empty, returns "shot".
std::string sanitizeFileToken(std::string_view s);

// Builds a unique output file path for a screenshot and ensures the output
// directory exists (creates it if needed).
//
// Returns an empty string on failure and optionally writes a human-readable
// error message to outErr.
std::string buildScreenshotPath(const ScreenshotRequest& req, std::string* outErr = nullptr);

// Captures the current backbuffer (GL_BACK) and writes it as a PNG.
//
// `width/height` are framebuffer pixel dimensions.
// Returns true on success.
bool captureBackbufferToPng(const std::string& path, int width, int height, std::string* outErr = nullptr);

} // namespace stellar::game
