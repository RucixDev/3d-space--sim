#include "Screenshot.h"

#include <SDL_opengl.h>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <vector>

// stb_image_write is a tiny single-header image writer (PNG, BMP, TGA, ...).
// We vendor it under apps/stellar_game/third_party to keep the core library
// headless and dependency-free.
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STBIW_ASSERT(x) ((void)0)
#include "third_party/stb_image_write.h"

namespace stellar::game {

std::string sanitizeFileToken(std::string_view s) {
  std::string out;
  out.reserve(s.size());

  for (unsigned char ch : s) {
    if (std::isalnum(ch)) {
      out.push_back((char)ch);
      continue;
    }
    if (ch == ' ' || ch == '\t') {
      out.push_back('_');
      continue;
    }
    if (ch == '_' || ch == '-') {
      out.push_back((char)ch);
      continue;
    }
    // Skip other characters (slashes, colons, unicode, etc).
  }

  // Collapse repeated underscores.
  {
    std::string tmp;
    tmp.reserve(out.size());
    char prev = '\0';
    for (char c : out) {
      if (c == '_' && prev == '_') continue;
      tmp.push_back(c);
      prev = c;
    }
    out.swap(tmp);
  }

  // Trim underscores.
  while (!out.empty() && out.front() == '_') out.erase(out.begin());
  while (!out.empty() && out.back() == '_') out.pop_back();

  if (out.empty()) out = "shot";
  return out;
}

static std::string fileTimestampNow() {
  using clock = std::chrono::system_clock;
  const auto now = clock::now();
  const auto t = clock::to_time_t(now);

  std::tm tm{};
#if defined(_WIN32)
  localtime_s(&tm, &t);
#else
  localtime_r(&t, &tm);
#endif

  const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y%m%d_%H%M%S")
      << '_' << std::setw(3) << std::setfill('0') << ms.count();
  return oss.str();
}

std::string buildScreenshotPath(const ScreenshotRequest& req, std::string* outErr) {
  namespace fs = std::filesystem;

  const std::string dir = req.outDir.empty() ? std::string("screenshots") : req.outDir;
  std::error_code ec;
  if (!fs::exists(dir, ec)) {
    if (!fs::create_directories(dir, ec)) {
      if (outErr) *outErr = "Failed to create directory: " + dir;
      return {};
    }
  }

  const std::string base = sanitizeFileToken(req.baseName);
  std::string file = base;
  if (req.timestamp) {
    file += '_';
    file += fileTimestampNow();
  }

  // Ensure uniqueness by appending _N if needed.
  fs::path p = fs::path(dir) / (file + ".png");
  if (fs::exists(p, ec)) {
    for (int i = 2; i < 10000; ++i) {
      fs::path cand = fs::path(dir) / (file + "_" + std::to_string(i) + ".png");
      if (!fs::exists(cand, ec)) {
        p = cand;
        break;
      }
    }
  }

  return p.string();
}

bool captureBackbufferToPng(const std::string& path, int width, int height, std::string* outErr) {
  if (width <= 0 || height <= 0) {
    if (outErr) *outErr = "Invalid framebuffer size.";
    return false;
  }
  if (path.empty()) {
    if (outErr) *outErr = "Empty output path.";
    return false;
  }

  // Read RGBA8 from the backbuffer.
  std::vector<unsigned char> pixels;
  pixels.resize((std::size_t)width * (std::size_t)height * 4u);

  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadBuffer(GL_BACK);
  glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());

  // Flip vertically: OpenGL origin is bottom-left, images are typically top-left.
  const int stride = width * 4;
  std::vector<unsigned char> flipped;
  flipped.resize(pixels.size());

  for (int y = 0; y < height; ++y) {
    const unsigned char* src = pixels.data() + (std::size_t)(height - 1 - y) * (std::size_t)stride;
    unsigned char* dst = flipped.data() + (std::size_t)y * (std::size_t)stride;
    std::copy(src, src + stride, dst);
  }

  const int ok = stbi_write_png(path.c_str(), width, height, 4, flipped.data(), stride);
  if (!ok) {
    if (outErr) *outErr = "Failed to write PNG.";
    return false;
  }
  return true;
}

} // namespace stellar::game
