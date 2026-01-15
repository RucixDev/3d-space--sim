#pragma once

#include "stellar/render/RenderTarget.h"
#include "stellar/render/ShaderToy.h"

#include <array>
#include <string>
#include <string_view>

namespace stellar::render {

// A tiny multi-pass wrapper around ShaderToy that mimics the common ShaderToy
// workflow:
//  - 4 feedback buffers (A..D) rendered in order each frame
//  - an Image pass rendered last
//
// Each pass can sample up to 4 input channels (iChannel0..3). A channel can be
// unbound (None) or routed to one of the buffers.
//
// Feedback is handled naturally:
//  - When a buffer samples itself, it receives the previous frame (ping buffer)
//  - When a buffer samples an earlier buffer (A before B, etc.), it sees the
//    freshly rendered current-frame result
//  - When a buffer samples a later buffer, it sees the previous frame

enum class ShaderToyPass : int {
  BufferA = 0,
  BufferB = 1,
  BufferC = 2,
  BufferD = 3,
  Image   = 4,
};

enum class ShaderToyChannelSource : int {
  None    = 0,
  BufferA = 1,
  BufferB = 2,
  BufferC = 3,
  BufferD = 4,
};

struct ShaderToyGraphPassSettings {
  bool enabled{true};

  // Resolution divisor for this pass (1 = full, 2 = half, 4 = quarter).
  int resolutionScale{1};

  std::array<ShaderToyChannelSource, 4> channels{
    ShaderToyChannelSource::None,
    ShaderToyChannelSource::None,
    ShaderToyChannelSource::None,
    ShaderToyChannelSource::None,
  };
};

class ShaderToyGraph {
public:
  ShaderToyGraph() = default;
  ~ShaderToyGraph() = default;

  ShaderToyGraph(const ShaderToyGraph&) = delete;
  ShaderToyGraph& operator=(const ShaderToyGraph&) = delete;

  bool init(std::string* outError = nullptr);

  // Ensure internal buffer render targets are allocated for the given output size.
  bool ensureSize(int outWidth, int outHeight, std::string* outError = nullptr);

  // Build / compile the shader for a specific pass.
  bool buildPass(ShaderToyPass pass, std::string_view userCode, std::string* outError = nullptr);

  // Optional wrapper header inserted into every pass right before the snippet.
  // Used by the Procedural Shader Lab to inject parsed user parameter uniforms.
  void setUserHeader(std::string_view extraHeader);
  const std::string& userHeader() const { return userHeader_; }

  // Configure per-pass routing / resolution.
  void setPassSettings(ShaderToyPass pass, const ShaderToyGraphPassSettings& s);
  const ShaderToyGraphPassSettings& passSettings(ShaderToyPass pass) const;

  bool passValid(ShaderToyPass pass) const;
  const std::string& passError(ShaderToyPass pass) const;

  // Clears all feedback buffers (both ping/pong) to black.
  void resetBuffers();

  // Render the full graph for a frame.
  // The Image pass is rendered into outImageTarget.
  void render(const ShaderToyInputs& baseInputs, RenderTarget2D& outImageTarget);

  // Access the current output for a buffer pass (A..D).
  const RenderTarget2D* bufferTarget(ShaderToyPass bufferPass) const;
  const Texture2D* bufferTexture(ShaderToyPass bufferPass) const;

  const ShaderToy& shader(ShaderToyPass pass) const;

private:
  static constexpr int kBufferCount = 4;
  static constexpr int kPassCount = 5;

  static int idx(ShaderToyPass p) { return static_cast<int>(p); }
  static bool isBuffer(ShaderToyPass p) { return idx(p) >= 0 && idx(p) < kBufferCount; }

  struct PassState {
    ShaderToyGraphPassSettings settings{};
    ShaderToy toy{};
    bool valid{false};
    std::string error{};

    // Double-buffered RTs for feedback buffers (A..D). Unused for Image.
    RenderTarget2D rt[2]{};
    bool rtInited[2]{false, false};
    int ping{0}; // index of the currently readable RT
  };

  PassState passes_[kPassCount]{};

  bool inited_{false};

  // Shared extra header inserted into every pass wrapper.
  // (e.g. auto-injected user parameter uniforms)
  std::string userHeader_{};

  int outW_{1};
  int outH_{1};

  // Extra wrapper header appended before each pass' user code.
  std::string userHeader_{};

  void clearRtIfInited_(PassState& p);
  const Texture2D* resolveChannelTexture_(ShaderToyPass currentPass, ShaderToyChannelSource src) const;
};

} // namespace stellar::render
