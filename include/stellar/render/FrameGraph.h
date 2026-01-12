#pragma once

// FrameGraph: a tiny, hackable render-graph style pass scheduler for OpenGL.
//
// Goals:
//  - Make it easy to express post-processing as a graph of passes with explicit
//    input/output textures.
//  - Provide simple compile-time dependency resolution (topo sort).
//  - Provide optional transient texture aliasing based on resource lifetimes
//    (classic "frame graph" memory reuse).
//
// This header intentionally avoids including any OpenGL headers. All GL calls
// live in the .cpp, and GL object handles are represented as unsigned int.

#include <array>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::render {

class FrameGraph {
public:
  struct TextureDesc {
    // If width/height are 0, they are derived from the graph viewport
    // (viewportW * scale, viewportH * scale).
    float scale{1.0f};
    int width{0};
    int height{0};

    // OpenGL allocation params (passed straight to glTexImage2D).
    int internalFormat{0};
    unsigned int format{0};
    unsigned int type{0};

    // Sampling params.
    bool linearFilter{true};
    bool clampToEdge{true};
  };

  struct TextureHandle {
    int id{-1};
    constexpr bool valid() const { return id >= 0; }
    constexpr bool isBackbuffer() const { return id == kBackbufferId; }
  };

  struct PassInfo {
    std::string name{};
    std::vector<TextureHandle> reads{};
    TextureHandle write{}; // Backbuffer() allowed.
  };

  struct TextureInfo {
    std::string name{};
    bool external{false};
    TextureDesc desc{};

    // Resolved size (after viewport scaling). For external textures, this is
    // the imported size.
    int resolvedW{0};
    int resolvedH{0};

    // Lifetime in scheduled pass order (inclusive).
    int firstUse{-1};
    int lastUse{-1};

    // Transient aliasing result.
    int physicalId{-1};

    // For external textures, the imported GL handle. For internal textures,
    // this is populated at execute() time and may be 0 until then.
    unsigned int glHandle{0};
  };

  struct PassContext {
    const FrameGraph& fg;
    int baseW{0};
    int baseH{0};
    int outW{0};
    int outH{0};

    unsigned int texture(TextureHandle h) const;
    int textureWidth(TextureHandle h) const;
    int textureHeight(TextureHandle h) const;
  };

  using PassCallback = std::function<void(const PassContext&)>;

  FrameGraph();
  ~FrameGraph();

  FrameGraph(const FrameGraph&) = delete;
  FrameGraph& operator=(const FrameGraph&) = delete;

  // Clears passes + textures for building a new graph, but keeps any allocated
  // GL objects for reuse.
  void reset();

  void setViewport(int w, int h);
  int viewportWidth() const { return viewW_; }
  int viewportHeight() const { return viewH_; }

  // Enable/disable transient aliasing. When disabled, each internal texture
  // receives a unique physical allocation (useful for debugging pass outputs).
  void setAliasingEnabled(bool v) { aliasingEnabled_ = v; }
  bool aliasingEnabled() const { return aliasingEnabled_; }

  // GPU profiling (timer queries). When enabled and supported by the GL context,
  // FrameGraph measures GPU time per pass using GL_TIME_ELAPSED queries.
  void setProfilingEnabled(bool v) { profilingEnabled_ = v; }
  bool profilingEnabled() const { return profilingEnabled_; }
  bool profilingSupported() const { return profilingSupported_; }

  // Last read-back GPU timings in milliseconds, indexed by pass id.
  const std::vector<double>& lastPassTimesMs() const { return lastPassTimesMs_; }
  double lastFrameTimeMs() const { return lastFrameTimeMs_; }

  // Create an internal transient texture.
  TextureHandle createTexture(std::string name, const TextureDesc& desc);

  // Import an external GL texture into the graph.
  TextureHandle importTexture(std::string name, unsigned int glHandle, int w, int h);

  // Add a pass.
  //
  // - reads:  input textures (may include external textures)
  // - write:  output texture (Backbuffer() allowed)
  //
  // The callback will be called with the pass output bound as the current
  // framebuffer and the viewport set to the output size.
  int addPass(std::string name,
              std::vector<TextureHandle> reads,
              TextureHandle write,
              PassCallback cb);

  // Compile the graph (topo sort + lifetime analysis + optional aliasing).
  bool compile(std::string* outError = nullptr);

  // Execute the compiled graph (requires a valid OpenGL context and that
  // stellar::render::gl::load() has been called).
  void execute();

  // Debug / introspection
  const std::vector<TextureInfo>& textures() const { return textures_; }
  const std::vector<PassInfo>& passes() const { return passes_; }
  const std::vector<int>& schedule() const { return schedule_; }

  // Number of physical transient allocations required for the last compile.
  int physicalTextureCount() const { return physRequired_; }

  // Returns the GL texture handle for a graph texture (0 if not available).
  unsigned int glTexture(TextureHandle h) const;

  static constexpr TextureHandle Backbuffer() { return TextureHandle{kBackbufferId}; }

private:
  static constexpr int kBackbufferId = -2;

  struct ResolvedDesc {
    int w{0};
    int h{0};
    int internalFormat{0};
    unsigned int format{0};
    unsigned int type{0};
    bool linearFilter{true};
    bool clampToEdge{true};
  };

  struct PhysicalTexture {
    ResolvedDesc desc{};
    unsigned int tex{0};
    unsigned int fbo{0};
  };

  ResolvedDesc resolveDesc(const TextureDesc& d) const;
  static bool sameDesc(const ResolvedDesc& a, const ResolvedDesc& b);

  void destroyGL();
  void ensurePhysicalAllocated();

  unsigned int ensureGlTexture(int physId, const ResolvedDesc& d);
  unsigned int ensureFbo(int physId);

  int viewW_{0};
  int viewH_{0};

  bool aliasingEnabled_{true};

  // GPU profiling (timer query ring-buffer). This is optional and defaults off.
  bool profilingEnabled_{false};
  bool profilingSupported_{false};
  static constexpr int kQueryFrames = 3;
  int queryFrame_{0};
  std::array<std::vector<unsigned int>, kQueryFrames> passQueries_{};
  std::vector<double> lastPassTimesMs_{};
  double lastFrameTimeMs_{0.0};

  std::vector<TextureInfo> textures_;
  std::vector<PassInfo> passes_;
  std::vector<PassCallback> callbacks_;
  std::vector<int> schedule_;

  // Compile output: how many physical textures are required and which desc each
  // physical allocation should have.
  int physRequired_{0};
  std::vector<ResolvedDesc> physPlan_;

  // GL objects (lazy allocated on execute()).
  std::vector<PhysicalTexture> phys_;
  bool glAllocated_{false};
};

} // namespace stellar::render
