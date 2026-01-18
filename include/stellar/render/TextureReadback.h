#pragma once

// Async texture/framebuffer readback helper for debug tooling.
//
// This is intentionally small and OpenGL-centric (desktop GL). It uses
// GL_PIXEL_PACK_BUFFER to stage GPU->CPU transfers into a PBO and, when
// available, GL sync objects (glFenceSync/glClientWaitSync) to avoid stalling
// when mapping.
//
// Supported sources:
//  - GL_TEXTURE_2D via glGetTexImage
//  - framebuffer/backbuffer via glReadPixels
//
// Notes:
//  - The API uses unsigned int for GL object/enums to keep this header light.

#include <cstddef>
#include <string>
#include <vector>

namespace stellar::render {

enum class ReadbackSource : unsigned char {
  Texture2D = 0,
  Framebuffer = 1,
};

// Describes a single readback.
//
// `format` and `type` correspond to the `format` and `type` parameters you would
// pass to glGetTexImage / glReadPixels.
struct TextureReadbackDesc {
  ReadbackSource source{ReadbackSource::Texture2D};

  // --- Texture2D source ---
  unsigned int texture{0};  // GL texture name (assumed GL_TEXTURE_2D)
  int level{0};             // mip level (usually 0)

  // --- Common ---
  int width{0};
  int height{0};
  unsigned int format{0};   // GLenum, e.g. GL_RGBA / GL_RG
  unsigned int type{0};     // GLenum, e.g. GL_UNSIGNED_BYTE / GL_FLOAT

  // --- Framebuffer source ---
  // `framebuffer` is the read framebuffer (0 = default/backbuffer). If
  // `readBuffer` is 0, it defaults to GL_BACK for the default framebuffer, and
  // GL_COLOR_ATTACHMENT0 for non-zero FBOs.
  unsigned int framebuffer{0};
  unsigned int readBuffer{0};
  int x{0};
  int y{0};
};

// Readback result.
//
// Exactly one of `bytes` or `floats` will be populated depending on `desc.type`.
struct TextureReadbackResult {
  TextureReadbackDesc desc{};
  std::string tag{}; // user-defined (e.g. a file path)
  std::vector<unsigned char> bytes{};
  std::vector<float> floats{};
};

// A small, self-contained async readback queue.
//
// This is designed for debug capture workflows (FrameGraph inspector,
// screenshots, etc.). When GL sync objects are available, we use a fence per job
// and poll it with a zero-timeout wait to avoid blocking. When unavailable,
// jobs are delayed by N frames before mapping the PBO (best-effort).
class AsyncTextureReadback {
public:
  AsyncTextureReadback() = default;
  ~AsyncTextureReadback();

  AsyncTextureReadback(const AsyncTextureReadback&) = delete;
  AsyncTextureReadback& operator=(const AsyncTextureReadback&) = delete;

  // Must be called after the OpenGL function loader (stellar::render::gl::load).
  // Returns false if required entry points for mapping are missing.
  bool init(std::string* outError = nullptr);

  void shutdown();

  bool supported() const { return supported_; }
  bool fencesSupported() const { return fencesSupported_; }

  // Enqueue a readback job. Returns false if unsupported or if params are invalid.
  //
  // `delayFrames` controls how many polls must occur before the job is eligible
  // to map. Values in the 1-3 range tend to avoid stalls on most drivers.
  bool enqueue(const TextureReadbackDesc& desc, std::string tag, int delayFrames = 2);

  // Poll the queue for a completed job.
  // Returns true and fills `out` when a job is available.
  bool poll(TextureReadbackResult& out);

  std::size_t pendingJobs() const { return jobs_.size(); }

  // Clear any queued jobs (deleting their PBOs and fences).
  void clear();

private:
  struct Job {
    TextureReadbackDesc desc{};
    std::string tag{};
    unsigned int pbo{0};
    std::size_t sizeBytes{0};
    int delay{0};

    // GLsync when available (stored as void* to keep the header GL-light).
    void* fence{nullptr};
  };

  std::vector<Job> jobs_{};
  bool supported_{false};
  bool fencesSupported_{false};
};

} // namespace stellar::render
