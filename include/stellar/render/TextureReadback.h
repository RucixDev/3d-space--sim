#pragma once

// Async texture readback helper for debug tooling.
//
// This is intentionally small and OpenGL-centric (desktop GL). It uses
// GL_PIXEL_PACK_BUFFER + glGetTexImage to stage GPU->CPU transfers into a PBO,
// then maps the buffer a few frames later to retrieve pixels.

#include <cstddef>
#include <string>
#include <vector>

namespace stellar::render {

// Describes a single texture readback.
//
// `format` and `type` correspond to the `format` and `type` parameters you would
// pass to `glGetTexImage` (or `glReadPixels`).
struct TextureReadbackDesc {
  unsigned int texture{0};  // GL texture name (assumed GL_TEXTURE_2D)
  int width{0};
  int height{0};
  unsigned int format{0};   // GLenum, e.g. GL_RGBA / GL_RG
  unsigned int type{0};     // GLenum, e.g. GL_UNSIGNED_BYTE / GL_FLOAT
  int level{0};             // mip level (usually 0)
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
// This is designed for *debug capture* workflows (FrameGraph inspector, etc.).
// It intentionally trades correctness guarantees for simplicity (no fences):
// jobs are delayed by N frames before mapping the PBO to reduce the chance of
// stalling on the GPU.
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

  // Enqueue a readback job. Returns false if unsupported or if params are invalid.
  //
  // `delayFrames` controls how many polls must occur before the job is eligible
  // to map. Values in the 2-4 range tend to avoid big stalls on most drivers.
  bool enqueue(const TextureReadbackDesc& desc, std::string tag, int delayFrames = 3);

  // Poll the queue for a completed job.
  // Returns true and fills `out` when a job is available.
  bool poll(TextureReadbackResult& out);

  std::size_t pendingJobs() const { return jobs_.size(); }

  // Clear any queued jobs (deleting their PBOs).
  void clear();

private:
  struct Job {
    TextureReadbackDesc desc{};
    std::string tag{};
    unsigned int pbo{0};
    std::size_t sizeBytes{0};
    int delay{0};
  };

  std::vector<Job> jobs_{};
  bool supported_{false};
};

} // namespace stellar::render
