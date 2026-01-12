#include "stellar/render/TextureReadback.h"

#include "stellar/render/Gl.h"

#include <SDL_opengl.h>

#include <algorithm>
#include <cstring>

namespace stellar::render {

static int channelsForFormat(unsigned int fmt) {
  switch (fmt) {
    case GL_RGBA:
    case GL_BGRA:
      return 4;
    case GL_RGB:
    case GL_BGR:
      return 3;
    case GL_RG:
      return 2;
    case GL_RED:
    case GL_LUMINANCE:
      return 1;
    default:
      return 4;
  }
}

static std::size_t bytesPerComponent(unsigned int type) {
  switch (type) {
    case GL_UNSIGNED_BYTE:
    case GL_BYTE:
      return 1;
    case GL_UNSIGNED_SHORT:
    case GL_SHORT:
    case GL_HALF_FLOAT:
      return 2;
    case GL_UNSIGNED_INT:
    case GL_INT:
    case GL_FLOAT:
      return 4;
    default:
      return 4;
  }
}

AsyncTextureReadback::~AsyncTextureReadback() {
  shutdown();
}

bool AsyncTextureReadback::init(std::string* outError) {
  // We require the PBO mapping entry points.
  supported_ = (gl::MapBufferRange != nullptr) && (gl::UnmapBuffer != nullptr);
  if (!supported_ && outError) {
    *outError = "Async readback unsupported (missing glMapBufferRange/glUnmapBuffer).";
  }
  return supported_;
}

void AsyncTextureReadback::shutdown() {
  clear();
  supported_ = false;
}

void AsyncTextureReadback::clear() {
  if (!jobs_.empty()) {
    for (auto& j : jobs_) {
      if (j.pbo) {
        unsigned int pbo = j.pbo;
        gl::DeleteBuffers(1, &pbo);
        j.pbo = 0;
      }
    }
    jobs_.clear();
  }
}

bool AsyncTextureReadback::enqueue(const TextureReadbackDesc& desc, std::string tag, int delayFrames) {
  if (!supported_) return false;
  if (desc.texture == 0) return false;
  if (desc.width <= 0 || desc.height <= 0) return false;
  if (desc.format == 0 || desc.type == 0) return false;

  delayFrames = std::max(0, delayFrames);

  const int comp = channelsForFormat(desc.format);
  const std::size_t bpc = bytesPerComponent(desc.type);
  const std::size_t sizeBytes =
      (std::size_t)desc.width * (std::size_t)desc.height * (std::size_t)comp * bpc;
  if (sizeBytes == 0) return false;

  unsigned int pbo = 0;
  gl::GenBuffers(1, &pbo);
  if (!pbo) return false;

  // Preserve pack state (debug tools shouldn't subtly break other readbacks).
  GLint prevAlign = 4;
  GLint prevRowLength = 0;
  GLint prevSkipRows = 0;
  GLint prevSkipPixels = 0;
  glGetIntegerv(GL_PACK_ALIGNMENT, &prevAlign);
#ifdef GL_PACK_ROW_LENGTH
  glGetIntegerv(GL_PACK_ROW_LENGTH, &prevRowLength);
  glGetIntegerv(GL_PACK_SKIP_ROWS, &prevSkipRows);
  glGetIntegerv(GL_PACK_SKIP_PIXELS, &prevSkipPixels);
#endif

  glPixelStorei(GL_PACK_ALIGNMENT, 1);
#ifdef GL_PACK_ROW_LENGTH
  glPixelStorei(GL_PACK_ROW_LENGTH, 0);
  glPixelStorei(GL_PACK_SKIP_ROWS, 0);
  glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
#endif

  gl::BindBuffer(GL_PIXEL_PACK_BUFFER, pbo);
  gl::BufferData(GL_PIXEL_PACK_BUFFER, (GLsizeiptr)sizeBytes, nullptr, GL_STREAM_READ);

  gl::BindTexture(GL_TEXTURE_2D, desc.texture);
  glGetTexImage(GL_TEXTURE_2D, desc.level, (GLenum)desc.format, (GLenum)desc.type, (void*)0);
  gl::BindTexture(GL_TEXTURE_2D, 0);
  gl::BindBuffer(GL_PIXEL_PACK_BUFFER, 0);

  // Restore pack state.
  glPixelStorei(GL_PACK_ALIGNMENT, prevAlign);
#ifdef GL_PACK_ROW_LENGTH
  glPixelStorei(GL_PACK_ROW_LENGTH, prevRowLength);
  glPixelStorei(GL_PACK_SKIP_ROWS, prevSkipRows);
  glPixelStorei(GL_PACK_SKIP_PIXELS, prevSkipPixels);
#endif

  Job j;
  j.desc = desc;
  j.tag = std::move(tag);
  j.pbo = pbo;
  j.sizeBytes = sizeBytes;
  j.delay = delayFrames;

  jobs_.push_back(std::move(j));
  return true;
}

bool AsyncTextureReadback::poll(TextureReadbackResult& out) {
  out = {};
  if (!supported_) return false;
  if (jobs_.empty()) return false;

  for (std::size_t i = 0; i < jobs_.size(); ++i) {
    Job& j = jobs_[i];

    if (j.delay > 0) {
      --j.delay;
      continue;
    }

    if (!j.pbo || j.sizeBytes == 0) {
      // Corrupt job; drop it.
      if (j.pbo) {
        unsigned int pbo = j.pbo;
        gl::DeleteBuffers(1, &pbo);
      }
      jobs_.erase(jobs_.begin() + (std::ptrdiff_t)i);
      return false;
    }

    gl::BindBuffer(GL_PIXEL_PACK_BUFFER, j.pbo);
    void* ptr = gl::MapBufferRange(GL_PIXEL_PACK_BUFFER, 0, (GLsizeiptr)j.sizeBytes, GL_MAP_READ_BIT);
    if (!ptr) {
      // Not ready (or mapping failed). Try again later.
      gl::BindBuffer(GL_PIXEL_PACK_BUFFER, 0);
      continue;
    }

    out.desc = j.desc;
    out.tag = j.tag;

    // Only GL_FLOAT is returned as a float buffer. Everything else is returned as raw bytes.
    if (j.desc.type == GL_FLOAT) {
      const std::size_t nFloats = j.sizeBytes / sizeof(float);
      out.floats.resize(nFloats);
      std::memcpy(out.floats.data(), ptr, nFloats * sizeof(float));
    } else {
      out.bytes.resize(j.sizeBytes);
      std::memcpy(out.bytes.data(), ptr, j.sizeBytes);
    }

    gl::UnmapBuffer(GL_PIXEL_PACK_BUFFER);
    gl::BindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    unsigned int pbo = j.pbo;
    gl::DeleteBuffers(1, &pbo);
    jobs_.erase(jobs_.begin() + (std::ptrdiff_t)i);
    return true;
  }

  return false;
}

} // namespace stellar::render
