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

static void deleteFence(void* fence) {
  if (!fence) return;
  if (gl::DeleteSync) {
    gl::DeleteSync(reinterpret_cast<GLsync>(fence));
  }
}

AsyncTextureReadback::~AsyncTextureReadback() {
  shutdown();
}

bool AsyncTextureReadback::init(std::string* outError) {
  // We require the PBO mapping entry points.
  supported_ = (gl::MapBufferRange != nullptr) && (gl::UnmapBuffer != nullptr);

  // Fences are optional, but when available they allow us to avoid blocking on MapBufferRange.
  fencesSupported_ = (gl::FenceSync != nullptr) && (gl::ClientWaitSync != nullptr) && (gl::DeleteSync != nullptr);

  if (!supported_ && outError) {
    *outError = "Async readback unsupported (missing glMapBufferRange/glUnmapBuffer).";
  }
  return supported_;
}

void AsyncTextureReadback::shutdown() {
  clear();
  supported_ = false;
  fencesSupported_ = false;
}

void AsyncTextureReadback::clear() {
  if (jobs_.empty()) return;

  for (auto& j : jobs_) {
    if (j.fence) {
      deleteFence(j.fence);
      j.fence = nullptr;
    }
    if (j.pbo) {
      unsigned int pbo = j.pbo;
      gl::DeleteBuffers(1, &pbo);
      j.pbo = 0;
    }
  }
  jobs_.clear();
}

bool AsyncTextureReadback::enqueue(const TextureReadbackDesc& descIn, std::string tag, int delayFrames) {
  if (!supported_) return false;
  if (descIn.width <= 0 || descIn.height <= 0) return false;
  if (descIn.format == 0 || descIn.type == 0) return false;

  TextureReadbackDesc desc = descIn;

  switch (desc.source) {
    case ReadbackSource::Texture2D:
      if (desc.texture == 0) return false;
      break;
    case ReadbackSource::Framebuffer:
      // framebuffer may be 0 (default framebuffer/backbuffer). readBuffer can be 0 (auto).
      if (desc.readBuffer == 0) {
        desc.readBuffer = (desc.framebuffer == 0) ? (unsigned int)GL_BACK : (unsigned int)GL_COLOR_ATTACHMENT0;
      }
      break;
    default:
      return false;
  }

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

  GLint prevPackPbo = 0;
#ifdef GL_PIXEL_PACK_BUFFER_BINDING
  glGetIntegerv(GL_PIXEL_PACK_BUFFER_BINDING, &prevPackPbo);
#endif

  glPixelStorei(GL_PACK_ALIGNMENT, 1);
#ifdef GL_PACK_ROW_LENGTH
  glPixelStorei(GL_PACK_ROW_LENGTH, 0);
  glPixelStorei(GL_PACK_SKIP_ROWS, 0);
  glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
#endif

  gl::BindBuffer(GL_PIXEL_PACK_BUFFER, pbo);
  gl::BufferData(GL_PIXEL_PACK_BUFFER, (GLsizeiptr)sizeBytes, nullptr, GL_STREAM_READ);

  if (desc.source == ReadbackSource::Texture2D) {
    GLint prevTex = 0;
#ifdef GL_TEXTURE_BINDING_2D
    glGetIntegerv(GL_TEXTURE_BINDING_2D, &prevTex);
#endif

    gl::BindTexture(GL_TEXTURE_2D, desc.texture);
    glGetTexImage(GL_TEXTURE_2D, desc.level, (GLenum)desc.format, (GLenum)desc.type, (void*)0);
    gl::BindTexture(GL_TEXTURE_2D, (GLuint)prevTex);
  } else {
    // Framebuffer readback (default framebuffer/backbuffer or an FBO).
    GLint prevReadBuf = (GLint)GL_BACK;
#ifdef GL_READ_BUFFER
    glGetIntegerv(GL_READ_BUFFER, &prevReadBuf);
#endif

    GLint prevReadFb = 0;
#ifdef GL_READ_FRAMEBUFFER_BINDING
    glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING, &prevReadFb);
#elif defined(GL_FRAMEBUFFER_BINDING)
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevReadFb);
#endif

#ifdef GL_READ_FRAMEBUFFER
    gl::BindFramebuffer(GL_READ_FRAMEBUFFER, desc.framebuffer);
#else
    gl::BindFramebuffer(GL_FRAMEBUFFER, desc.framebuffer);
#endif
    glReadBuffer((GLenum)desc.readBuffer);
    glReadPixels(desc.x, desc.y, desc.width, desc.height,
                 (GLenum)desc.format, (GLenum)desc.type, (void*)0);

    // Restore framebuffer + read buffer state.
    glReadBuffer((GLenum)prevReadBuf);
#ifdef GL_READ_FRAMEBUFFER
    gl::BindFramebuffer(GL_READ_FRAMEBUFFER, (GLuint)prevReadFb);
#else
    gl::BindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevReadFb);
#endif
  }

  void* fence = nullptr;
  if (fencesSupported_) {
    GLsync f = gl::FenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
    fence = reinterpret_cast<void*>(f);
  }

  // Restore pack buffer binding and pack state.
  gl::BindBuffer(GL_PIXEL_PACK_BUFFER, (GLuint)prevPackPbo);

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
  j.fence = fence;

  jobs_.push_back(std::move(j));
  return true;
}

bool AsyncTextureReadback::poll(TextureReadbackResult& out) {
  out = {};
  if (!supported_) return false;
  if (jobs_.empty()) return false;

  GLint prevPackPbo = 0;
#ifdef GL_PIXEL_PACK_BUFFER_BINDING
  glGetIntegerv(GL_PIXEL_PACK_BUFFER_BINDING, &prevPackPbo);
#endif

  for (std::size_t i = 0; i < jobs_.size(); ++i) {
    Job& j = jobs_[i];

    if (j.delay > 0) {
      --j.delay;
      continue;
    }

    if (!j.pbo || j.sizeBytes == 0) {
      // Corrupt job; drop it.
      if (j.fence) {
        deleteFence(j.fence);
        j.fence = nullptr;
      }
      if (j.pbo) {
        unsigned int pbo = j.pbo;
        gl::DeleteBuffers(1, &pbo);
        j.pbo = 0;
      }
      jobs_.erase(jobs_.begin() + (std::ptrdiff_t)i);
      --i;
      continue;
    }

    // If we have a fence, only map once the GPU has completed the readback.
    if (j.fence && gl::ClientWaitSync) {
      const GLenum r = gl::ClientWaitSync(reinterpret_cast<GLsync>(j.fence),
                                          GL_SYNC_FLUSH_COMMANDS_BIT,
                                          0);
      if (r == GL_TIMEOUT_EXPIRED) {
        continue; // not ready yet
      }
      if (r == GL_WAIT_FAILED) {
        // Fence is broken; drop the job.
        deleteFence(j.fence);
        j.fence = nullptr;

        unsigned int pbo = j.pbo;
        gl::DeleteBuffers(1, &pbo);
        jobs_.erase(jobs_.begin() + (std::ptrdiff_t)i);
        --i;
        continue;
      }
      // GL_ALREADY_SIGNALED or GL_CONDITION_SATISFIED => ready.
    }

    gl::BindBuffer(GL_PIXEL_PACK_BUFFER, j.pbo);
    void* ptr = gl::MapBufferRange(GL_PIXEL_PACK_BUFFER, 0, (GLsizeiptr)j.sizeBytes, GL_MAP_READ_BIT);
    if (!ptr) {
      // Not ready (or mapping failed). Try again later.
      gl::BindBuffer(GL_PIXEL_PACK_BUFFER, (GLuint)prevPackPbo);
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
    gl::BindBuffer(GL_PIXEL_PACK_BUFFER, (GLuint)prevPackPbo);

    if (j.fence) {
      deleteFence(j.fence);
      j.fence = nullptr;
    }

    unsigned int pbo = j.pbo;
    gl::DeleteBuffers(1, &pbo);

    jobs_.erase(jobs_.begin() + (std::ptrdiff_t)i);
    return true;
  }

  gl::BindBuffer(GL_PIXEL_PACK_BUFFER, (GLuint)prevPackPbo);
  return false;
}

} // namespace stellar::render
