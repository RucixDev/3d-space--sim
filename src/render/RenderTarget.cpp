#include "stellar/render/RenderTarget.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <string>

namespace stellar::render {

namespace {

bool sameDesc(const RenderTarget2DDesc& a, const RenderTarget2DDesc& b) {
  return a.hasDepth == b.hasDepth &&
         a.generateMips == b.generateMips &&
         a.nearestFilter == b.nearestFilter &&
         a.clampToEdge == b.clampToEdge;
}

} // namespace

RenderTarget2D::~RenderTarget2D() { destroy(); }

void RenderTarget2D::destroy() {
  if (depthRbo_) {
    gl::DeleteRenderbuffers(1, &depthRbo_);
    depthRbo_ = 0;
  }
  if (fbo_) {
    gl::DeleteFramebuffers(1, &fbo_);
    fbo_ = 0;
  }
  color_ = Texture2D{};
  w_ = h_ = 0;
  inited_ = false;
  desc_ = RenderTarget2DDesc{};
}

bool RenderTarget2D::init(int w, int h, std::string* outError) {
  return init(w, h, RenderTarget2DDesc{}, outError);
}

bool RenderTarget2D::init(int w, int h, const RenderTarget2DDesc& desc, std::string* outError) {
  return createOrResize(w, h, desc, outError);
}

bool RenderTarget2D::createOrResize(int w, int h, const RenderTarget2DDesc& desc, std::string* outError) {
  w = std::clamp(w, 1, 4096);
  h = std::clamp(h, 1, 4096);

  const RenderTarget2DDesc requested = desc;

  // Recreate only if size or descriptor changes.
  if (inited_ && w == w_ && h == h_ && sameDesc(desc_, requested)) return true;

  destroy();

  // Store descriptor (so ensureSize() can preserve it).
  desc_ = requested;

  // Allocate color texture.
  color_.allocateRGBA(w, h,
                      /*generateMips=*/desc_.generateMips,
                      /*nearestFilter=*/desc_.nearestFilter,
                      /*clampToEdge=*/desc_.clampToEdge);

  gl::GenFramebuffers(1, &fbo_);
  gl::BindFramebuffer(GL_FRAMEBUFFER, fbo_);
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, color_.handle(), 0);

  if (desc_.hasDepth) {
    // Depth buffer (no stencil needed, but DEPTH24_STENCIL8 is widely supported and convenient).
    gl::GenRenderbuffers(1, &depthRbo_);
    gl::BindRenderbuffer(GL_RENDERBUFFER, depthRbo_);
    gl::RenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, w, h);
    gl::FramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, depthRbo_);
  }

  const GLenum status = gl::CheckFramebufferStatus(GL_FRAMEBUFFER);
  gl::BindFramebuffer(GL_FRAMEBUFFER, 0);

  if (status != GL_FRAMEBUFFER_COMPLETE) {
    if (outError) {
      *outError = "RenderTarget2D: framebuffer incomplete (status=" + std::to_string((unsigned int)status) + ")";
    }
    destroy();
    return false;
  }

  w_ = w;
  h_ = h;
  inited_ = true;
  return true;
}

void RenderTarget2D::ensureSize(int w, int h) {
  (void)createOrResize(w, h, desc_, nullptr);
}

void RenderTarget2D::begin() const {
  if (!inited_) return;
  gl::BindFramebuffer(GL_FRAMEBUFFER, fbo_);
}

void RenderTarget2D::end() { gl::BindFramebuffer(GL_FRAMEBUFFER, 0); }

void RenderTarget2D::generateMips() const {
  if (!inited_ || !desc_.generateMips || color_.handle() == 0) return;

  // Preserve the caller's texture binding (ImGui and other render code can be stateful).
  GLint prevTex = 0;
  glGetIntegerv(GL_TEXTURE_BINDING_2D, &prevTex);

  gl::BindTexture(GL_TEXTURE_2D, color_.handle());
  gl::GenerateMipmap(GL_TEXTURE_2D);
  gl::BindTexture(GL_TEXTURE_2D, (unsigned int)prevTex);
}

} // namespace stellar::render
