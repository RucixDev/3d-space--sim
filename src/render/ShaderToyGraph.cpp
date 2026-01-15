#include "stellar/render/ShaderToyGraph.h"

#include "stellar/render/Gl.h"

#include <algorithm>

namespace stellar::render {

namespace {

static int clampScale(int s) {
  if (s <= 1) return 1;
  if (s <= 2) return 2;
  if (s <= 4) return 4;
  return 4;
}

} // namespace

bool ShaderToyGraph::init(std::string* outError) {
  if (inited_) return true;

  for (int i = 0; i < kPassCount; ++i) {
    passes_[i].valid = false;
    passes_[i].error.clear();
    passes_[i].ping = 0;
    passes_[i].rtInited[0] = false;
    passes_[i].rtInited[1] = false;
    // Default: buffers enabled, image enabled.
    passes_[i].settings = ShaderToyGraphPassSettings{};
  }

  // Image defaults to enabled even if the caller forgets.
  passes_[idx(ShaderToyPass::Image)].settings.enabled = true;

  for (int i = 0; i < kPassCount; ++i) {
    if (!passes_[i].toy.init(outError)) {
      return false;
    }
    passes_[i].valid = true;
  }

  inited_ = true;
  return true;
}

bool ShaderToyGraph::ensureSize(int outWidth, int outHeight, std::string* outError) {
  if (!init(outError)) return false;

  outW_ = std::max(1, outWidth);
  outH_ = std::max(1, outHeight);

  // Allocate/resize feedback buffers.
  for (int i = 0; i < kBufferCount; ++i) {
    PassState& p = passes_[i];
    const int scale = clampScale(p.settings.resolutionScale);
    p.settings.resolutionScale = scale;

    const int w = std::max(1, outW_ / scale);
    const int h = std::max(1, outH_ / scale);

    for (int b = 0; b < 2; ++b) {
      if (!p.rtInited[b]) {
        if (!p.rt[b].init(w, h, outError)) return false;
        p.rtInited[b] = true;
      } else {
        p.rt[b].ensureSize(w, h);
      }
    }
  }

  return true;
}

bool ShaderToyGraph::buildPass(ShaderToyPass pass, std::string_view userCode, std::string* outError) {
  if (!init(outError)) return false;

  PassState& p = passes_[idx(pass)];
  p.error.clear();
  p.valid = p.toy.buildWithHeader(userCode, userHeader_, &p.error);
  if (!p.valid && outError) *outError = p.error;
  return p.valid;
}

void ShaderToyGraph::setUserHeader(std::string_view extraHeader) {
  userHeader_.assign(extraHeader.data(), extraHeader.size());
}

void ShaderToyGraph::setUserHeader(std::string_view extraHeader) {
  userHeader_.assign(extraHeader.data(), extraHeader.size());
}

void ShaderToyGraph::setPassSettings(ShaderToyPass pass, const ShaderToyGraphPassSettings& s) {
  PassState& p = passes_[idx(pass)];
  p.settings = s;
  p.settings.resolutionScale = clampScale(p.settings.resolutionScale);

  // The Image pass is always logically enabled.
  if (pass == ShaderToyPass::Image) {
    p.settings.enabled = true;
  }
}

const ShaderToyGraphPassSettings& ShaderToyGraph::passSettings(ShaderToyPass pass) const {
  return passes_[idx(pass)].settings;
}

bool ShaderToyGraph::passValid(ShaderToyPass pass) const {
  return passes_[idx(pass)].valid;
}

const std::string& ShaderToyGraph::passError(ShaderToyPass pass) const {
  return passes_[idx(pass)].error;
}

void ShaderToyGraph::clearRtIfInited_(PassState& p) {
  for (int b = 0; b < 2; ++b) {
    if (!p.rtInited[b]) continue;

    p.rt[b].begin();
    glViewport(0, 0, p.rt[b].width(), p.rt[b].height());
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    RenderTarget2D::end();
  }
}

void ShaderToyGraph::resetBuffers() {
  for (int i = 0; i < kBufferCount; ++i) {
    PassState& p = passes_[i];
    p.ping = 0;
    clearRtIfInited_(p);
  }
}

const Texture2D* ShaderToyGraph::resolveChannelTexture_(ShaderToyPass /*currentPass*/, ShaderToyChannelSource src) const {
  if (src == ShaderToyChannelSource::None) return nullptr;

  const int srcIdx = static_cast<int>(src) - 1;
  if (srcIdx < 0 || srcIdx >= kBufferCount) return nullptr;

  const PassState& p = passes_[srcIdx];
  const int ping = (p.ping < 0 || p.ping > 1) ? 0 : p.ping;
  if (!p.rtInited[ping]) return nullptr;

  return &p.rt[ping].color();
}

void ShaderToyGraph::render(const ShaderToyInputs& baseInputs, RenderTarget2D& outImageTarget) {
  if (!inited_) return;

  // The graph's target size is driven by the output image target.
  // (Call ensureSize before render; this keeps the function void and cheap.)
  outW_ = std::max(1, outImageTarget.width());
  outH_ = std::max(1, outImageTarget.height());

  // Normalize mouse from baseInputs so per-pass resolution scaling stays coherent.
  const float baseW = (float)std::max(1, baseInputs.width);
  const float baseH = (float)std::max(1, baseInputs.height);

  const float mxN = baseInputs.mouse[0] / baseW;
  const float myN = baseInputs.mouse[1] / baseH;
  const float mdxN = baseInputs.mouse[2] / baseW;
  const float mdyN = baseInputs.mouse[3] / baseH;

  // ---- Buffers A..D ----
  for (int i = 0; i < kBufferCount; ++i) {
    PassState& p = passes_[i];
    if (!p.settings.enabled) continue;
    if (!p.valid) continue;

    const int writeIdx = 1 - p.ping;
    if (!p.rtInited[writeIdx]) continue;

    RenderTarget2D& dst = p.rt[writeIdx];
    dst.begin();
    glViewport(0, 0, dst.width(), dst.height());
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    ShaderToyInputs in = baseInputs;
    in.width = dst.width();
    in.height = dst.height();
    in.passIndex = i;

    // Re-scale mouse to this pass' resolution.
    in.mouse[0] = mxN * (float)in.width;
    in.mouse[1] = myN * (float)in.height;
    in.mouse[2] = mdxN * (float)in.width;
    in.mouse[3] = mdyN * (float)in.height;

    for (int c = 0; c < 4; ++c) {
      in.channels[c] = resolveChannelTexture_(static_cast<ShaderToyPass>(i), p.settings.channels[c]);
    }

    p.toy.draw(in);
    RenderTarget2D::end();

    p.ping = writeIdx;
  }

  // ---- Image pass ----
  {
    PassState& img = passes_[idx(ShaderToyPass::Image)];

    outImageTarget.begin();
    glViewport(0, 0, outImageTarget.width(), outImageTarget.height());
    glClearColor(0.02f, 0.02f, 0.03f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (img.valid) {
      ShaderToyInputs in = baseInputs;
      in.width = outImageTarget.width();
      in.height = outImageTarget.height();
      in.passIndex = idx(ShaderToyPass::Image);

      // (mouse already in output resolution, but keep it coherent if baseInputs differs.)
      in.mouse[0] = mxN * (float)in.width;
      in.mouse[1] = myN * (float)in.height;
      in.mouse[2] = mdxN * (float)in.width;
      in.mouse[3] = mdyN * (float)in.height;

      for (int c = 0; c < 4; ++c) {
        in.channels[c] = resolveChannelTexture_(ShaderToyPass::Image, img.settings.channels[c]);
      }

      img.toy.draw(in);
    }

    RenderTarget2D::end();
  }
}

const RenderTarget2D* ShaderToyGraph::bufferTarget(ShaderToyPass bufferPass) const {
  if (!isBuffer(bufferPass)) return nullptr;
  const PassState& p = passes_[idx(bufferPass)];
  const int ping = (p.ping < 0 || p.ping > 1) ? 0 : p.ping;
  if (!p.rtInited[ping]) return nullptr;
  return &p.rt[ping];
}

const Texture2D* ShaderToyGraph::bufferTexture(ShaderToyPass bufferPass) const {
  const RenderTarget2D* t = bufferTarget(bufferPass);
  if (!t) return nullptr;
  return &t->color();
}

const ShaderToy& ShaderToyGraph::shader(ShaderToyPass pass) const {
  return passes_[idx(pass)].toy;
}

} // namespace stellar::render
