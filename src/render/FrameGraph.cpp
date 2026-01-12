#include "stellar/render/FrameGraph.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <sstream>
#include <unordered_map>

namespace stellar::render {

FrameGraph::FrameGraph() = default;

FrameGraph::~FrameGraph() {
  destroyGL();
}

void FrameGraph::reset() {
  textures_.clear();
  passes_.clear();
  callbacks_.clear();
  schedule_.clear();
  physPlan_.clear();
  physRequired_ = 0;
}

void FrameGraph::setViewport(int w, int h) {
  viewW_ = std::max(0, w);
  viewH_ = std::max(0, h);
}

FrameGraph::TextureHandle FrameGraph::createTexture(std::string name, const TextureDesc& desc) {
  TextureInfo info;
  info.name = std::move(name);
  info.external = false;
  info.desc = desc;
  info.resolvedW = 0;
  info.resolvedH = 0;
  info.firstUse = -1;
  info.lastUse = -1;
  info.physicalId = -1;
  info.glHandle = 0;

  const int id = (int)textures_.size();
  textures_.push_back(std::move(info));
  return TextureHandle{id};
}

FrameGraph::TextureHandle FrameGraph::importTexture(std::string name, unsigned int glHandle, int w, int h) {
  TextureInfo info;
  info.name = std::move(name);
  info.external = true;
  info.desc = {};
  info.resolvedW = std::max(0, w);
  info.resolvedH = std::max(0, h);
  info.firstUse = -1;
  info.lastUse = -1;
  info.physicalId = -1;
  info.glHandle = glHandle;

  const int id = (int)textures_.size();
  textures_.push_back(std::move(info));
  return TextureHandle{id};
}

int FrameGraph::addPass(std::string name,
                        std::vector<TextureHandle> reads,
                        TextureHandle write,
                        PassCallback cb) {
  PassInfo p;
  p.name = std::move(name);
  p.reads = std::move(reads);
  p.write = write;
  const int id = (int)passes_.size();
  passes_.push_back(std::move(p));
  callbacks_.push_back(std::move(cb));
  return id;
}

FrameGraph::ResolvedDesc FrameGraph::resolveDesc(const TextureDesc& d) const {
  ResolvedDesc r;
  const float sc = (d.scale <= 0.0f) ? 1.0f : d.scale;

  const int wFromScale = (viewW_ > 0) ? (int)std::lround((float)viewW_ * sc) : 0;
  const int hFromScale = (viewH_ > 0) ? (int)std::lround((float)viewH_ * sc) : 0;

  r.w = std::max(1, (d.width > 0) ? d.width : wFromScale);
  r.h = std::max(1, (d.height > 0) ? d.height : hFromScale);
  r.internalFormat = d.internalFormat;
  r.format = d.format;
  r.type = d.type;
  r.linearFilter = d.linearFilter;
  r.clampToEdge = d.clampToEdge;
  return r;
}

bool FrameGraph::sameDesc(const ResolvedDesc& a, const ResolvedDesc& b) {
  return a.w == b.w && a.h == b.h &&
         a.internalFormat == b.internalFormat &&
         a.format == b.format &&
         a.type == b.type &&
         a.linearFilter == b.linearFilter &&
         a.clampToEdge == b.clampToEdge;
}

bool FrameGraph::compile(std::string* outError) {
  schedule_.clear();
  physPlan_.clear();
  physRequired_ = 0;

  // Keep profiling arrays sized to the pass count without resetting values unless needed.
  if (lastPassTimesMs_.size() != passes_.size()) {
    lastPassTimesMs_.assign(passes_.size(), 0.0);
  }

  // Resolve sizes for internal textures.
  for (auto& t : textures_) {
    if (t.external) continue;
    const ResolvedDesc rd = resolveDesc(t.desc);
    t.resolvedW = rd.w;
    t.resolvedH = rd.h;
  }

  const int passCount = (int)passes_.size();
  if (passCount == 0) {
    // Nothing to do.
    return true;
  }

  // Build dependencies using a simple last-writer map per texture.
  std::vector<std::vector<int>> adj(passCount);
  std::vector<int> indeg(passCount, 0);
  std::unordered_map<int, int> lastWriter; // textureId -> passId

  auto addEdge = [&](int a, int b) {
    if (a < 0 || b < 0 || a == b) return;
    adj[a].push_back(b);
  };

  for (int p = 0; p < passCount; ++p) {
    // Reads depend on last writer.
    for (const auto& r : passes_[p].reads) {
      if (!r.valid()) continue;
      auto it = lastWriter.find(r.id);
      if (it != lastWriter.end()) {
        addEdge(it->second, p);
      }
    }

    // Writes also depend on last writer (write-after-write ordering).
    if (passes_[p].write.valid()) {
      auto it = lastWriter.find(passes_[p].write.id);
      if (it != lastWriter.end()) {
        addEdge(it->second, p);
      }
      lastWriter[passes_[p].write.id] = p;
    }
  }

  // Compute indegrees (dedupe edges per node).
  for (int a = 0; a < passCount; ++a) {
    auto& v = adj[a];
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    for (int b : v) {
      if (b < 0 || b >= passCount) continue;
      indeg[b]++;
    }
  }

  // Kahn topo sort. Use pass insertion order as a tie-breaker.
  std::queue<int> q;
  for (int i = 0; i < passCount; ++i) {
    if (indeg[i] == 0) q.push(i);
  }

  while (!q.empty()) {
    const int n = q.front();
    q.pop();
    schedule_.push_back(n);
    for (int m : adj[n]) {
      if (--indeg[m] == 0) q.push(m);
    }
  }

  if ((int)schedule_.size() != passCount) {
    if (outError) {
      *outError = "FrameGraph compile failed: dependency cycle detected.";
    }
    return false;
  }

  // Lifetime analysis (inclusive intervals) in scheduled order.
  for (auto& t : textures_) {
    t.firstUse = -1;
    t.lastUse = -1;
    t.physicalId = -1;
  }

  auto markUse = [&](TextureHandle h, int passIdx) {
    if (!h.valid()) return;
    auto& t = textures_[(std::size_t)h.id];
    if (t.firstUse < 0) t.firstUse = passIdx;
    t.lastUse = std::max(t.lastUse, passIdx);
  };

  for (int si = 0; si < (int)schedule_.size(); ++si) {
    const int p = schedule_[si];
    for (const auto& r : passes_[p].reads) {
      markUse(r, si);
    }
    if (passes_[p].write.valid()) {
      markUse(passes_[p].write, si);
    }
  }

  // Build a list of internal textures that actually participate.
  struct TexAlloc {
    int texId;
    ResolvedDesc desc;
    int firstUse;
    int lastUse;
  };

  std::vector<TexAlloc> allocs;
  allocs.reserve(textures_.size());

  for (int i = 0; i < (int)textures_.size(); ++i) {
    const auto& t = textures_[(std::size_t)i];
    if (t.external) continue;
    if (t.firstUse < 0 || t.lastUse < 0) continue;

    TexAlloc a;
    a.texId = i;
    a.desc = resolveDesc(t.desc);
    a.firstUse = t.firstUse;
    a.lastUse = t.lastUse;
    allocs.push_back(a);
  }

  std::sort(allocs.begin(), allocs.end(), [](const TexAlloc& a, const TexAlloc& b) {
    if (a.firstUse != b.firstUse) return a.firstUse < b.firstUse;
    return a.texId < b.texId;
  });

  if (!aliasingEnabled_) {
    // No aliasing: each internal texture gets a unique physical allocation.
    physRequired_ = (int)allocs.size();
    physPlan_.resize((std::size_t)physRequired_);
    for (int i = 0; i < (int)allocs.size(); ++i) {
      const TexAlloc& a = allocs[(std::size_t)i];
      textures_[(std::size_t)a.texId].physicalId = i;
      physPlan_[(std::size_t)i] = a.desc;
    }
    return true;
  }

  // Aliasing: linear-scan reuse by lifetime + matching desc.
  struct PhysSlot {
    ResolvedDesc desc;
    int lastUse{-1};
  };

  std::vector<PhysSlot> slots;

  for (const TexAlloc& a : allocs) {
    int reuse = -1;
    for (int i = 0; i < (int)slots.size(); ++i) {
      if (slots[(std::size_t)i].lastUse < a.firstUse && sameDesc(slots[(std::size_t)i].desc, a.desc)) {
        reuse = i;
        break;
      }
    }

    if (reuse < 0) {
      PhysSlot s;
      s.desc = a.desc;
      s.lastUse = a.lastUse;
      reuse = (int)slots.size();
      slots.push_back(s);
    } else {
      slots[(std::size_t)reuse].lastUse = a.lastUse;
    }

    textures_[(std::size_t)a.texId].physicalId = reuse;
  }

  physRequired_ = (int)slots.size();
  physPlan_.resize((std::size_t)physRequired_);
  for (int i = 0; i < physRequired_; ++i) {
    physPlan_[(std::size_t)i] = slots[(std::size_t)i].desc;
  }

  return true;
}

unsigned int FrameGraph::glTexture(TextureHandle h) const {
  if (!h.valid()) return 0;
  const auto& t = textures_[(std::size_t)h.id];
  if (t.external) return t.glHandle;
  const int pid = t.physicalId;
  if (pid < 0 || pid >= (int)phys_.size()) return 0;
  return phys_[(std::size_t)pid].tex;
}

unsigned int FrameGraph::PassContext::texture(TextureHandle h) const {
  return fg.glTexture(h);
}

int FrameGraph::PassContext::textureWidth(TextureHandle h) const {
  if (!h.valid()) return 0;
  return fg.textures()[(std::size_t)h.id].resolvedW;
}

int FrameGraph::PassContext::textureHeight(TextureHandle h) const {
  if (!h.valid()) return 0;
  return fg.textures()[(std::size_t)h.id].resolvedH;
}

void FrameGraph::destroyGL() {
  // Delete timer queries used for GPU profiling (if supported).
  if (gl::DeleteQueries) {
    for (auto& frameQs : passQueries_) {
      if (!frameQs.empty()) {
        gl::DeleteQueries((GLsizei)frameQs.size(), frameQs.data());
      }
      frameQs.clear();
    }
  } else {
    for (auto& frameQs : passQueries_) frameQs.clear();
  }
  queryFrame_ = 0;
  profilingSupported_ = false;

  if (!glAllocated_) return;
  for (auto& p : phys_) {
    if (p.fbo) {
      gl::DeleteFramebuffers(1, &p.fbo);
      p.fbo = 0;
    }
    if (p.tex) {
      gl::DeleteTextures(1, &p.tex);
      p.tex = 0;
    }
  }
  phys_.clear();
  glAllocated_ = false;
}

unsigned int FrameGraph::ensureGlTexture(int physId, const ResolvedDesc& d) {
  if (physId < 0) return 0;
  if (physId >= (int)phys_.size()) phys_.resize((std::size_t)physId + 1);

  auto& p = phys_[(std::size_t)physId];

  const bool needAlloc = (p.tex == 0) || !sameDesc(p.desc, d);
  if (needAlloc) {
    if (!p.tex) gl::GenTextures(1, &p.tex);
    gl::BindTexture(GL_TEXTURE_2D, p.tex);
    gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, d.linearFilter ? GL_LINEAR : GL_NEAREST);
    gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, d.linearFilter ? GL_LINEAR : GL_NEAREST);
    gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, d.clampToEdge ? GL_CLAMP_TO_EDGE : GL_REPEAT);
    gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, d.clampToEdge ? GL_CLAMP_TO_EDGE : GL_REPEAT);
    gl::TexImage2D(GL_TEXTURE_2D, 0, d.internalFormat, d.w, d.h, 0, d.format, d.type, nullptr);
    gl::BindTexture(GL_TEXTURE_2D, 0);
    p.desc = d;
  }

  return p.tex;
}

unsigned int FrameGraph::ensureFbo(int physId) {
  if (physId < 0 || physId >= (int)phys_.size()) return 0;
  auto& p = phys_[(std::size_t)physId];
  if (!p.fbo) gl::GenFramebuffers(1, &p.fbo);
  gl::BindFramebuffer(GL_FRAMEBUFFER, p.fbo);
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, p.tex, 0);
  glDrawBuffer(GL_COLOR_ATTACHMENT0);
  gl::BindFramebuffer(GL_FRAMEBUFFER, 0);
  return p.fbo;
}

void FrameGraph::ensurePhysicalAllocated() {
  if (physRequired_ <= 0) return;
  if ((int)phys_.size() < physRequired_) phys_.resize((std::size_t)physRequired_);

  // Create/resize textures and attach FBOs.
  for (int i = 0; i < physRequired_; ++i) {
    const ResolvedDesc& d = physPlan_[(std::size_t)i];
    ensureGlTexture(i, d);
    ensureFbo(i);
  }

  glAllocated_ = true;
}


void FrameGraph::execute() {
  if (schedule_.empty()) return;
  ensurePhysicalAllocated();

  // Update internal texture GL handles for debug.
  for (auto& t : textures_) {
    if (t.external) continue;
    const int pid = t.physicalId;
    t.glHandle = (pid >= 0 && pid < (int)phys_.size()) ? phys_[(std::size_t)pid].tex : 0;
  }

  // Determine whether timer queries are supported by the current GL context.
  profilingSupported_ = (gl::GenQueries && gl::DeleteQueries && gl::BeginQuery && gl::EndQuery && gl::GetQueryObjectui64v);

  const int passCount = (int)passes_.size();
  const bool doProfile = profilingEnabled_ && profilingSupported_ && passCount > 0;

  int readFrame = 0;
  int writeFrame = 0;

  if (doProfile) {
    // Ensure query objects exist for each pass in a small ring-buffer to avoid stalls.
    for (int fi = 0; fi < kQueryFrames; ++fi) {
      auto& qs = passQueries_[(std::size_t)fi];
      if ((int)qs.size() != passCount) {
        if (!qs.empty()) {
          gl::DeleteQueries((GLsizei)qs.size(), qs.data());
        }
        qs.assign((std::size_t)passCount, 0u);
        gl::GenQueries((GLsizei)passCount, qs.data());
      }
    }

    readFrame = (queryFrame_ + kQueryFrames - 1) % kQueryFrames;
    writeFrame = queryFrame_;

    // Read back results from the previous frame (if available).
    // We check availability to avoid stalling the CPU.
    double sumMs = 0.0;
    bool anyReady = false;

    auto& rqs = passQueries_[(std::size_t)readFrame];
    const int n = std::min(passCount, (int)rqs.size());
    for (int i = 0; i < n; ++i) {
      const unsigned int q = rqs[(std::size_t)i];
      if (!q) {
        sumMs += (i < (int)lastPassTimesMs_.size()) ? lastPassTimesMs_[(std::size_t)i] : 0.0;
        continue;
      }

      GLuint64 avail = 0;
      gl::GetQueryObjectui64v(q, GL_QUERY_RESULT_AVAILABLE, &avail);

      if (avail) {
        GLuint64 ns = 0;
        gl::GetQueryObjectui64v(q, GL_QUERY_RESULT, &ns);
        const double ms = (double)ns / 1000000.0;
        if (i < (int)lastPassTimesMs_.size()) lastPassTimesMs_[(std::size_t)i] = ms;
        sumMs += ms;
        anyReady = true;
      } else {
        sumMs += (i < (int)lastPassTimesMs_.size()) ? lastPassTimesMs_[(std::size_t)i] : 0.0;
      }
    }

    if (anyReady) {
      lastFrameTimeMs_ = sumMs;
    }
  }

  // Execute passes.
  for (int si = 0; si < (int)schedule_.size(); ++si) {
    const int pIdx = schedule_[(std::size_t)si];
    if (pIdx < 0 || pIdx >= (int)passes_.size()) continue;
    const PassInfo& p = passes_[(std::size_t)pIdx];

    int outW = viewW_;
    int outH = viewH_;

    if (p.write.isBackbuffer() || !p.write.valid()) {
      gl::BindFramebuffer(GL_FRAMEBUFFER, 0);
    } else {
      const auto& tex = textures_[(std::size_t)p.write.id];
      outW = tex.resolvedW;
      outH = tex.resolvedH;
      const int pid = tex.physicalId;
      if (pid >= 0 && pid < (int)phys_.size()) {
        gl::BindFramebuffer(GL_FRAMEBUFFER, phys_[(std::size_t)pid].fbo);
      } else {
        gl::BindFramebuffer(GL_FRAMEBUFFER, 0);
      }
    }

    glViewport(0, 0, std::max(1, outW), std::max(1, outH));

    bool beganQuery = false;
    if (doProfile) {
      const auto& wqs = passQueries_[(std::size_t)writeFrame];
      if (pIdx >= 0 && pIdx < (int)wqs.size() && wqs[(std::size_t)pIdx] != 0u) {
        gl::BeginQuery(GL_TIME_ELAPSED, wqs[(std::size_t)pIdx]);
        beganQuery = true;
      }
    }

    if (pIdx < (int)callbacks_.size()) {
      PassContext ctx{*this, viewW_, viewH_, outW, outH};
      callbacks_[(std::size_t)pIdx](ctx);
    }

    if (beganQuery) {
      gl::EndQuery(GL_TIME_ELAPSED);
    }
  }

  if (doProfile) {
    queryFrame_ = (queryFrame_ + 1) % kQueryFrames;
  }

  // Leave whatever the last pass bound (often the backbuffer).
}

} // namespace stellar::render
