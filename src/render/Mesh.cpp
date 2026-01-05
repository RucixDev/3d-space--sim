#include "stellar/render/Mesh.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <cmath>

namespace stellar::render {

Mesh::~Mesh() { destroy(); }

Mesh::Mesh(Mesh&& o) noexcept {
  vao_ = o.vao_; vbo_ = o.vbo_; ebo_ = o.ebo_; indexCount_ = o.indexCount_;
  o.vao_ = o.vbo_ = o.ebo_ = 0;
  o.indexCount_ = 0;
}

Mesh& Mesh::operator=(Mesh&& o) noexcept {
  if (this == &o) return *this;
  destroy();
  vao_ = o.vao_; vbo_ = o.vbo_; ebo_ = o.ebo_; indexCount_ = o.indexCount_;
  o.vao_ = o.vbo_ = o.ebo_ = 0;
  o.indexCount_ = 0;
  return *this;
}

void Mesh::destroy() {
  if (ebo_) gl::DeleteBuffers(1, &ebo_);
  if (vbo_) gl::DeleteBuffers(1, &vbo_);
  if (vao_) gl::DeleteVertexArrays(1, &vao_);
  vao_ = vbo_ = ebo_ = 0;
  indexCount_ = 0;
}

void Mesh::upload(const std::vector<VertexPNUT>& vertices, const std::vector<std::uint32_t>& indices) {
  destroy();

  gl::GenVertexArrays(1, &vao_);
  gl::GenBuffers(1, &vbo_);
  gl::GenBuffers(1, &ebo_);

  gl::BindVertexArray(vao_);

  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);
  gl::BufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(vertices.size() * sizeof(VertexPNUT)),
                 vertices.data(),
                 GL_STATIC_DRAW);

  gl::BindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);
  gl::BufferData(GL_ELEMENT_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(indices.size() * sizeof(std::uint32_t)),
                 indices.data(),
                 GL_STATIC_DRAW);

  // layout (location=0) position
  gl::EnableVertexAttribArray(0);
  gl::VertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPNUT), (void*)offsetof(VertexPNUT, px));

  // normal (location=1)
  gl::EnableVertexAttribArray(1);
  gl::VertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPNUT), (void*)offsetof(VertexPNUT, nx));

  // uv (location=2)
  gl::EnableVertexAttribArray(2);
  gl::VertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(VertexPNUT), (void*)offsetof(VertexPNUT, u));

  gl::BindVertexArray(0);

  indexCount_ = static_cast<std::uint32_t>(indices.size());
}

void Mesh::bind() const {
  gl::BindVertexArray(vao_);
}

void Mesh::draw() const {
  gl::BindVertexArray(vao_);
  glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(indexCount_), GL_UNSIGNED_INT, nullptr);
  gl::BindVertexArray(0);
}

void Mesh::drawInstanced(std::uint32_t instanceCount) const {
  gl::BindVertexArray(vao_);
  gl::DrawElementsInstanced(GL_TRIANGLES, static_cast<GLsizei>(indexCount_), GL_UNSIGNED_INT, nullptr, static_cast<GLsizei>(instanceCount));
  gl::BindVertexArray(0);
}

Mesh Mesh::makeCube() {
  // Unit cube centered at origin, each face has its own vertices for correct normals.
  std::vector<VertexPNUT> v;
  std::vector<std::uint32_t> i;

  const float p = 0.5f;

  auto addFace = [&](float nx, float ny, float nz,
                     float x0,float y0,float z0,
                     float x1,float y1,float z1,
                     float x2,float y2,float z2,
                     float x3,float y3,float z3) {
    const std::uint32_t base = static_cast<std::uint32_t>(v.size());
    v.push_back({x0,y0,z0,nx,ny,nz,0,0});
    v.push_back({x1,y1,z1,nx,ny,nz,1,0});
    v.push_back({x2,y2,z2,nx,ny,nz,1,1});
    v.push_back({x3,y3,z3,nx,ny,nz,0,1});
    i.insert(i.end(), {base+0, base+1, base+2, base+0, base+2, base+3});
  };

  // +Z
  addFace(0,0,1, -p,-p, p,  p,-p, p,  p, p, p, -p, p, p);
  // -Z
  addFace(0,0,-1,  p,-p,-p, -p,-p,-p, -p, p,-p,  p, p,-p);
  // +X
  addFace(1,0,0,  p,-p, p,  p,-p,-p,  p, p,-p,  p, p, p);
  // -X
  addFace(-1,0,0, -p,-p,-p, -p,-p, p, -p, p, p, -p, p,-p);
  // +Y
  addFace(0,1,0, -p, p, p,  p, p, p,  p, p,-p, -p, p,-p);
  // -Y
  addFace(0,-1,0, -p,-p,-p,  p,-p,-p,  p,-p, p, -p,-p, p);

  Mesh m;
  m.upload(v, i);
  return m;
}

Mesh Mesh::makeUvSphere(int slices, int stacks) {
  slices = std::max(3, slices);
  stacks = std::max(2, stacks);

  std::vector<VertexPNUT> v;
  std::vector<std::uint32_t> idx;

  for (int y = 0; y <= stacks; ++y) {
    const float vty = static_cast<float>(y) / static_cast<float>(stacks);
    const float phi = vty * 3.14159265358979323846f; // 0..pi

    for (int x = 0; x <= slices; ++x) {
      const float vtx = static_cast<float>(x) / static_cast<float>(slices);
      const float theta = vtx * 2.0f * 3.14159265358979323846f; // 0..2pi

      const float sx = std::sin(phi) * std::cos(theta);
      const float sy = std::cos(phi);
      const float sz = std::sin(phi) * std::sin(theta);

      VertexPNUT vert{};
      vert.px = sx * 0.5f;
      vert.py = sy * 0.5f;
      vert.pz = sz * 0.5f;
      vert.nx = sx;
      vert.ny = sy;
      vert.nz = sz;
      vert.u = vtx;
      vert.v = 1.0f - vty;
      v.push_back(vert);
    }
  }

  const int stride = slices + 1;
  for (int y = 0; y < stacks; ++y) {
    for (int x = 0; x < slices; ++x) {
      const std::uint32_t i0 = static_cast<std::uint32_t>(y * stride + x);
      const std::uint32_t i1 = static_cast<std::uint32_t>((y+1) * stride + x);
      const std::uint32_t i2 = static_cast<std::uint32_t>((y+1) * stride + (x+1));
      const std::uint32_t i3 = static_cast<std::uint32_t>(y * stride + (x+1));

      idx.push_back(i0); idx.push_back(i1); idx.push_back(i2);
      idx.push_back(i0); idx.push_back(i2); idx.push_back(i3);
    }
  }

  Mesh m;
  m.upload(v, idx);
  return m;
}

Mesh Mesh::makeRing(int segments, float innerRadius, float outerRadius, bool doubleSided) {
  segments = std::max(3, segments);
  innerRadius = std::max(0.0f, innerRadius);
  outerRadius = std::max(innerRadius + 1e-4f, outerRadius);

  // Two vertices per segment (inner + outer). Duplicate the first segment at the end
  // so u can reach 1.0 for a seamless seam.
  const int ringVerts = (segments + 1) * 2;

  std::vector<VertexPNUT> v;
  v.reserve(doubleSided ? ringVerts * 2 : ringVerts);
  std::vector<std::uint32_t> idx;
  idx.reserve(static_cast<std::size_t>(segments) * 6u * (doubleSided ? 2u : 1u));

  const float twoPi = 2.0f * 3.14159265358979323846f;

  // Top face (+Y)
  for (int s = 0; s <= segments; ++s) {
    const float u = static_cast<float>(s) / static_cast<float>(segments);
    const float th = u * twoPi;
    const float c = std::cos(th);
    const float si = std::sin(th);

    // inner
    v.push_back(VertexPNUT{innerRadius * c, 0.0f, innerRadius * si,
                           0.0f, 1.0f, 0.0f,
                           u, 0.0f});
    // outer
    v.push_back(VertexPNUT{outerRadius * c, 0.0f, outerRadius * si,
                           0.0f, 1.0f, 0.0f,
                           u, 1.0f});
  }

  for (int s = 0; s < segments; ++s) {
    const std::uint32_t i0 = static_cast<std::uint32_t>(s * 2 + 0);
    const std::uint32_t i1 = static_cast<std::uint32_t>(s * 2 + 1);
    const std::uint32_t i2 = static_cast<std::uint32_t>((s + 1) * 2 + 0);
    const std::uint32_t i3 = static_cast<std::uint32_t>((s + 1) * 2 + 1);

    // Two triangles for the quad strip.
    // Winding is CCW when viewed from +Y.
    idx.insert(idx.end(), {i0, i1, i3, i0, i3, i2});
  }

  if (doubleSided) {
    // Bottom face (-Y): duplicate vertices with flipped normals and reversed winding.
    const std::uint32_t base = static_cast<std::uint32_t>(v.size());
    for (int k = 0; k < ringVerts; ++k) {
      VertexPNUT vv = v[static_cast<std::size_t>(k)];
      vv.nx = 0.0f;
      vv.ny = -1.0f;
      vv.nz = 0.0f;
      v.push_back(vv);
    }

    const std::size_t topTriCount = static_cast<std::size_t>(segments) * 2u;
    const std::size_t topIndexCount = topTriCount * 3u;
    // Mirror the top-face indices in reverse order to flip winding.
    for (std::size_t t = 0; t < topIndexCount; t += 3) {
      const std::uint32_t a = idx[t + 0];
      const std::uint32_t b = idx[t + 1];
      const std::uint32_t c = idx[t + 2];
      idx.push_back(base + a);
      idx.push_back(base + c);
      idx.push_back(base + b);
    }
  }

  Mesh m;
  m.upload(v, idx);
  return m;
}

} // namespace stellar::render
