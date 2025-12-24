#pragma once

#include "stellar/render/Mesh.h"
#include "stellar/render/Shader.h"
#include "stellar/render/Texture.h"

#include <vector>

namespace stellar::render {

struct InstanceData {
  float px, py, pz;
  float scale;

  // Per-instance orientation quaternion (w,x,y,z).
  // Identity = (1,0,0,0).
  float qw, qx, qy, qz;

  float cr, cg, cb;
};

class MeshRenderer {
public:
  MeshRenderer() = default;
  ~MeshRenderer();

  MeshRenderer(const MeshRenderer&) = delete;
  MeshRenderer& operator=(const MeshRenderer&) = delete;

  bool init(std::string* outError = nullptr);

  void setMesh(const Mesh* mesh) { mesh_ = mesh; }
  void setTexture(const Texture2D* tex) { tex_ = tex; }

  void setViewProj(const float* view, const float* proj);
  void drawInstances(const std::vector<InstanceData>& instances);

private:
  const Mesh* mesh_{nullptr};
  const Texture2D* tex_{nullptr};

  ShaderProgram shader_{};

  unsigned int instanceVbo_{0};

  float view_[16]{};
  float proj_[16]{};
};

} // namespace stellar::render
