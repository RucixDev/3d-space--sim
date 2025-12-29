#pragma once

#include "stellar/render/Mesh.h"
#include "stellar/render/Shader.h"
#include "stellar/render/Texture.h"

#include <vector>

namespace stellar::render {

// NOTE: This struct is tightly packed as floats and uploaded to the GPU.
// Layout must match the vertex shader instance attributes.
struct InstanceData {
  // World position (render units)
  float px, py, pz;

  // Non-uniform scale (render units)
  float sx, sy, sz;

  // Rotation quaternion (x,y,z,w)
  float qx, qy, qz, qw;

  // Base color multiplier
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
  // When enabled, skips directional lighting and renders the mesh as emissive/unlit.
  // Useful for stars, UI 3D previews, and debug visualizations.
  void setUnlit(bool unlit) { unlit_ = unlit; }

  // Point-light position in world/render units.
  //
  // The prototype treats the main star as a point light. By default the light is at
  // the origin, but UI preview scenes (hangar) can set this to a nicer angle.
  void setLightPos(float x, float y, float z) {
    lightPos_[0] = x;
    lightPos_[1] = y;
    lightPos_[2] = z;
  }

  // When enabled, the fragment shader outputs alpha from the bound texture's A channel.
  // This is primarily used for translucent shells like planet cloud layers.
  void setAlphaFromTexture(bool enabled) { alphaFromTexture_ = enabled; }
  // Global alpha multiplier applied when alpha-from-texture is enabled.
  void setAlphaMul(float a) { alphaMul_ = a; }

  void setViewProj(const float* view, const float* proj);
  void drawInstances(const std::vector<InstanceData>& instances);

private:
  const Mesh* mesh_{nullptr};
  const Texture2D* tex_{nullptr};
  bool unlit_{false};

	float lightPos_[3]{0.0f, 0.0f, 0.0f};

  bool alphaFromTexture_{false};
  float alphaMul_{1.0f};

  ShaderProgram shader_{};

  unsigned int instanceVbo_{0};

  float view_[16]{};
  float proj_[16]{};
};

} // namespace stellar::render
