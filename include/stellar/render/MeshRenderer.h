#pragma once

#include "stellar/render/Mesh.h"
#include "stellar/render/Shader.h"
#include "stellar/render/Texture.h"

#include <vector>
#include <cstddef>

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

  void setMesh(const Mesh* mesh) {
    if (mesh_ != mesh) {
      mesh_ = mesh;
      // Instance attributes live on the mesh VAO. If the mesh changes we must
      // re-bind our instance attribute layout once.
      instanceLayoutBoundForMesh_ = nullptr;
    }
  }
  void setTexture(const Texture2D* tex) { tex_ = tex; }
  // Optional tangent-space normal map (RGBA8 normal, encoded in RGB).
  // When null, normal mapping is disabled.
  void setNormalTexture(const Texture2D* tex) { normalTex_ = tex; }
  void setNormalStrength(float strength) { normalStrength_ = strength; }

  // Optional emissive map (added on top of lighting).
  // When set, the fragment shader can reveal the emissive texture on the night side
  // based on the dot(N,L) term.
  void setEmissiveTexture(const Texture2D* tex) { emissiveTex_ = tex; }
  void setEmissiveStrength(float strength) { emissiveStrength_ = strength; }
  // Controls how emissive fades in across the terminator.
  // start/end are applied as smoothstep(start,end,-dot(N,L)).
  void setEmissiveNightFade(float start, float end) {
    emissiveNightFadeStart_ = start;
    emissiveNightFadeEnd_ = end;
  }

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

  // Camera position in world/render units. Used for specular highlights.
  void setCameraPos(float x, float y, float z) {
    cameraPos_[0] = x;
    cameraPos_[1] = y;
    cameraPos_[2] = z;
  }

  // Simple Blinn-Phong specular controls.
  // strength: 0..1
  // shininess: typical range 8..128
  void setSpecular(float strength, float shininess) {
    specularStrength_ = strength;
    shininess_ = shininess;
  }

  // When enabled, the fragment shader outputs alpha from the bound texture's A channel.
  // This is primarily used for translucent shells like planet cloud layers.
  void setAlphaFromTexture(bool enabled) { alphaFromTexture_ = enabled; }
  // Global alpha multiplier applied when alpha-from-texture is enabled.
  void setAlphaMul(float a) { alphaMul_ = a; }

  void setViewProj(const float* view, const float* proj);
  void drawInstances(const std::vector<InstanceData>& instances);

private:
  void ensureInstanceAttribLayoutBound_();

  const Mesh* mesh_{nullptr};
  const Texture2D* tex_{nullptr};
  const Texture2D* normalTex_{nullptr};
  const Texture2D* emissiveTex_{nullptr};
  bool unlit_{false};

	float lightPos_[3]{0.0f, 0.0f, 0.0f};
	float cameraPos_[3]{0.0f, 0.0f, 0.0f};

  float normalStrength_{1.0f};

  float emissiveStrength_{1.0f};
  float emissiveNightFadeStart_{0.02f};
  float emissiveNightFadeEnd_{0.25f};
  float specularStrength_{0.08f};
  float shininess_{48.0f};

  bool alphaFromTexture_{false};
  float alphaMul_{1.0f};

  ShaderProgram shader_{};

  unsigned int instanceVbo_{0};

  // The instance attribute pointers are stored on the mesh VAO. Re-binding
  // them every draw call is surprisingly expensive on some drivers. We
  // therefore bind the layout once per mesh.
  const Mesh* instanceLayoutBoundForMesh_{nullptr};

  // Capacity of the instance buffer allocation (bytes). Avoids reallocating
  // the data store every draw call.
  std::size_t instanceVboCapacityBytes_{0};

  float view_[16]{};
  float proj_[16]{};
};

} // namespace stellar::render
