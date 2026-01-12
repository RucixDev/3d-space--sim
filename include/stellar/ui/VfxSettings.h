#pragma once

#include "stellar/core/Types.h"

#include "stellar/render/ProceduralArtifact.h"
#include "stellar/render/ProceduralAsteroid.h"
#include "stellar/render/ProceduralPlanet.h"
#include "stellar/render/ProceduralSky.h"

#include <string>

namespace stellar::ui {

// Persistent configuration for the game's procedural VFX stack.
//
// The VFX Lab started as an experimentation panel; VfxSettings makes those tweaks
// part of the base game by persisting them to disk.
struct VfxSettings {
  int version{2};
  bool autoSaveOnExit{true};

  // --- Background ---
  bool proceduralSkyEnabled{false};
  render::ProceduralSkySettings proceduralSky{};

  bool starfieldEnabled{true};
  bool starfieldTextured{true};
  int starCount{5200};
  double starRadiusU{18000.0};

  bool nebulaEnabled{true};
  int nebulaPuffCount{1400};
  int nebulaVariant{0};
  double nebulaInnerRadiusU{9000.0};
  double nebulaOuterRadiusU{22000.0};
  double nebulaParallax{0.25};
  float nebulaIntensity{1.4f};
  float nebulaOpacity{0.18f};
  float nebulaSizeMinPx{90.0f};
  float nebulaSizeMaxPx{320.0f};
  float nebulaBandPower{1.8f};
  float nebulaTurbulence{0.35f};
  float nebulaTurbulenceSpeed{0.10f};

  // --- Particles ---
  bool particlesEnabled{true};
  bool particlesTextured{true};
  bool thrustersEnabled{true};
  bool impactsEnabled{true};
  bool explosionsEnabled{true};
  float particleIntensity{1.0f};

  // --- GPU dust field (GPGPU) ---
  // Separate high-count particle field simulated on the GPU via ping-pong float textures.
  bool gpuDustEnabled{false};
  bool gpuDustTextured{true};
  int gpuDustDim{256};
  int gpuDustSeed{1337};
  float gpuDustBoundsU{45.0f};
  float gpuDustDrag{0.18f};
  float gpuDustNoiseFreq{0.07f};
  float gpuDustNoiseStrength{1.10f};
  float gpuDustAttract{0.0f};
  float gpuDustMaxSpeed{18.0f};
  float gpuDustPointSizePx{2.25f};
  float gpuDustIntensity{1.0f};
  float gpuDustAlpha{0.65f};

  // --- Procedural asteroid meshes ---
  bool asteroidProcMeshEnabled{true};
  bool asteroidUseRockyTexture{true};
  int asteroidVariantCount{8};
  render::AsteroidParams asteroidParams{};
  core::u64 asteroidStyleNonce{1};

  // --- Procedural SDF artifacts ---
  bool artifactProcEnabled{false};
  bool artifactUseProceduralTexture{true};
  render::SurfaceKind artifactSurfaceKind{render::SurfaceKind::Ice};
  int artifactVariantCount{3};
  int artifactInstanceCount{18};
  double artifactRingRadiusU{1.25};
  double artifactRingHeightU{0.0};
  double artifactScaleU{0.55};
  render::ArtifactParams artifactParams{};
  core::u64 artifactStyleNonce{1};
};

// Default location for the persisted file.
std::string defaultVfxSettingsPath();

// Default settings (used when no file exists).
VfxSettings makeDefaultVfxSettings();

// Save/load in a human-editable text format.
bool saveVfxSettingsToFile(const VfxSettings& settings, const std::string& path);
bool loadVfxSettingsFromFile(const std::string& path, VfxSettings& outSettings);

} // namespace stellar::ui
