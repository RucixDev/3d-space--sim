#include "test_harness.h"

#include "stellar/ui/VfxSettings.h"

#include <cmath>
#include <cstdio>

using stellar::ui::VfxSettings;

static bool feq(double a, double b, double eps = 1e-6) {
  return std::fabs(a - b) <= eps;
}

static bool feqf(float a, float b, float eps = 1e-4f) {
  return std::fabs(a - b) <= eps;
}

int test_vfx_settings() {
  int failures = 0;

  VfxSettings s = stellar::ui::makeDefaultVfxSettings();
  s.autoSaveOnExit = false;

  // Background
  s.proceduralSkyEnabled = true;
  s.proceduralSky.seed = 123;
  s.proceduralSky.intensity = 2.5f;
  s.proceduralSky.starIntensity = 1.25f;
  s.proceduralSky.starDensity = 777.0f;
  s.proceduralSky.starProbability = 0.123f;
  s.proceduralSky.nebulaEnabled = true;
  s.proceduralSky.nebulaSteps = 13;
  s.proceduralSky.nebulaWorldScale = 7.5e-5f;
  s.proceduralSky.nebulaDrift = 0.07f;

  s.starfieldEnabled = false;
  s.starfieldTextured = true;
  s.starCount = 1337;
  s.starRadiusU = 15000.0;

  s.nebulaEnabled = true;
  s.nebulaPuffCount = 2048;
  s.nebulaVariant = 42;
  s.nebulaInnerRadiusU = 7777.0;
  s.nebulaOuterRadiusU = 33333.0;
  s.nebulaParallax = 0.3;
  s.nebulaIntensity = 2.2f;
  s.nebulaOpacity = 0.21f;
  s.nebulaSizeMinPx = 12.0f;
  s.nebulaSizeMaxPx = 456.0f;
  s.nebulaBandPower = 2.3f;
  s.nebulaTurbulence = 0.4f;
  s.nebulaTurbulenceSpeed = 0.09f;

  // Particles
  s.particlesEnabled = true;
  s.particlesTextured = false;
  s.thrustersEnabled = true;
  s.impactsEnabled = false;
  s.explosionsEnabled = true;
  s.particleIntensity = 0.33f;

  // Asteroids
  s.asteroidProcMeshEnabled = true;
  s.asteroidUseRockyTexture = false;
  s.asteroidVariantCount = 9;
  s.asteroidStyleNonce = 99;
  s.asteroidParams.slices = 33;
  s.asteroidParams.stacks = 17;
  s.asteroidParams.noiseFrequency = 3.14f;
  s.asteroidParams.noiseAmplitude = 0.25f;
  s.asteroidParams.noiseOctaves = 5;
  s.asteroidParams.craterCount = 12;
  s.asteroidParams.craterRadiusMinDeg = 3.0f;
  s.asteroidParams.craterRadiusMaxDeg = 25.0f;
  s.asteroidParams.craterDepth = 0.12f;
  s.asteroidParams.craterRim = 0.22f;

  // Artifacts
  s.artifactProcEnabled = true;
  s.artifactUseProceduralTexture = true;
  s.artifactSurfaceKind = stellar::render::SurfaceKind::Desert;
  s.artifactVariantCount = 7;
  s.artifactInstanceCount = 11;
  s.artifactRingRadiusU = 2.25;
  s.artifactRingHeightU = -0.5;
  s.artifactScaleU = 0.75;
  s.artifactStyleNonce = 1234;
  s.artifactParams.resolution = 64;
  s.artifactParams.bounds = 1.6f;
  s.artifactParams.iso = 0.1f;
  s.artifactParams.baseRadius = 0.9f;
  s.artifactParams.noiseFrequency = 2.25f;
  s.artifactParams.noiseAmplitude = 0.5f;
  s.artifactParams.noiseOctaves = 4;
  s.artifactParams.blobCount = 8;
  s.artifactParams.cutCount = 3;
  s.artifactParams.normalEps = 0.009f;

  const std::string tmpPath = "vfx_settings_test_tmp.txt";
  CHECK(stellar::ui::saveVfxSettingsToFile(s, tmpPath));

  VfxSettings out = stellar::ui::makeDefaultVfxSettings();
  CHECK(stellar::ui::loadVfxSettingsFromFile(tmpPath, out));

  CHECK(out.autoSaveOnExit == s.autoSaveOnExit);

  CHECK(out.proceduralSkyEnabled == s.proceduralSkyEnabled);
  CHECK(out.proceduralSky.seed == s.proceduralSky.seed);
  CHECK(feqf(out.proceduralSky.intensity, s.proceduralSky.intensity));
  CHECK(feqf(out.proceduralSky.starIntensity, s.proceduralSky.starIntensity));
  CHECK(feqf(out.proceduralSky.starDensity, s.proceduralSky.starDensity));
  CHECK(feqf(out.proceduralSky.starProbability, s.proceduralSky.starProbability));
  CHECK(out.proceduralSky.nebulaEnabled == s.proceduralSky.nebulaEnabled);
  CHECK(out.proceduralSky.nebulaSteps == s.proceduralSky.nebulaSteps);
  CHECK(feqf(out.proceduralSky.nebulaWorldScale, s.proceduralSky.nebulaWorldScale));
  CHECK(feqf(out.proceduralSky.nebulaDrift, s.proceduralSky.nebulaDrift));

  CHECK(out.starfieldEnabled == s.starfieldEnabled);
  CHECK(out.starfieldTextured == s.starfieldTextured);
  CHECK(out.starCount == s.starCount);
  CHECK(feq(out.starRadiusU, s.starRadiusU));

  CHECK(out.nebulaEnabled == s.nebulaEnabled);
  CHECK(out.nebulaPuffCount == s.nebulaPuffCount);
  CHECK(out.nebulaVariant == s.nebulaVariant);
  CHECK(feq(out.nebulaInnerRadiusU, s.nebulaInnerRadiusU));
  CHECK(feq(out.nebulaOuterRadiusU, s.nebulaOuterRadiusU));
  CHECK(feq(out.nebulaParallax, s.nebulaParallax));
  CHECK(feqf(out.nebulaIntensity, s.nebulaIntensity));
  CHECK(feqf(out.nebulaOpacity, s.nebulaOpacity));

  CHECK(out.particlesEnabled == s.particlesEnabled);
  CHECK(out.particlesTextured == s.particlesTextured);
  CHECK(out.thrustersEnabled == s.thrustersEnabled);
  CHECK(out.impactsEnabled == s.impactsEnabled);
  CHECK(out.explosionsEnabled == s.explosionsEnabled);
  CHECK(feqf(out.particleIntensity, s.particleIntensity));

  CHECK(out.asteroidProcMeshEnabled == s.asteroidProcMeshEnabled);
  CHECK(out.asteroidUseRockyTexture == s.asteroidUseRockyTexture);
  CHECK(out.asteroidVariantCount == s.asteroidVariantCount);
  CHECK(out.asteroidStyleNonce == s.asteroidStyleNonce);
  CHECK(out.asteroidParams.slices == s.asteroidParams.slices);
  CHECK(out.asteroidParams.stacks == s.asteroidParams.stacks);
  CHECK(feqf(out.asteroidParams.noiseFrequency, s.asteroidParams.noiseFrequency));
  CHECK(feqf(out.asteroidParams.noiseAmplitude, s.asteroidParams.noiseAmplitude));
  CHECK(out.asteroidParams.noiseOctaves == s.asteroidParams.noiseOctaves);

  CHECK(out.artifactProcEnabled == s.artifactProcEnabled);
  CHECK(out.artifactUseProceduralTexture == s.artifactUseProceduralTexture);
  CHECK(out.artifactSurfaceKind == s.artifactSurfaceKind);
  CHECK(out.artifactVariantCount == s.artifactVariantCount);
  CHECK(out.artifactInstanceCount == s.artifactInstanceCount);
  CHECK(feq(out.artifactRingRadiusU, s.artifactRingRadiusU));
  CHECK(feq(out.artifactRingHeightU, s.artifactRingHeightU));
  CHECK(feq(out.artifactScaleU, s.artifactScaleU));
  CHECK(out.artifactStyleNonce == s.artifactStyleNonce);
  CHECK(out.artifactParams.resolution == s.artifactParams.resolution);
  CHECK(feqf(out.artifactParams.bounds, s.artifactParams.bounds));
  CHECK(feqf(out.artifactParams.iso, s.artifactParams.iso));
  CHECK(feqf(out.artifactParams.baseRadius, s.artifactParams.baseRadius));
  CHECK(feqf(out.artifactParams.noiseFrequency, s.artifactParams.noiseFrequency));
  CHECK(feqf(out.artifactParams.noiseAmplitude, s.artifactParams.noiseAmplitude));
  CHECK(out.artifactParams.noiseOctaves == s.artifactParams.noiseOctaves);
  CHECK(out.artifactParams.blobCount == s.artifactParams.blobCount);
  CHECK(out.artifactParams.cutCount == s.artifactParams.cutCount);
  CHECK(feqf(out.artifactParams.normalEps, s.artifactParams.normalEps));

  (void)std::remove(tmpPath.c_str());
  return failures;
}
