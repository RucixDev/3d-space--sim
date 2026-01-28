#include "stellar/ui/VfxSettings.h"

#include "stellar/core/Clamp.h"
#include "stellar/core/Log.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <string>

namespace stellar::ui {
namespace {

static std::string lowerCopy(std::string s) {
  for (char& c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}

static bool parseBool(const std::string& token, bool* out) {
  const std::string t = lowerCopy(token);
  if (t == "1" || t == "true" || t == "yes" || t == "y" || t == "on") {
    *out = true;
    return true;
  }
  if (t == "0" || t == "false" || t == "no" || t == "n" || t == "off") {
    *out = false;
    return true;
  }
  return false;
}

static const char* surfaceKindToString(render::SurfaceKind k) {
  switch (k) {
    case render::SurfaceKind::Rocky: return "rocky";
    case render::SurfaceKind::Desert: return "desert";
    case render::SurfaceKind::Ocean: return "ocean";
    case render::SurfaceKind::Ice: return "ice";
    case render::SurfaceKind::GasGiant: return "gas";
    case render::SurfaceKind::Star: return "star";
    case render::SurfaceKind::Clouds: return "clouds";
    case render::SurfaceKind::CityLights: return "city_lights";
    default: return "rocky";
  }
}

static bool surfaceKindFromString(const std::string& s, render::SurfaceKind* out) {
  const std::string t = lowerCopy(s);
  if (t == "rocky" || t == "rock") { *out = render::SurfaceKind::Rocky; return true; }
  if (t == "desert" || t == "sand") { *out = render::SurfaceKind::Desert; return true; }
  if (t == "ocean" || t == "water") { *out = render::SurfaceKind::Ocean; return true; }
  if (t == "ice" || t == "icy") { *out = render::SurfaceKind::Ice; return true; }
  if (t == "gas" || t == "gasgiant" || t == "gas_giant") { *out = render::SurfaceKind::GasGiant; return true; }
  if (t == "star" || t == "sun") { *out = render::SurfaceKind::Star; return true; }
  if (t == "clouds" || t == "cloud") { *out = render::SurfaceKind::Clouds; return true; }
  if (t == "city_lights" || t == "citylights" || t == "lights" || t == "city") { *out = render::SurfaceKind::CityLights; return true; }
  return false;
}

} // namespace

std::string defaultVfxSettingsPath() {
  return "vfx_settings.txt";
}

VfxSettings makeDefaultVfxSettings() {
  return VfxSettings{};
}

bool saveVfxSettingsToFile(const VfxSettings& s, const std::string& path) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f.good()) return false;

  f.setf(std::ios::fixed);
  f.precision(6);

  f << "StellarForgeVfxSettings " << s.version << "\n";
  f << "autoSaveOnExit " << (s.autoSaveOnExit ? 1 : 0) << "\n";

  // --- Background ---
  f << "proceduralSkyEnabled " << (s.proceduralSkyEnabled ? 1 : 0) << "\n";
  f << "procSkySeed " << s.proceduralSky.seed << "\n";
  f << "procSkyIntensity " << s.proceduralSky.intensity << "\n";
  f << "procSkyStarIntensity " << s.proceduralSky.starIntensity << "\n";
  f << "procSkyStarDensity " << s.proceduralSky.starDensity << "\n";
  f << "procSkyStarProbability " << s.proceduralSky.starProbability << "\n";
  f << "procSkyStarSize " << s.proceduralSky.starSize << "\n";
  f << "procSkyStarTwinkle " << s.proceduralSky.starTwinkle << "\n";
  f << "procSkyStarClusterStrength " << s.proceduralSky.starClusterStrength << "\n";
  f << "procSkyStarClusterFrequency " << s.proceduralSky.starClusterFrequency << "\n";
  f << "procSkyStarSpikeStrength " << s.proceduralSky.starSpikeStrength << "\n";
  f << "procSkyMilkyWayEnabled " << (s.proceduralSky.milkyWayEnabled ? 1 : 0) << "\n";
  f << "procSkyMilkyWayIntensity " << s.proceduralSky.milkyWayIntensity << "\n";
  f << "procSkyMilkyWayFrequency " << s.proceduralSky.milkyWayFrequency << "\n";
  f << "procSkyMilkyWayBandPower " << s.proceduralSky.milkyWayBandPower << "\n";
  f << "procSkyMilkyWayDust " << s.proceduralSky.milkyWayDust << "\n";
  f << "procSkyNebulaEnabled " << (s.proceduralSky.nebulaEnabled ? 1 : 0) << "\n";
  f << "procSkyNebulaIntensity " << s.proceduralSky.nebulaIntensity << "\n";
  f << "procSkyNebulaFrequency " << s.proceduralSky.nebulaFrequency << "\n";
  f << "procSkyNebulaThreshold " << s.proceduralSky.nebulaThreshold << "\n";
  f << "procSkyNebulaSoftness " << s.proceduralSky.nebulaSoftness << "\n";
  f << "procSkyNebulaBandPower " << s.proceduralSky.nebulaBandPower << "\n";
  f << "procSkyNebulaParallax " << s.proceduralSky.nebulaParallax << "\n";
  f << "procSkyNebulaSteps " << s.proceduralSky.nebulaSteps << "\n";
  f << "procSkyNebulaWorldScale " << s.proceduralSky.nebulaWorldScale << "\n";
  f << "procSkyNebulaDrift " << s.proceduralSky.nebulaDrift << "\n";

  f << "starfieldEnabled " << (s.starfieldEnabled ? 1 : 0) << "\n";
  f << "starfieldTextured " << (s.starfieldTextured ? 1 : 0) << "\n";
  f << "starfieldGpuEnabled " << (s.starfieldGpuEnabled ? 1 : 0) << "\n";
  f << "starCount " << s.starCount << "\n";
  f << "starRadiusU " << s.starRadiusU << "\n";

  f << "nebulaEnabled " << (s.nebulaEnabled ? 1 : 0) << "\n";
  f << "nebulaPuffCount " << s.nebulaPuffCount << "\n";
  f << "nebulaVariant " << s.nebulaVariant << "\n";
  f << "nebulaInnerRadiusU " << s.nebulaInnerRadiusU << "\n";
  f << "nebulaOuterRadiusU " << s.nebulaOuterRadiusU << "\n";
  f << "nebulaParallax " << s.nebulaParallax << "\n";
  f << "nebulaIntensity " << s.nebulaIntensity << "\n";
  f << "nebulaOpacity " << s.nebulaOpacity << "\n";
  f << "nebulaSizeMinPx " << s.nebulaSizeMinPx << "\n";
  f << "nebulaSizeMaxPx " << s.nebulaSizeMaxPx << "\n";
  f << "nebulaBandPower " << s.nebulaBandPower << "\n";
  f << "nebulaTurbulence " << s.nebulaTurbulence << "\n";
  f << "nebulaTurbulenceSpeed " << s.nebulaTurbulenceSpeed << "\n";

  // --- Particles ---
  f << "particlesEnabled " << (s.particlesEnabled ? 1 : 0) << "\n";
  f << "particlesTextured " << (s.particlesTextured ? 1 : 0) << "\n";
  f << "thrustersEnabled " << (s.thrustersEnabled ? 1 : 0) << "\n";
  f << "impactsEnabled " << (s.impactsEnabled ? 1 : 0) << "\n";
  f << "explosionsEnabled " << (s.explosionsEnabled ? 1 : 0) << "\n";
  f << "particleIntensity " << s.particleIntensity << "\n";

  // --- GPU dust field (GPGPU) ---
  f << "gpuDustEnabled " << (s.gpuDustEnabled ? 1 : 0) << "\n";
  f << "gpuDustTextured " << (s.gpuDustTextured ? 1 : 0) << "\n";
  f << "gpuDustDim " << s.gpuDustDim << "\n";
  f << "gpuDustSeed " << s.gpuDustSeed << "\n";
  f << "gpuDustBoundsU " << s.gpuDustBoundsU << "\n";
  f << "gpuDustDrag " << s.gpuDustDrag << "\n";
  f << "gpuDustNoiseFreq " << s.gpuDustNoiseFreq << "\n";
  f << "gpuDustNoiseStrength " << s.gpuDustNoiseStrength << "\n";
  f << "gpuDustAttract " << s.gpuDustAttract << "\n";
  f << "gpuDustMaxSpeed " << s.gpuDustMaxSpeed << "\n";
  f << "gpuDustPointSizePx " << s.gpuDustPointSizePx << "\n";
  f << "gpuDustIntensity " << s.gpuDustIntensity << "\n";
  f << "gpuDustAlpha " << s.gpuDustAlpha << "\n";

  // --- Procedural asteroid meshes ---
  f << "asteroidProcMeshEnabled " << (s.asteroidProcMeshEnabled ? 1 : 0) << "\n";
  f << "asteroidUseSdfMesher " << (s.asteroidUseSdfMesher ? 1 : 0) << "\n";
  f << "asteroidUseRockyTexture " << (s.asteroidUseRockyTexture ? 1 : 0) << "\n";
  f << "asteroidVariantCount " << s.asteroidVariantCount << "\n";
  f << "asteroidStyleNonce " << (unsigned long long)s.asteroidStyleNonce << "\n";

  const auto& ap = s.asteroidParams;
  f << "asteroidSlices " << ap.slices << "\n";
  f << "asteroidStacks " << ap.stacks << "\n";
  f << "asteroidBaseRadius " << ap.baseRadius << "\n";
  f << "asteroidNoiseFrequency " << ap.noiseFrequency << "\n";
  f << "asteroidNoiseAmplitude " << ap.noiseAmplitude << "\n";
  f << "asteroidNoiseOctaves " << ap.noiseOctaves << "\n";
  f << "asteroidNoiseLacunarity " << ap.noiseLacunarity << "\n";
  f << "asteroidNoiseGain " << ap.noiseGain << "\n";
  f << "asteroidCraterCount " << ap.craterCount << "\n";
  f << "asteroidCraterRadiusMinDeg " << ap.craterRadiusMinDeg << "\n";
  f << "asteroidCraterRadiusMaxDeg " << ap.craterRadiusMaxDeg << "\n";
  f << "asteroidCraterDepth " << ap.craterDepth << "\n";
  f << "asteroidCraterRim " << ap.craterRim << "\n";
  f << "asteroidMinRadius " << ap.minRadius << "\n";
  f << "asteroidMaxRadius " << ap.maxRadius << "\n";

  const auto& asp = s.asteroidSdfParams;
  f << "asteroidSdfResolution " << asp.resolution << "\n";
  f << "asteroidSdfBounds " << asp.bounds << "\n";
  f << "asteroidSdfIso " << asp.iso << "\n";
  f << "asteroidSdfBaseRadius " << asp.baseRadius << "\n";
  f << "asteroidSdfAxisScaleX " << asp.axisScaleX << "\n";
  f << "asteroidSdfAxisScaleY " << asp.axisScaleY << "\n";
  f << "asteroidSdfAxisScaleZ " << asp.axisScaleZ << "\n";
  f << "asteroidSdfNoise1Frequency " << asp.noise1Frequency << "\n";
  f << "asteroidSdfNoise1Amplitude " << asp.noise1Amplitude << "\n";
  f << "asteroidSdfNoise1Octaves " << asp.noise1Octaves << "\n";
  f << "asteroidSdfNoise1Lacunarity " << asp.noise1Lacunarity << "\n";
  f << "asteroidSdfNoise1Gain " << asp.noise1Gain << "\n";
  f << "asteroidSdfNoise2Frequency " << asp.noise2Frequency << "\n";
  f << "asteroidSdfNoise2Amplitude " << asp.noise2Amplitude << "\n";
  f << "asteroidSdfNoise2Octaves " << asp.noise2Octaves << "\n";
  f << "asteroidSdfNoise2Lacunarity " << asp.noise2Lacunarity << "\n";
  f << "asteroidSdfNoise2Gain " << asp.noise2Gain << "\n";
  f << "asteroidSdfCraterCount " << asp.craterCount << "\n";
  f << "asteroidSdfCraterRadiusMinDeg " << asp.craterRadiusMinDeg << "\n";
  f << "asteroidSdfCraterRadiusMaxDeg " << asp.craterRadiusMaxDeg << "\n";
  f << "asteroidSdfCraterDepth " << asp.craterDepth << "\n";
  f << "asteroidSdfCraterSmoothK " << asp.craterSmoothK << "\n";
  f << "asteroidSdfCutCount " << asp.cutCount << "\n";
  f << "asteroidSdfCutOffsetMin " << asp.cutOffsetMin << "\n";
  f << "asteroidSdfCutOffsetMax " << asp.cutOffsetMax << "\n";
  f << "asteroidSdfGrooveStrength " << asp.grooveStrength << "\n";
  f << "asteroidSdfGrooveFrequency " << asp.grooveFrequency << "\n";
  f << "asteroidSdfNormalEps " << asp.normalEps << "\n";

  // --- Procedural SDF artifacts ---
  f << "artifactProcEnabled " << (s.artifactProcEnabled ? 1 : 0) << "\n";
  f << "artifactUseProceduralTexture " << (s.artifactUseProceduralTexture ? 1 : 0) << "\n";
  f << "artifactSurfaceKind " << surfaceKindToString(s.artifactSurfaceKind) << "\n";
  f << "artifactVariantCount " << s.artifactVariantCount << "\n";
  f << "artifactInstanceCount " << s.artifactInstanceCount << "\n";
  f << "artifactRingRadiusU " << s.artifactRingRadiusU << "\n";
  f << "artifactRingHeightU " << s.artifactRingHeightU << "\n";
  f << "artifactScaleU " << s.artifactScaleU << "\n";
  f << "artifactStyleNonce " << (unsigned long long)s.artifactStyleNonce << "\n";

  const auto& af = s.artifactParams;
  f << "artifactResolution " << af.resolution << "\n";
  f << "artifactBounds " << af.bounds << "\n";
  f << "artifactIso " << af.iso << "\n";
  f << "artifactBaseRadius " << af.baseRadius << "\n";
  f << "artifactNoiseFrequency " << af.noiseFrequency << "\n";
  f << "artifactNoiseAmplitude " << af.noiseAmplitude << "\n";
  f << "artifactNoiseOctaves " << af.noiseOctaves << "\n";
  f << "artifactNoiseLacunarity " << af.noiseLacunarity << "\n";
  f << "artifactNoiseGain " << af.noiseGain << "\n";
  f << "artifactBlobCount " << af.blobCount << "\n";
  f << "artifactBlobRadiusMin " << af.blobRadiusMin << "\n";
  f << "artifactBlobRadiusMax " << af.blobRadiusMax << "\n";
  f << "artifactBlobCenterRadius " << af.blobCenterRadius << "\n";
  f << "artifactBlobSmoothK " << af.blobSmoothK << "\n";
  f << "artifactCutCount " << af.cutCount << "\n";
  f << "artifactCutOffsetMin " << af.cutOffsetMin << "\n";
  f << "artifactCutOffsetMax " << af.cutOffsetMax << "\n";
  f << "artifactGrooveStrength " << af.grooveStrength << "\n";
  f << "artifactGrooveFrequency " << af.grooveFrequency << "\n";
  f << "artifactNormalEps " << af.normalEps << "\n";

  return true;
}

bool loadVfxSettingsFromFile(const std::string& path, VfxSettings& outSettings) {
  std::ifstream f(path);
  if (!f.good()) return false;

  std::string header;
  int version = 0;
  {
    std::getline(f, header);
    std::stringstream ss(header);
    std::string magic;
    ss >> magic >> version;
    if (magic != "StellarForgeVfxSettings") {
      core::log(core::LogLevel::Warn, "VFX settings: bad header in " + path);
      return false;
    }
  }

  // Start from defaults (current version) and overlay any keys from disk.
  // This keeps newly added fields working even when loading older files.
  VfxSettings s = makeDefaultVfxSettings();

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::stringstream ss(line);
    std::string key;
    ss >> key;
    if (key.empty()) continue;

    // bool
    if (key == "autoSaveOnExit") {
      std::string v; ss >> v; bool b = s.autoSaveOnExit; if (parseBool(v, &b)) s.autoSaveOnExit = b;
      continue;
    }

    if (key == "proceduralSkyEnabled") {
      std::string v; ss >> v; bool b = s.proceduralSkyEnabled; if (parseBool(v, &b)) s.proceduralSkyEnabled = b;
      continue;
    }
    if (key == "procSkyMilkyWayEnabled") {
      std::string v; ss >> v;
      bool b = s.proceduralSky.milkyWayEnabled;
      if (parseBool(v, &b)) s.proceduralSky.milkyWayEnabled = b;
      continue;
    }
    if (key == "procSkyNebulaEnabled") {
      std::string v; ss >> v; bool b = s.proceduralSky.nebulaEnabled; if (parseBool(v, &b)) s.proceduralSky.nebulaEnabled = b;
      continue;
    }

    if (key == "starfieldEnabled") {
      std::string v; ss >> v; bool b = s.starfieldEnabled; if (parseBool(v, &b)) s.starfieldEnabled = b;
      continue;
    }
    if (key == "starfieldTextured") {
      std::string v; ss >> v; bool b = s.starfieldTextured; if (parseBool(v, &b)) s.starfieldTextured = b;
      continue;
    }
    if (key == "starfieldGpuEnabled") {
      std::string v; ss >> v; bool b = s.starfieldGpuEnabled; if (parseBool(v, &b)) s.starfieldGpuEnabled = b;
      continue;
    }

    if (key == "nebulaEnabled") {
      std::string v; ss >> v; bool b = s.nebulaEnabled; if (parseBool(v, &b)) s.nebulaEnabled = b;
      continue;
    }

    if (key == "particlesEnabled") {
      std::string v; ss >> v; bool b = s.particlesEnabled; if (parseBool(v, &b)) s.particlesEnabled = b;
      continue;
    }
    if (key == "particlesTextured") {
      std::string v; ss >> v; bool b = s.particlesTextured; if (parseBool(v, &b)) s.particlesTextured = b;
      continue;
    }
    if (key == "thrustersEnabled") {
      std::string v; ss >> v; bool b = s.thrustersEnabled; if (parseBool(v, &b)) s.thrustersEnabled = b;
      continue;
    }
    if (key == "impactsEnabled") {
      std::string v; ss >> v; bool b = s.impactsEnabled; if (parseBool(v, &b)) s.impactsEnabled = b;
      continue;
    }
    if (key == "explosionsEnabled") {
      std::string v; ss >> v; bool b = s.explosionsEnabled; if (parseBool(v, &b)) s.explosionsEnabled = b;
      continue;
    }

    if (key == "gpuDustEnabled") {
      std::string v; ss >> v; bool b = s.gpuDustEnabled; if (parseBool(v, &b)) s.gpuDustEnabled = b;
      continue;
    }
    if (key == "gpuDustTextured") {
      std::string v; ss >> v; bool b = s.gpuDustTextured; if (parseBool(v, &b)) s.gpuDustTextured = b;
      continue;
    }

    if (key == "asteroidProcMeshEnabled") {
      std::string v; ss >> v; bool b = s.asteroidProcMeshEnabled; if (parseBool(v, &b)) s.asteroidProcMeshEnabled = b;
      continue;
    }
    if (key == "asteroidUseSdfMesher") {
      std::string v; ss >> v; bool b = s.asteroidUseSdfMesher; if (parseBool(v, &b)) s.asteroidUseSdfMesher = b;
      continue;
    }
    if (key == "asteroidUseRockyTexture") {
      std::string v; ss >> v; bool b = s.asteroidUseRockyTexture; if (parseBool(v, &b)) s.asteroidUseRockyTexture = b;
      continue;
    }

    if (key == "artifactProcEnabled") {
      std::string v; ss >> v; bool b = s.artifactProcEnabled; if (parseBool(v, &b)) s.artifactProcEnabled = b;
      continue;
    }
    if (key == "artifactUseProceduralTexture") {
      std::string v; ss >> v; bool b = s.artifactUseProceduralTexture; if (parseBool(v, &b)) s.artifactUseProceduralTexture = b;
      continue;
    }

    if (key == "artifactSurfaceKind") {
      std::string v; ss >> v;
      render::SurfaceKind k = s.artifactSurfaceKind;
      if (surfaceKindFromString(v, &k)) s.artifactSurfaceKind = k;
      continue;
    }

    // ints
    if (key == "procSkySeed") { ss >> s.proceduralSky.seed; continue; }
    if (key == "procSkyNebulaSteps") { ss >> s.proceduralSky.nebulaSteps; continue; }

    if (key == "starCount") { ss >> s.starCount; continue; }
    if (key == "nebulaPuffCount") { ss >> s.nebulaPuffCount; continue; }
    if (key == "nebulaVariant") { ss >> s.nebulaVariant; continue; }

    if (key == "gpuDustDim") { ss >> s.gpuDustDim; continue; }
    if (key == "gpuDustSeed") { ss >> s.gpuDustSeed; continue; }

    if (key == "artifactVariantCount") { ss >> s.artifactVariantCount; continue; }
    if (key == "artifactInstanceCount") { ss >> s.artifactInstanceCount; continue; }

    if (key == "asteroidVariantCount") { ss >> s.asteroidVariantCount; continue; }
    if (key == "asteroidSlices") { ss >> s.asteroidParams.slices; continue; }
    if (key == "asteroidStacks") { ss >> s.asteroidParams.stacks; continue; }
    if (key == "asteroidNoiseOctaves") { ss >> s.asteroidParams.noiseOctaves; continue; }
    if (key == "asteroidCraterCount") { ss >> s.asteroidParams.craterCount; continue; }

    if (key == "asteroidSdfResolution") { ss >> s.asteroidSdfParams.resolution; continue; }
    if (key == "asteroidSdfNoise1Octaves") { ss >> s.asteroidSdfParams.noise1Octaves; continue; }
    if (key == "asteroidSdfNoise2Octaves") { ss >> s.asteroidSdfParams.noise2Octaves; continue; }
    if (key == "asteroidSdfCraterCount") { ss >> s.asteroidSdfParams.craterCount; continue; }
    if (key == "asteroidSdfCutCount") { ss >> s.asteroidSdfParams.cutCount; continue; }

    if (key == "artifactResolution") { ss >> s.artifactParams.resolution; continue; }
    if (key == "artifactNoiseOctaves") { ss >> s.artifactParams.noiseOctaves; continue; }
    if (key == "artifactBlobCount") { ss >> s.artifactParams.blobCount; continue; }
    if (key == "artifactCutCount") { ss >> s.artifactParams.cutCount; continue; }

    // u64
    if (key == "asteroidStyleNonce") {
      unsigned long long v = 0; ss >> v; s.asteroidStyleNonce = (core::u64)v; continue;
    }
    if (key == "artifactStyleNonce") {
      unsigned long long v = 0; ss >> v; s.artifactStyleNonce = (core::u64)v; continue;
    }

    // doubles
    if (key == "starRadiusU") { ss >> s.starRadiusU; continue; }
    if (key == "nebulaInnerRadiusU") { ss >> s.nebulaInnerRadiusU; continue; }
    if (key == "nebulaOuterRadiusU") { ss >> s.nebulaOuterRadiusU; continue; }
    if (key == "nebulaParallax") { ss >> s.nebulaParallax; continue; }
    if (key == "artifactRingRadiusU") { ss >> s.artifactRingRadiusU; continue; }
    if (key == "artifactRingHeightU") { ss >> s.artifactRingHeightU; continue; }
    if (key == "artifactScaleU") { ss >> s.artifactScaleU; continue; }

    // floats
    if (key == "procSkyIntensity") { ss >> s.proceduralSky.intensity; continue; }
    if (key == "procSkyStarIntensity") { ss >> s.proceduralSky.starIntensity; continue; }
    if (key == "procSkyStarDensity") { ss >> s.proceduralSky.starDensity; continue; }
    if (key == "procSkyStarProbability") { ss >> s.proceduralSky.starProbability; continue; }
    if (key == "procSkyStarSize") { ss >> s.proceduralSky.starSize; continue; }
    if (key == "procSkyStarTwinkle") { ss >> s.proceduralSky.starTwinkle; continue; }
    if (key == "procSkyStarClusterStrength") { ss >> s.proceduralSky.starClusterStrength; continue; }
    if (key == "procSkyStarClusterFrequency") { ss >> s.proceduralSky.starClusterFrequency; continue; }
    if (key == "procSkyStarSpikeStrength") { ss >> s.proceduralSky.starSpikeStrength; continue; }
    if (key == "procSkyMilkyWayIntensity") { ss >> s.proceduralSky.milkyWayIntensity; continue; }
    if (key == "procSkyMilkyWayFrequency") { ss >> s.proceduralSky.milkyWayFrequency; continue; }
    if (key == "procSkyMilkyWayBandPower") { ss >> s.proceduralSky.milkyWayBandPower; continue; }
    if (key == "procSkyMilkyWayDust") { ss >> s.proceduralSky.milkyWayDust; continue; }
    if (key == "procSkyNebulaIntensity") { ss >> s.proceduralSky.nebulaIntensity; continue; }
    if (key == "procSkyNebulaFrequency") { ss >> s.proceduralSky.nebulaFrequency; continue; }
    if (key == "procSkyNebulaThreshold") { ss >> s.proceduralSky.nebulaThreshold; continue; }
    if (key == "procSkyNebulaSoftness") { ss >> s.proceduralSky.nebulaSoftness; continue; }
    if (key == "procSkyNebulaBandPower") { ss >> s.proceduralSky.nebulaBandPower; continue; }
    if (key == "procSkyNebulaParallax") { ss >> s.proceduralSky.nebulaParallax; continue; }
    if (key == "procSkyNebulaWorldScale") { ss >> s.proceduralSky.nebulaWorldScale; continue; }
    if (key == "procSkyNebulaDrift") { ss >> s.proceduralSky.nebulaDrift; continue; }

    if (key == "nebulaIntensity") { ss >> s.nebulaIntensity; continue; }
    if (key == "nebulaOpacity") { ss >> s.nebulaOpacity; continue; }
    if (key == "nebulaSizeMinPx") { ss >> s.nebulaSizeMinPx; continue; }
    if (key == "nebulaSizeMaxPx") { ss >> s.nebulaSizeMaxPx; continue; }
    if (key == "nebulaBandPower") { ss >> s.nebulaBandPower; continue; }
    if (key == "nebulaTurbulence") { ss >> s.nebulaTurbulence; continue; }
    if (key == "nebulaTurbulenceSpeed") { ss >> s.nebulaTurbulenceSpeed; continue; }

    if (key == "particleIntensity") { ss >> s.particleIntensity; continue; }

    if (key == "gpuDustBoundsU") { ss >> s.gpuDustBoundsU; continue; }
    if (key == "gpuDustDrag") { ss >> s.gpuDustDrag; continue; }
    if (key == "gpuDustNoiseFreq") { ss >> s.gpuDustNoiseFreq; continue; }
    if (key == "gpuDustNoiseStrength") { ss >> s.gpuDustNoiseStrength; continue; }
    if (key == "gpuDustAttract") { ss >> s.gpuDustAttract; continue; }
    if (key == "gpuDustMaxSpeed") { ss >> s.gpuDustMaxSpeed; continue; }
    if (key == "gpuDustPointSizePx") { ss >> s.gpuDustPointSizePx; continue; }
    if (key == "gpuDustIntensity") { ss >> s.gpuDustIntensity; continue; }
    if (key == "gpuDustAlpha") { ss >> s.gpuDustAlpha; continue; }

    if (key == "asteroidBaseRadius") { ss >> s.asteroidParams.baseRadius; continue; }
    if (key == "asteroidNoiseFrequency") { ss >> s.asteroidParams.noiseFrequency; continue; }
    if (key == "asteroidNoiseAmplitude") { ss >> s.asteroidParams.noiseAmplitude; continue; }
    if (key == "asteroidNoiseLacunarity") { ss >> s.asteroidParams.noiseLacunarity; continue; }
    if (key == "asteroidNoiseGain") { ss >> s.asteroidParams.noiseGain; continue; }
    if (key == "asteroidCraterRadiusMinDeg") { ss >> s.asteroidParams.craterRadiusMinDeg; continue; }
    if (key == "asteroidCraterRadiusMaxDeg") { ss >> s.asteroidParams.craterRadiusMaxDeg; continue; }
    if (key == "asteroidCraterDepth") { ss >> s.asteroidParams.craterDepth; continue; }
    if (key == "asteroidCraterRim") { ss >> s.asteroidParams.craterRim; continue; }
    if (key == "asteroidMinRadius") { ss >> s.asteroidParams.minRadius; continue; }
    if (key == "asteroidMaxRadius") { ss >> s.asteroidParams.maxRadius; continue; }

    // SDF asteroid params (isosurface mesher)
    if (key == "asteroidSdfBounds") { ss >> s.asteroidSdfParams.bounds; continue; }
    if (key == "asteroidSdfIso") { ss >> s.asteroidSdfParams.iso; continue; }
    if (key == "asteroidSdfBaseRadius") { ss >> s.asteroidSdfParams.baseRadius; continue; }
    if (key == "asteroidSdfAxisScaleX") { ss >> s.asteroidSdfParams.axisScaleX; continue; }
    if (key == "asteroidSdfAxisScaleY") { ss >> s.asteroidSdfParams.axisScaleY; continue; }
    if (key == "asteroidSdfAxisScaleZ") { ss >> s.asteroidSdfParams.axisScaleZ; continue; }

    if (key == "asteroidSdfNoise1Frequency") { ss >> s.asteroidSdfParams.noise1Frequency; continue; }
    if (key == "asteroidSdfNoise1Amplitude") { ss >> s.asteroidSdfParams.noise1Amplitude; continue; }
    if (key == "asteroidSdfNoise1Lacunarity") { ss >> s.asteroidSdfParams.noise1Lacunarity; continue; }
    if (key == "asteroidSdfNoise1Gain") { ss >> s.asteroidSdfParams.noise1Gain; continue; }

    if (key == "asteroidSdfNoise2Frequency") { ss >> s.asteroidSdfParams.noise2Frequency; continue; }
    if (key == "asteroidSdfNoise2Amplitude") { ss >> s.asteroidSdfParams.noise2Amplitude; continue; }
    if (key == "asteroidSdfNoise2Lacunarity") { ss >> s.asteroidSdfParams.noise2Lacunarity; continue; }
    if (key == "asteroidSdfNoise2Gain") { ss >> s.asteroidSdfParams.noise2Gain; continue; }

    if (key == "asteroidSdfCraterRadiusMinDeg") { ss >> s.asteroidSdfParams.craterRadiusMinDeg; continue; }
    if (key == "asteroidSdfCraterRadiusMaxDeg") { ss >> s.asteroidSdfParams.craterRadiusMaxDeg; continue; }
    if (key == "asteroidSdfCraterDepth") { ss >> s.asteroidSdfParams.craterDepth; continue; }
    if (key == "asteroidSdfCraterSmoothK") { ss >> s.asteroidSdfParams.craterSmoothK; continue; }

    if (key == "asteroidSdfCutOffsetMin") { ss >> s.asteroidSdfParams.cutOffsetMin; continue; }
    if (key == "asteroidSdfCutOffsetMax") { ss >> s.asteroidSdfParams.cutOffsetMax; continue; }
    if (key == "asteroidSdfGrooveStrength") { ss >> s.asteroidSdfParams.grooveStrength; continue; }
    if (key == "asteroidSdfGrooveFrequency") { ss >> s.asteroidSdfParams.grooveFrequency; continue; }
    if (key == "asteroidSdfNormalEps") { ss >> s.asteroidSdfParams.normalEps; continue; }

    if (key == "artifactBounds") { ss >> s.artifactParams.bounds; continue; }
    if (key == "artifactIso") { ss >> s.artifactParams.iso; continue; }
    if (key == "artifactBaseRadius") { ss >> s.artifactParams.baseRadius; continue; }
    if (key == "artifactNoiseFrequency") { ss >> s.artifactParams.noiseFrequency; continue; }
    if (key == "artifactNoiseAmplitude") { ss >> s.artifactParams.noiseAmplitude; continue; }
    if (key == "artifactNoiseLacunarity") { ss >> s.artifactParams.noiseLacunarity; continue; }
    if (key == "artifactNoiseGain") { ss >> s.artifactParams.noiseGain; continue; }
    if (key == "artifactBlobRadiusMin") { ss >> s.artifactParams.blobRadiusMin; continue; }
    if (key == "artifactBlobRadiusMax") { ss >> s.artifactParams.blobRadiusMax; continue; }
    if (key == "artifactBlobCenterRadius") { ss >> s.artifactParams.blobCenterRadius; continue; }
    if (key == "artifactBlobSmoothK") { ss >> s.artifactParams.blobSmoothK; continue; }
    if (key == "artifactCutOffsetMin") { ss >> s.artifactParams.cutOffsetMin; continue; }
    if (key == "artifactCutOffsetMax") { ss >> s.artifactParams.cutOffsetMax; continue; }
    if (key == "artifactGrooveStrength") { ss >> s.artifactParams.grooveStrength; continue; }
    if (key == "artifactGrooveFrequency") { ss >> s.artifactParams.grooveFrequency; continue; }
    if (key == "artifactNormalEps") { ss >> s.artifactParams.normalEps; continue; }

    // Unknown keys are ignored for forward-compat.
  }

  // --- Clamp unsafe values ---
  s.starCount = std::clamp(s.starCount, 0, 20000);
  s.starRadiusU = std::clamp(s.starRadiusU, 2500.0, 50000.0);

  s.nebulaPuffCount = std::clamp(s.nebulaPuffCount, 0, 8000);
  s.nebulaVariant = std::clamp(s.nebulaVariant, 0, 31);
  s.nebulaInnerRadiusU = std::clamp(s.nebulaInnerRadiusU, 1000.0, 60000.0);
  s.nebulaOuterRadiusU = std::clamp(s.nebulaOuterRadiusU, s.nebulaInnerRadiusU + 1.0, 80000.0);
  s.nebulaParallax = std::clamp(s.nebulaParallax, 0.0, 1.0);
  s.nebulaIntensity = std::clamp(s.nebulaIntensity, 0.0f, 10.0f);
  s.nebulaOpacity = std::clamp(s.nebulaOpacity, 0.0f, 1.0f);
  s.nebulaSizeMinPx = std::clamp(s.nebulaSizeMinPx, 2.0f, 1200.0f);
  s.nebulaSizeMaxPx = std::clamp(s.nebulaSizeMaxPx, s.nebulaSizeMinPx, 2500.0f);
  s.nebulaBandPower = std::clamp(s.nebulaBandPower, 0.2f, 6.0f);
  s.nebulaTurbulence = std::clamp(s.nebulaTurbulence, 0.0f, 3.0f);
  s.nebulaTurbulenceSpeed = std::clamp(s.nebulaTurbulenceSpeed, 0.0f, 5.0f);

  s.particleIntensity = std::clamp(s.particleIntensity, 0.0f, 10.0f);

  s.gpuDustDim = std::clamp(s.gpuDustDim, 64, 1024);
  s.gpuDustSeed = std::max(s.gpuDustSeed, 0);
  s.gpuDustBoundsU = std::clamp(s.gpuDustBoundsU, 2.0f, 800.0f);
  s.gpuDustDrag = std::clamp(s.gpuDustDrag, 0.0f, 25.0f);
  s.gpuDustNoiseFreq = std::clamp(s.gpuDustNoiseFreq, 0.0f, 2.0f);
  s.gpuDustNoiseStrength = std::clamp(s.gpuDustNoiseStrength, 0.0f, 50.0f);
  s.gpuDustAttract = std::clamp(s.gpuDustAttract, 0.0f, 5.0f);
  s.gpuDustMaxSpeed = std::clamp(s.gpuDustMaxSpeed, 0.05f, 500.0f);
  s.gpuDustPointSizePx = std::clamp(s.gpuDustPointSizePx, 0.05f, 64.0f);
  s.gpuDustIntensity = std::clamp(s.gpuDustIntensity, 0.0f, 50.0f);
  s.gpuDustAlpha = std::clamp(s.gpuDustAlpha, 0.0f, 1.0f);

  s.asteroidVariantCount = std::clamp(s.asteroidVariantCount, 1, 32);
  s.asteroidParams.slices = std::clamp(s.asteroidParams.slices, 8, 128);
  s.asteroidParams.stacks = std::clamp(s.asteroidParams.stacks, 6, 96);
  s.asteroidParams.noiseOctaves = std::clamp(s.asteroidParams.noiseOctaves, 1, 12);
  s.asteroidParams.craterCount = std::clamp(s.asteroidParams.craterCount, 0, 96);
  s.asteroidParams.baseRadius = std::clamp(s.asteroidParams.baseRadius, 0.25f, 4.0f);
  s.asteroidParams.noiseFrequency = std::clamp(s.asteroidParams.noiseFrequency, 0.10f, 12.0f);
  s.asteroidParams.noiseAmplitude = std::clamp(s.asteroidParams.noiseAmplitude, 0.0f, 2.0f);
  s.asteroidParams.noiseLacunarity = std::clamp(s.asteroidParams.noiseLacunarity, 1.0f, 4.0f);
  s.asteroidParams.noiseGain = std::clamp(s.asteroidParams.noiseGain, 0.1f, 0.95f);
  s.asteroidParams.craterRadiusMinDeg = std::clamp(s.asteroidParams.craterRadiusMinDeg, 0.1f, 60.0f);
  s.asteroidParams.craterRadiusMaxDeg = std::clamp(s.asteroidParams.craterRadiusMaxDeg, s.asteroidParams.craterRadiusMinDeg, 85.0f);
  s.asteroidParams.craterDepth = std::clamp(s.asteroidParams.craterDepth, 0.0f, 0.65f);
  s.asteroidParams.craterRim = std::clamp(s.asteroidParams.craterRim, 0.0f, 1.0f);
  s.asteroidParams.minRadius = std::clamp(s.asteroidParams.minRadius, 0.05f, 10.0f);
  s.asteroidParams.maxRadius = std::max(s.asteroidParams.maxRadius, s.asteroidParams.minRadius);

  s.asteroidSdfParams.resolution = std::clamp(s.asteroidSdfParams.resolution, 12, 192);
  s.asteroidSdfParams.bounds = std::clamp(s.asteroidSdfParams.bounds, 0.25f, 6.0f);
  s.asteroidSdfParams.baseRadius = std::clamp(s.asteroidSdfParams.baseRadius, 0.25f, 4.0f);
  s.asteroidSdfParams.axisScaleX = std::clamp(s.asteroidSdfParams.axisScaleX, 0.15f, 4.0f);
  s.asteroidSdfParams.axisScaleY = std::clamp(s.asteroidSdfParams.axisScaleY, 0.15f, 4.0f);
  s.asteroidSdfParams.axisScaleZ = std::clamp(s.asteroidSdfParams.axisScaleZ, 0.15f, 4.0f);

  s.asteroidSdfParams.noise1Octaves = std::clamp(s.asteroidSdfParams.noise1Octaves, 1, 12);
  s.asteroidSdfParams.noise2Octaves = std::clamp(s.asteroidSdfParams.noise2Octaves, 1, 12);
  s.asteroidSdfParams.noise1Frequency = std::clamp(s.asteroidSdfParams.noise1Frequency, 0.0f, 24.0f);
  s.asteroidSdfParams.noise2Frequency = std::clamp(s.asteroidSdfParams.noise2Frequency, 0.0f, 64.0f);
  s.asteroidSdfParams.noise1Amplitude = std::clamp(s.asteroidSdfParams.noise1Amplitude, 0.0f, 2.0f);
  s.asteroidSdfParams.noise2Amplitude = std::clamp(s.asteroidSdfParams.noise2Amplitude, 0.0f, 2.0f);
  s.asteroidSdfParams.noise1Lacunarity = std::clamp(s.asteroidSdfParams.noise1Lacunarity, 1.01f, 4.0f);
  s.asteroidSdfParams.noise2Lacunarity = std::clamp(s.asteroidSdfParams.noise2Lacunarity, 1.01f, 4.0f);
  s.asteroidSdfParams.noise1Gain = std::clamp(s.asteroidSdfParams.noise1Gain, 0.05f, 0.95f);
  s.asteroidSdfParams.noise2Gain = std::clamp(s.asteroidSdfParams.noise2Gain, 0.05f, 0.95f);

  s.asteroidSdfParams.craterCount = std::clamp(s.asteroidSdfParams.craterCount, 0, 48);
  s.asteroidSdfParams.craterRadiusMinDeg = std::clamp(s.asteroidSdfParams.craterRadiusMinDeg, 0.1f, 60.0f);
  s.asteroidSdfParams.craterRadiusMaxDeg = std::clamp(s.asteroidSdfParams.craterRadiusMaxDeg, s.asteroidSdfParams.craterRadiusMinDeg, 85.0f);
  s.asteroidSdfParams.craterDepth = std::clamp(s.asteroidSdfParams.craterDepth, 0.0f, 1.0f);
  s.asteroidSdfParams.craterSmoothK = std::clamp(s.asteroidSdfParams.craterSmoothK, 0.0f, 0.75f);
  s.asteroidSdfParams.cutCount = std::clamp(s.asteroidSdfParams.cutCount, 0, 16);
  s.asteroidSdfParams.cutOffsetMin = std::clamp(s.asteroidSdfParams.cutOffsetMin, 0.0f, 2.5f);
  s.asteroidSdfParams.cutOffsetMax = std::clamp(s.asteroidSdfParams.cutOffsetMax, s.asteroidSdfParams.cutOffsetMin, 3.5f);
  s.asteroidSdfParams.grooveStrength = std::clamp(s.asteroidSdfParams.grooveStrength, 0.0f, 1.0f);
  s.asteroidSdfParams.grooveFrequency = std::clamp(s.asteroidSdfParams.grooveFrequency, 0.0f, 64.0f);
  s.asteroidSdfParams.normalEps = std::clamp(s.asteroidSdfParams.normalEps, 0.0001f, 0.05f);

  s.artifactVariantCount = std::clamp(s.artifactVariantCount, 1, 16);
  s.artifactInstanceCount = std::clamp(s.artifactInstanceCount, 1, 256);
  s.artifactRingRadiusU = std::clamp(s.artifactRingRadiusU, 0.0, 25.0);
  s.artifactRingHeightU = std::clamp(s.artifactRingHeightU, -25.0, 25.0);
  s.artifactScaleU = std::clamp(s.artifactScaleU, 0.05, 6.0);
  s.artifactParams.resolution = std::clamp(s.artifactParams.resolution, 12, 128);
  s.artifactParams.bounds = std::clamp(s.artifactParams.bounds, 0.25f, 5.0f);
  s.artifactParams.normalEps = std::clamp(s.artifactParams.normalEps, 0.0001f, 0.05f);

  s.proceduralSky.intensity = std::clamp(s.proceduralSky.intensity, 0.0f, 10.0f);
  s.proceduralSky.starIntensity = std::clamp(s.proceduralSky.starIntensity, 0.0f, 10.0f);
  s.proceduralSky.starDensity = std::clamp(s.proceduralSky.starDensity, 10.0f, 5000.0f);
  s.proceduralSky.starProbability = std::clamp(s.proceduralSky.starProbability, 0.0f, 1.0f);
  s.proceduralSky.starSize = std::clamp(s.proceduralSky.starSize, 0.05f, 10.0f);
  s.proceduralSky.starTwinkle = std::clamp(s.proceduralSky.starTwinkle, 0.0f, 5.0f);
  s.proceduralSky.nebulaIntensity = std::clamp(s.proceduralSky.nebulaIntensity, 0.0f, 10.0f);
  s.proceduralSky.nebulaFrequency = std::clamp(s.proceduralSky.nebulaFrequency, 0.05f, 12.0f);
  s.proceduralSky.nebulaThreshold = std::clamp(s.proceduralSky.nebulaThreshold, 0.0f, 1.0f);
  s.proceduralSky.nebulaSoftness = std::clamp(s.proceduralSky.nebulaSoftness, 0.0f, 1.0f);
  s.proceduralSky.nebulaBandPower = std::clamp(s.proceduralSky.nebulaBandPower, 0.2f, 6.0f);
  s.proceduralSky.nebulaParallax = std::clamp(s.proceduralSky.nebulaParallax, 0.0f, 1.0f);
  s.proceduralSky.nebulaSteps = std::clamp(s.proceduralSky.nebulaSteps, 2, 64);
  s.proceduralSky.nebulaWorldScale = std::clamp(s.proceduralSky.nebulaWorldScale, 1e-7f, 0.05f);
  s.proceduralSky.nebulaDrift = std::clamp(s.proceduralSky.nebulaDrift, 0.0f, 3.0f);

  outSettings = s;
  return true;
}

} // namespace stellar::ui
