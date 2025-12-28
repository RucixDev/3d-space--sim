#include "stellar/render/ParticleSystem.h"

#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>

namespace stellar::render {

static math::Vec3d orthogonal(const math::Vec3d& n) {
  // Pick a vector guaranteed not parallel to n, then cross.
  const math::Vec3d a = (std::abs(n.z) < 0.9) ? math::Vec3d{0,0,1} : math::Vec3d{0,1,0};
  const math::Vec3d t = math::cross(n, a);
  const double l2 = t.lengthSq();
  if (l2 < 1e-12) return {1,0,0};
  return t / std::sqrt(l2);
}

static math::Vec3d sampleCone(core::SplitMix64& rng,
                              const math::Vec3d& axis,
                              double angleRad) {
  // Uniform-ish sample inside a cone around axis.
  // angleRad is the max deviation.
  if (angleRad <= 1e-6) return axis;

  const math::Vec3d w = axis.normalized();
  const math::Vec3d u = orthogonal(w);
  const math::Vec3d v = math::cross(w, u).normalized();

  const double u1 = rng.nextUnit();
  const double u2 = rng.nextUnit();

  const double cosMax = std::cos(angleRad);
  const double cosTheta = (1.0 - u1) + u1 * cosMax; // lerp(1, cosMax, u1)
  const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
  const double phi = 2.0 * math::kPi * u2;

  const double x = std::cos(phi) * sinTheta;
  const double y = std::sin(phi) * sinTheta;
  const double z = cosTheta;

  // Local (x,y,z) where z aligns with axis.
  return (u * x + v * y + w * z).normalized();
}

static math::Vec3d sampleHemisphere(core::SplitMix64& rng,
                                    const math::Vec3d& normal) {
  // Cosine-ish hemisphere.
  const math::Vec3d n = normal.normalized();
  const math::Vec3d u = orthogonal(n);
  const math::Vec3d v = math::cross(n, u).normalized();

  const double u1 = rng.nextUnit();
  const double u2 = rng.nextUnit();

  const double r = std::sqrt(u1);
  const double theta = 2.0 * math::kPi * u2;
  const double x = r * std::cos(theta);
  const double y = r * std::sin(theta);
  const double z = std::sqrt(std::max(0.0, 1.0 - u1));

  return (u * x + v * y + n * z).normalized();
}

void ParticleSystem::push(const Particle& p) {
  if (particles_.size() >= maxParticles_) {
    // Drop the oldest particles first (simple and stable).
    const std::size_t drop = std::max<std::size_t>(1, maxParticles_ / 8);
    if (drop >= particles_.size()) {
      particles_.clear();
    } else {
      particles_.erase(particles_.begin(), particles_.begin() + (std::ptrdiff_t)drop);
    }
  }
  particles_.push_back(p);
}

void ParticleSystem::update(double dtSeconds) {
  if (dtSeconds <= 0.0) return;

  for (auto& p : particles_) {
    p.life -= (float)dtSeconds;
    if (p.life <= 0.0f) continue;

    // Drag is implemented as a simple exponential decay.
    if (p.drag > 0.0f) {
      const double k = std::exp(-(double)p.drag * dtSeconds);
      p.velU = p.velU * k;
    }
    p.posU += p.velU * dtSeconds;
  }

  particles_.erase(std::remove_if(particles_.begin(), particles_.end(),
                                  [](const Particle& p) { return p.life <= 0.0f; }),
                   particles_.end());
}

void ParticleSystem::buildPoints(std::vector<PointVertex>& out) const {
  out.clear();
  out.reserve(particles_.size());

  for (const auto& p : particles_) {
    const float t = (p.maxLife > 1e-6f) ? std::clamp(p.life / p.maxLife, 0.0f, 1.0f) : 0.0f;

    // t=1 at spawn, t=0 at death.
    const float age = 1.0f - t;

    auto smoothstep = [](float e0, float e1, float x) {
      const float u = std::clamp((x - e0) / (e1 - e0), 0.0f, 1.0f);
      return u * u * (3.0f - 2.0f * u);
    };

    // Fast fade-in, then hold, then fade-out near end.
    const float fadeIn = smoothstep(0.0f, 0.08f, age);
    const float fadeOut = smoothstep(0.0f, 0.22f, t);
    const float a = p.baseAlpha * fadeIn * fadeOut;

    // Slight shrink as particles die.
    const float size = std::max(0.5f, p.sizePx * (0.65f + 0.35f * t));

    PointVertex v;
    v.px = (float)p.posU.x;
    v.py = (float)p.posU.y;
    v.pz = (float)p.posU.z;
    v.cr = p.r;
    v.cg = p.g;
    v.cb = p.b;
    v.a = a;
    v.size = size;
    out.push_back(v);
  }
}

void ParticleSystem::emitThruster(const math::Vec3d& posU,
                                  const math::Vec3d& dirWorld,
                                  double intensity,
                                  double dtSeconds,
                                  bool boost) {
  intensity = std::clamp(intensity, 0.0, 1.0);
  if (intensity <= 1e-4) return;
  if (dtSeconds <= 0.0) return;

  const math::Vec3d axis = dirWorld.lengthSq() > 1e-12 ? dirWorld.normalized() : math::Vec3d{0,0,-1};

  // Particles per second.
  double rate = 260.0 * intensity;
  if (boost) rate *= 1.55;

  emitCarry_ += rate * dtSeconds;
  const int count = (int)std::floor(emitCarry_);
  emitCarry_ -= (double)count;
  if (count <= 0) return;

  for (int i = 0; i < count; ++i) {
    Particle p;

    // Spawn in a small disk so it doesn't look like a single pixel.
    const double jitterR = rng_.range(0.0, 0.12) * std::sqrt(rng_.nextUnit());
    const double jitterA = rng_.range(0.0, 2.0 * math::kPi);
    const math::Vec3d t = orthogonal(axis);
    const math::Vec3d b = math::cross(axis, t).normalized();
    const math::Vec3d jitter = t * (std::cos(jitterA) * jitterR) + b * (std::sin(jitterA) * jitterR);

    p.posU = posU + jitter;

    const math::Vec3d d = sampleCone(rng_, axis, (boost ? 0.22 : 0.30));
    const double speed = (boost ? rng_.range(2.2, 4.2) : rng_.range(1.4, 3.2)) * intensity;
    p.velU = d * speed;

    // Color: orange-white to blue-white based on intensity.
    const float tCol = (float)std::clamp(intensity, 0.0, 1.0);
    p.r = 1.00f;
    p.g = 0.55f + 0.35f * tCol;
    p.b = 0.25f + 0.65f * tCol;

    p.baseAlpha = (float)std::clamp(0.35 + 0.65 * intensity, 0.0, 1.0);
    p.sizePx = (float)(boost ? rng_.range(10.0, 18.0) : rng_.range(8.0, 14.0));
    p.life = p.maxLife = (float)(boost ? rng_.range(0.20, 0.34) : rng_.range(0.18, 0.28));
    p.drag = (float)rng_.range(1.2, 2.8);

    push(p);
  }
}

void ParticleSystem::spawnExplosion(const math::Vec3d& posU,
                                    const math::Vec3d& inheritVelU,
                                    double energy,
                                    int count) {
  energy = std::max(0.0, energy);
  if (energy <= 1e-6) return;

  if (count < 0) {
    count = (int)std::clamp(energy * 220.0, 20.0, 900.0);
  }

  for (int i = 0; i < count; ++i) {
    Particle p;

    // Random unit vector on sphere.
    const double z = rng_.range(-1.0, 1.0);
    const double a = rng_.range(0.0, 2.0 * math::kPi);
    const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
    const math::Vec3d dir{r * std::cos(a), r * std::sin(a), z};

    const double shell = rng_.range(0.0, 0.35) * std::pow(rng_.nextUnit(), 0.35);
    p.posU = posU + dir * shell;

    const double sp = rng_.range(1.2, 6.0) * std::pow(rng_.nextUnit(), 0.25) * std::sqrt(energy);
    p.velU = inheritVelU + dir * sp;

    // Explosion palette: hot core + embers.
    const double hot = rng_.nextUnit();
    p.r = 1.0f;
    p.g = (float)std::clamp(0.25 + 0.75 * hot, 0.0, 1.0);
    p.b = (float)std::clamp(0.08 + 0.45 * hot, 0.0, 1.0);

    p.baseAlpha = (float)std::clamp(0.55 + 0.45 * std::sqrt(energy), 0.0, 1.0);
    p.sizePx = (float)rng_.range(10.0, 26.0);
    p.life = p.maxLife = (float)rng_.range(0.35, 0.90);
    p.drag = (float)rng_.range(0.6, 1.8);

    push(p);
  }
}

void ParticleSystem::spawnSparks(const math::Vec3d& posU,
                                 const math::Vec3d& normalWorld,
                                 const math::Vec3d& inheritVelU,
                                 double energy,
                                 int count) {
  energy = std::max(0.0, energy);
  if (energy <= 1e-6) return;

  if (count < 0) {
    count = (int)std::clamp(energy * 120.0, 8.0, 220.0);
  }

  const math::Vec3d n = (normalWorld.lengthSq() > 1e-12) ? normalWorld.normalized() : math::Vec3d{0,1,0};

  for (int i = 0; i < count; ++i) {
    Particle p;
    p.posU = posU + sampleCone(rng_, n, 1.15) * rng_.range(0.0, 0.12);

    const math::Vec3d dir = sampleHemisphere(rng_, n);
    const double sp = rng_.range(1.8, 5.0) * std::pow(rng_.nextUnit(), 0.25) * std::sqrt(energy);
    p.velU = inheritVelU + dir * sp;

    // White-ish sparks.
    const float w = (float)rng_.range(0.85, 1.0);
    p.r = w;
    p.g = w;
    p.b = (float)std::clamp(w * rng_.range(0.65, 1.0), 0.0, 1.0);

    p.baseAlpha = (float)std::clamp(0.25 + 0.75 * std::sqrt(energy), 0.0, 1.0);
    p.sizePx = (float)rng_.range(5.0, 11.0);
    p.life = p.maxLife = (float)rng_.range(0.12, 0.35);
    p.drag = (float)rng_.range(1.6, 3.4);

    push(p);
  }
}

} // namespace stellar::render
