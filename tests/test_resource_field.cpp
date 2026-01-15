#include "stellar/sim/ResourceField.h"

#include <cmath>
#include <iostream>
#include <unordered_map>

static bool near(double a, double b, double eps = 1e-6) {
  return std::abs(a - b) <= eps;
}

static bool near3(const stellar::math::Vec3d& a, const stellar::math::Vec3d& b, double eps = 1e-6) {
  return near(a.x, b.x, eps) && near(a.y, b.y, eps) && near(a.z, b.z, eps);
}

int test_resource_field() {
  int fails = 0;

  const stellar::core::u64 seed = 123456789ull;
  const stellar::sim::SystemId systemId = 42ull;
  const stellar::math::Vec3d anchorPosKm{1000.0, -2000.0, 500.0};
  const double commsKm = 120000.0;

  const stellar::math::Vec3d planeN{0.0, 1.0, 0.0};

  const auto p0 = stellar::sim::generateResourceFields(seed, systemId, anchorPosKm, commsKm, 3, planeN);
  const auto p1 = stellar::sim::generateResourceFields(seed, systemId, anchorPosKm, commsKm, 3, planeN);

  if (p0.fields.size() != p1.fields.size() || p0.asteroids.size() != p1.asteroids.size() ||
      p0.features.size() != p1.features.size()) {
    std::cerr << "[test_resource_field] determinism: size mismatch\n";
    ++fails;
  }

  // Basic determinism and id tagging.
  const std::size_t nFields = std::min(p0.fields.size(), p1.fields.size());
  for (std::size_t i = 0; i < nFields; ++i) {
    const auto& a = p0.fields[i];
    const auto& b = p1.fields[i];
    if (a.id != b.id || a.kind != b.kind) {
      std::cerr << "[test_resource_field] determinism: field mismatch at i=" << i << "\n";
      ++fails;
      break;
    }

    if (a.layout != b.layout) {
      std::cerr << "[test_resource_field] determinism: field layout mismatch at i=" << i << "\n";
      ++fails;
      break;
    }
    if (!stellar::sim::isDeterministicWorldId(a.id)) {
      std::cerr << "[test_resource_field] expected deterministic id bit set for field id=" << a.id << "\n";
      ++fails;
    }
    if (!near3(a.posKm, b.posKm, 1e-9)) {
      std::cerr << "[test_resource_field] determinism: field position mismatch at i=" << i << "\n";
      ++fails;
    }

    // Basis + params should be deterministic too.
    if (!near3(a.basisX, b.basisX, 1e-9) || !near3(a.basisY, b.basisY, 1e-9) || !near3(a.basisZ, b.basisZ, 1e-9)) {
      std::cerr << "[test_resource_field] determinism: field basis mismatch at i=" << i << "\n";
      ++fails;
    }
    if (!near(a.majorRadiusKm, b.majorRadiusKm, 1e-6) || !near(a.minorRadiusKm, b.minorRadiusKm, 1e-6) ||
        !near(a.arcRad, b.arcRad, 1e-6) || !near(a.arcCenterRad, b.arcCenterRad, 1e-6)) {
      std::cerr << "[test_resource_field] determinism: field layout params mismatch at i=" << i << "\n";
      ++fails;
    }

    // Orthonormal-ish basis.
    const double lx = stellar::math::dot(a.basisX, a.basisX);
    const double ly = stellar::math::dot(a.basisY, a.basisY);
    const double lz = stellar::math::dot(a.basisZ, a.basisZ);
    const double xy = stellar::math::dot(a.basisX, a.basisY);
    const double xz = stellar::math::dot(a.basisX, a.basisZ);
    const double yz = stellar::math::dot(a.basisY, a.basisZ);

    if (!near(lx, 1.0, 1e-6) || !near(ly, 1.0, 1e-6) || !near(lz, 1.0, 1e-6) ||
        !near(xy, 0.0, 1e-5) || !near(xz, 0.0, 1e-5) || !near(yz, 0.0, 1e-5)) {
      std::cerr << "[test_resource_field] basis not orthonormal-ish at i=" << i << "\n";
      ++fails;
    }

    // Preferred plane normal should be honored (within tolerance).
    const double align = stellar::math::dot(a.basisY, planeN.normalized());
    if (align < 0.98) {
      std::cerr << "[test_resource_field] basisY not aligned with preferred plane normal at i=" << i
                << " (align=" << align << ")\n";
      ++fails;
    }
  }

  // Feature determinism.
  const std::size_t nFeat = std::min(p0.features.size(), p1.features.size());
  for (std::size_t i = 0; i < nFeat; ++i) {
    const auto& a = p0.features[i];
    const auto& b = p1.features[i];
    if (a.fieldId != b.fieldId || a.kind != b.kind) {
      std::cerr << "[test_resource_field] determinism: feature mismatch at i=" << i << "\n";
      ++fails;
      break;
    }
    if (!near(a.angleRad, b.angleRad, 1e-12) || !near(a.width, b.width, 1e-12) || !near(a.strength01, b.strength01, 1e-12) ||
        !near(a.param, b.param, 1e-12) || !near3(a.localPos, b.localPos, 1e-12)) {
      std::cerr << "[test_resource_field] determinism: feature params mismatch at i=" << i << "\n";
      ++fails;
      break;
    }
  }

  // Build a lookup for field metadata.
  std::unordered_map<stellar::core::u64, const stellar::sim::ResourceFieldSite*> fieldById;
  fieldById.reserve(p0.fields.size() * 2 + 1);
  for (const auto& f : p0.fields) fieldById[f.id] = &f;

  // Every field should have at least one structural feature.
  for (const auto& f : p0.fields) {
    int n = 0;
    for (const auto& ft : p0.features) {
      if (ft.fieldId == f.id) ++n;
    }
    if (n <= 0) {
      std::cerr << "[test_resource_field] expected features for field id=" << f.id << "\n";
      ++fails;
      break;
    }
  }
  if (fails != 0) return fails;

  // Asteroid ids should be deterministic and belong to known fields.
  std::size_t checked = 0;
  for (std::size_t i = 0; i < p0.asteroids.size() && checked < 64; ++i, ++checked) {
    const auto& a = p0.asteroids[i];
    if (!stellar::sim::isDeterministicWorldId(a.id)) {
      std::cerr << "[test_resource_field] expected deterministic id bit set for asteroid id=" << a.id << "\n";
      ++fails;
      break;
    }
    bool foundField = false;
    for (const auto& f : p0.fields) {
      if (f.id == a.fieldId) { foundField = true; break; }
    }
    if (!foundField) {
      std::cerr << "[test_resource_field] asteroid references unknown field id=" << a.fieldId << "\n";
      ++fails;
      break;
    }

    // Layout bounds sanity for a sample of rocks.
    const auto* f = fieldById[a.fieldId];
    if (f) {
      const stellar::math::Vec3d d = a.posKm - f->posKm;
      const double x = stellar::math::dot(d, f->basisX);
      const double y = stellar::math::dot(d, f->basisY);
      const double z = stellar::math::dot(d, f->basisZ);

      const double major = std::max(1.0, f->majorRadiusKm);
      const double minor = std::max(1.0, f->minorRadiusKm);

      if (f->layout == stellar::sim::ResourceFieldLayout::Cluster) {
        const double nx = x / major;
        const double ny = y / minor;
        const double nz = z / major;
        const double ell = nx * nx + ny * ny + nz * nz;
        if (ell > 1.55) {
          std::cerr << "[test_resource_field] asteroid outside cluster bounds\n";
          ++fails;
          break;
        }
      } else if (f->layout == stellar::sim::ResourceFieldLayout::Sheet) {
        const double r = std::sqrt(x * x + z * z);
        if (r > major * 1.08 || std::abs(y) > minor * 1.25) {
          std::cerr << "[test_resource_field] asteroid outside sheet bounds\n";
          ++fails;
          break;
        }
      } else if (f->layout == stellar::sim::ResourceFieldLayout::Torus) {
        const double r = std::sqrt(x * x + z * z);
        const double dr = std::abs(r - major);
        const double tube = std::sqrt(dr * dr + y * y);
        if (tube > minor * 1.30) {
          std::cerr << "[test_resource_field] asteroid outside torus bounds\n";
          ++fails;
          break;
        }
      }

      // Density should be well-formed and match the density function at the asteroid position.
      if (a.density01 < -1e-9 || a.density01 > 1.0 + 1e-9) {
        std::cerr << "[test_resource_field] asteroid density01 out of range\n";
        ++fails;
        break;
      }
      const double d01 = stellar::sim::resourceFieldDensity01(*f, p0.features, a.posKm);
      if (!near(a.density01, d01, 1e-9)) {
        std::cerr << "[test_resource_field] asteroid density01 mismatch vs density function\n";
        ++fails;
        break;
      }
    }

    if (a.baseUnits <= 0.0) {
      std::cerr << "[test_resource_field] asteroid baseUnits should be > 0\n";
      ++fails;
      break;
    }
    if (a.radiusKm <= 0.0) {
      std::cerr << "[test_resource_field] asteroid radiusKm should be > 0\n";
      ++fails;
      break;
    }

  }

  // Ensure asteroids don't overlap within each field (simple check for a small sample).
  for (const auto& f : p0.fields) {
    int checkedPairs = 0;
    for (std::size_t i = 0; i < p0.asteroids.size(); ++i) {
      const auto& a = p0.asteroids[i];
      if (a.fieldId != f.id) continue;
      for (std::size_t j = i + 1; j < p0.asteroids.size(); ++j) {
        const auto& b = p0.asteroids[j];
        if (b.fieldId != f.id) continue;
        const double d = (a.posKm - b.posKm).length();
        const double min = (a.radiusKm + b.radiusKm) * 0.75; // should not intersect
        if (d < min) {
          std::cerr << "[test_resource_field] asteroid overlap detected\n";
          ++fails;
          checkedPairs = 999999;
          break;
        }
        if (++checkedPairs > 256) break;
      }
      if (checkedPairs > 256) break;
    }
    if (fails != 0) break;
  }

  if (fails == 0) std::cout << "[test_resource_field] pass\n";
  return fails;
}
