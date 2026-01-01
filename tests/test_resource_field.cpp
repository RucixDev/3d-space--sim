#include "stellar/sim/ResourceField.h"

#include <cmath>
#include <iostream>

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

  const auto p0 = stellar::sim::generateResourceFields(seed, systemId, anchorPosKm, commsKm, 3);
  const auto p1 = stellar::sim::generateResourceFields(seed, systemId, anchorPosKm, commsKm, 3);

  if (p0.fields.size() != p1.fields.size() || p0.asteroids.size() != p1.asteroids.size()) {
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
    if (!stellar::sim::isDeterministicWorldId(a.id)) {
      std::cerr << "[test_resource_field] expected deterministic id bit set for field id=" << a.id << "\n";
      ++fails;
    }
    if (!near3(a.posKm, b.posKm, 1e-9)) {
      std::cerr << "[test_resource_field] determinism: field position mismatch at i=" << i << "\n";
      ++fails;
    }
  }

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

  if (fails == 0) std::cout << "[test_resource_field] pass\n";
  return fails;
}
