#include "stellar/sim/Mining.h"

#include "test_harness.h"

#include <cmath>

using namespace stellar;

int test_mining() {
  int failures = 0;

  // Efficiency curve basic sanity.
  {
    const double r = 100.0;
    const double e0 = sim::miningEfficiency(0.0, r);
    const double eMid = sim::miningEfficiency(0.5 * r, r);
    const double eMax = sim::miningEfficiency(r, r);

    CHECK(std::abs(e0 - 1.0) < 1e-12);
    CHECK(eMid <= 1.0 + 1e-12);
    CHECK(eMid >= 0.25 - 1e-12);
    CHECK(std::abs(eMax - 0.25) < 1e-12);
  }

  // Prospecting should boost yield by 20% when not clamped by remaining.
  {
    sim::MiningHitInput a{};
    a.universeSeed = 1234u;
    a.asteroidId = 55u;
    a.fieldKind = sim::ResourceFieldKind::OreBelt;
    a.distKm = 0.0;
    a.rangeKm = 100.0;
    a.baseUnitsPerHit = 10.0;
    a.prospected = false;
    a.baseUnits = 200.0;
    a.remainingUnits = 200.0;
    a.fractureAlreadyTriggered = false;

    sim::MiningHitInput b = a;
    b.prospected = true;

    const auto ra = sim::computeMiningHit(a);
    const auto rb = sim::computeMiningHit(b);

    CHECK(std::abs(ra.extractedUnits - 10.0) < 1e-9);
    CHECK(std::abs(rb.extractedUnits - 12.0) < 1e-9);
  }


  // Seam direction should be deterministic and seam-aware mining should modulate yield.
  {
    const core::u64 seed = 1234u;
    const core::u64 id = 55u;
    const auto seam = sim::miningSeamDir(seed, id, sim::ResourceFieldKind::OreBelt);
    CHECK(std::abs(seam.length() - 1.0) < 1e-9);

    sim::MiningHitInput in{};
    in.universeSeed = seed;
    in.asteroidId = id;
    in.fieldKind = sim::ResourceFieldKind::OreBelt;
    in.distKm = 0.0;
    in.rangeKm = 100.0;
    in.baseUnitsPerHit = 10.0;
    in.prospected = true;
    in.baseUnits = 200.0;
    in.remainingUnits = 200.0;

    // Best-case: hit directly on the seam pole.
    in.hitDirUnit = seam;
    const auto rPole = sim::computeMiningHit(in);
    CHECK(rPole.seamMultiplier > 1.05);
    CHECK(rPole.extractedUnits > 12.0);

    // Worst-case: hit near the seam equator (orthogonal direction).
    // Build an orthogonal unit vector robustly.
    math::Vec3d basis{1.0, 0.0, 0.0};
    if (std::abs(math::dot(basis, seam)) > 0.85) basis = {0.0, 1.0, 0.0};
    const math::Vec3d ortho = math::cross(seam, basis).normalized();
    in.hitDirUnit = ortho;
    const auto rEq = sim::computeMiningHit(in);
    CHECK(rEq.seamMultiplier < 0.98);
    CHECK(rEq.extractedUnits < rPole.extractedUnits);
  }

  // Traits should be deterministic and the fracture band should be sane.
  {
    const core::u64 seed = 9999u;

    // Find a volatile asteroid deterministically (bounded scan).
    core::u64 volatileId = 0;
    sim::MiningAsteroidTraits traits{};
    for (core::u64 id = 1; id < 10000; ++id) {
      traits = sim::miningAsteroidTraits(seed, id, sim::ResourceFieldKind::MetalPocket);
      if (traits.volatilePocket) {
        volatileId = id;
        break;
      }
    }

    CHECK(volatileId != 0);
    CHECK(traits.volatilePocket);
    CHECK(traits.fractureFrac >= 0.22 - 1e-12);
    CHECK(traits.fractureFrac <= 0.48 + 1e-12);

    // Determinism: repeated calls match.
    const auto traits2 = sim::miningAsteroidTraits(seed, volatileId, sim::ResourceFieldKind::MetalPocket);
    CHECK(traits2.volatilePocket == traits.volatilePocket);
    CHECK(std::abs(traits2.fractureFrac - traits.fractureFrac) < 1e-12);

    // Fracture triggers when we cross the threshold.
    {
      const double base = 100.0;
      const double thresh = base * traits.fractureFrac;

      sim::MiningHitInput in{};
      in.universeSeed = seed;
      in.asteroidId = volatileId;
      in.fieldKind = sim::ResourceFieldKind::MetalPocket;
      in.distKm = 0.0;
      in.rangeKm = 100.0;
      in.baseUnitsPerHit = 10.0;
      in.prospected = false;
      in.baseUnits = base;
      in.remainingUnits = thresh + 1.0; // just above
      in.fractureAlreadyTriggered = false;

      const auto r0 = sim::computeMiningHit(in);
      CHECK(r0.extractedUnits > 0.0);
      CHECK(r0.fractureTriggered);

      // If we're already below the threshold, it should not trigger.
      sim::MiningHitInput in2 = in;
      in2.remainingUnits = thresh * 0.95;
      const auto r1 = sim::computeMiningHit(in2);
      CHECK(!r1.fractureTriggered);

      // If we've already triggered once, it should not trigger again.
      sim::MiningHitInput in3 = in;
      in3.fractureAlreadyTriggered = true;
      const auto r2 = sim::computeMiningHit(in3);
      CHECK(!r2.fractureTriggered);
    }
  }

  return failures;
}
