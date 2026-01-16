#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"

#include <unordered_map>
#include <vector>

namespace stellar::game {

// Lightweight sensor contact tracking for the game app.
//
// The goal is to provide a "fog-of-war" layer over world entities (ships, signals,
// stations) without forcing that policy into the simulation layer.
//
// NOTE: This module is intentionally self-contained and does not depend on main.cpp
// contact structs. Callers feed it a list of targets + occluders each frame.

struct SensorTarget {
  core::u64 id{0};
  math::Vec3d positionKm{0, 0, 0};

  // 0..1-ish (but values >1 are allowed) multiplier that describes how "loud" a
  // target is to sensors. Larger values extend effective detection range.
  double signature{1.0};

  // Optional physical radius used when callers want size-aware occlusion.
  // (Not required by the current LOS test but kept for future expansion.)
  double radiusKm{0.0};
};

struct SensorOccluder {
  math::Vec3d positionKm{0, 0, 0};
  double radiusKm{1.0};
};

struct SensorContactsSettings {
  // Passive sensor range (km) at signature == 1.
  double baseRangeKm{6500.0};

  // Track "memory" (seconds). Confidence decays toward 0 over this window.
  double trackMemorySec{10.0};

  // Identification accumulates while a contact is detected, decays slowly when not.
  double identifyTimeSec{2.5};
  double identifyDecaySec{18.0};
  double identifyThreshold01{0.65};

  // Line-of-sight attenuation.
  bool occlusionEnabled{true};
  // 0..1 fraction of detection removed when occluded (0 = none, 1 = fully blocked).
  double occlusionPenalty01{0.70};

  // Active ping ("radar pulse") that temporarily boosts detection.
  double pingBoost01{0.35};
  double pingDurationSec{1.25};
  double pingCooldownSec{9.0};
};

struct SensorTrack {
  core::u64 id{0};
  math::Vec3d lastPosKm{0, 0, 0};
  double lastRangeKm{0.0};

  // 0..1 detection confidence (roughly: blip strength).
  double confidence01{0.0};

  // 0..1 identification progress.
  double identify01{0.0};

  bool identified{false};
  bool detectedNow{false};

  double lastSeenDays{0.0};
};

struct SensorView {
  core::u64 id{0};
  math::Vec3d positionKm{0, 0, 0};
  double rangeKm{0.0};
  double confidence01{0.0};
  double identify01{0.0};
  bool identified{false};
  bool detectedNow{false};
  double ageSec{0.0};
};

class SensorContacts {
public:
  void reset();

  // Returns true if a ping was started (i.e., not on cooldown).
  bool tryStartPing(double nowDays, const SensorContactsSettings& settings);

  // External hooks can boost identification (e.g., after a scan completes).
  void boostIdentify(core::u64 id, double amount01 = 1.0);

  // Compute a 0..1 detection value for a single target.
  // This does NOT mutate internal track state.
  double computeDetection01(const math::Vec3d& sensorPosKm,
                            const SensorTarget& target,
                            const std::vector<SensorOccluder>& occluders,
                            double nowDays,
                            const SensorContactsSettings& settings) const;

  // Update internal tracks and return the current views.
  //
  // - nowDays: simulation time in days (same convention as main.cpp).
  // - dtSec:  time step in seconds.
  std::vector<SensorView> update(double nowDays,
                                 double dtSec,
                                 const math::Vec3d& sensorPosKm,
                                 const std::vector<SensorTarget>& targets,
                                 const std::vector<SensorOccluder>& occluders,
                                 const SensorContactsSettings& settings);

  const std::vector<SensorView>& views() const { return views_; }

  bool pingActive(double nowDays) const { return nowDays < pingUntilDays_; }
  double pingCooldownRemainingSec(double nowDays) const;

private:
  std::unordered_map<core::u64, SensorTrack> tracks_;
  std::vector<SensorView> views_;

  double pingUntilDays_{0.0};
  double nextPingAllowedDays_{0.0};
};

} // namespace stellar::game
