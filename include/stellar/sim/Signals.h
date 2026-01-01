#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Distress.h"
#include "stellar/sim/ResourceField.h"
#include "stellar/sim/System.h"

#include <vector>

namespace stellar::sim {

struct Mission; // fwd

enum class SignalKind : core::u8 {
  ResourceField = 0,
  Derelict = 1,
  Distress = 2,
  MissionSalvage = 3,
};

const char* signalKindName(SignalKind k);

struct SignalSite {
  core::u64 id{0};
  SignalKind kind{SignalKind::Derelict};
  math::Vec3d posKm{0,0,0};

  // When expireDay > 0 and timeDays >= expireDay, the site should be treated as expired.
  // Persistent sites can leave this at 0.
  double expireDay{0.0};

  // Optional: resolved/consumed by the player. This is typically tracked by SaveGame.resolvedSignalIds.
  bool resolved{false};

  // Optional payloads (only valid for specific kinds)
  ResourceFieldKind fieldKind{ResourceFieldKind::OreBelt};

  bool hasDistressPlan{false};
  DistressPlan distress{};

  // For mission-related signals, helps correlate UI -> mission list.
  core::u64 missionId{0};
};

struct SignalGenParams {
  // Persistent mining sites.
  int resourceFieldCount{1};

  // A single "daily" derelict per system (expires at the next integer day).
  bool includeDailyDerelict{true};

  // Distress calls: deterministic encounters with a short TTL.
  bool includeDistress{true};
  int distressPerDay{1};
  double distressTtlDays{1.0};
};

struct SystemSignalPlan {
  // Mining sites + asteroids.
  ResourceFieldPlan resourceFields{};

  // Flat list of all sites (includes resource fields as SignalKind::ResourceField entries).
  std::vector<SignalSite> sites{};
};

// Deterministic in-system "signal" generation (resource fields, derelicts, distress, and mission salvage sites).
SystemSignalPlan generateSystemSignals(core::u64 universeSeed,
                                      const StarSystem& system,
                                      double timeDays,
                                      const std::vector<Mission>& activeMissions,
                                      const std::vector<core::u64>& resolvedSignalIds,
                                      const SignalGenParams& params = {});

// O(n) helper (resolvedSignalIds is usually small). Returns true if signalId is contained.
bool isSignalResolved(const std::vector<core::u64>& resolvedSignalIds, core::u64 signalId);

} // namespace stellar::sim
