#pragma once

#include "stellar/sim/SaveGame.h"
#include "stellar/sim/SystemEvents.h"
#include "stellar/sim/SystemSecurityDynamics.h"

#include <span>
#include <string>
#include <vector>

namespace stellar::sim {

class Universe;
struct StarSystem;
struct Station;

// -----------------------------------------------------------------------------
// Mission Briefings (procedural narrative + risk model)
// -----------------------------------------------------------------------------
//
// The game prototype has a rich mission system, but missions in the save file
// are intentionally compact (no long strings). This module generates an
// immersive contract briefing text deterministically from:
//   - universe seed
//   - mission fields (id/type/faction/etc.)
//   - system conditions (security/piracy/law)
//
// The sim layer does not interpret any markup. If useMarkup is true, the
// generated strings may contain optional ui::textfx tags so the renderer/UI can
// add flair (gradient titles, subtle typewriter, etc.). Unknown tags remain
// literal text.

struct MissionBriefingParams {
  SystemSecurityDynamicsParams dynamicsParams{};
  SystemEventParams eventParams{};

  // If true, the risk model uses per-system security delta states when provided.
  bool applySecurityDeltas{true};

  // If true, embed optional ui::textfx markup tags in generated text.
  bool useMarkup{true};

  // If true, include detailed risk breakdowns (component values, customs/fine
  // estimates, etc.). If false, the briefing stays more "in-universe" and less
  // numeric.
  bool includeRiskHints{true};

  // If true, include a lightweight cue about the player's standing with the
  // issuing faction (useful for smuggling/black-market flavored missions).
  bool includeReputationCues{false};
};

struct MissionRisk {
  // Straight-line distance from origin system to the mission's *next objective*
  // system. (Multi-leg contracts use the via stop when leg==0.)
  double distanceLy{0.0};

  // Overall summary risk in [0,1]. Used for tier labels in UI.
  double overall01{0.0};

  // Component risks (all in [0,1]).
  double danger01{0.0};   // piracy + contestedness + travel/delivery profile
  double lawRisk01{0.0};  // customs/enforcement (mostly relevant for smuggling)

  // Raw system condition signals (effective).
  double security01{0.5};
  double piracy01{0.5};
  double traffic01{0.5};
  double contest01{0.0};

  // Smuggling-only estimates.
  bool contrabandIllegalAtDestination{false};
  double blackMarketAccess01{0.0};
  double stingChance01{0.0};
  double expectedFineCr{0.0};
};

struct MissionBriefing {
  // Short, stable identifiers useful in UI/logging.
  std::string contractCode;
  std::string contactName;

  // UI-facing strings. If params.useMarkup is true, these may contain optional
  // markup tags compatible with stellar::ui::textfx.
  std::string titleMarkup;
  std::string synopsisMarkup;
  std::vector<std::string> bulletsMarkup;

  MissionRisk risk{};
};

// Convenience label for UI.
const char* riskTierName(double risk01);

// Compute a compact mission risk model without generating narrative text.
MissionRisk computeMissionRisk(Universe& universe,
                               const StarSystem& originSystem,
                               const Station& originStation,
                               double timeDays,
                               double playerRepWithIssuerFaction,
                               const Mission& mission,
                               std::span<const SystemSecurityDeltaState> securityDeltas = {},
                               const MissionBriefingParams& params = {});

// Generate a deterministic mission contract briefing text + risk model.
MissionBriefing generateMissionBriefing(Universe& universe,
                                        const StarSystem& originSystem,
                                        const Station& originStation,
                                        double timeDays,
                                        double playerRepWithIssuerFaction,
                                        const Mission& mission,
                                        std::span<const SystemSecurityDeltaState> securityDeltas = {},
                                        const MissionBriefingParams& params = {});

} // namespace stellar::sim
