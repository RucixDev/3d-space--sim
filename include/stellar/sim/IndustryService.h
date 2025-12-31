#pragma once

#include "stellar/sim/Industry.h"
#include "stellar/sim/System.h"
#include "stellar/sim/Warehouse.h"

#include <array>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Industry services (headless)
// -----------------------------------------------------------------------------
//
// apps/stellar_game historically implemented industry order submission / queueing
// and claiming directly in the giant main.cpp. This header extracts the
// deterministic *business logic* into the core library so:
//   - the SDL prototype and headless tools/tests share the same rules
//   - save/load semantics remain stable
//   - UI code focuses on presentation
//
// Key behavior differences vs the early prototype:
//   - Orders are queued per-station (a simple FIFO using readyDay as "queue end").
//     This makes processing feel like an actual station service instead of
//     infinitely parallel instant factories.

struct IndustrySubmitResult {
  bool ok{false};
  // Short machine-friendly failure reason (nullptr on success).
  // Common values:
  //  - "invalid"        (bad inputs)
  //  - "no_service"     (station can't run this recipe)
  //  - "no_batches"     (batches <= 0)
  //  - "need_inputs"    (not enough cargo inputs)
  //  - "need_credits"   (not enough credits for service fee)
  const char* reason{nullptr};

  // Filled on success.
  IndustryOrder order{};
  IndustryQuote quote{};
  double creditsPaid{0.0};
};

// Submit a new industry order.
//
// Mutates:
//  - shipCargo (consumes inputs)
//  - credits (pays the station service fee)
//  - orders / nextOrderId (adds a new order)
IndustrySubmitResult submitIndustryOrder(core::u64& nextOrderId,
                                        std::vector<IndustryOrder>& orders,
                                        std::array<double, econ::kCommodityCount>& shipCargo,
                                        double& credits,
                                        const Station& station,
                                        double timeDays,
                                        double rep,
                                        const IndustryRecipeDef& recipe,
                                        int batches,
                                        double effectiveStationFeeRate);

struct IndustryClaimResult {
  bool ok{false};
  const char* reason{nullptr};

  // Output movement.
  double unitsMoved{0.0};

  // True if the order became fully claimed as a result of this operation.
  bool completed{false};
};

// Claim ready output into the player's cargo hold.
//
// If cargo is full, this will claim as much as fits (partial claim) and leave
// the remainder in the order.
IndustryClaimResult claimIndustryOrderToCargo(IndustryOrder& order,
                                             const Station& station,
                                             double timeDays,
                                             std::array<double, econ::kCommodityCount>& shipCargo,
                                             double cargoCapacityKg);

// Move all remaining output into station warehouse storage.
//
// This is a convenience for "cargo full" situations; it accrues storage fees
// first so fee projection remains consistent.
IndustryClaimResult moveIndustryOrderOutputToWarehouse(IndustryOrder& order,
                                                       const Station& station,
                                                       double timeDays,
                                                       double rep,
                                                       std::vector<StationStorage>& stationStorage);

// Remove orders that are fully claimed.
void pruneClaimedIndustryOrders(std::vector<IndustryOrder>& orders);

} // namespace stellar::sim
