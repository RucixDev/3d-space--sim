#include "test_harness.h"

#include "stellar/sim/LawLedger.h"

#include <algorithm>
#include <cmath>

static bool nearly(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

int test_law_ledger() {
  int failures = 0;

  using namespace stellar;
  using namespace stellar::sim;

  // --- addFine merges amounts and keeps earliest dueDay ----------------------
  {
    LawLedger l;
    l.addFine(1, 100.0, 10.0);
    l.addFine(1, 50.0, 12.0);

    CHECK(nearly(l.fine(1), 150.0));
    CHECK(nearly(l.fineDueDay(1), 10.0));

    // Earlier due date should win.
    l.addFine(1, 25.0, 9.0);
    CHECK(nearly(l.fine(1), 175.0));
    CHECK(nearly(l.fineDueDay(1), 9.0));
  }

  // --- payFine supports partial pay and clears when paid --------------------
  {
    LawLedger l;
    l.addFine(2, 100.0, 10.0);

    double credits = 60.0;
    const double paid1 = l.payFine(2, credits, 100.0);
    CHECK(nearly(paid1, 60.0));
    CHECK(nearly(credits, 0.0));
    CHECK(nearly(l.fine(2), 40.0));

    credits = 100.0;
    const double paid2 = l.payFine(2, credits, 100.0);
    CHECK(nearly(paid2, 40.0));
    CHECK(nearly(credits, 60.0));
    CHECK(nearly(l.fine(2), 0.0));
    CHECK(l.fines().find(2) == l.fines().end());
  }

  // --- overdue fines convert into bounties ---------------------------------
  {
    LawLedger l;
    l.addFine(3, 100.0, 5.0);
    const auto conv = l.processOverdueFines(6.0, 0.25);

    CHECK(conv.size() == 1);
    CHECK(conv[0].factionId == 3);
    CHECK(nearly(conv[0].fineCr, 100.0));
    CHECK(nearly(conv[0].lateFeeCr, 25.0));
    CHECK(nearly(conv[0].bountyAddedCr, 125.0));

    CHECK(nearly(l.fine(3), 0.0));
    CHECK(nearly(l.bounty(3), 125.0));
  }

  // --- Interstellar Factors pays all other-jurisdiction fines ---------------
  {
    LawLedger l;
    l.addFine(1, 50.0, 10.0);
    l.addFine(2, 100.0, 11.0);
    l.addFine(3, 30.0, 12.0);

    double credits = 500.0;
    const auto res = l.payAllOtherFines(2, credits, 0.20);

    CHECK(nearly(res.baseFinesCr, 80.0));
    CHECK(nearly(res.serviceFeeCr, 16.0));
    CHECK(nearly(res.totalCostCr, 96.0));
    CHECK(nearly(res.paidCr, 96.0));
    CHECK(res.clearedFactions == 2);

    CHECK(nearly(credits, 404.0));
    CHECK(nearly(l.fine(1), 0.0));
    CHECK(nearly(l.fine(3), 0.0));
    CHECK(nearly(l.fine(2), 100.0));
  }

  // --- save round-trip aggregates + sorts ----------------------------------
  {
    SaveGame s;
    s.bounties = { {7, 100.0}, {7, 5.0}, {2, 10.0} };
    s.bountyVouchers = { {9, 1.0}, {9, 2.0}, {1, 5.0} };
    s.fines = { {4, 50.0, 10.0}, {4, 25.0, 8.0}, {3, 10.0, 2.0} };

    LawLedger l;
    l.importFromSave(s);

    SaveGame out;
    l.exportToSave(out);

    CHECK(out.bounties.size() == 2);
    CHECK(out.bounties[0].factionId == 2);
    CHECK(nearly(out.bounties[0].bountyCr, 10.0));
    CHECK(out.bounties[1].factionId == 7);
    CHECK(nearly(out.bounties[1].bountyCr, 105.0));

    CHECK(out.bountyVouchers.size() == 2);
    CHECK(out.bountyVouchers[0].factionId == 1);
    CHECK(nearly(out.bountyVouchers[0].bountyCr, 5.0));
    CHECK(out.bountyVouchers[1].factionId == 9);
    CHECK(nearly(out.bountyVouchers[1].bountyCr, 3.0));

    CHECK(out.fines.size() == 2);
    CHECK(out.fines[0].factionId == 3);
    CHECK(nearly(out.fines[0].fineCr, 10.0));
    CHECK(nearly(out.fines[0].dueDay, 2.0));

    CHECK(out.fines[1].factionId == 4);
    CHECK(nearly(out.fines[1].fineCr, 75.0));
    CHECK(nearly(out.fines[1].dueDay, 8.0));
  }

  return failures;
}
