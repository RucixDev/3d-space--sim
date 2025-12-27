#include "stellar/core/StableHash.h"
#include "stellar/sim/Signature.h"
#include "stellar/sim/Universe.h"

#include <iostream>

using namespace stellar;

int test_signature() {
  // A small procedural-determinism regression check.
  //
  // If you intentionally change galaxy/system generation, update the expected
  // constants below.
  sim::Universe u(1337);
  const math::Vec3d posLy{0, 0, 0};

  const auto stubs = u.queryNearby(posLy, 120.0, 64);
  if (stubs.empty()) {
    std::cerr << "[test_signature] expected non-empty stub list\n";
    return 1;
  }

  core::StableHash64 h;
  for (const auto& s : stubs) {
    h.addU64(sim::signatureSystemStub(s));
  }

  const core::u64 stubListSig = h.value();
  const core::u64 expectedStubListSig = 6216983436330206766ull;

  const auto& sys0 = u.getSystem(stubs.front().id, &stubs.front());
  const core::u64 sys0Sig = sim::signatureStarSystem(sys0);
  const core::u64 expectedSys0Sig = 5099655091060190414ull;

  int fails = 0;
  if (stubListSig != expectedStubListSig) {
    std::cerr << "[test_signature] stubListSig mismatch. got=" << (unsigned long long)stubListSig
              << " expected=" << (unsigned long long)expectedStubListSig << "\n";
    ++fails;
  }
  if (sys0Sig != expectedSys0Sig) {
    std::cerr << "[test_signature] sys0Sig mismatch. got=" << (unsigned long long)sys0Sig
              << " expected=" << (unsigned long long)expectedSys0Sig << "\n";
    ++fails;
  }

  if (fails == 0) {
    std::cout << "[test_signature] PASS\n";
  }
  return fails;
}
