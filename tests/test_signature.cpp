#include "stellar/sim/Signature.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <iostream>

using namespace stellar;

int test_signature() {
  // A small procedural-determinism regression check.
  //
  // If you intentionally change galaxy/system generation, update the expected
  // constants below.
  sim::Universe u(1337);
  const math::Vec3d posLy{0, 0, 0};

  auto stubs = u.queryNearby(posLy, 120.0, 64);
  if (stubs.empty()) {
    std::cerr << "[test_signature] expected non-empty stub list\n";
    return 1;
  }

  // Ensure the signature is based on content, not incidental query ordering.
  std::sort(stubs.begin(), stubs.end(),
            [](const sim::SystemStub& a, const sim::SystemStub& b) { return a.id < b.id; });

  const core::u64 stubListSig = sim::signatureSystemStubList(stubs);
  const core::u64 expectedStubListSig = 7978128231034067079ull;

  const auto& sys0 = u.getSystem(stubs.front().id, &stubs.front());
  const core::u64 sys0Sig = sim::signatureStarSystem(sys0, u.factions());
  const core::u64 expectedSys0Sig = 7801071842472395155ull;

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
