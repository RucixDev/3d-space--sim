#include "test_harness.h"

#include "stellar/core/Hash.h"
#include "stellar/proc/Hyperlanes.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

namespace {

struct DSU {
  std::vector<int> p;
  std::vector<int> r;

  explicit DSU(int n) : p(n), r(n, 0) {
    for (int i = 0; i < n; ++i) p[i] = i;
  }

  int find(int a) {
    while (p[a] != a) {
      p[a] = p[p[a]];
      a = p[a];
    }
    return a;
  }

  void unite(int a, int b) {
    a = find(a);
    b = find(b);
    if (a == b) return;
    if (r[a] < r[b]) std::swap(a, b);
    p[b] = a;
    if (r[a] == r[b]) r[a] += 1;
  }
};

} // namespace

int test_hyperlanes() {
  int failures = 0;

  constexpr stellar::core::u64 kSeed = 1337ull;
  stellar::sim::Universe u(kSeed);

  const auto stubs = u.queryNearby(stellar::math::Vec3d{0.0, 0.0, 0.0}, 120.0, 256);
  CHECK(stubs.size() > 8);

  stellar::proc::HyperlaneParams p;
  p.maxNeighborDistanceLy = 16.0;
  p.neighborK = 4;
  p.forceConnected = true;
  p.minDegree = 2;
  p.extraEdgeChance = 0.22;
  p.regionCellSizeLy = 900.0;
  p.maxEdges = 0;

  const auto net = stellar::proc::generateHyperlaneNetwork(u.seed(), stubs, p);
  CHECK(!net.edges.empty());

  // Invariants: undirected edges, no duplicates.
  std::vector<std::pair<stellar::sim::SystemId, stellar::sim::SystemId>> pairs;
  pairs.reserve(net.edges.size());

  for (const auto& e : net.edges) {
    CHECK(e.a < e.b);
    CHECK(e.distanceLy > 0.0);
    CHECK(e.bandwidth01 >= 0.0 && e.bandwidth01 <= 1.0);
    CHECK(e.risk01 >= 0.0 && e.risk01 <= 1.0);
    pairs.emplace_back(e.a, e.b);
  }

  std::sort(pairs.begin(), pairs.end());
  CHECK(std::unique(pairs.begin(), pairs.end()) == pairs.end());

  // Connectivity (best-effort): with forceConnected, all nodes in the stub set should be connected.
  std::unordered_map<stellar::sim::SystemId, int> idToIdx;
  idToIdx.reserve(stubs.size() * 2);
  for (int i = 0; i < static_cast<int>(stubs.size()); ++i) {
    idToIdx.emplace(stubs[static_cast<std::size_t>(i)].id, i);
  }

  DSU dsu(static_cast<int>(stubs.size()));
  for (const auto& e : net.edges) {
    const auto ita = idToIdx.find(e.a);
    const auto itb = idToIdx.find(e.b);
    if (ita != idToIdx.end() && itb != idToIdx.end()) {
      dsu.unite(ita->second, itb->second);
    }
  }

  int root0 = dsu.find(0);
  for (int i = 1; i < static_cast<int>(stubs.size()); ++i) {
    CHECK(dsu.find(i) == root0);
  }

  // Determinism: identical inputs must produce identical outputs.
  auto signature = [](const stellar::proc::HyperlaneNetwork& n) {
    stellar::core::u64 sig = 0;
    for (const auto& e : n.edges) {
      sig = stellar::core::hashCombine(sig, e.a);
      sig = stellar::core::hashCombine(sig, e.b);
      sig = stellar::core::hashCombine(sig, (stellar::core::u64)std::llround(e.distanceLy * 1000.0));
      sig = stellar::core::hashCombine(sig, (stellar::core::u64)std::llround(e.bandwidth01 * 1000000.0));
      sig = stellar::core::hashCombine(sig, (stellar::core::u64)std::llround(e.risk01 * 1000000.0));
    }
    return sig;
  };

  const stellar::core::u64 sig1 = signature(net);
  CHECK(sig1 != 0ull);

  const auto net2 = stellar::proc::generateHyperlaneNetwork(u.seed(), stubs, p);
  const stellar::core::u64 sig2 = signature(net2);
  CHECK(sig2 == sig1);

  const auto net3 = stellar::proc::generateHyperlaneNetwork(u.seed() + 1ull, stubs, p);
  const stellar::core::u64 sig3 = signature(net3);
  CHECK(sig3 != sig1);


  return failures;
}
