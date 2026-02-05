#pragma once

#include "stellar/math/Frustum.h"
#include "stellar/math/Geometry.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace stellar::math {

// BVH build mode.
//
// Median:
//   - classic recursive median split along largest centroid extent.
//   - good build quality, deterministic.
//
// Morton:
//   - sort primitives by Morton (Z-order) code of their centroids.
//   - build a radix-tree style BVH by splitting on the highest differing Morton bit.
//   - deterministic and tends to preserve spatial locality.
//
// Notes:
// - This is a CPU-side, single-threaded approximation of LBVH ideas.
// - We intentionally keep this small, header-only, and dependency-free.
enum class BvhBuildMode : std::uint8_t {
  Median = 0,
  Morton = 1,
};

// A lightweight, header-only BVH (bounding volume hierarchy) for AABB primitives.
//
// Goals:
// - Fast broad-phase queries for gameplay + rendering.
// - Deterministic builds (stable split rule).
// - No dynamic allocations during queries (except a small stack vector).
//
// This BVH is intentionally simple:
// - Median split along the largest centroid extent, or a Morton-code build mode.
// - Leaves store a small range of primitives.
// - Queries are conservative: frustum culling uses AABB-frustum intersection.

struct BvhItem3d {
  Aabb3d bounds{};
  int id{0};
};

class Bvh3d {
public:
  Bvh3d() = default;

  void clear() {
    items_.clear();
    prims_.clear();
    nodes_.clear();
    leafSize_ = 4;
    buildMode_ = BvhBuildMode::Median;
  }

  bool empty() const { return nodes_.empty(); }
  std::size_t itemCount() const { return items_.size(); }
  std::size_t nodeCount() const { return nodes_.size(); }
  BvhBuildMode buildMode() const { return buildMode_; }

  // Build a BVH over the provided items.
  //
  // Notes:
  // - The BVH copies `items` (callers can std::move).
  // - Non-finite bounds are ignored.
  // - leafSize is clamped to [1, 64].
  // - mode selects the build strategy.
  void build(std::vector<BvhItem3d> items,
             std::size_t leafSize = 4,
             BvhBuildMode mode = BvhBuildMode::Median) {
    clear();

    buildMode_ = mode;
    leafSize_ = std::clamp<std::size_t>(leafSize, 1, 64);

    // Filter invalid bounds, and precompute centroids.
    items_.reserve(items.size());
    for (auto& it : items) {
      if (!it.bounds.isFinite()) continue;
      Item in{};
      in.bounds = it.bounds;
      in.center = it.bounds.center();
      in.id = it.id;
      items_.push_back(in);
    }

    if (items_.empty()) return;

    prims_.resize(items_.size());
    for (std::size_t i = 0; i < prims_.size(); ++i) prims_[i] = static_cast<std::uint32_t>(i);

    // A binary BVH over N primitives has at most (2N-1) nodes.
    nodes_.reserve(items_.size() * 2);

    if (buildMode_ == BvhBuildMode::Morton) {
      computeMortonCodes_();

      // Deterministic sort by (mortonCode, originalIndex).
      std::stable_sort(prims_.begin(), prims_.end(), [&](std::uint32_t a, std::uint32_t b) {
        const std::uint32_t ca = items_[a].morton;
        const std::uint32_t cb = items_[b].morton;
        if (ca != cb) return ca < cb;
        return a < b;
      });

      buildNodeMorton_(0, static_cast<std::uint32_t>(items_.size()));
    } else {
      buildNodeMedian_(0, static_cast<std::uint32_t>(items_.size()));
    }
  }

  // AABB overlap query.
  // Calls `onHit(id)` for each intersecting primitive.
  template <class F>
  void queryAabb(const Aabb3d& q, F&& onHit) const {
    if (nodes_.empty() || !q.isFinite()) return;

    std::vector<int> stack;
    stack.reserve(64);
    stack.push_back(0);

    while (!stack.empty()) {
      const int ni = stack.back();
      stack.pop_back();
      const Node& n = nodes_[static_cast<std::size_t>(ni)];
      if (!n.bounds.intersectsAabb(q)) continue;

      if (n.isLeaf()) {
        for (std::uint32_t k = 0; k < n.count; ++k) {
          const Item& it = items_[prims_[static_cast<std::size_t>(n.start + k)]];
          if (it.bounds.intersectsAabb(q)) {
            onHit(it.id);
          }
        }
      } else {
        if (n.left >= 0) stack.push_back(n.left);
        if (n.right >= 0) stack.push_back(n.right);
      }
    }
  }

  // Frustum query.
  // Calls `onHit(id)` for each primitive whose AABB intersects the frustum.
  template <class F>
  void queryFrustum(const Frustumd& fr, F&& onHit) const {
    if (nodes_.empty() || !fr.isFinite()) return;

    std::vector<int> stack;
    stack.reserve(64);
    stack.push_back(0);

    while (!stack.empty()) {
      const int ni = stack.back();
      stack.pop_back();
      const Node& n = nodes_[static_cast<std::size_t>(ni)];
      if (!fr.intersectsAabb(n.bounds)) continue;

      if (n.isLeaf()) {
        for (std::uint32_t k = 0; k < n.count; ++k) {
          const Item& it = items_[prims_[static_cast<std::size_t>(n.start + k)]];
          if (fr.intersectsAabb(it.bounds)) {
            onHit(it.id);
          }
        }
      } else {
        if (n.left >= 0) stack.push_back(n.left);
        if (n.right >= 0) stack.push_back(n.right);
      }
    }
  }

  // Segment query.
  //
  // Calls `onHit(id)` for each primitive whose AABB intersects the segment [a,b].
  //
  // Return/early-out behavior:
  // - If `onHit` returns void, all hits are reported.
  // - If `onHit` returns bool, returning false stops traversal early.
  //
  // Returns true if the segment intersects at least one primitive AABB.
  template <class F>
  bool querySegment(const Vec3d& a, const Vec3d& b, F&& onHit) const {
    if (nodes_.empty()) return false;
    if (!stellar::math::isFinite(a) || !stellar::math::isFinite(b)) return false;

    std::vector<int> stack;
    stack.reserve(64);
    stack.push_back(0);

    bool anyHit = false;

    while (!stack.empty()) {
      const int ni = stack.back();
      stack.pop_back();

      const Node& n = nodes_[static_cast<std::size_t>(ni)];
      double t0 = 0.0, t1 = 0.0;
      if (!n.bounds.segmentIntersectionT(a, b, t0, t1)) continue;

      if (n.isLeaf()) {
        for (std::uint32_t k = 0; k < n.count; ++k) {
          const Item& it = items_[prims_[static_cast<std::size_t>(n.start + k)]];
          double tp0 = 0.0, tp1 = 0.0;
          if (!it.bounds.segmentIntersectionT(a, b, tp0, tp1)) continue;

          anyHit = true;
          if constexpr (std::is_same_v<std::invoke_result_t<F, int>, void>) {
            onHit(it.id);
          } else {
            if (!onHit(it.id)) return true;
          }
        }
      } else {
        if (n.left >= 0) stack.push_back(n.left);
        if (n.right >= 0) stack.push_back(n.right);
      }
    }

    return anyHit;
  }

  // Raycast against primitive AABBs.
  //
  // Returns true when a hit is found, writing:
  //   outT  - ray entry distance (>= 0)
  //   outId - primitive id
  //
  // The direction does not need to be normalized.
  bool raycast(const Vec3d& origin, const Vec3d& dir, double* outT, int* outId) const {
    if (nodes_.empty()) return false;
    if (!stellar::math::isFinite(origin) || !stellar::math::isFinite(dir)) return false;

    const double dirLenSq = dir.lengthSq();
    if (!(dirLenSq > 1e-24)) return false;

    const Vec3d dUnit = dir * (1.0 / std::sqrt(dirLenSq));
    Ray r{origin, dUnit};

    double bestT = std::numeric_limits<double>::infinity();
    int bestId = -1;

    std::vector<int> stack;
    stack.reserve(64);
    stack.push_back(0);

    while (!stack.empty()) {
      const int ni = stack.back();
      stack.pop_back();

      const Node& n = nodes_[static_cast<std::size_t>(ni)];
      double tEnter = 0.0, tExit = 0.0;
      if (!rayAabb_(n.bounds, r, tEnter, tExit)) continue;
      if (tEnter > bestT) continue;

      if (n.isLeaf()) {
        for (std::uint32_t k = 0; k < n.count; ++k) {
          const Item& it = items_[prims_[static_cast<std::size_t>(n.start + k)]];
          double tp0 = 0.0, tp1 = 0.0;
          if (!rayAabb_(it.bounds, r, tp0, tp1)) continue;
          if (tp0 < bestT) {
            bestT = tp0;
            bestId = it.id;
          }
        }
      } else {
        // Traverse near-first by sorting children by entry distance.
        struct ChildHit {
          int idx{-1};
          double t{0.0};
        } hit[2];
        int hitCount = 0;

        if (n.left >= 0) {
          double tl0 = 0.0, tl1 = 0.0;
          if (rayAabb_(nodes_[static_cast<std::size_t>(n.left)].bounds, r, tl0, tl1) && tl0 <= bestT) {
            hit[hitCount++] = ChildHit{n.left, tl0};
          }
        }
        if (n.right >= 0) {
          double tr0 = 0.0, tr1 = 0.0;
          if (rayAabb_(nodes_[static_cast<std::size_t>(n.right)].bounds, r, tr0, tr1) && tr0 <= bestT) {
            hit[hitCount++] = ChildHit{n.right, tr0};
          }
        }

        if (hitCount == 2 && hit[0].t > hit[1].t) {
          std::swap(hit[0], hit[1]);
        }

        // Push far child first so the near child is processed next.
        for (int i = hitCount - 1; i >= 0; --i) {
          stack.push_back(hit[i].idx);
        }
      }
    }

    if (bestId < 0 || !std::isfinite(bestT)) return false;
    if (outT) *outT = bestT;
    if (outId) *outId = bestId;
    return true;
  }

private:
  struct Item {
    Aabb3d bounds{};
    Vec3d center{0.0, 0.0, 0.0};
    int id{0};
    std::uint32_t morton{0};
  };

  struct Node {
    Aabb3d bounds{};
    int left{-1};
    int right{-1};
    std::uint32_t start{0};
    std::uint32_t count{0};

    bool isLeaf() const { return left < 0 && right < 0; }
  };

  struct Ray {
    Vec3d o;
    Vec3d d; // unit length
  };

  static bool rayAabb_(const Aabb3d& box,
                       const Ray& ray,
                       double& outTEnter,
                       double& outTExit,
                       double eps = 1e-12) {
    outTEnter = 0.0;
    outTExit = 0.0;
    if (!box.isFinite()) return false;

    double tMin = 0.0;
    double tMax = std::numeric_limits<double>::infinity();

    const auto slab = [&](double o, double d, double lo, double hi) {
      if (std::abs(d) <= eps) {
        // Parallel: must be inside slab.
        return (o >= lo && o <= hi);
      }
      const double inv = 1.0 / d;
      double t1 = (lo - o) * inv;
      double t2 = (hi - o) * inv;
      if (t1 > t2) std::swap(t1, t2);
      tMin = std::max(tMin, t1);
      tMax = std::min(tMax, t2);
      return tMax >= tMin;
    };

    if (!slab(ray.o.x, ray.d.x, box.min.x, box.max.x)) return false;
    if (!slab(ray.o.y, ray.d.y, box.min.y, box.max.y)) return false;
    if (!slab(ray.o.z, ray.d.z, box.min.z, box.max.z)) return false;

    if (tMax < 0.0) return false;
    outTEnter = std::max(0.0, tMin);
    outTExit = tMax;
    return true;
  }

  static std::uint32_t morton3D10_(std::uint32_t x, std::uint32_t y, std::uint32_t z) {
    // Interleave 10 bits from each coordinate into a 30-bit Morton code:
    //   X0 Y0 Z0 X1 Y1 Z1 ... X9 Y9 Z9
    x &= 1023u;
    y &= 1023u;
    z &= 1023u;

    std::uint32_t code = 0u;
    for (int i = 0; i < 10; ++i) {
      code |= ((x >> i) & 1u) << (3 * i);
      code |= ((y >> i) & 1u) << (3 * i + 1);
      code |= ((z >> i) & 1u) << (3 * i + 2);
    }
    return code;
  }

  static std::uint32_t quantize10_(double v, double lo, double size) {
    if (!std::isfinite(v) || !std::isfinite(lo) || !std::isfinite(size)) return 0u;
    if (!(size > 1e-18)) return 0u;

    double u = (v - lo) / size;
    if (!std::isfinite(u)) return 0u;
    u = std::clamp(u, 0.0, 1.0);

    // Map [0,1] -> [0,1023] with deterministic truncation.
    int q = static_cast<int>(u * 1024.0);
    q = std::clamp(q, 0, 1023);
    return static_cast<std::uint32_t>(q);
  }

  void computeMortonCodes_() {
    // Build centroid bounds.
    Aabb3d cb;
    for (const auto& it : items_) {
      cb.expand(it.center);
    }

    const Vec3d diag = cb.size();

    // Compute Morton code for each item centroid.
    for (auto& it : items_) {
      const std::uint32_t x = quantize10_(it.center.x, cb.min.x, diag.x);
      const std::uint32_t y = quantize10_(it.center.y, cb.min.y, diag.y);
      const std::uint32_t z = quantize10_(it.center.z, cb.min.z, diag.z);
      it.morton = morton3D10_(x, y, z);
    }
  }

  static int highestBitIndex_(std::uint32_t v) {
    if (v == 0u) return -1;
    return 31 - static_cast<int>(std::countl_zero(v));
  }

  std::uint32_t codeAt_(std::uint32_t primIndexInPrims) const {
    const std::uint32_t prim = prims_[static_cast<std::size_t>(primIndexInPrims)];
    return items_[prim].morton;
  }

  // Median-split builder (legacy).
  int buildNodeMedian_(std::uint32_t start, std::uint32_t end) {
    const std::uint32_t count = end - start;
    Node n{};
    n.start = start;
    n.count = count;

    // Compute bounds for this node.
    Aabb3d b;
    for (std::uint32_t i = start; i < end; ++i) {
      b.expand(items_[prims_[static_cast<std::size_t>(i)]].bounds);
    }
    n.bounds = b;

    const int nodeIndex = static_cast<int>(nodes_.size());
    nodes_.push_back(n);

    // Leaf criteria.
    if (count <= static_cast<std::uint32_t>(leafSize_)) {
      return nodeIndex;
    }

    // Compute centroid bounds to choose split axis.
    Aabb3d cb;
    for (std::uint32_t i = start; i < end; ++i) {
      cb.expand(items_[prims_[static_cast<std::size_t>(i)]].center);
    }

    const Vec3d diag = cb.size();
    int axis = 0;
    double best = diag.x;
    if (diag.y > best) {
      axis = 1;
      best = diag.y;
    }
    if (diag.z > best) {
      axis = 2;
      best = diag.z;
    }

    // Degenerate: all centers are (nearly) equal.
    if (!(best > 1e-18)) {
      return nodeIndex;
    }

    const std::uint32_t mid = start + count / 2;

    auto coord = [&](std::uint32_t primIndex) {
      const Vec3d c = items_[primIndex].center;
      switch (axis) {
        case 0: return c.x;
        case 1: return c.y;
        default: return c.z;
      }
    };

    std::nth_element(prims_.begin() + start,
                     prims_.begin() + mid,
                     prims_.begin() + end,
                     [&](std::uint32_t a, std::uint32_t b2) {
                       return coord(a) < coord(b2);
                     });

    const int left = buildNodeMedian_(start, mid);
    const int right = buildNodeMedian_(mid, end);

    nodes_[static_cast<std::size_t>(nodeIndex)].left = left;
    nodes_[static_cast<std::size_t>(nodeIndex)].right = right;

    // Mark internal nodes as non-leaf ranges.
    nodes_[static_cast<std::size_t>(nodeIndex)].count = 0;

    return nodeIndex;
  }

  // Morton-code builder (LBVH-inspired).
  int buildNodeMorton_(std::uint32_t start, std::uint32_t end) {
    const std::uint32_t count = end - start;

    Node n{};
    n.start = start;
    n.count = count;

    // Compute bounds for this node.
    Aabb3d b;
    for (std::uint32_t i = start; i < end; ++i) {
      b.expand(items_[prims_[static_cast<std::size_t>(i)]].bounds);
    }
    n.bounds = b;

    const int nodeIndex = static_cast<int>(nodes_.size());
    nodes_.push_back(n);

    // Leaf criteria.
    if (count <= static_cast<std::uint32_t>(leafSize_)) {
      return nodeIndex;
    }

    // Split by the highest differing Morton bit among the range.
    const std::uint32_t first = codeAt_(start);
    const std::uint32_t last = codeAt_(end - 1);
    const std::uint32_t diff = first ^ last;

    std::uint32_t split = start + count / 2;

    if (diff != 0u) {
      const int bit = highestBitIndex_(diff);
      if (bit >= 0) {
        const std::uint32_t mask = 1u << static_cast<std::uint32_t>(bit);

        // Since the array is sorted by Morton code, all entries with mask==0
        // appear before mask==1 at this highest-differing bit.
        std::uint32_t lo = start;
        std::uint32_t hi = end;
        while (lo < hi) {
          const std::uint32_t mid = lo + (hi - lo) / 2;
          const std::uint32_t c = codeAt_(mid);
          if ((c & mask) == 0u) {
            lo = mid + 1;
          } else {
            hi = mid;
          }
        }
        split = lo;

        // Safety fallback: if all codes ended up on one side, split evenly.
        if (split <= start || split >= end) {
          split = start + count / 2;
        }
      }
    }

    const int left = buildNodeMorton_(start, split);
    const int right = buildNodeMorton_(split, end);

    nodes_[static_cast<std::size_t>(nodeIndex)].left = left;
    nodes_[static_cast<std::size_t>(nodeIndex)].right = right;

    // Mark internal nodes as non-leaf ranges.
    nodes_[static_cast<std::size_t>(nodeIndex)].count = 0;

    return nodeIndex;
  }

  std::vector<Item> items_;
  std::vector<std::uint32_t> prims_;
  std::vector<Node> nodes_;
  std::size_t leafSize_{4};
  BvhBuildMode buildMode_{BvhBuildMode::Median};
};

} // namespace stellar::math
