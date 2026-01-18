#include "stellar/proc/FluidSim2D.h"

#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace stellar::proc {

namespace {

inline float safeInv(float x) {
  return (std::fabs(x) > 1.0e-12f) ? (1.0f / x) : 0.0f;
}

inline int idxS(int x, int y, int stride) {
  return x + stride * y;
}

static void setBoundaryGeneric(int b, int n, int stride, std::vector<float>& x) {
  // b: 0 scalar, 1 horizontal vel, 2 vertical vel
  for (int i = 1; i <= n; ++i) {
    x[(std::size_t)idxS(0, i, stride)]      = (b == 1) ? -x[(std::size_t)idxS(1, i, stride)] : x[(std::size_t)idxS(1, i, stride)];
    x[(std::size_t)idxS(n + 1, i, stride)]  = (b == 1) ? -x[(std::size_t)idxS(n, i, stride)] : x[(std::size_t)idxS(n, i, stride)];
    x[(std::size_t)idxS(i, 0, stride)]      = (b == 2) ? -x[(std::size_t)idxS(i, 1, stride)] : x[(std::size_t)idxS(i, 1, stride)];
    x[(std::size_t)idxS(i, n + 1, stride)]  = (b == 2) ? -x[(std::size_t)idxS(i, n, stride)] : x[(std::size_t)idxS(i, n, stride)];
  }

  x[(std::size_t)idxS(0, 0, stride)]           = 0.5f * (x[(std::size_t)idxS(1, 0, stride)] + x[(std::size_t)idxS(0, 1, stride)]);
  x[(std::size_t)idxS(0, n + 1, stride)]       = 0.5f * (x[(std::size_t)idxS(1, n + 1, stride)] + x[(std::size_t)idxS(0, n, stride)]);
  x[(std::size_t)idxS(n + 1, 0, stride)]       = 0.5f * (x[(std::size_t)idxS(n, 0, stride)] + x[(std::size_t)idxS(n + 1, 1, stride)]);
  x[(std::size_t)idxS(n + 1, n + 1, stride)]   = 0.5f * (x[(std::size_t)idxS(n, n + 1, stride)] + x[(std::size_t)idxS(n + 1, n, stride)]);
}

// --- Multigrid helpers for the pressure Poisson solve ------------------------

static void relaxPressureGS(int n, int stride,
                            std::vector<float>& p,
                            const std::vector<float>& rhs,
                            int iterations) {
  iterations = std::max(0, iterations);
  for (int k = 0; k < iterations; ++k) {
    for (int j = 1; j <= n; ++j) {
      for (int i = 1; i <= n; ++i) {
        const int id = idxS(i, j, stride);
        p[(std::size_t)id] = 0.25f * (rhs[(std::size_t)id] +
                                     p[(std::size_t)idxS(i - 1, j, stride)] +
                                     p[(std::size_t)idxS(i + 1, j, stride)] +
                                     p[(std::size_t)idxS(i, j - 1, stride)] +
                                     p[(std::size_t)idxS(i, j + 1, stride)]);
      }
    }
    setBoundaryGeneric(0, n, stride, p);
  }
}

static void computePressureResidual(int n, int stride,
                                    const std::vector<float>& p,
                                    const std::vector<float>& rhs,
                                    std::vector<float>& outRes) {
  // residual r = rhs - A p, where A p = 4p - neighbors.
  for (int j = 1; j <= n; ++j) {
    for (int i = 1; i <= n; ++i) {
      const int id = idxS(i, j, stride);
      const float Ap = 4.0f * p[(std::size_t)id] -
                       p[(std::size_t)idxS(i - 1, j, stride)] -
                       p[(std::size_t)idxS(i + 1, j, stride)] -
                       p[(std::size_t)idxS(i, j - 1, stride)] -
                       p[(std::size_t)idxS(i, j + 1, stride)];
      outRes[(std::size_t)id] = rhs[(std::size_t)id] - Ap;
    }
  }
  setBoundaryGeneric(0, n, stride, outRes);
}

static void restrictResidualFullWeighting(int nFine, int strideFine,
                                          const std::vector<float>& fineRes,
                                          int nCoarse, int strideCoarse,
                                          std::vector<float>& coarseRhs) {
  // Scale rhs by (h_c^2 / h_f^2) ~= (nFine / nCoarse)^2 because the Stable-Fluids
  // formulation embeds h^2 into the rhs.
  const float ratio = (nCoarse > 0) ? ((float)nFine / (float)nCoarse) : 1.0f;
  const float scale = ratio * ratio;

  for (int J = 1; J <= nCoarse; ++J) {
    for (int I = 1; I <= nCoarse; ++I) {
      const int i = 2 * I - 1;
      const int j = 2 * J - 1;

      const float c = fineRes[(std::size_t)idxS(i, j, strideFine)] * 0.25f +
                      (fineRes[(std::size_t)idxS(i - 1, j, strideFine)] +
                       fineRes[(std::size_t)idxS(i + 1, j, strideFine)] +
                       fineRes[(std::size_t)idxS(i, j - 1, strideFine)] +
                       fineRes[(std::size_t)idxS(i, j + 1, strideFine)]) * 0.125f +
                      (fineRes[(std::size_t)idxS(i - 1, j - 1, strideFine)] +
                       fineRes[(std::size_t)idxS(i + 1, j - 1, strideFine)] +
                       fineRes[(std::size_t)idxS(i - 1, j + 1, strideFine)] +
                       fineRes[(std::size_t)idxS(i + 1, j + 1, strideFine)]) * 0.0625f;

      coarseRhs[(std::size_t)idxS(I, J, strideCoarse)] = c * scale;
    }
  }

  setBoundaryGeneric(0, nCoarse, strideCoarse, coarseRhs);
}

static void prolongateAndAdd(int nFine, int strideFine,
                             std::vector<float>& fineP,
                             int nCoarse, int strideCoarse,
                             const std::vector<float>& coarseE) {
  // Cell-centered bilinear prolongation. Coarse cell centers align with odd
  // fine indices (2I-1).
  for (int j = 1; j <= nFine; ++j) {
    const int J = (j + 1) / 2;
    const float wy = (j % 2 == 0) ? 0.5f : 0.0f;

    for (int i = 1; i <= nFine; ++i) {
      const int I = (i + 1) / 2;
      const float wx = (i % 2 == 0) ? 0.5f : 0.0f;

      const float c00 = coarseE[(std::size_t)idxS(I, J, strideCoarse)];
      const float c10 = coarseE[(std::size_t)idxS(I + 1, J, strideCoarse)];
      const float c01 = coarseE[(std::size_t)idxS(I, J + 1, strideCoarse)];
      const float c11 = coarseE[(std::size_t)idxS(I + 1, J + 1, strideCoarse)];

      const float c0 = c00 * (1.0f - wx) + c10 * wx;
      const float c1 = c01 * (1.0f - wx) + c11 * wx;
      const float e = c0 * (1.0f - wy) + c1 * wy;

      fineP[(std::size_t)idxS(i, j, strideFine)] += e;
    }
  }

  setBoundaryGeneric(0, nFine, strideFine, fineP);
}

static float rmsInterior(int n, int stride, const std::vector<float>& x) {
  if (n <= 0) return 0.0f;
  double sum2 = 0.0;
  int count = 0;
  for (int j = 1; j <= n; ++j) {
    for (int i = 1; i <= n; ++i) {
      const float v = x[(std::size_t)idxS(i, j, stride)];
      sum2 += (double)v * (double)v;
      ++count;
    }
  }
  if (count <= 0) return 0.0f;
  return (float)std::sqrt(sum2 / (double)count);
}

static void computeVelocityDivergenceStats(int n, int stride,
                                           const std::vector<float>& u,
                                           const std::vector<float>& v,
                                           FluidSim2DStats& stats) {
  stats.maxAbsDivergence = 0.0f;
  stats.rmsDivergence = 0.0f;

  if (n < 3) return;

  double sum2 = 0.0;
  int count = 0;
  float maxAbs = 0.0f;

  for (int j = 2; j <= n - 1; ++j) {
    for (int i = 2; i <= n - 1; ++i) {
      const int id = idxS(i, j, stride);

      // Physical-ish divergence (dudx + dvdy) on the unit square domain.
      const float div = 0.5f * (float)n *
                        (u[(std::size_t)(id + 1)] - u[(std::size_t)(id - 1)] +
                         v[(std::size_t)(id + stride)] - v[(std::size_t)(id - stride)]);

      const float a = std::fabs(div);
      maxAbs = std::max(maxAbs, a);
      sum2 += (double)div * (double)div;
      ++count;
    }
  }

  stats.maxAbsDivergence = maxAbs;
  stats.rmsDivergence = (count > 0) ? (float)std::sqrt(sum2 / (double)count) : 0.0f;
}

} // namespace

void FluidSim2D::resize(int gridSize) {
  n_ = std::max(1, gridSize);
  stride_ = n_ + 2;
  size_ = stride_ * stride_;

  u_.assign(size_, 0.0f);
  v_.assign(size_, 0.0f);
  uTmp_.assign(size_, 0.0f);
  vTmp_.assign(size_, 0.0f);
  p_.assign(size_, 0.0f);
  div_.assign(size_, 0.0f);
  curl_.assign(size_, 0.0f);

  dyeR_.assign(size_, 0.0f);
  dyeG_.assign(size_, 0.0f);
  dyeB_.assign(size_, 0.0f);
  dyeTmp_.assign(size_, 0.0f);
  dyeTmp2_.assign(size_, 0.0f);

  rebuildMultigridLevels();
}

void FluidSim2D::rebuildMultigridLevels() {
  mg_.clear();
  mgRes0_.assign((std::size_t)size_, 0.0f);

  int n = n_;
  while (n > 2) {
    n = (n + 1) / 2;

    MgLevel lvl;
    lvl.n = n;
    lvl.stride = n + 2;
    lvl.size = lvl.stride * lvl.stride;

    lvl.p.assign((std::size_t)lvl.size, 0.0f);
    lvl.rhs.assign((std::size_t)lvl.size, 0.0f);
    lvl.res.assign((std::size_t)lvl.size, 0.0f);

    mg_.push_back(std::move(lvl));
  }
}

void FluidSim2D::clear() {
  std::fill(u_.begin(), u_.end(), 0.0f);
  std::fill(v_.begin(), v_.end(), 0.0f);
  std::fill(uTmp_.begin(), uTmp_.end(), 0.0f);
  std::fill(vTmp_.begin(), vTmp_.end(), 0.0f);
  std::fill(p_.begin(), p_.end(), 0.0f);
  std::fill(div_.begin(), div_.end(), 0.0f);
  std::fill(curl_.begin(), curl_.end(), 0.0f);
  std::fill(dyeR_.begin(), dyeR_.end(), 0.0f);
  std::fill(dyeG_.begin(), dyeG_.end(), 0.0f);
  std::fill(dyeB_.begin(), dyeB_.end(), 0.0f);
  std::fill(dyeTmp_.begin(), dyeTmp_.end(), 0.0f);
  std::fill(dyeTmp2_.begin(), dyeTmp2_.end(), 0.0f);

  for (auto& lvl : mg_) {
    std::fill(lvl.p.begin(), lvl.p.end(), 0.0f);
    std::fill(lvl.rhs.begin(), lvl.rhs.end(), 0.0f);
    std::fill(lvl.res.begin(), lvl.res.end(), 0.0f);
  }
  std::fill(mgRes0_.begin(), mgRes0_.end(), 0.0f);

  stats_ = FluidSim2DStats{};
}

void FluidSim2D::addVelocity(int x, int y, float vx, float vy) {
  if (x < 1 || y < 1 || x > n_ || y > n_) return;
  const int id = idx(x, y);
  u_[(std::size_t)id] += vx;
  v_[(std::size_t)id] += vy;
}

void FluidSim2D::addDye(int x, int y, float r, float g, float b) {
  if (x < 1 || y < 1 || x > n_ || y > n_) return;
  const int id = idx(x, y);
  dyeR_[(std::size_t)id] += r;
  dyeG_[(std::size_t)id] += g;
  dyeB_[(std::size_t)id] += b;
}

void FluidSim2D::dyeAt(int x, int y, float* outR, float* outG, float* outB) const {
  if (outR) *outR = 0.0f;
  if (outG) *outG = 0.0f;
  if (outB) *outB = 0.0f;
  if (x < 1 || y < 1 || x > n_ || y > n_) return;
  const int id = idx(x, y);
  if (outR) *outR = dyeR_[(std::size_t)id];
  if (outG) *outG = dyeG_[(std::size_t)id];
  if (outB) *outB = dyeB_[(std::size_t)id];
}

void FluidSim2D::splat(float x01, float y01, float radius01,
                       float velX, float velY,
                       float dyeR, float dyeG, float dyeB) {
  if (n_ <= 0) return;

  x01 = std::clamp(x01, 0.0f, 1.0f);
  y01 = std::clamp(y01, 0.0f, 1.0f);
  radius01 = std::max(0.0f, radius01);

  const float cx = 1.0f + x01 * (float)n_;
  const float cy = 1.0f + y01 * (float)n_;
  const float rCells = std::max(1.0f, radius01 * (float)n_);
  const float r2 = rCells * rCells;

  const int x0 = std::max(1, (int)std::floor(cx - rCells));
  const int x1 = std::min(n_, (int)std::ceil(cx + rCells));
  const int y0 = std::max(1, (int)std::floor(cy - rCells));
  const int y1 = std::min(n_, (int)std::ceil(cy + rCells));

  // Gaussian-ish kernel for a smooth splat.
  const float invSigma2 = 1.0f / std::max(1.0e-6f, 0.35f * r2);

  for (int j = y0; j <= y1; ++j) {
    for (int i = x0; i <= x1; ++i) {
      const float dx = (float)i - cx;
      const float dy = (float)j - cy;
      const float d2 = dx * dx + dy * dy;
      if (d2 > r2) continue;
      const float w = std::exp(-d2 * invSigma2);
      const int id = idx(i, j);
      u_[(std::size_t)id] += velX * w;
      v_[(std::size_t)id] += velY * w;
      dyeR_[(std::size_t)id] += dyeR * w;
      dyeG_[(std::size_t)id] += dyeG * w;
      dyeB_[(std::size_t)id] += dyeB * w;
    }
  }
}

void FluidSim2D::setBoundary(int b, std::vector<float>& x) {
  setBoundaryGeneric(b, n_, stride_, x);
}

void FluidSim2D::linSolve(int b, std::vector<float>& x, const std::vector<float>& x0,
                          float a, float c, int iterations) {
  const float invC = safeInv(c);
  for (int k = 0; k < iterations; ++k) {
    for (int j = 1; j <= n_; ++j) {
      for (int i = 1; i <= n_; ++i) {
        x[(std::size_t)idx(i, j)] = (x0[(std::size_t)idx(i, j)] + a * (x[(std::size_t)idx(i - 1, j)] + x[(std::size_t)idx(i + 1, j)] + x[(std::size_t)idx(i, j - 1)] + x[(std::size_t)idx(i, j + 1)])) * invC;
      }
    }
    setBoundary(b, x);
  }
}

void FluidSim2D::diffuse(int b, std::vector<float>& x, const std::vector<float>& x0,
                         float diff, float dt, int iterations) {
  const float a = dt * diff * (float)(n_ * n_);
  linSolve(b, x, x0, a, 1.0f + 4.0f * a, iterations);
}

void FluidSim2D::advect(int b, std::vector<float>& d, const std::vector<float>& d0,
                        const std::vector<float>& u, const std::vector<float>& v, float dt) {
  const float dt0 = dt * (float)n_;

  for (int j = 1; j <= n_; ++j) {
    for (int i = 1; i <= n_; ++i) {
      float x = (float)i - dt0 * u[(std::size_t)idx(i, j)];
      float y = (float)j - dt0 * v[(std::size_t)idx(i, j)];

      x = std::clamp(x, 0.5f, (float)n_ + 0.5f);
      y = std::clamp(y, 0.5f, (float)n_ + 0.5f);

      const int i0 = (int)std::floor(x);
      const int i1 = i0 + 1;
      const int j0 = (int)std::floor(y);
      const int j1 = j0 + 1;

      const float s1 = x - (float)i0;
      const float s0 = 1.0f - s1;
      const float t1 = y - (float)j0;
      const float t0 = 1.0f - t1;

      d[(std::size_t)idx(i, j)] =
          s0 * (t0 * d0[(std::size_t)idx(i0, j0)] + t1 * d0[(std::size_t)idx(i0, j1)]) +
          s1 * (t0 * d0[(std::size_t)idx(i1, j0)] + t1 * d0[(std::size_t)idx(i1, j1)]);
    }
  }

  setBoundary(b, d);
}

void FluidSim2D::advectBFECC(int b, std::vector<float>& d, const std::vector<float>& d0,
                             const std::vector<float>& u, const std::vector<float>& v, float dt,
                             float correction, bool clampExtrema,
                             std::vector<float>& tmpBack) {
  correction = std::clamp(correction, 0.0f, 1.0f);
  if (correction <= 0.0f) {
    advect(b, d, d0, u, v, dt);
    return;
  }

  // Ensure scratch has the right size.
  if ((int)tmpBack.size() != size_) tmpBack.assign((std::size_t)size_, 0.0f);

  // 1) Forward advect (semi-Lagrangian) directly into 'd'.
  advect(b, d, d0, u, v, dt);

  // 2) Backward advect the result to estimate the local truncation error.
  advect(b, tmpBack, d, u, v, -dt);

  // 3) Error compensation/correction.
  //    d_corrected = d_forward + 0.5 * (d0 - d_back)
  const float dt0 = dt * (float)n_;

  for (int j = 1; j <= n_; ++j) {
    for (int i = 1; i <= n_; ++i) {
      const int id = idx(i, j);

      float val = d[(std::size_t)id] + 0.5f * correction * (d0[(std::size_t)id] - tmpBack[(std::size_t)id]);

      if (clampExtrema) {
        // Clamp to the min/max of the *source* interpolation stencil to prevent
        // introducing new extrema (reduces ringing / overshoot).
        float x = (float)i - dt0 * u[(std::size_t)id];
        float y = (float)j - dt0 * v[(std::size_t)id];

        x = std::clamp(x, 0.5f, (float)n_ + 0.5f);
        y = std::clamp(y, 0.5f, (float)n_ + 0.5f);

        const int i0 = (int)std::floor(x);
        const int i1 = i0 + 1;
        const int j0 = (int)std::floor(y);
        const int j1 = j0 + 1;

        const float s00 = d0[(std::size_t)idx(i0, j0)];
        const float s01 = d0[(std::size_t)idx(i0, j1)];
        const float s10 = d0[(std::size_t)idx(i1, j0)];
        const float s11 = d0[(std::size_t)idx(i1, j1)];

        const float mn = std::min(std::min(s00, s01), std::min(s10, s11));
        const float mx = std::max(std::max(s00, s01), std::max(s10, s11));
        val = std::clamp(val, mn, mx);
      }

      d[(std::size_t)id] = val;
    }
  }

  setBoundary(b, d);
}

void FluidSim2D::solvePressureMultigrid(int vCycles) {
  if (mg_.empty()) {
    // Too small to build a useful hierarchy.
    return;
  }

  vCycles = std::clamp(vCycles, 1, 16);

  const int pre = std::clamp(params_.multigridPreSmooth, 1, 8);
  const int post = std::clamp(params_.multigridPostSmooth, 1, 8);
  const int coarseIters = std::clamp(params_.multigridCoarseIters, 4, 400);

  // Make sure scratch is valid (resize can be called with N=1).
  if ((int)mgRes0_.size() != size_) mgRes0_.assign((std::size_t)size_, 0.0f);

  for (int c = 0; c < vCycles; ++c) {
    mgVCycle(0, pre, post, coarseIters);
  }
}

void FluidSim2D::mgVCycle(int level, int preSmooth, int postSmooth, int coarseIters) {
  const int maxLevel = (int)mg_.size();
  if (level < 0 || level > maxLevel) return;

  int n = 0;
  int stride = 0;
  std::vector<float>* p = nullptr;
  const std::vector<float>* rhs = nullptr;
  std::vector<float>* res = nullptr;

  if (level == 0) {
    n = n_;
    stride = stride_;
    p = &p_;
    rhs = &div_;
    res = &mgRes0_;
  } else {
    MgLevel& L = mg_[(std::size_t)(level - 1)];
    n = L.n;
    stride = L.stride;
    p = &L.p;
    rhs = &L.rhs;
    res = &L.res;
  }

  if (!p || !rhs || !res) return;

  if (level == maxLevel) {
    // Coarsest solve.
    relaxPressureGS(n, stride, *p, *rhs, coarseIters);
    return;
  }

  // Pre-smooth.
  relaxPressureGS(n, stride, *p, *rhs, preSmooth);

  // Compute residual.
  if ((int)res->size() != (stride * stride)) res->assign((std::size_t)(stride * stride), 0.0f);
  computePressureResidual(n, stride, *p, *rhs, *res);

  // Restrict residual to the next level RHS.
  MgLevel& coarse = mg_[(std::size_t)level]; // next level
  restrictResidualFullWeighting(n, stride, *res, coarse.n, coarse.stride, coarse.rhs);

  // Solve A e = r on coarse grid, starting from 0.
  std::fill(coarse.p.begin(), coarse.p.end(), 0.0f);
  setBoundaryGeneric(0, coarse.n, coarse.stride, coarse.p);

  mgVCycle(level + 1, preSmooth, postSmooth, coarseIters);

  // Prolongate and correct.
  prolongateAndAdd(n, stride, *p, coarse.n, coarse.stride, coarse.p);

  // Post-smooth.
  relaxPressureGS(n, stride, *p, *rhs, postSmooth);
}

void FluidSim2D::project(std::vector<float>& u, std::vector<float>& v, int iterations) {
  // Compute divergence into div_ and clear pressure p_.
  for (int j = 1; j <= n_; ++j) {
    for (int i = 1; i <= n_; ++i) {
      div_[(std::size_t)idx(i, j)] =
          -0.5f * (u[(std::size_t)idx(i + 1, j)] - u[(std::size_t)idx(i - 1, j)] + v[(std::size_t)idx(i, j + 1)] - v[(std::size_t)idx(i, j - 1)]) / (float)n_;
      p_[(std::size_t)idx(i, j)] = 0.0f;
    }
  }

  setBoundary(0, div_);
  setBoundary(0, p_);

  // Choose projection solver.
  const bool canMG = params_.multigridProjection && !mg_.empty() && n_ >= 8;

  stats_ = FluidSim2DStats{};
  stats_.usedMultigrid = false;
  stats_.multigridLevels = 0;
  stats_.multigridVCycles = 0;
  stats_.gaussSeidelIterations = 0;

  if (canMG) {
    // Map the UI iteration budget to a small number of V-cycles.
    const int it = std::clamp(iterations, 1, 120);
    const int vCycles = std::clamp((it + 9) / 10, 1, 12);

    stats_.usedMultigrid = true;
    stats_.multigridLevels = (int)mg_.size() + 1;
    stats_.multigridVCycles = vCycles;

    solvePressureMultigrid(vCycles);
  } else {
    iterations = std::clamp(iterations, 1, 200);
    stats_.gaussSeidelIterations = iterations;

    linSolve(0, p_, div_, 1.0f, 4.0f, iterations);
  }

  // Pressure solve residual RMS (diagnostic).
  if ((int)mgRes0_.size() != size_) mgRes0_.assign((std::size_t)size_, 0.0f);
  computePressureResidual(n_, stride_, p_, div_, mgRes0_);
  stats_.pressureResidualRms = rmsInterior(n_, stride_, mgRes0_);

  // Subtract pressure gradient.
  for (int j = 1; j <= n_; ++j) {
    for (int i = 1; i <= n_; ++i) {
      u[(std::size_t)idx(i, j)] -= 0.5f * (float)n_ * (p_[(std::size_t)idx(i + 1, j)] - p_[(std::size_t)idx(i - 1, j)]);
      v[(std::size_t)idx(i, j)] -= 0.5f * (float)n_ * (p_[(std::size_t)idx(i, j + 1)] - p_[(std::size_t)idx(i, j - 1)]);
    }
  }

  setBoundary(1, u);
  setBoundary(2, v);

  computeVelocityDivergenceStats(n_, stride_, u, v, stats_);
}

void FluidSim2D::applyVorticityConfinement(float dt) {
  if (params_.vorticityConfinement <= 0.0f) return;

  // Compute scalar curl (vorticity) for each cell.
  for (int j = 1; j <= n_; ++j) {
    for (int i = 1; i <= n_; ++i) {
      const float dw = (v_[(std::size_t)idx(i + 1, j)] - v_[(std::size_t)idx(i - 1, j)]) * 0.5f -
                       (u_[(std::size_t)idx(i, j + 1)] - u_[(std::size_t)idx(i, j - 1)]) * 0.5f;
      curl_[(std::size_t)idx(i, j)] = dw;
    }
  }

  const float eps = params_.vorticityConfinement;

  for (int j = 2; j <= n_ - 1; ++j) {
    for (int i = 2; i <= n_ - 1; ++i) {
      const float cL = std::fabs(curl_[(std::size_t)idx(i - 1, j)]);
      const float cR = std::fabs(curl_[(std::size_t)idx(i + 1, j)]);
      const float cB = std::fabs(curl_[(std::size_t)idx(i, j - 1)]);
      const float cT = std::fabs(curl_[(std::size_t)idx(i, j + 1)]);

      float nx = cR - cL;
      float ny = cT - cB;

      const float len = std::sqrt(nx * nx + ny * ny) + 1.0e-6f;
      nx /= len;
      ny /= len;

      const float w = curl_[(std::size_t)idx(i, j)];
      const float f = eps * w;

      // Force is perpendicular to grad(|w|) and scaled by vorticity.
      u_[(std::size_t)idx(i, j)] += dt * ny * f;
      v_[(std::size_t)idx(i, j)] -= dt * nx * f;
    }
  }
}

void FluidSim2D::applyCurlNoise(core::u64 seed, float timeSec, float dt) {
  if (params_.curlNoiseStrength <= 0.0f) return;

  const float freq = std::max(0.0001f, params_.curlNoiseFrequency);
  const float t = timeSec * params_.curlNoiseTimeScale;
  const float strength = params_.curlNoiseStrength;

  // Finite difference step in *domain* units.
  const float h = 1.0f / (float)n_;

  auto psi = [&](float x, float y) -> float {
    // perlin2D returns [0,1]; map to [-1,1].
    const double n = perlin2D(seed, (double)x * (double)freq + (double)t, (double)y * (double)freq + (double)t);
    return (float)(2.0 * (n - 0.5));
  };

  for (int j = 1; j <= n_; ++j) {
    const float y = ((float)j - 0.5f) / (float)n_;
    for (int i = 1; i <= n_; ++i) {
      const float x = ((float)i - 0.5f) / (float)n_;

      // Velocity from stream function psi:
      //   u = dpsi/dy
      //   v = -dpsi/dx
      const float dpsi_dx = (psi(x + h, y) - psi(x - h, y)) * (0.5f / h);
      const float dpsi_dy = (psi(x, y + h) - psi(x, y - h)) * (0.5f / h);

      const float fx = dpsi_dy;
      const float fy = -dpsi_dx;

      u_[(std::size_t)idx(i, j)] += strength * fx * dt;
      v_[(std::size_t)idx(i, j)] += strength * fy * dt;
    }
  }
}

void FluidSim2D::applyDissipation(std::vector<float>& x, float decayPerSec, float dt) {
  if (decayPerSec <= 0.0f) return;
  const float k = std::exp(-decayPerSec * dt);
  for (int j = 1; j <= n_; ++j) {
    for (int i = 1; i <= n_; ++i) {
      x[(std::size_t)idx(i, j)] *= k;
    }
  }
}

void FluidSim2D::clampFields() {
  const float maxD = std::max(0.0f, params_.maxDye);
  const float maxV = std::max(0.0f, params_.maxSpeed);

  for (int j = 1; j <= n_; ++j) {
    for (int i = 1; i <= n_; ++i) {
      const int id = idx(i, j);
      u_[(std::size_t)id] = std::clamp(u_[(std::size_t)id], -maxV, maxV);
      v_[(std::size_t)id] = std::clamp(v_[(std::size_t)id], -maxV, maxV);

      dyeR_[(std::size_t)id] = std::clamp(dyeR_[(std::size_t)id], 0.0f, maxD);
      dyeG_[(std::size_t)id] = std::clamp(dyeG_[(std::size_t)id], 0.0f, maxD);
      dyeB_[(std::size_t)id] = std::clamp(dyeB_[(std::size_t)id], 0.0f, maxD);
    }
  }
}

void FluidSim2D::step(float dtSec, int iterations, core::u64 noiseSeed, float timeSec) {
  if (n_ <= 0) return;
  if (dtSec <= 0.0f) return;
  iterations = std::clamp(iterations, 1, 80);

  // Procedural forcing is applied as an external force before the viscous solve.
  applyCurlNoise(noiseSeed, timeSec, dtSec);

  // Velocity damping (simple exponential decay).
  if (params_.velocityDamping > 0.0f) {
    const float k = std::exp(-params_.velocityDamping * dtSec);
    for (int j = 1; j <= n_; ++j) {
      for (int i = 1; i <= n_; ++i) {
        const int id = idx(i, j);
        u_[(std::size_t)id] *= k;
        v_[(std::size_t)id] *= k;
      }
    }
  }

  // Viscosity (diffuse velocities).
  diffuse(1, uTmp_, u_, params_.viscosity, dtSec, iterations);
  diffuse(2, vTmp_, v_, params_.viscosity, dtSec, iterations);
  project(uTmp_, vTmp_, iterations);

  // Advect velocity.
  advect(1, u_, uTmp_, uTmp_, vTmp_, dtSec);
  advect(2, v_, vTmp_, uTmp_, vTmp_, dtSec);

  // Artist-friendly detail booster.
  applyVorticityConfinement(dtSec);

  // Enforce incompressibility.
  project(u_, v_, iterations);

  // Dye: diffuse + advect each channel.
  diffuse(0, dyeTmp_, dyeR_, params_.diffusion, dtSec, iterations);
  advectBFECC(0, dyeR_, dyeTmp_, u_, v_, dtSec,
              params_.dyeAdvectionCorrection, params_.dyeAdvectionClamp,
              dyeTmp2_);

  diffuse(0, dyeTmp_, dyeG_, params_.diffusion, dtSec, iterations);
  advectBFECC(0, dyeG_, dyeTmp_, u_, v_, dtSec,
              params_.dyeAdvectionCorrection, params_.dyeAdvectionClamp,
              dyeTmp2_);

  diffuse(0, dyeTmp_, dyeB_, params_.diffusion, dtSec, iterations);
  advectBFECC(0, dyeB_, dyeTmp_, u_, v_, dtSec,
              params_.dyeAdvectionCorrection, params_.dyeAdvectionClamp,
              dyeTmp2_);

  applyDissipation(dyeR_, params_.dyeDissipation, dtSec);
  applyDissipation(dyeG_, params_.dyeDissipation, dtSec);
  applyDissipation(dyeB_, params_.dyeDissipation, dtSec);

  clampFields();
}

} // namespace stellar::proc
