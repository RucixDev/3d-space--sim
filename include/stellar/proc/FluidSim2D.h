#pragma once

#include "stellar/core/Types.h"

#include <vector>

namespace stellar::proc {

// A small 2D incompressible “stable fluids” style solver suitable for
// procedural VFX / texture generation and lightweight gameplay sims.
//
// Design goals:
//  - dependency-light (CPU only; safe for headless builds)
//  - deterministic for a given sequence of inputs
//  - tuned for artist-friendly procedural forcing (curl noise + vorticity)
//
// Grid layout:
//  - The internal storage uses a padded (N+2)x(N+2) grid.
//  - Valid simulation cells are (1..N, 1..N).
//  - Cells at 0 and N+1 are boundary “ghost” cells.

struct FluidSim2DParams {
  // Kinematic viscosity (velocity diffusion). Larger values produce smoother
  // / more laminar flow.
  float viscosity{0.0001f};

  // Diffusion for dye channels.
  float diffusion{0.0f};

  // Exponential dye dissipation (per second). 0 disables.
  float dyeDissipation{0.15f};

  // Exponential velocity damping (per second). 0 disables.
  float velocityDamping{0.0f};

  // Dye advection quality: Back-and-Forth Error Compensation and Correction (BFECC / MacCormack).
  //
  //  0.0 = plain semi-Lagrangian (very stable, but noticeably diffusive)
  //  1.0 = full BFECC correction (sharper dye, less numerical diffusion)
  //
  // Values in between blend the correction.
  float dyeAdvectionCorrection{1.0f};

  // Clamp corrected dye to the local min/max of the source stencil to avoid
  // introducing new extrema (prevents ringing/overshoot).
  bool dyeAdvectionClamp{true};

  // Vorticity confinement strength. 0 disables.
  // Adds back small-scale “rolling” motion that semi-Lagrangian advection can
  // otherwise damp out.
  float vorticityConfinement{0.0f};

  // Procedural divergence-free forcing generated from a time-varying stream
  // function (2D curl noise). 0 disables.
  float curlNoiseStrength{0.0f};
  float curlNoiseFrequency{2.5f};
  float curlNoiseTimeScale{0.15f};

  // Projection solver (pressure Poisson solve).
  //
  // The classic Stable Fluids approach uses Gauss-Seidel relaxation for the
  // pressure solve, which is robust but can require many iterations for large
  // grids. When multigridProjection is enabled, the pressure solve is
  // accelerated with a geometric multigrid V-cycle while keeping the same
  // boundary conditions and determinism.
  bool multigridProjection{true};

  // Smoothing iterations per V-cycle (only used when multigridProjection==true).
  // Keep these small (1–3); the V-cycle provides the heavy lifting.
  int multigridPreSmooth{2};
  int multigridPostSmooth{2};

  // Extra relax iterations on the coarsest grid (only used when multigridProjection==true).
  int multigridCoarseIters{40};

  // Safety clamps (keeps the sim numerically well-behaved under aggressive
  // input splats).
  float maxSpeed{50.0f};
  float maxDye{25.0f};
};

struct FluidSim2DStats {
  // Divergence metrics after the last projection step.
  //
  // These are reported in *physical-ish* units (~ dudx + dvdy) so values are
  // comparable across grid resolutions.
  float maxAbsDivergence{0.0f};
  float rmsDivergence{0.0f};

  // RMS of the pressure solve residual b - A p after the last projection.
  // (In the same units as the internal rhs 'div_' array.)
  float pressureResidualRms{0.0f};

  // Projection method bookkeeping.
  bool usedMultigrid{false};
  int multigridLevels{0};
  int multigridVCycles{0};
  int gaussSeidelIterations{0};
};

class FluidSim2D {
public:
  FluidSim2D() = default;
  explicit FluidSim2D(int gridSize) { resize(gridSize); }

  void resize(int gridSize);
  void clear();

  // Advance simulation by dt seconds.
  //
  // `iterations` controls the relaxation steps in diffusion and pressure solve
  // passes. Typical values: 10–40.
  //
  // Note: when multigridProjection is enabled, `iterations` is interpreted as
  // a *work budget* and is mapped to a small number of V-cycles.
  //
  // `noiseSeed` and `timeSec` are used only when curl noise forcing is enabled.
  void step(float dtSec, int iterations = 20, core::u64 noiseSeed = 0, float timeSec = 0.0f);

  // Inject dye + velocity using normalized coordinates (0..1) in sim space.
  //
  // `radius01` is relative to the unit domain; e.g. 0.05 means 5% of width.
  void splat(float x01, float y01,
             float radius01,
             float velX, float velY,
             float dyeR, float dyeG, float dyeB);

  // Low-level injection in cell coordinates (1..N).
  void addVelocity(int x, int y, float vx, float vy);
  void addDye(int x, int y, float r, float g, float b);

  int gridSize() const { return n_; }
  int paddedStride() const { return stride_; }

  FluidSim2DParams& params() { return params_; }
  const FluidSim2DParams& params() const { return params_; }

  const FluidSim2DStats& stats() const { return stats_; }

  const std::vector<float>& u() const { return u_; }
  const std::vector<float>& v() const { return v_; }
  const std::vector<float>& dyeR() const { return dyeR_; }
  const std::vector<float>& dyeG() const { return dyeG_; }
  const std::vector<float>& dyeB() const { return dyeB_; }

  // Dye sample at cell coordinates (1..N). Returns 0 outside.
  void dyeAt(int x, int y, float* outR, float* outG, float* outB) const;

private:
  int n_{0};
  int stride_{0};
  int size_{0};
  FluidSim2DParams params_{};
  FluidSim2DStats stats_{};

  std::vector<float> u_, v_;
  std::vector<float> uTmp_, vTmp_;
  std::vector<float> p_, div_;
  std::vector<float> curl_;

  std::vector<float> dyeR_, dyeG_, dyeB_;
  std::vector<float> dyeTmp_;
  std::vector<float> dyeTmp2_;

  // Multigrid acceleration state (allocated on resize).
  struct MgLevel {
    int n{0};
    int stride{0};
    int size{0};
    std::vector<float> p{};
    std::vector<float> rhs{};
    std::vector<float> res{};
  };
  std::vector<MgLevel> mg_{};
  std::vector<float> mgRes0_{};

  int idx(int x, int y) const { return x + stride_ * y; }

  void setBoundary(int b, std::vector<float>& x);
  void linSolve(int b, std::vector<float>& x, const std::vector<float>& x0,
                float a, float c, int iterations);
  void diffuse(int b, std::vector<float>& x, const std::vector<float>& x0,
               float diff, float dt, int iterations);
  void advect(int b, std::vector<float>& d, const std::vector<float>& d0,
              const std::vector<float>& u, const std::vector<float>& v, float dt);

  // BFECC (MacCormack-style) corrected advection for scalar fields.
  //
  // Implementation note: the forward advect is written directly into 'd', and
  // tmpBack stores the backward-advected estimate used for error compensation.
  void advectBFECC(int b, std::vector<float>& d, const std::vector<float>& d0,
                   const std::vector<float>& u, const std::vector<float>& v, float dt,
                   float correction, bool clampExtrema,
                   std::vector<float>& tmpBack);

  void project(std::vector<float>& u, std::vector<float>& v, int iterations);

  // Multigrid pressure solve (used by project()).
  void rebuildMultigridLevels();
  void solvePressureMultigrid(int vCycles);
  void mgVCycle(int level, int preSmooth, int postSmooth, int coarseIters);

  void applyVorticityConfinement(float dt);
  void applyCurlNoise(core::u64 seed, float timeSec, float dt);

  void applyDissipation(std::vector<float>& x, float decayPerSec, float dt);
  void clampFields();
};

} // namespace stellar::proc
