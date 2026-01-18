#include "test_harness.h"

#include "stellar/proc/FluidSim2D.h"

#include <cmath>
#include <cstddef>
#include <limits>

// A lightweight deterministic smoke test for the Stable-Fluids style solver.
int test_fluid_sim2d() {
  int failures = 0;

  using stellar::proc::FluidSim2D;

  FluidSim2D sim(48);

  auto& p = sim.params();
  p.viscosity = 0.0001f;
  p.diffusion = 0.0f;
  p.dyeDissipation = 0.0f;
  p.velocityDamping = 0.0f;
  p.dyeAdvectionCorrection = 1.0f;
  p.dyeAdvectionClamp = true;
  p.vorticityConfinement = 0.0f;
  p.curlNoiseStrength = 0.0f;
  p.maxSpeed = 200.0f;
  p.maxDye = 200.0f;

  // Inject a deterministic blob and a lateral velocity impulse.
  sim.splat(0.50f, 0.50f, 0.07f, 35.0f, 0.0f, 12.0f, 0.0f, 0.0f);

  const int n = sim.gridSize();
  const int stride = sim.paddedStride();

  auto maxAbsDiv = [&]() -> float {
    float m = 0.0f;
    const auto& u = sim.u();
    const auto& v = sim.v();
    for (int j = 2; j <= n - 1; ++j) {
      for (int i = 2; i <= n - 1; ++i) {
        const int id = i + stride * j;
        const float div = 0.5f * (u[(std::size_t)(id + 1)] - u[(std::size_t)(id - 1)] +
                                  v[(std::size_t)(id + stride)] - v[(std::size_t)(id - stride)]) / (float)n;
        m = std::max(m, std::fabs(div));
      }
    }
    return m;
  };

  // Run a handful of steps and ensure:
  //  - no NaNs/Infs
  //  - divergence stays reasonably small
  for (int k = 0; k < 24; ++k) {
    sim.step(1.0f / 60.0f, 25, 1234u, (float)k * (1.0f / 60.0f));

    const auto& u = sim.u();
    const auto& v = sim.v();
    const auto& r = sim.dyeR();
    const auto& g = sim.dyeG();
    const auto& b = sim.dyeB();

    for (int j = 1; j <= n; ++j) {
      for (int i = 1; i <= n; ++i) {
        const std::size_t id = (std::size_t)(i + stride * j);
        CHECK(std::isfinite(u[id]));
        CHECK(std::isfinite(v[id]));
        CHECK(std::isfinite(r[id]));
        CHECK(std::isfinite(g[id]));
        CHECK(std::isfinite(b[id]));

        // Dye fields are clamped every step; ensure advection never produces runaway values.
        const float eps = 1.0e-3f;
        CHECK(r[id] >= -eps);
        CHECK(g[id] >= -eps);
        CHECK(b[id] >= -eps);
        CHECK(r[id] <= p.maxDye + eps);
        CHECK(g[id] <= p.maxDye + eps);
        CHECK(b[id] <= p.maxDye + eps);
      }
    }

    const float d = maxAbsDiv();
    CHECK(d < 5.0e-2f);
  }

  // Ensure dye exists somewhere.
  float total = 0.0f;
  {
    const auto& r = sim.dyeR();
    for (int j = 1; j <= n; ++j) {
      for (int i = 1; i <= n; ++i) {
        total += r[(std::size_t)(i + stride * j)];
      }
    }
  }
  CHECK(total > 0.0f);

  return failures;
}
