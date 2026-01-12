#include "stellar/render/FrameGraph.h"

#include "test_harness.h"

#include <iostream>
#include <string>

// NOTE:
// This file intentionally avoids third-party test frameworks so headless builds
// (and Windows users building from scratch) do not need extra dependencies.

int test_frame_graph() {
  using namespace stellar::render;

  int failures = 0;

  // ---------------------------------------------------------------------------
  // Topological sort respects last-writer dependencies.
  // ---------------------------------------------------------------------------
  {
    FrameGraph fg;
    fg.setViewport(1280, 720);

    FrameGraph::TextureDesc d;
    d.scale = 1.0f;
    d.internalFormat = 1;
    d.format = 1;
    d.type = 1;

    const auto a = fg.createTexture("A", d);
    const auto b = fg.createTexture("B", d);
    const auto c = fg.createTexture("C", d);

    // Add in intentionally "wrong" order; compile should topo-sort.
    fg.addPass("P2", {b}, c, [](const FrameGraph::PassContext&) {});
    fg.addPass("P1", {a}, b, [](const FrameGraph::PassContext&) {});
    fg.addPass("P0", {}, a, [](const FrameGraph::PassContext&) {});

    std::string err;
    const bool ok = fg.compile(&err);
    if (!ok) {
      std::cerr << "[test_frame_graph] compile failed: " << err << "\n";
      ++failures;
    } else {
      CHECK(err.empty());
      CHECK(fg.schedule().size() == 3);
      const auto& order = fg.schedule();
      // Expect P0 -> P1 -> P2
      if (order.size() == 3) {
        CHECK(fg.passes()[(std::size_t)order[0]].name == "P0");
        CHECK(fg.passes()[(std::size_t)order[1]].name == "P1");
        CHECK(fg.passes()[(std::size_t)order[2]].name == "P2");
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Transient aliasing: reuses physical textures when lifetimes do not overlap.
  // ---------------------------------------------------------------------------
  {
    FrameGraph fg;
    fg.setViewport(800, 600);

    FrameGraph::TextureDesc d;
    d.scale = 0.5f;
    d.internalFormat = 7;
    d.format = 7;
    d.type = 7;

    // Chain of full-screen passes: each writes a new intermediate.
    // This should need only ~2 physical textures due to ping-pong constraints.
    const int N = 10;
    FrameGraph::TextureHandle prev = fg.createTexture("T0", d);
    fg.addPass("Write0", {}, prev, [](const FrameGraph::PassContext&) {});

    for (int i = 1; i < N; ++i) {
      const auto out = fg.createTexture("T" + std::to_string(i), d);
      fg.addPass("Pass" + std::to_string(i), {prev}, out, [](const FrameGraph::PassContext&) {});
      prev = out;
    }

    const bool ok = fg.compile(nullptr);
    if (!ok) {
      std::cerr << "[test_frame_graph] compile failed in aliasing test\n";
      ++failures;
    } else {
      // With inclusive lifetimes at a pass boundary (read+write), the chain cannot
      // collapse to 1 physical texture, but it should be much less than N.
      CHECK(fg.physicalTextureCount() <= 3);
    }
  }

  // ---------------------------------------------------------------------------
  // Debug: disabling aliasing forces unique physical textures.
  // ---------------------------------------------------------------------------
  {
    FrameGraph fg;
    fg.setViewport(640, 360);
    fg.setAliasingEnabled(false);

    FrameGraph::TextureDesc d;
    d.scale = 1.0f;
    d.internalFormat = 2;
    d.format = 2;
    d.type = 2;

    const auto t0 = fg.createTexture("t0", d);
    const auto t1 = fg.createTexture("t1", d);
    const auto t2 = fg.createTexture("t2", d);

    fg.addPass("p0", {}, t0, [](const FrameGraph::PassContext&) {});
    fg.addPass("p1", {t0}, t1, [](const FrameGraph::PassContext&) {});
    fg.addPass("p2", {t1}, t2, [](const FrameGraph::PassContext&) {});

    const bool ok = fg.compile(nullptr);
    if (!ok) {
      std::cerr << "[test_frame_graph] compile failed in no-aliasing test\n";
      ++failures;
    } else {
      CHECK(fg.physicalTextureCount() == 3);
    }
  }

  if (failures == 0) {
    std::cout << "[test_frame_graph] PASS\n";
  }
  return failures;
}
