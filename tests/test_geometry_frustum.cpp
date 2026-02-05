#include "stellar/math/Frustum.h"

#include <cmath>
#include <iostream>

using namespace stellar;

int test_geometry_frustum() {
  int fails = 0;

  const double pi = 3.14159265358979323846;

  // 90 degree fov, symmetric frustum, standard OpenGL-style projection.
  const math::Mat4d proj = math::Mat4d::perspective(0.5 * pi, 1.0, 1.0, 10.0);
  const math::Frustumd fr = math::frustumFromClipMatrix(proj);

  if (!fr.isFinite()) {
    std::cerr << "[test_geometry_frustum] expected finite frustum\n";
    ++fails;
  }

  // Point containment sanity.
  if (!fr.containsPoint({0.0, 0.0, -5.0})) {
    std::cerr << "[test_geometry_frustum] expected point inside\n";
    ++fails;
  }
  if (fr.containsPoint({0.0, 0.0, -0.5})) {
    std::cerr << "[test_geometry_frustum] expected point in front of near plane to be outside\n";
    ++fails;
  }
  if (fr.containsPoint({0.0, 0.0, 5.0})) {
    std::cerr << "[test_geometry_frustum] expected point behind the camera to be outside\n";
    ++fails;
  }

  // Sphere/frustum tests.
  if (!fr.intersectsSphere({0.0, 0.0, -5.0}, 0.25)) {
    std::cerr << "[test_geometry_frustum] expected sphere inside\n";
    ++fails;
  }
  if (fr.intersectsSphere({0.0, 0.0, -0.5}, 0.1)) {
    std::cerr << "[test_geometry_frustum] expected small sphere before near plane to be outside\n";
    ++fails;
  }
  if (!fr.intersectsSphere({0.0, 0.0, -0.5}, 0.6)) {
    std::cerr << "[test_geometry_frustum] expected sphere straddling near plane to intersect\n";
    ++fails;
  }
  if (fr.intersectsSphere({0.0, 0.0, -11.0}, 0.25)) {
    std::cerr << "[test_geometry_frustum] expected sphere beyond far plane to be outside\n";
    ++fails;
  }

  // AABB/frustum tests.
  {
    const math::Aabb3d inside = math::Aabb3d::fromCenterExtents({0.0, 0.0, -5.0}, {0.5, 0.5, 0.5});
    if (!fr.intersectsAabb(inside)) {
      std::cerr << "[test_geometry_frustum] expected AABB inside\n";
      ++fails;
    }

    const math::Aabb3d outsideX = math::Aabb3d::fromCenterExtents({10.0, 0.0, -5.0}, {0.5, 0.5, 0.5});
    if (fr.intersectsAabb(outsideX)) {
      std::cerr << "[test_geometry_frustum] expected far-X AABB outside\n";
      ++fails;
    }

    const math::Aabb3d straddleNear = math::Aabb3d::fromMinMax({-0.5, -0.5, -1.5}, {0.5, 0.5, -0.2});
    if (!fr.intersectsAabb(straddleNear)) {
      std::cerr << "[test_geometry_frustum] expected AABB straddling near plane to intersect\n";
      ++fails;
    }
  }

  return fails;
}
