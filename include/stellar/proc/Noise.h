#pragma once

#include "stellar/core/Types.h"

namespace stellar::proc {

// 2D/3D value noise helpers for procedural fields.
// This is intentionally tiny; replace with Perlin/Simplex later if needed.

double valueNoise2D(core::u64 seed, int x, int y);
double valueNoise3D(core::u64 seed, int x, int y, int z);

// Smooth value noise sampled at fractional coordinates.
double smoothNoise2D(core::u64 seed, double x, double y);
double smoothNoise3D(core::u64 seed, double x, double y, double z);

// Fractional Brownian motion (fBm) using smoothNoise2D.
double fbm2D(core::u64 seed, double x, double y, int octaves = 5, double lacunarity = 2.0, double gain = 0.5);

// Fractional Brownian motion (fBm) using smoothNoise3D.
double fbm3D(core::u64 seed, double x, double y, double z, int octaves = 5, double lacunarity = 2.0, double gain = 0.5);

// -----------------------------------------------------------------------------
// Gradient noise (Perlin-style) + derived fractals.
//
// These helpers exist alongside the legacy value/smooth noise so existing
// procedural content remains stable. Newer generators can opt into gradient
// noise for fewer grid artifacts and better “flow”.
//
// Notes:
//  - perlin*() returns a value in [0,1].
//  - fbmPerlin*() returns a sum in [0, sumAmp] with amp0=0.5.
//  - ridgedFbmPerlin*() returns a sum in [0, sumAmp] where each octave is a
//    ridge transform of the underlying Perlin noise.

double perlin2D(core::u64 seed, double x, double y);
double perlin3D(core::u64 seed, double x, double y, double z);

double fbmPerlin2D(core::u64 seed,
                   double x, double y,
                   int octaves = 5,
                   double lacunarity = 2.0,
                   double gain = 0.5);

double fbmPerlin3D(core::u64 seed,
                   double x, double y, double z,
                   int octaves = 5,
                   double lacunarity = 2.0,
                   double gain = 0.5);

double ridgedFbmPerlin2D(core::u64 seed,
                         double x, double y,
                         int octaves = 5,
                         double lacunarity = 2.0,
                         double gain = 0.5);

double ridgedFbmPerlin3D(core::u64 seed,
                         double x, double y, double z,
                         int octaves = 5,
                         double lacunarity = 2.0,
                         double gain = 0.5);

} // namespace stellar::proc
