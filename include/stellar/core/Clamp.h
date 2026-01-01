#pragma once

#include <type_traits>

namespace stellar::core {

// Simple, dependency-light clamp helpers.
//
// We intentionally avoid std::clamp here because a few call sites clamp values across
// different integer widths (u8/u32/int) and some compilers can produce confusing
// overload-resolution diagnostics when the types don't match exactly.

template <class T>
constexpr T clamp(T v, T lo, T hi) {
  return (v < lo) ? lo : ((v > hi) ? hi : v);
}

// Clamp a value in a common type and cast to the requested output type.
// Useful for safely clamping ints into small integer fields.
template <class Out, class In, class Lo, class Hi>
constexpr Out clampCast(In v, Lo lo, Hi hi) {
  using C = std::common_type_t<In, Lo, Hi>;
  C cv = static_cast<C>(v);
  C clo = static_cast<C>(lo);
  C chi = static_cast<C>(hi);
  if (cv < clo) cv = clo;
  if (cv > chi) cv = chi;
  return static_cast<Out>(cv);
}

} // namespace stellar::core
