#include "CameraRigWindow.h"

#include "stellar/math/Math.h"
#include "stellar/math/Fov.h"
#include "stellar/sim/Gravity.h"
#include "stellar/sim/Ship.h"
#include "stellar/sim/System.h"

#include <imgui.h>

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <string>

namespace stellar::game {

void applyCameraRigPreset(CameraRigWindowState& st, CameraRigPreset preset, const ToastFn& toast) {
  st.lastPresetId = (int)preset;

  // Presets intentionally only touch a subset of knobs, so advanced users
  // can still tweak in the UI without constantly fighting the system.
  switch (preset) {
    case CameraRigPreset::DefaultOrbit: {
      st.mode = CameraRigMode::Orbit;
      st.orbitPitchDeg = 15.0f;
      st.orbitDistanceU = 6.0f;
      st.orbitPivotHeightU = 1.2f;
      st.orbitPivotShoulderU = 0.0f;
      st.orbitPanU = {0.0, 0.0, 0.0};
      st.orbitUseReferenceUp = true;
      st.orbitLookAhead = true;
      st.orbitLookAheadSec = 1.15f;
      st.orbitLookAheadMaxU = 2.0f;

      st.baseFovDeg = 60.0f;
      st.speedFov = true;
      st.speedFovGainDeg = 14.0f;
      st.speedFovMaxAddDeg = 18.0f;

      st.dollyZoomLock = false;
      st.dollyAutoCapture = true;
      st.dollyApplySpeedFov = false;
      break;
    }
    case CameraRigPreset::Travel: {
      // Wide, stable chase view.
      st.mode = CameraRigMode::Chase;
      st.chaseOffsetLocalU = {0.0, 2.4, -9.5};
      st.chaseLookAtShip = false;
      st.chaseUseReferenceUp = false;

      st.baseFovDeg = 72.0f;
      st.speedFov = true;
      st.speedFovGainDeg = 18.0f;
      st.speedFovMaxAddDeg = 24.0f;

      st.dollyZoomLock = false;
      break;
    }
    case CameraRigPreset::Docking: {
      // Orbit view with horizon lock and reduced look-ahead so the station approach feels steadier.
      st.mode = CameraRigMode::Orbit;
      st.orbitUseReferenceUp = true;
      st.orbitLookAhead = false;
      st.orbitDistanceU = 4.0f;
      st.orbitPivotHeightU = 1.1f;
      st.baseFovDeg = 52.0f;
      st.speedFov = false;
      st.dollyZoomLock = false;
      break;
    }
    case CameraRigPreset::Combat: {
      // Tight chase for situational awareness.
      st.mode = CameraRigMode::Chase;
      st.chaseOffsetLocalU = {0.65, 1.55, -5.0};
      st.chaseLookAtShip = false;
      st.chaseUseReferenceUp = false;

      st.baseFovDeg = 62.0f;
      st.speedFov = true;
      st.speedFovGainDeg = 10.0f;
      st.speedFovMaxAddDeg = 14.0f;
      st.dollyZoomLock = false;
      break;
    }
    case CameraRigPreset::Cinematic: {
      // Orbit + dolly zoom lock for "cinema" feeling while still allowing the player to fly.
      st.mode = CameraRigMode::Orbit;
      st.orbitUseReferenceUp = true;
      st.orbitLookAhead = true;
      st.orbitLookAheadSec = 0.75f;
      st.orbitLookAheadMaxU = 1.5f;
      st.orbitDistanceU = 7.5f;
      st.baseFovDeg = 50.0f;
      st.speedFov = false;

      st.dollyZoomLock = true;
      st.dollyAutoCapture = true;
      st.dollyApplySpeedFov = false;
      st.dollyMinFovDeg = 18.0f;
      st.dollyMaxFovDeg = 100.0f;
      break;
    }
    default:
      break;
  }

  // Reset smoothing so the new camera configuration doesn't "ease" from stale state.
  st.smoothValid = false;

  if (toast) {
    toast(std::string("Camera preset: ") + cameraRigPresetName(preset), 1.2);
  }
}

namespace {

constexpr double kLN2 = 0.693147180559945309417232121458176568;

static double halfLifeAlpha(double halfLifeSec, double dtSec) {
  if (halfLifeSec <= 1e-6 || dtSec <= 0.0) return 1.0;
  return 1.0 - std::exp(-kLN2 * dtSec / halfLifeSec);
}

static double quatDot(const math::Quatd& a, const math::Quatd& b) {
  return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
}

static double clampDeg(double v, double lo, double hi) {
  return std::max(lo, std::min(hi, v));
}

static math::Vec3d safeNormalized(const math::Vec3d& v, const math::Vec3d& fallback) {
  const double lsq = v.lengthSq();
  if (lsq <= 1e-18) return fallback;
  return v * (1.0 / std::sqrt(lsq));
}

static math::Quatd quatFromBasis(const math::Vec3d& right,
                                 const math::Vec3d& up,
                                 const math::Vec3d& forward) {
  // Columns are the world-space axes of the local basis.
  const double m00 = right.x;
  const double m01 = up.x;
  const double m02 = forward.x;

  const double m10 = right.y;
  const double m11 = up.y;
  const double m12 = forward.y;

  const double m20 = right.z;
  const double m21 = up.z;
  const double m22 = forward.z;

  const double trace = m00 + m11 + m22;
  double w = 1.0, x = 0.0, y = 0.0, z = 0.0;

  if (trace > 0.0) {
    const double s = std::sqrt(trace + 1.0) * 2.0;
    w = 0.25 * s;
    x = (m21 - m12) / s;
    y = (m02 - m20) / s;
    z = (m10 - m01) / s;
  } else if (m00 > m11 && m00 > m22) {
    const double s = std::sqrt(1.0 + m00 - m11 - m22) * 2.0;
    w = (m21 - m12) / s;
    x = 0.25 * s;
    y = (m01 + m10) / s;
    z = (m02 + m20) / s;
  } else if (m11 > m22) {
    const double s = std::sqrt(1.0 + m11 - m00 - m22) * 2.0;
    w = (m02 - m20) / s;
    x = (m01 + m10) / s;
    y = 0.25 * s;
    z = (m12 + m21) / s;
  } else {
    const double s = std::sqrt(1.0 + m22 - m00 - m11) * 2.0;
    w = (m10 - m01) / s;
    x = (m02 + m20) / s;
    y = (m12 + m21) / s;
    z = 0.25 * s;
  }

  return math::Quatd{w, x, y, z}.normalized();
}

static math::Quatd quatLookRotation(const math::Vec3d& forwardDir,
                                    const math::Vec3d& upHint) {
  const math::Vec3d f = safeNormalized(forwardDir, {0.0, 0.0, 1.0});
  math::Vec3d u = safeNormalized(upHint, {0.0, 1.0, 0.0});

  // If up is nearly parallel to forward, pick a fallback up.
  if (std::abs(math::dot(u, f)) > 0.999) {
    u = (std::abs(f.y) < 0.999) ? math::Vec3d{0.0, 1.0, 0.0} : math::Vec3d{1.0, 0.0, 0.0};
  }

  math::Vec3d r = math::cross(u, f);
  r = safeNormalized(r, {1.0, 0.0, 0.0});
  const math::Vec3d u2 = math::cross(f, r);

  return quatFromBasis(r, u2, f);
}

static math::Quatd expSmoothQuat(const math::Quatd& current,
                                 const math::Quatd& target,
                                 double alpha) {
  // Blend in axis-angle space so we get consistent angular response.
  if (alpha >= 1.0) return target;

  math::Quatd t = target;
  math::Quatd c = current;

  // Ensure shortest arc (avoid sudden flips when w crosses 0).
  if (quatDot(t, c) < 0.0) {
    t.w = -t.w;
    t.x = -t.x;
    t.y = -t.y;
    t.z = -t.z;
  }

  const math::Quatd dq = t * c.conjugate();

  double angle = 2.0 * std::acos(std::max(-1.0, std::min(1.0, dq.w)));
  if (angle < 1e-8) return t;

  // Map angle into [-pi, pi] for stability.
  if (angle > math::kPi) angle -= 2.0 * math::kPi;

  const double s = std::sqrt(std::max(0.0, 1.0 - dq.w * dq.w));
  math::Vec3d axis{1.0, 0.0, 0.0};
  if (s > 1e-8) {
    axis = math::Vec3d{dq.x / s, dq.y / s, dq.z / s};
  }

  const math::Quatd step = math::Quatd::fromAxisAngle(axis, angle * alpha);
  return (step * c).normalized();
}

static math::Vec3d computeReferenceUp(const sim::Ship& ship,
                                      const sim::StarSystem* system,
                                      double timeDays,
                                      const sim::GravityParams& gravityParams) {
  // Mirror the ship attitude instruments: prefer gravity up when meaningful,
  // otherwise fall back to orbit plane normal, then ship up.
  math::Vec3d refUp = ship.up();

  if (system) {
    const math::Vec3d gWorld = sim::systemGravityAccelKmS2(*system, timeDays, ship.positionKm(), gravityParams);
    const double gMag = gWorld.length();
    // Convert to ~g's for an intuitive threshold. 9.80665 m/s^2 = 0.00980665 km/s^2.
    const double gG = (gMag * 1000.0) / 9.80665;

    if (gMag > 1e-12 && gG > 0.0025) {
      // "Up" is opposite gravity.
      refUp = (-gWorld).normalized();
    } else {
      // Use orbit plane normal.
      const math::Vec3d h = math::cross(ship.positionKm(), ship.velocityKmS());
      if (h.lengthSq() > 1e-18) {
        refUp = h.normalized();
        // Flip so up generally points +Y to reduce sudden upside-down transitions.
        if (refUp.y < 0.0) refUp = -refUp;
      }
    }
  }

  return safeNormalized(refUp, {0.0, 1.0, 0.0});
}

static void applyBodyCollision(math::Vec3d& camPosU,
                               const math::Vec3d& shipPosU,
                               const sim::Ship& ship,
                               const sim::StarSystem* system,
                               double timeDays,
                               const sim::GravityParams& gravityParams,
                               double paddingU) {
  if (!system) return;

  const sim::DominantGravityBody body = sim::dominantGravityBody(*system, timeDays, ship.positionKm(), gravityParams);
	  if (!body.valid) return;
	  const auto& gb = body.body;

  // Convert body center to render-units using the same scale as the visual instancing.
  // NOTE: This code mirrors the instancing scales in main.cpp:
  //   star radius scale factor ~= 3.0
  //   planet/moon radius scale factor ~= 200.0
  constexpr double kRENDER_UNIT_KM = 1e6;
  auto toRenderPosU = [&](const math::Vec3d& pKm) {
    // Camera collision is local; "origin" shifting doesn't matter as long as
    // ship and bodies are expressed in the same render frame.
    // The main loop subtracts gRenderOriginKm before rendering; here we only
    // need relative positions, so shipPosU must already be in the same frame.
    // body.posKm is in the same universe frame as ship.positionKm.
    // We reconstruct the same shifting by anchoring on ship position.
    const math::Vec3d relKm = pKm - ship.positionKm();
    return shipPosU + (relKm * (1.0 / kRENDER_UNIT_KM));
  };

	  const math::Vec3d bodyPosU = toRenderPosU(gb.posKm);
	  const double scaleFactor = (gb.kind == sim::GravityBody::Kind::Star) ? 3.0 : 200.0;
	  const double bodyRadiusU = std::max(0.1, (gb.radiusKm / kRENDER_UNIT_KM) * scaleFactor);
  const double minDist = bodyRadiusU + paddingU;

  math::Vec3d d = camPosU - bodyPosU;
  const double dist = d.length();
  if (dist < minDist) {
    if (dist < 1e-9) {
      // Degenerate; push away from ship.
      d = safeNormalized(camPosU - shipPosU, {0.0, 1.0, 0.0});
    } else {
      d = d * (1.0 / dist);
    }
    camPosU = bodyPosU + d * minDist;
  }
}

static double dynamicNearPlaneU(double camDistU) {
  // Improve depth precision by pushing near plane out when the camera is far.
  // Clamp so we don't clip close HUD/cockpit work.
  const double nearBase = 0.01;
  const double nearByDist = camDistU * 0.02;
  return std::max(nearBase, std::min(0.5, nearByDist));
}

} // namespace

bool handleCameraRigEvent(CameraRigWindowState& st,
                          const SDL_Event& event,
                          const ImGuiIO& io,
                          bool mouseSteerActive,
                          bool actionWheelOpen,
                          const ToastFn& toast) {
  // Don't steal mouse if UI wants it, mouse-steer is active, or the action wheel is up.
  const bool allowMouse = !io.WantCaptureMouse && !mouseSteerActive && !actionWheelOpen;

  const SDL_Keymod mods = SDL_GetModState();
  const bool altDown = (mods & KMOD_ALT) != 0;
  const bool shiftDown = (mods & KMOD_SHIFT) != 0;

  // Hotkeys
  if (event.type == SDL_KEYDOWN && !io.WantCaptureKeyboard) {
    if (event.key.keysym.sym == SDLK_F6) {
      // Cycle between chase and orbit.
      st.mode = (st.mode == CameraRigMode::Chase) ? CameraRigMode::Orbit : CameraRigMode::Chase;
      st.smoothValid = false;
      if (toast) {
        toast(std::string("Camera mode: ") + ((st.mode == CameraRigMode::Orbit) ? "Orbit" : "Chase"), 1.4);
      }
      return true;
    }
    if (event.key.keysym.sym == SDLK_F7) {
      // Quick reset.
      st.orbitYawDeg = 0.0f;
      st.orbitPitchDeg = 15.0f;
      st.orbitDistanceU = 6.0f;
      st.orbitPanU = {0.0, 0.0, 0.0};
      st.smoothValid = false;
      if (toast) toast("Camera reset.", 1.2);
      return true;
    }
  }

  if (!st.enableAltMouseOrbit || !altDown || !allowMouse) {
    // If orbit controls aren't active, make sure we stop dragging.
    if (event.type == SDL_MOUSEBUTTONUP) {
      if (event.button.button == SDL_BUTTON_LEFT) st.draggingOrbit = false;
      if (event.button.button == SDL_BUTTON_MIDDLE) st.draggingPan = false;
    }
    return false;
  }

  // Orbit controls
  switch (event.type) {
    case SDL_MOUSEBUTTONDOWN: {
      if (event.button.button == SDL_BUTTON_LEFT) {
        st.draggingOrbit = true;
        st.lastMouseX = event.button.x;
        st.lastMouseY = event.button.y;
        return true;
      }
      if (event.button.button == SDL_BUTTON_MIDDLE) {
        st.draggingPan = true;
        st.lastMouseX = event.button.x;
        st.lastMouseY = event.button.y;
        return true;
      }
    } break;

    case SDL_MOUSEBUTTONUP: {
      if (event.button.button == SDL_BUTTON_LEFT) {
        st.draggingOrbit = false;
        return true;
      }
      if (event.button.button == SDL_BUTTON_MIDDLE) {
        st.draggingPan = false;
        return true;
      }
    } break;

    case SDL_MOUSEMOTION: {
      const int dx = event.motion.xrel;
      const int dy = event.motion.yrel;
      if (st.draggingOrbit) {
        const float sx = st.mouseSensitivityDegPerPx;
        st.orbitYawDeg += (float)dx * sx;
        const float dySigned = st.invertY ? (float)dy : (float)-dy;
        st.orbitPitchDeg += dySigned * sx;
        st.orbitPitchDeg = (float)clampDeg(st.orbitPitchDeg, -89.0, 89.0);
        return true;
      }
      if (st.draggingPan) {
        const float panScale = 0.01f * st.orbitDistanceU;
        // Pan in reference-space right/up.
        st.orbitPanU.x += (double)(-dx) * (double)panScale;
        st.orbitPanU.y += (double)(dy) * (double)panScale;
        return true;
      }
    } break;

    case SDL_MOUSEWHEEL: {
      if (shiftDown) {
        // Fine zoom.
        st.orbitDistanceU *= (float)std::exp(-(double)event.wheel.y * (double)(st.scrollZoomSpeed * 0.25f));
      } else {
        st.orbitDistanceU *= (float)std::exp(-(double)event.wheel.y * (double)st.scrollZoomSpeed);
      }
      st.orbitDistanceU = (float)std::max((double)st.minCamDistU, std::min(800.0, (double)st.orbitDistanceU));
      return true;
    } break;

    default: break;
  }

  return false;
}

void tickCameraRig(CameraRigWindowState& st, double /*dtRealSec*/, const Uint8* /*keys*/, const ImGuiIO& /*io*/) {
  // Placeholder for future continuous behaviors (e.g., camera shake, auto-centering).
  // Kept as a function to preserve a clean integration point.
  (void)st;
}

CameraRigFrame computeCameraRigFrame(CameraRigWindowState& st,
                                     const sim::Ship& ship,
                                     const math::Vec3d& shipPosU,
                                     const sim::StarSystem* system,
                                     double timeDays,
                                     const sim::GravityParams& gravityParams,
                                     double dtRealSec,
                                     bool cinematicEnabled,
                                     const math::Vec3d& cinematicOffsetLocalU,
                                     bool cinematicLookAtShip,
                                     double cinematicFovDeg,
                                     bool cinematicOverrideFov) {
  CameraRigFrame out{};

  // Reference up is used for horizon-locking in orbit (and optionally chase).
  const math::Vec3d refUp = computeReferenceUp(ship, system, timeDays, gravityParams);
  out.refUp = refUp;

  // Compute ship basis in render space (unit vectors).
  const math::Vec3d shipF = safeNormalized(ship.forward(), {0.0, 0.0, 1.0});
  const math::Vec3d shipR = safeNormalized(ship.right(), {1.0, 0.0, 0.0});
  const math::Vec3d shipU = safeNormalized(ship.up(), {0.0, 1.0, 0.0});

  // Determine the driver: either cinematic camera (if enabled) or the active mode.
  const bool useCinematic = cinematicEnabled && st.followCinematicWhenEnabled;

  // Desired pose (unsmoothed).
  math::Vec3d desiredPosU{0.0, 0.0, 0.0};
  math::Quatd desiredOrient = math::Quatd::identity();
  double desiredFovDeg = (double)st.baseFovDeg;

  // Look-ahead target point (in render units).
  math::Vec3d lookTargetU = shipPosU;

  if (useCinematic) {
    // Ship-local offset.
    const math::Vec3d off = cinematicOffsetLocalU;
    desiredPosU = shipPosU + shipR * off.x + shipU * off.y + shipF * off.z;

    const bool lookAt = cinematicLookAtShip;
    if (lookAt) {
      desiredOrient = quatLookRotation(shipPosU - desiredPosU, refUp);
    } else {
      desiredOrient = ship.orientation();
    }

    if (cinematicOverrideFov) desiredFovDeg = cinematicFovDeg;
  } else if (st.mode == CameraRigMode::Chase) {
    const math::Vec3d off = st.chaseOffsetLocalU;
    desiredPosU = shipPosU + shipR * off.x + shipU * off.y + shipF * off.z;
    if (st.chaseLookAtShip) {
      desiredOrient = quatLookRotation(shipPosU - desiredPosU, st.chaseUseReferenceUp ? refUp : shipU);
    } else {
      desiredOrient = ship.orientation();
    }
    desiredFovDeg = (double)st.baseFovDeg;
  } else {
    // Orbit mode (horizon locked).
    // Forward reference is ship forward projected onto the ref-up plane.
    math::Vec3d fFlat = shipF - refUp * math::dot(shipF, refUp);
    fFlat = safeNormalized(fFlat, shipF);
    math::Vec3d rRef = safeNormalized(math::cross(refUp, fFlat), shipR);

    // Pivot point.
    const math::Vec3d pivotU = shipPosU + refUp * (double)st.orbitPivotHeightU + rRef * (double)st.orbitPivotShoulderU +
                               rRef * st.orbitPanU.x + refUp * st.orbitPanU.y + fFlat * st.orbitPanU.z;

    const double yawRad = math::degToRad((double)st.orbitYawDeg);
    const double pitchRad = math::degToRad((double)st.orbitPitchDeg);

    // Start behind the ship (negative forward).
    math::Vec3d dir = -fFlat;
    dir = math::Quatd::fromAxisAngle(refUp, yawRad).rotate(dir);
    const math::Vec3d pitchAxis = safeNormalized(math::cross(refUp, dir), rRef);
    dir = math::Quatd::fromAxisAngle(pitchAxis, pitchRad).rotate(dir);
    dir = safeNormalized(dir, -fFlat);

    desiredPosU = pivotU + dir * (double)st.orbitDistanceU;

    // Look target with optional look-ahead.
    lookTargetU = pivotU;
    if (st.orbitLookAhead) {
      // Convert km/s -> render units/s. kRENDER_UNIT_KM = 1e6.
      constexpr double kINV_RENDER_UNIT_KM = 1.0 / 1e6;
      const math::Vec3d velU = ship.velocityKmS() * kINV_RENDER_UNIT_KM;
      math::Vec3d lead = velU * (double)st.orbitLookAheadSec;
      const double leadLen = lead.length();
      if (leadLen > (double)st.orbitLookAheadMaxU) {
        lead = lead * ((double)st.orbitLookAheadMaxU / leadLen);
      }
      lookTargetU = lookTargetU + lead;
    }

    desiredOrient = quatLookRotation(lookTargetU - desiredPosU, st.orbitUseReferenceUp ? refUp : shipU);
    desiredFovDeg = (double)st.baseFovDeg;
  }

  // Estimate the camera distance after applying hard constraints (body avoidance + minimum distance).
  // This makes distance-dependent FOV tools stable near large bodies and docking corridors.
  double constrainedDistU = (shipPosU - desiredPosU).length();
  {
    math::Vec3d cposU = desiredPosU;
    if (st.avoidBodies) {
      applyBodyCollision(cposU, shipPosU, ship, system, timeDays, gravityParams, (double)st.bodyPaddingU);
    }

    const math::Vec3d toShipC = shipPosU - cposU;
    const double d = toShipC.length();
    if (d < (double)st.minCamDistU) {
      const math::Vec3d dir = safeNormalized(cposU - shipPosU, {0.0, 0.0, 1.0});
      cposU = shipPosU + dir * (double)st.minCamDistU;
    }
    constrainedDistU = (shipPosU - cposU).length();
  }

  const auto addSpeedFovKick = [&]() {
    constexpr double kINV_RENDER_UNIT_KM = 1.0 / 1e6;
    const double speedU = ship.velocityKmS().length() * kINV_RENDER_UNIT_KM;
    const double add = std::min((double)st.speedFovMaxAddDeg, speedU * (double)st.speedFovGainDeg);
    desiredFovDeg += add;
  };

  // ---- Distance-aware FOV tools ----
  // Dolly zoom lock keeps framing roughly constant while changing camera distance.
  if (st.dollyZoomLock) {
    desiredFovDeg = math::fov::dollyZoomFovDeg((double)st.dollyRefFovDeg,
                                              (double)st.dollyRefDistanceU,
                                              constrainedDistU);
    desiredFovDeg = clampDeg(desiredFovDeg, (double)st.dollyMinFovDeg, (double)st.dollyMaxFovDeg);

    if (st.dollyApplySpeedFov && st.speedFov) {
      addSpeedFovKick();
    }
  } else if (st.speedFov) {
    addSpeedFovKick();
  }

  // Smooth
  if (!st.smoothValid) {
    st.smoothPosU = desiredPosU;
    st.smoothOrient = desiredOrient;
    st.smoothFovDeg = desiredFovDeg;
    st.smoothValid = true;
  } else {
    // Clamp dt to avoid huge camera jumps after a breakpoint / tab-out.
    const double dt = std::clamp(dtRealSec, 0.0, 0.25);

    const double alphaPos = halfLifeAlpha((double)st.posHalfLifeSec, dt);
    const double alphaRot = halfLifeAlpha((double)st.rotHalfLifeSec, dt);
    const double alphaFov = halfLifeAlpha((double)st.fovHalfLifeSec, dt);

    st.smoothPosU = st.smoothPosU * (1.0 - alphaPos) + desiredPosU * alphaPos;
    st.smoothOrient = expSmoothQuat(st.smoothOrient, desiredOrient, alphaRot);
    st.smoothFovDeg = st.smoothFovDeg * (1.0 - alphaFov) + desiredFovDeg * alphaFov;
  }

  math::Vec3d finalPosU = st.smoothPosU;
  if (st.avoidBodies) {
    applyBodyCollision(finalPosU, shipPosU, ship, system, timeDays, gravityParams, (double)st.bodyPaddingU);
  }

  // Ensure minimum distance to the ship.
  const math::Vec3d toShip = shipPosU - finalPosU;
  const double distToShip = toShip.length();
  if (distToShip < (double)st.minCamDistU) {
    const math::Vec3d dir = safeNormalized(finalPosU - shipPosU, {0.0, 0.0, 1.0});
    finalPosU = shipPosU + dir * (double)st.minCamDistU;
  }

  // Re-aim after collision/clamp for orbit/look-at modes.
  math::Quatd finalOrient = st.smoothOrient;
  if (useCinematic) {
    if (cinematicLookAtShip) {
      finalOrient = quatLookRotation(shipPosU - finalPosU, refUp);
    }
  } else if (st.mode == CameraRigMode::Chase) {
    if (st.chaseLookAtShip) {
      finalOrient = quatLookRotation(shipPosU - finalPosU, st.chaseUseReferenceUp ? refUp : shipU);
    }
  } else {
    // Orbit
    finalOrient = quatLookRotation(lookTargetU - finalPosU, st.orbitUseReferenceUp ? refUp : shipU);
  }

  out.posU = finalPosU;
  out.orient = finalOrient;
  out.fovDeg = st.smoothFovDeg;

  const double camDist = (shipPosU - finalPosU).length();
  out.nearPlane = dynamicNearPlaneU(camDist);
  out.farPlane = 20000.0;

  // ---- Tooling telemetry ----
  st.lastCamDistU = camDist;
  st.lastCamDistKm = camDist * 1.0e6;
  st.lastFinalFovDeg = out.fovDeg;
  st.lastCamOffsetU = finalPosU - shipPosU;

  // Write back final for continuity.
  st.smoothPosU = finalPosU;
  st.smoothOrient = finalOrient;
  st.smoothFovDeg = out.fovDeg;

  return out;
}

void cameraRigFocusOrbitPivot(CameraRigWindowState& st,
                              const sim::Ship& ship,
                              const sim::StarSystem* system,
                              double timeDays,
                              const sim::GravityParams& gravityParams,
                              const math::Vec3d& shipPosU,
                              const math::Vec3d& targetPosU) {
  // Ensure orbit mode so the pan offset has meaning.
  st.mode = CameraRigMode::Orbit;

  // Recreate the same basis used by computeCameraRigFrame so the pan axes feel
  // consistent with the active horizon lock.
  const math::Vec3d refUp = computeReferenceUp(ship, system, timeDays, gravityParams);

  const math::Vec3d shipF = safeNormalized(ship.forward(), {0.0, 0.0, 1.0});
  const math::Vec3d shipR = safeNormalized(ship.right(), {1.0, 0.0, 0.0});

  math::Vec3d fFlat = shipF - refUp * math::dot(shipF, refUp);
  fFlat = safeNormalized(fFlat, shipF);
  const math::Vec3d rRef = safeNormalized(math::cross(refUp, fFlat), shipR);

  // Base pivot before panning.
  const math::Vec3d basePivotU = shipPosU + refUp * (double)st.orbitPivotHeightU +
                                rRef * (double)st.orbitPivotShoulderU;

  // Solve for orbitPanU in the (rRef, refUp, fFlat) basis.
  const math::Vec3d delta = targetPosU - basePivotU;
  st.orbitPanU = {
      math::dot(delta, rRef),
      math::dot(delta, refUp),
      math::dot(delta, fFlat),
  };
}

void drawCameraRigWindow(CameraRigWindowState& st,
                         const sim::Ship& ship,
                         const sim::StarSystem* system,
                         double timeDays,
                         const sim::GravityParams& gravityParams,
                         const ToastFn& toast) {
  if (!st.open) return;

  if (ImGui::Begin("3D Camera Rig", &st.open)) {
    ImGui::TextUnformatted("Alt + LMB: orbit\nAlt + MMB: pan\nAlt + wheel: zoom\nF6: toggle Orbit/Chase\nF7: reset orbit");
    ImGui::Separator();

    int mode = (int)st.mode;
    const char* modeLabels[] = {"Chase", "Orbit"};
    if (ImGui::Combo("Mode", &mode, modeLabels, IM_ARRAYSIZE(modeLabels))) {
      st.mode = (mode == 0) ? CameraRigMode::Chase : CameraRigMode::Orbit;
      st.smoothValid = false;
    }
    ImGui::Checkbox("Follow cinematic rig when enabled", &st.followCinematicWhenEnabled);

    ImGui::Separator();
    ImGui::TextUnformatted("Presets");
    ImGui::TextDisabled("Quick camera configurations (also usable from Integration Hub automations).");

    auto presetBtn = [&](const char* label, CameraRigPreset preset) {
      const bool active = (st.lastPresetId == (int)preset);
      if (active) ImGui::PushStyleVar(ImGuiStyleVar_Alpha, 0.85f);
      if (ImGui::SmallButton(label)) {
        applyCameraRigPreset(st, preset, toast);
      }
      if (active) ImGui::PopStyleVar();
    };

    presetBtn("Default",  CameraRigPreset::DefaultOrbit);
    ImGui::SameLine(); presetBtn("Travel",   CameraRigPreset::Travel);
    ImGui::SameLine(); presetBtn("Docking",  CameraRigPreset::Docking);
    ImGui::SameLine(); presetBtn("Combat",   CameraRigPreset::Combat);
    ImGui::SameLine(); presetBtn("Cinematic", CameraRigPreset::Cinematic);

    if (st.lastPresetId >= 0) {
      const int pid = std::clamp(st.lastPresetId, 0, 4);
      ImGui::TextDisabled("Last preset: %s", cameraRigPresetName((CameraRigPreset)pid));
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Input");
    ImGui::Checkbox("Enable Alt-mouse orbit", &st.enableAltMouseOrbit);
    ImGui::SliderFloat("Mouse sensitivity (deg/px)", &st.mouseSensitivityDegPerPx, 0.02f, 0.60f, "%.3f");
    ImGui::Checkbox("Invert Y", &st.invertY);
    ImGui::SliderFloat("Zoom speed", &st.scrollZoomSpeed, 0.05f, 1.2f, "%.2f");

    ImGui::Separator();
    if (st.mode == CameraRigMode::Orbit) {
      ImGui::TextUnformatted("Orbit");
      ImGui::SliderFloat("Yaw (deg)", &st.orbitYawDeg, -180.0f, 180.0f, "%.1f");
      ImGui::SliderFloat("Pitch (deg)", &st.orbitPitchDeg, -89.0f, 89.0f, "%.1f");
      ImGui::SliderFloat("Distance (U)", &st.orbitDistanceU, 0.75f, 200.0f, "%.2f");
      ImGui::SliderFloat("Pivot height (U)", &st.orbitPivotHeightU, -10.0f, 10.0f, "%.2f");
      ImGui::SliderFloat("Pivot shoulder (U)", &st.orbitPivotShoulderU, -10.0f, 10.0f, "%.2f");
      ImGui::Checkbox("Use reference up (horizon lock)", &st.orbitUseReferenceUp);
      ImGui::Checkbox("Look-ahead", &st.orbitLookAhead);
      if (st.orbitLookAhead) {
        ImGui::SliderFloat("Look-ahead sec", &st.orbitLookAheadSec, 0.0f, 5.0f, "%.2f");
        ImGui::SliderFloat("Look-ahead max (U)", &st.orbitLookAheadMaxU, 0.0f, 25.0f, "%.2f");
      }
      if (ImGui::Button("Reset orbit")) {
        st.orbitYawDeg = 0.0f;
        st.orbitPitchDeg = 15.0f;
        st.orbitDistanceU = 6.0f;
        st.orbitPanU = {0.0, 0.0, 0.0};
        st.smoothValid = false;
      }
    } else {
      ImGui::TextUnformatted("Chase");
      ImGui::DragScalarN("Offset (ship-local U)", ImGuiDataType_Double, &st.chaseOffsetLocalU.x, 3, 0.05, nullptr, nullptr, "%.3f");
      ImGui::Checkbox("Look-at ship", &st.chaseLookAtShip);
      ImGui::Checkbox("Use reference up", &st.chaseUseReferenceUp);
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Smoothing");
    ImGui::SliderFloat("Pos half-life (s)", &st.posHalfLifeSec, 0.0f, 2.0f, "%.2f");
    ImGui::SliderFloat("Rot half-life (s)", &st.rotHalfLifeSec, 0.0f, 2.0f, "%.2f");
    ImGui::SliderFloat("FOV half-life (s)", &st.fovHalfLifeSec, 0.0f, 2.0f, "%.2f");
    if (ImGui::Button("Reset smoothing")) {
      st.smoothValid = false;
    }

    ImGui::Separator();
    ImGui::TextUnformatted("FOV");
    ImGui::SliderFloat("Base FOV", &st.baseFovDeg, 30.0f, 110.0f, "%.1f deg");
    ImGui::Checkbox("Speed FOV", &st.speedFov);
    if (st.speedFov) {
      ImGui::SliderFloat("Speed gain", &st.speedFovGainDeg, 0.0f, 60.0f, "%.1f deg/(U/s)");
      ImGui::SliderFloat("Max add", &st.speedFovMaxAddDeg, 0.0f, 60.0f, "%.1f deg");
    }

    ImGui::Separator();
    if (ImGui::CollapsingHeader("Distance & FOV tools", ImGuiTreeNodeFlags_DefaultOpen)) {
      constexpr double kRENDER_UNIT_KM = 1.0e6;

      ImGui::Checkbox("Show distance telemetry", &st.showDistanceTelemetry);
      if (st.showDistanceTelemetry) {
        const ImGuiIO& io = ImGui::GetIO();
        const double vpW = (double)io.DisplaySize.x;
        const double vpH = (double)io.DisplaySize.y;
        const double aspect = (vpH > 1.0) ? (vpW / vpH) : 1.0;

        const double fovV = st.lastFinalFovDeg;
        const double fovH = math::fov::verticalToHorizontalDeg(fovV, aspect);
        ImGui::Text("Viewport: %.0f x %.0f  (aspect %.3f)", vpW, vpH, aspect);
        ImGui::Text("Final FOV: %.2f deg V / %.2f deg H", fovV, fovH);
        ImGui::Text("Camera distance: %.3f U  (%.0f km)", st.lastCamDistU, st.lastCamDistKm);

        const double fovRad = math::degToRad(fovV);
        const double viewH = math::fov::viewHeightAtDistance(st.lastCamDistU, fovRad);
        const double upp = math::fov::unitsPerPixelAtDistance(st.lastCamDistU, fovRad, vpH);
        ImGui::Text("View height @ depth: %.6g U", viewH);
        ImGui::Text("Scale @ depth: %.6g U/px (center)", upp);
      }

      ImGui::Separator();
      ImGui::TextUnformatted("Dolly zoom (framing lock)");
      bool dollyChanged = ImGui::Checkbox("Enable dolly zoom lock", &st.dollyZoomLock);
      ImGui::SameLine();
      ImGui::Checkbox("Auto-capture", &st.dollyAutoCapture);

      if (dollyChanged && st.dollyZoomLock && st.dollyAutoCapture) {
        st.dollyRefDistanceU = (float)std::max(0.001, st.lastCamDistU);
        st.dollyRefFovDeg = (float)st.lastFinalFovDeg;
        if (toast) toast("Dolly zoom reference captured (distance + FOV).", 2.0);
      }

      if (st.dollyZoomLock) {
        if (ImGui::Button("Capture current framing")) {
          st.dollyRefDistanceU = (float)std::max(0.001, st.lastCamDistU);
          st.dollyRefFovDeg = (float)st.lastFinalFovDeg;
          if (toast) toast("Dolly zoom reference captured.", 2.0);
        }
        ImGui::SameLine();
        ImGui::Checkbox("Apply speed FOV on top", &st.dollyApplySpeedFov);

        ImGui::DragFloat("Ref distance (U)", &st.dollyRefDistanceU, 0.05f, 0.01f, 2000.0f, "%.3f");
        ImGui::DragFloat("Ref FOV (deg)", &st.dollyRefFovDeg, 0.10f, 5.0f, 140.0f, "%.2f");
        ImGui::SliderFloat("Min FOV (deg)", &st.dollyMinFovDeg, 5.0f, 120.0f, "%.1f");
        ImGui::SliderFloat("Max FOV (deg)", &st.dollyMaxFovDeg, 10.0f, 160.0f, "%.1f");

        const double predicted = math::fov::dollyZoomFovDeg((double)st.dollyRefFovDeg,
                                                          (double)st.dollyRefDistanceU,
                                                          st.lastCamDistU);
        ImGui::Text("Locked FOV @ current distance: %.2f deg", predicted);
      }

      ImGui::Separator();
      ImGui::TextUnformatted("Physical lens equivalence");
      ImGui::DragFloat("Sensor width (mm)", &st.sensorWidthMm, 0.1f, 1.0f, 100.0f, "%.1f");
      ImGui::DragFloat("Sensor height (mm)", &st.sensorHeightMm, 0.1f, 1.0f, 100.0f, "%.1f");
      ImGui::DragFloat("Focal length (mm)", &st.focalLengthMm, 0.1f, 1.0f, 500.0f, "%.1f");

      const double lensVFovDeg = math::radToDeg(
          math::fov::fovRadFromFocalLengthMm((double)st.focalLengthMm, (double)st.sensorHeightMm));
      const double lensHFovDeg = math::radToDeg(
          math::fov::fovRadFromFocalLengthMm((double)st.focalLengthMm, (double)st.sensorWidthMm));
      ImGui::Text("Lens FOV: %.2f deg V / %.2f deg H", lensVFovDeg, lensHFovDeg);

      const double eqFocalMm = math::fov::focalLengthMmFromFovRad(
          math::degToRad(st.lastFinalFovDeg), (double)st.sensorHeightMm);
      ImGui::Text("Current FOV ~= %.1fmm focal (vertical, sensor %.1fmm)", eqFocalMm, st.sensorHeightMm);

      if (ImGui::Button("Apply lens to base FOV (vertical)")) {
        st.baseFovDeg = (float)math::clamp(lensVFovDeg, 5.0, 140.0);
        st.smoothValid = false;
        if (toast) toast("Base FOV set from lens vertical FOV.", 2.0);
      }

      // --- Dominant body framing helper ---
      if (system) {
        const math::Vec3d camPosKm = ship.positionKm() + st.lastCamOffsetU * kRENDER_UNIT_KM;
        const sim::DominantGravityBody dom = sim::dominantGravityBody(*system, timeDays, camPosKm, gravityParams);
        if (dom.valid) {
          const auto& b = dom.body;
          const char* kind = (b.kind == sim::GravityBody::Kind::Star) ? "Star" : "Planet";

          const double distCamKm = (camPosKm - b.posKm).length();
          const double distShipKm = (ship.positionKm() - b.posKm).length();
          const double altKm = distShipKm - b.radiusKm;

          const double angDeg = math::fov::angularDiameterDeg(b.radiusKm, distCamKm);

          ImGui::Separator();
          ImGui::TextUnformatted("Dominant body framing");
          if (!b.name.empty()) {
            ImGui::Text("Body: %s (%s)", b.name.c_str(), kind);
          } else {
            ImGui::Text("Body: [%s]", kind);
          }
          ImGui::Text("Ship altitude: %.0f km", altKm);
          ImGui::Text("Cam distance: %.0f km  |  Angular diameter: %.2f deg", distCamKm, angDeg);

          ImGui::SliderFloat("Frame fill", &st.frameBodyFill, 0.10f, 1.00f, "%.2f");

          if (ImGui::Button("Set base FOV to frame body")) {
            const double fill = std::clamp((double)st.frameBodyFill, 0.10, 1.0);
            const double newFov = angDeg / fill;
            st.baseFovDeg = (float)math::clamp(newFov, 5.0, 140.0);
            st.smoothValid = false;
            if (toast) toast("Base FOV set to frame dominant body.", 2.0);
          }
        }
      }
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Collision");
    ImGui::Checkbox("Avoid bodies", &st.avoidBodies);
    ImGui::SliderFloat("Body padding (U)", &st.bodyPaddingU, 0.0f, 2.0f, "%.3f");
    ImGui::SliderFloat("Min cam dist (U)", &st.minCamDistU, 0.1f, 50.0f, "%.2f");
  }
  ImGui::End();
}

} // namespace stellar::game
