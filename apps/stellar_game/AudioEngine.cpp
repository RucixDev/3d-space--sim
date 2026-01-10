#include "AudioEngine.h"

#include <SDL.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace stellar::game {
namespace {

static float clamp01(float v) { return std::clamp(v, 0.0f, 1.0f); }

static float clamp11(float v) { return std::clamp(v, -1.0f, 1.0f); }

static float lerp(float a, float b, float t) { return a + (b - a) * t; }

// Quick-and-cheap white noise in [-1, 1].
static float noise01(std::uint32_t& state) {
  state = state * 1664525u + 1013904223u;
  const std::uint32_t v = (state >> 9) | 0x3f800000u; // 1.x float bits
  float f;
  std::memcpy(&f, &v, sizeof(float));
  return (f - 1.5f) * 2.0f; // ~[-1,1]
}

// Constant-power panning (roughly preserves perceived loudness across pan).
static void panGains(float pan, float& gL, float& gR) {
  const float p = clamp11(pan);
  const float a = 0.5f * (p + 1.0f); // 0..1
  const float angle = a * 1.57079632679f; // 0..pi/2
  gL = std::cos(angle);
  gR = std::sin(angle);
}

static float twoPi() { return 6.28318530718f; }

} // namespace

AudioEngine::AudioEngine() = default;

AudioEngine::~AudioEngine() {
  shutdown();
}

bool AudioEngine::init() {
  if (active()) return true;

  if ((SDL_WasInit(SDL_INIT_AUDIO) & SDL_INIT_AUDIO) == 0) {
    // SDL audio subsystem must be initialized by the app (SDL_Init).
    return false;
  }

  SDL_AudioSpec desired{};
  SDL_zero(desired);
  desired.freq = 48000;
  desired.format = AUDIO_F32SYS;
  desired.channels = 2;
  desired.samples = 1024;
  desired.callback = &AudioEngine::sdlCallback;
  desired.userdata = this;

  SDL_AudioSpec obtained{};
  SDL_zero(obtained);

  SDL_AudioDeviceID dev = SDL_OpenAudioDevice(nullptr, 0, &desired, &obtained, 0);
  if (dev == 0) {
    SDL_Log("AudioEngine: SDL_OpenAudioDevice failed: %s", SDL_GetError());
    return false;
  }

  device_ = (void*)(std::uintptr_t)dev;
  obtained_ = new SDL_AudioSpec(obtained);
  sampleRate_ = (float)obtained.freq;

  // Reset synth state.
  enginePhase_ = 0.0f;
  lfoPhase_ = 0.0f;
  noiseState_ = 0x12345678u;
  thrusterNoiseLP_ = 0.0f;
  voiceCount_ = 0;

  active_.store(true, std::memory_order_release);

  // Start playback if enabled.
  SDL_PauseAudioDevice(dev, enabled_.load(std::memory_order_relaxed) ? 0 : 1);

  SDL_Log("AudioEngine: device opened (%d Hz, %d ch, samples=%d)",
          obtained.freq, (int)obtained.channels, (int)obtained.samples);
  return true;
}

void AudioEngine::shutdown() {
  if (!active()) return;

  const SDL_AudioDeviceID dev = (SDL_AudioDeviceID)(std::uintptr_t)device_;
  if (dev != 0) {
    // Pause before closing to avoid races with the callback.
    SDL_PauseAudioDevice(dev, 1);
    SDL_CloseAudioDevice(dev);
  }
  device_ = nullptr;

  delete obtained_;
  obtained_ = nullptr;

  active_.store(false, std::memory_order_release);
}

void AudioEngine::setEnabled(bool enabled) {
  enabled_.store(enabled, std::memory_order_relaxed);
  if (!active()) return;

  const SDL_AudioDeviceID dev = (SDL_AudioDeviceID)(std::uintptr_t)device_;
  if (dev != 0) {
    SDL_PauseAudioDevice(dev, enabled ? 0 : 1);
  }
}

void AudioEngine::setMaster(float v01) {
  master_.store(clamp01(v01), std::memory_order_relaxed);
}

void AudioEngine::setBus(float engine, float weapons, float ui, float world) {
  volEngine_.store(clamp01(engine), std::memory_order_relaxed);
  volWeapons_.store(clamp01(weapons), std::memory_order_relaxed);
  volUi_.store(clamp01(ui), std::memory_order_relaxed);
  volWorld_.store(clamp01(world), std::memory_order_relaxed);
}

void AudioEngine::setSpatialize(bool enabled) {
  spatialize_.store(enabled, std::memory_order_relaxed);
}

void AudioEngine::setShipParams(float engine01,
                                float boost01,
                                float thruster01,
                                float angular01,
                                float speed01,
                                bool docked,
                                bool inHyperspace) {
  shipEngine01_.store(clamp01(engine01), std::memory_order_relaxed);
  shipBoost01_.store(clamp01(boost01), std::memory_order_relaxed);
  shipThruster01_.store(clamp01(thruster01), std::memory_order_relaxed);
  shipAngular01_.store(clamp01(angular01), std::memory_order_relaxed);
  shipSpeed01_.store(clamp01(speed01), std::memory_order_relaxed);
  shipDocked_.store(docked, std::memory_order_relaxed);
  shipInHyperspace_.store(inHyperspace, std::memory_order_relaxed);
}

void AudioEngine::play(Sfx sfx, float gain, float pan) {
  if (!enabled_.load(std::memory_order_relaxed)) return;
  (void)enqueue(Event{sfx, gain, pan});
}

bool AudioEngine::enqueue(const Event& e) {
  const std::uint32_t w = qWrite_.load(std::memory_order_relaxed);
  const std::uint32_t r = qRead_.load(std::memory_order_acquire);
  if ((w - r) >= kQueueSize) {
    // Drop when saturated; callback must never block.
    return false;
  }
  q_[w & kQueueMask] = e;
  qWrite_.store(w + 1, std::memory_order_release);
  return true;
}

bool AudioEngine::dequeue(Event& out) {
  const std::uint32_t r = qRead_.load(std::memory_order_relaxed);
  const std::uint32_t w = qWrite_.load(std::memory_order_acquire);
  if (r == w) return false;
  out = q_[r & kQueueMask];
  qRead_.store(r + 1, std::memory_order_release);
  return true;
}

void AudioEngine::sdlCallback(void* userdata, std::uint8_t* stream, int lenBytes) {
  auto* self = reinterpret_cast<AudioEngine*>(userdata);
  if (!self) return;

  // We requested AUDIO_F32SYS stereo.
  const int frames = lenBytes / (int)(sizeof(float) * 2);
  float* out = reinterpret_cast<float*>(stream);
  self->mix(out, frames);
}

void AudioEngine::spawnVoice(const Event& e) {
  if (voiceCount_ >= kMaxVoices) {
    // Replace the oldest voice (simple priority: keep newer events).
    int oldest = 0;
    std::uint32_t bestAge = 0;
    for (int i = 0; i < voiceCount_; ++i) {
      if (voices_[i].age > bestAge) {
        bestAge = voices_[i].age;
        oldest = i;
      }
    }
    voices_[oldest] = voices_[voiceCount_ - 1];
    --voiceCount_;
  }

  Voice v{};
  v.sfx = e.sfx;
  v.gain = std::max(0.0f, e.gain);
  v.pan = clamp11(e.pan);
  v.phase = 0.0f;
  v.phase2 = 0.0f;
  v.rng = (std::uint32_t)(noiseState_ ^ (std::uint32_t)e.sfx * 2654435761u) + 0x9E3779B9u;

  // Per-sfx tuning (dur in seconds, f0/f1 are used as "pitch envelope" endpoints).
  float durSec = 0.10f;
  switch (e.sfx) {
    case Sfx::UiClick:    durSec = 0.045f; v.f0 = 1800.0f; v.f1 = 1200.0f; break;
    case Sfx::UiConfirm:  durSec = 0.085f; v.f0 = 1200.0f; v.f1 = 1800.0f; break;
    case Sfx::UiError:    durSec = 0.120f; v.f0 = 420.0f;  v.f1 = 180.0f;  break;

    case Sfx::WeaponLaser:  durSec = 0.090f; v.f0 = 900.0f;  v.f1 = 520.0f;  break;
    case Sfx::WeaponPulse:  durSec = 0.070f; v.f0 = 1100.0f; v.f1 = 700.0f;  break;
    case Sfx::WeaponCannon: durSec = 0.140f; v.f0 = 140.0f;  v.f1 = 70.0f;   break;
    case Sfx::WeaponRailgun: durSec = 0.110f; v.f0 = 2100.0f; v.f1 = 680.0f;  break;
    case Sfx::WeaponMining: durSec = 0.095f; v.f0 = 650.0f;  v.f1 = 350.0f;  break;
    case Sfx::WeaponMissile:durSec = 0.160f; v.f0 = 260.0f;  v.f1 = 90.0f;   break;

    case Sfx::Explosion: durSec = 0.85f; v.f0 = 120.0f; v.f1 = 38.0f; break;

    case Sfx::Docked:    durSec = 0.180f; v.f0 = 520.0f; v.f1 = 320.0f; break;
    case Sfx::Undocked:  durSec = 0.180f; v.f0 = 480.0f; v.f1 = 900.0f; break;

    case Sfx::FsdCharge: durSec = 0.40f;  v.f0 = 180.0f; v.f1 = 520.0f; break;
    case Sfx::FsdJump:   durSec = 0.80f;  v.f0 = 620.0f; v.f1 = 90.0f;  break;
    case Sfx::FsdArrive: durSec = 0.32f;  v.f0 = 260.0f; v.f1 = 140.0f; break;
  }

  v.age = 0;
  v.dur = std::max<std::uint32_t>(1u, (std::uint32_t)std::round(durSec * sampleRate_));

  voices_[voiceCount_++] = v;
}

float AudioEngine::voiceSample(Voice& v) {
  const float t = (v.dur > 1) ? (float)v.age / (float)(v.dur - 1) : 1.0f;

  // Generic envelope: quick attack, quadratic decay.
  float env = 1.0f - t;
  env *= env;
  const float attack = std::min(1.0f, t * 18.0f);
  env *= attack;

  // Pitch sweep / "chirp".
  const float f = lerp(v.f0, v.f1, t);
  const float dp = twoPi() * f / sampleRate_;
  v.phase += dp;
  if (v.phase > twoPi()) v.phase -= twoPi();

  float s = 0.0f;

  switch (v.sfx) {
    case Sfx::UiClick:
    case Sfx::UiConfirm:
    case Sfx::UiError: {
      // Bright sine/triangle-ish click.
      const float base = std::sin(v.phase);
      const float tri = (2.0f / 3.14159265f) * std::asin(base);
      s = 0.65f * tri + 0.35f * base;
      break;
    }

    case Sfx::WeaponLaser:
    case Sfx::WeaponPulse:
    case Sfx::WeaponMining: {
      // Zap: sine + noisy edge.
      const float base = std::sin(v.phase);
      const float nz = noise01(v.rng) * 0.25f;
      s = base + nz;
      break;
    }

    case Sfx::WeaponCannon:
    case Sfx::WeaponMissile: {
      // Thump: low tone + noise burst.
      const float base = std::sin(v.phase);
      const float nz = noise01(v.rng);
      s = 0.35f * base + 0.75f * nz;
      break;
    }

    case Sfx::Explosion: {
      // Explosion: decaying noise + sub rumble + a faint "ring".
      const float nz = noise01(v.rng);
      const float subF = lerp(v.f0, v.f1, t);
      const float dsub = twoPi() * subF / sampleRate_;
      v.phase2 += dsub;
      if (v.phase2 > twoPi()) v.phase2 -= twoPi();
      const float sub = std::sin(v.phase2);

      // Add a short "crack" early on.
      const float crack = (t < 0.07f) ? (nz * (1.0f - t / 0.07f)) : 0.0f;

      s = 0.85f * nz + 0.55f * sub + 0.25f * crack;
      break;
    }

    case Sfx::Docked:
    case Sfx::Undocked:
    case Sfx::FsdCharge:
    case Sfx::FsdJump:
    case Sfx::FsdArrive: {
      // Whoosh-ish: tone + filtered noise.
      const float base = std::sin(v.phase);
      const float nz = noise01(v.rng) * 0.6f;
      s = base + nz;
      break;
    }
  }

  v.age++;
  return s * env * v.gain;
}

void AudioEngine::mix(float* outInterleaved, int frames) {
  if (!outInterleaved || frames <= 0) return;

  // Pull queued SFX events.
  Event e{};
  while (dequeue(e)) {
    spawnVoice(e);
  }

  const bool enabled = enabled_.load(std::memory_order_relaxed);
  if (!enabled) {
    std::memset(outInterleaved, 0, (std::size_t)frames * sizeof(float) * 2);
    return;
  }

  const float master = master_.load(std::memory_order_relaxed);
  const float volEngine = volEngine_.load(std::memory_order_relaxed);
  const float volWeapons = volWeapons_.load(std::memory_order_relaxed);
  const float volUi = volUi_.load(std::memory_order_relaxed);
  const float volWorld = volWorld_.load(std::memory_order_relaxed);
  const bool spatial = spatialize_.load(std::memory_order_relaxed);

  const float eng01 = shipEngine01_.load(std::memory_order_relaxed);
  const float boost01 = shipBoost01_.load(std::memory_order_relaxed);
  const float thr01 = shipThruster01_.load(std::memory_order_relaxed);
  const float ang01 = shipAngular01_.load(std::memory_order_relaxed);
  const float spd01 = shipSpeed01_.load(std::memory_order_relaxed);
  const bool docked = shipDocked_.load(std::memory_order_relaxed);
  const bool hyper = shipInHyperspace_.load(std::memory_order_relaxed);

  // Continuous ship audio (engine hum + thrusters).
  // Keep it subtle; this is a prototype and we don't want fatigue.
  const float shipActive = (!docked && !hyper) ? 1.0f : 0.0f;
  const float engineAmp = (0.03f + 0.08f * eng01 + 0.10f * boost01) * shipActive;
  const float thrAmp = (0.01f + 0.10f * thr01 + 0.06f * ang01) * shipActive;

  // Slight vibrato to avoid "pure tone" fatigue.
  const float lfoF = 0.8f + 1.8f * spd01;
  const float dlfo = twoPi() * lfoF / sampleRate_;

  // Engine frequency influenced by speed + thrust.
  const float baseF = 45.0f + 120.0f * eng01 + 80.0f * spd01 + 65.0f * boost01;
  const float dphBase = twoPi() * baseF / sampleRate_;

  for (int i = 0; i < frames; ++i) {
    // LFO
    lfoPhase_ += dlfo;
    if (lfoPhase_ > twoPi()) lfoPhase_ -= twoPi();
    const float lfo = 0.008f * std::sin(lfoPhase_);

    // Engine
    enginePhase_ += dphBase * (1.0f + lfo);
    if (enginePhase_ > twoPi()) enginePhase_ -= twoPi();
    const float s = std::sin(enginePhase_);
    const float harm = s * s * s; // soft 3rd harmonic
    float mono = (s + 0.35f * harm) * engineAmp;

    // Thruster noise (low-pass filtered).
    const float nz = noise01(noiseState_);
    thrusterNoiseLP_ = thrusterNoiseLP_ + 0.10f * (nz - thrusterNoiseLP_);
    mono += thrusterNoiseLP_ * thrAmp;

    // Hyperspace ambience (optional): a low, soft drone.
    if (hyper) {
      const float hyperF = 58.0f;
      enginePhase_ += twoPi() * hyperF / sampleRate_;
      if (enginePhase_ > twoPi()) enginePhase_ -= twoPi();
      mono += 0.035f * std::sin(enginePhase_);
    }

    float left = mono * master * volEngine;
    float right = mono * master * volEngine;

    // One-shots
    for (int v = 0; v < voiceCount_;) {
      Voice& voice = voices_[v];
      const float sfx = voiceSample(voice);

      // Remove finished voices.
      if (voice.age >= voice.dur) {
        voices_[v] = voices_[voiceCount_ - 1];
        --voiceCount_;
        continue;
      }

      float busVol = 1.0f;
      switch (voice.sfx) {
        case Sfx::UiClick:
        case Sfx::UiConfirm:
        case Sfx::UiError:
        case Sfx::Docked:
        case Sfx::Undocked:
          busVol = volUi;
          break;

        case Sfx::WeaponLaser:
        case Sfx::WeaponPulse:
        case Sfx::WeaponRailgun:
        case Sfx::WeaponCannon:
        case Sfx::WeaponMining:
        case Sfx::WeaponMissile:
          busVol = volWeapons;
          break;

        case Sfx::Explosion:
        case Sfx::FsdCharge:
        case Sfx::FsdJump:
        case Sfx::FsdArrive:
          busVol = volWorld;
          break;
      }

      const float pan = spatial ? voice.pan : 0.0f;
      float gL, gR;
      panGains(pan, gL, gR);

      left += sfx * gL * master * busVol;
      right += sfx * gR * master * busVol;

      ++v;
    }

    // Clamp to avoid NaNs/clipping.
    outInterleaved[i * 2 + 0] = clamp11(left);
    outInterleaved[i * 2 + 1] = clamp11(right);
  }
}

} // namespace stellar::game
