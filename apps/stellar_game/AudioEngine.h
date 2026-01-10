#pragma once

#include <atomic>
#include <cstdint>

struct SDL_AudioSpec;

namespace stellar::game {

// Procedural audio synth + mixer built directly on SDL2's audio callback.
//
// Goals:
//  - zero external assets (everything synthesized)
//  - real-time safe callback (no allocations/locks inside callback)
//  - simple SFX hooks for combat/flight feedback
//
// This is intentionally "game-side" (apps/stellar_game) rather than in the core
// library so the simulation/tests remain SDL-free.
class AudioEngine {
public:
  enum class Sfx : std::uint8_t {
    UiClick,
    UiConfirm,
    UiError,

    WeaponLaser,
    WeaponPulse,
    WeaponCannon,
    WeaponRailgun,
    WeaponMining,
    WeaponMissile,

    Explosion,

    Docked,
    Undocked,

    FsdCharge,
    FsdJump,
    FsdArrive,
  };

  AudioEngine();
  ~AudioEngine();

  // Open the audio device and start playback.
  // Returns false if audio could not be initialized (game should continue silently).
  bool init();

  void shutdown();

  bool active() const { return active_.load(std::memory_order_relaxed); }

  // Settings
  void setEnabled(bool enabled);
  void setMaster(float v01);
  void setBus(float engine, float weapons, float ui, float world);
  void setSpatialize(bool enabled);

  // Continuous ship state (polled each frame from main thread).
  void setShipParams(float engine01,
                     float boost01,
                     float thruster01,
                     float angular01,
                     float speed01,
                     bool docked,
                     bool inHyperspace);

  // One-shot SFX.
  //  gain: nominal 0..1 (can exceed slightly for emphasis; mixer clamps)
  //  pan:  -1 (left) .. +1 (right)
  void play(Sfx sfx, float gain = 1.0f, float pan = 0.0f);

private:
  struct Event {
    Sfx sfx{Sfx::UiClick};
    float gain{1.0f};
    float pan{0.0f};
  };

  static void sdlCallback(void* userdata, std::uint8_t* stream, int lenBytes);

  void mix(float* outInterleaved, int frames);

  bool enqueue(const Event& e);
  bool dequeue(Event& out);

  // Device
  std::atomic<bool> active_{false};
  std::atomic<bool> enabled_{true};
  void* device_{nullptr}; // SDL_AudioDeviceID stored as void* to avoid SDL headers here.
  SDL_AudioSpec* obtained_{nullptr};

  // Mixer parameters (atomics: written from main thread, read by callback)
  std::atomic<float> master_{0.7f};
  std::atomic<float> volEngine_{0.75f};
  std::atomic<float> volWeapons_{0.85f};
  std::atomic<float> volUi_{0.65f};
  std::atomic<float> volWorld_{0.85f};
  std::atomic<bool> spatialize_{true};

  std::atomic<float> shipEngine01_{0.0f};
  std::atomic<float> shipBoost01_{0.0f};
  std::atomic<float> shipThruster01_{0.0f};
  std::atomic<float> shipAngular01_{0.0f};
  std::atomic<float> shipSpeed01_{0.0f};
  std::atomic<bool> shipDocked_{false};
  std::atomic<bool> shipInHyperspace_{false};

  // SPSC queue (main thread -> audio callback)
  static constexpr std::uint32_t kQueueSize = 512;
  static constexpr std::uint32_t kQueueMask = kQueueSize - 1;
  static_assert((kQueueSize & (kQueueSize - 1)) == 0, "kQueueSize must be power of two");
  std::atomic<std::uint32_t> qWrite_{0};
  std::atomic<std::uint32_t> qRead_{0};
  Event q_[kQueueSize];

  // Synth state (audio thread only)
  float sampleRate_{48000.0f};
  float enginePhase_{0.0f};
  float lfoPhase_{0.0f};
  std::uint32_t noiseState_{0x12345678u};
  float thrusterNoiseLP_{0.0f};

  struct Voice {
    Sfx sfx{Sfx::UiClick};
    float gain{1.0f};
    float pan{0.0f};
    float phase{0.0f};
    float phase2{0.0f};
    float f0{440.0f};
    float f1{440.0f};
    std::uint32_t rng{0xC0FFEEu};
    std::uint32_t age{0};
    std::uint32_t dur{1};
  };

  static constexpr int kMaxVoices = 48;
  Voice voices_[kMaxVoices];
  int voiceCount_{0};

  void spawnVoice(const Event& e);
  float voiceSample(Voice& v);
};

} // namespace stellar::game
