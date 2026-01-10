#pragma once

#include <string>

namespace stellar::ui {

// Persistent audio configuration.
//
// This is intentionally lightweight and human-editable (plain text key/value pairs).
// It is saved separately from:
//  - ui_settings.txt (UI theme/scale)
//  - hud_settings.txt (HUD + tactical overlay settings)
//  - imgui.ini (window positions)
//  - savegame.txt (player progression)
//
// The game prototype currently synthesizes audio procedurally at runtime, so these
// settings mostly tune mixing/volume rather than selecting external assets.
struct AudioSettings {
  bool enabled{true};

  // Master output gain (0..1).
  float master{0.70f};

  // Per-bus volumes (0..1). These are multiplied by master.
  float engine{0.75f};  // ship hum + thrusters (continuous)
  float weapons{0.85f}; // weapon shots
  float ui{0.65f};      // UI beeps / confirmations
  float world{0.85f};   // explosions + jump sounds

  // If enabled, some SFX are stereo panned based on relative direction.
  bool spatialize{true};
};

// Default config file path (relative to the working directory).
std::string defaultAudioSettingsPath();

AudioSettings makeDefaultAudioSettings();

bool saveAudioSettingsToFile(const AudioSettings& s, const std::string& path);
bool loadAudioSettingsFromFile(const std::string& path, AudioSettings& out);

} // namespace stellar::ui
