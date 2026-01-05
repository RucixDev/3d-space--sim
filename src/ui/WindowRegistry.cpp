#include "stellar/ui/WindowRegistry.h"

#include "stellar/ui/UiWorkspaces.h"

#include <algorithm>
#include <cctype>

namespace stellar::ui {

static std::string trimAscii(std::string s) {
  std::size_t a = 0;
  while (a < s.size() && std::isspace((unsigned char)s[a])) ++a;
  std::size_t b = s.size();
  while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
  return s.substr(a, b - a);
}

static std::string sanitizeKey(std::string_view key) {
  std::string s(key);
  // Keep keys file-friendly so they can round-trip through UiWorkspaces.
  for (char& c : s) {
    if (c == '\n' || c == '\r' || c == '\t') c = ' ';
  }
  s = trimAscii(std::move(s));
  for (char& c : s) {
    if (std::isspace((unsigned char)c)) c = '_';
  }
  return s;
}

WindowBinding& WindowRegistry::add(WindowBinding b) {
  b.desc.key = sanitizeKey(b.desc.key);
  // Ignore invalid entries (but keep API total).
  if (b.desc.key.empty()) {
    static WindowBinding dummy;
    return dummy;
  }

  for (auto& existing : items_) {
    if (existing.desc.key == b.desc.key) {
      existing = std::move(b);
      return existing;
    }
  }

  items_.push_back(std::move(b));
  return items_.back();
}

WindowBinding* WindowRegistry::find(std::string_view key) {
  const std::string k = sanitizeKey(key);
  for (auto& b : items_) {
    if (b.desc.key == k) return &b;
  }
  return nullptr;
}

const WindowBinding* WindowRegistry::find(std::string_view key) const {
  const std::string k = sanitizeKey(key);
  for (const auto& b : items_) {
    if (b.desc.key == k) return &b;
  }
  return nullptr;
}

bool WindowRegistry::setOpen(std::string_view key, bool open, bool fireCallbacks) {
  WindowBinding* b = find(key);
  if (!b || !b->open) return false;
  const bool was = *b->open;
  if (was == open) return false;
  *b->open = open;

  if (fireCallbacks) {
    if (open) {
      if (b->onOpened) b->onOpened();
    } else {
      if (b->onClosed) b->onClosed();
    }
  }
  return true;
}

void WindowRegistry::setAll(bool open, bool fireCallbacks) {
  for (auto& b : items_) {
    if (!b.open) continue;
    if (*b.open == open) continue;
    *b.open = open;
    if (fireCallbacks) {
      if (open) {
        if (b.onOpened) b.onOpened();
      } else {
        if (b.onClosed) b.onClosed();
      }
    }
  }
}

void WindowRegistry::resetToDefaults(bool fireCallbacks) {
  for (auto& b : items_) {
    if (!b.open) continue;
    const bool want = b.desc.defaultOpen;
    if (*b.open == want) continue;
    *b.open = want;
    if (fireCallbacks) {
      if (want) {
        if (b.onOpened) b.onOpened();
      } else {
        if (b.onClosed) b.onClosed();
      }
    }
  }
}

void captureWorkspaceWindows(const WindowRegistry& reg, UiWorkspace& ws) {
  ws.windows.clear();
  for (const auto& b : reg.items()) {
    if (!b.desc.persistInWorkspace) continue;
    if (!b.open) continue;
    if (b.desc.key.empty()) continue;
    ws.windows[b.desc.key] = *b.open;
  }
}

void applyWorkspaceWindows(WindowRegistry& reg, const UiWorkspace& ws, bool fireCallbacks) {
  for (const auto& b : reg.items()) {
    if (!b.desc.persistInWorkspace) continue;
    if (!b.open) continue;
    const auto it = ws.windows.find(b.desc.key);
    if (it == ws.windows.end()) continue;
    reg.setOpen(b.desc.key, it->second, fireCallbacks);
  }
}

} // namespace stellar::ui
