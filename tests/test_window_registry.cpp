#include "test_harness.h"

#include "stellar/ui/WindowRegistry.h"
#include "stellar/ui/UiWorkspaces.h"

using stellar::ui::UiWorkspace;
using stellar::ui::WindowBinding;
using stellar::ui::WindowDesc;
using stellar::ui::WindowRegistry;

int test_window_registry() {
  int failures = 0;

  bool aOpen = true;
  bool bOpen = false;
  int aOpened = 0;
  int aClosed = 0;

  WindowRegistry reg;
  reg.add(WindowBinding{WindowDesc{"A", "Alpha", "Main", 10, true, true}, &aOpen,
                        [&]() { ++aOpened; }, [&]() { ++aClosed; }, {}});
  reg.add(WindowBinding{WindowDesc{"B", "Beta", "Main", 0, false, true}, &bOpen, {}, {}, {}});

  // Capture.
  UiWorkspace ws;
  ws.name = "Test";
  stellar::ui::captureWorkspaceWindows(reg, ws);
  CHECK(ws.windows.size() == 2);
  CHECK(ws.windows["A"] == true);
  CHECK(ws.windows["B"] == false);

  // Apply (with callbacks).
  aOpen = false;
  bOpen = true;
  CHECK(aOpened == 0);
  CHECK(aClosed == 0);
  stellar::ui::applyWorkspaceWindows(reg, ws, /*fireCallbacks=*/true);
  CHECK(aOpen == true);
  CHECK(bOpen == false);
  CHECK(aOpened == 1);
  CHECK(aClosed == 0);

  // Bulk reset -> defaults.
  aOpen = false;
  bOpen = true;
  reg.resetToDefaults(/*fireCallbacks=*/false);
  CHECK(aOpen == true);
  CHECK(bOpen == false);

  // Replacement by key.
  bool a2Open = false;
  reg.add(WindowBinding{WindowDesc{"A", "Alpha2", "Main", 0, false, true}, &a2Open, {}, {}, {}});
  CHECK(reg.find("A") != nullptr);
  CHECK(reg.find("A")->desc.label == "Alpha2");
  CHECK(reg.find("A")->open == &a2Open);

  return failures;
}
