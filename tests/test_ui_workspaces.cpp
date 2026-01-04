#include "test_harness.h"

#include "stellar/ui/UiWorkspaces.h"

#include <fstream>

using stellar::ui::UiWorkspace;
using stellar::ui::UiWorkspaces;

int test_ui_workspaces() {
  int failures = 0;

  const std::string path = "test_ui_workspaces_tmp.txt";

  // Roundtrip save/load.
  {
    UiWorkspaces w = stellar::ui::makeDefaultUiWorkspaces();
    w.autoSaveOnExit = false;

    UiWorkspace a;
    a.name = "Trading";
    a.imguiIniFile = "ui_layouts/trade.ini";
    a.windows["Galaxy"] = true;
    a.windows["Trade"] = true;
    a.windows["Console"] = false;

    UiWorkspace b;
    b.name = "Combat";
    b.imguiIniFile = "ui_layouts/combat.ini";
    b.windows["Contacts"] = true;
    b.windows["Ship"] = true;

    w.items.push_back(a);
    w.items.push_back(b);
    w.active = "Combat";

    CHECK(stellar::ui::saveUiWorkspacesToFile(w, path));

    UiWorkspaces out;
    CHECK(stellar::ui::loadUiWorkspacesFromFile(path, out));
    CHECK(out.autoSaveOnExit == false);
    CHECK(out.active == "Combat");
    CHECK(out.items.size() == 2);
    {
      const UiWorkspace* c = stellar::ui::findWorkspace(out, "Combat");
      CHECK(c != nullptr);
      CHECK(c->imguiIniFile == "ui_layouts/combat.ini");
      CHECK(stellar::ui::getWindowOpen(*c, "Contacts", false) == true);
      CHECK(stellar::ui::getWindowOpen(*c, "MissingKey", true) == true);
    }
  }

  // Resilience: unknown tokens + duplicates (keep last).
  {
    std::ofstream f(path, std::ios::out | std::ios::trunc);
    f << "StellarForgeUiWorkspaces 1\n";
    f << "autoSaveOnExit 1\n";
    f << "active Alpha\n";
    f << "workspace Alpha\n";
    f << "ini imgui.ini\n";
    f << "win Galaxy 1\n";
    f << "garbageToken shouldBeIgnored\n";
    f << "workspace Alpha\n"; // duplicate name
    f << "ini ui_layouts/alpha.ini\n";
    f << "win Galaxy 0\n";
    f.close();

    UiWorkspaces out;
    CHECK(stellar::ui::loadUiWorkspacesFromFile(path, out));
    CHECK(out.items.size() == 1);
    const UiWorkspace* a = stellar::ui::findWorkspace(out, "Alpha");
    CHECK(a != nullptr);
    CHECK(a->imguiIniFile == "ui_layouts/alpha.ini");
    CHECK(stellar::ui::getWindowOpen(*a, "Galaxy", true) == false);
  }

  return failures;
}
