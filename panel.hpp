// SPDX-License-Identifier: MIT
// Simple ImGui control panel for the BlackholeSim engine
// Keeps dependencies minimal and avoids redundant includes.

#pragma once

#include <string>
#include <vector>
#include <functional>

#include "imgui.h"
#include "engine.hpp"

namespace BlackholeSim {

    using PanelList = std::vector<std::pair<std::string, std::function<void(Engine&, const std::string&)>>>;

    inline PanelList createPanels() {
        return {
            {"Photonen", [](Engine& eng, const std::string& name) {
                if (ImGui::SliderInt(name.c_str(), &eng.photonCount, 10, 2000)) {
                    eng.resetPhotons(eng.photonCount);
                }
            }},
            {"Spin a/M", [](Engine& eng, const std::string& name) {
                double amin = 0.0, amax = 1.0;
                if (ImGui::SliderScalar(name.c_str(), ImGuiDataType_Double, &eng.a_spin, &amin, &amax)) {
                    // Integrator reads a_spin each step
                }
            }},
            {"Trail Length", [](Engine& eng, const std::string& name) {
                size_t trmin = 10, trmax = 5000;
                uint64_t tmp = static_cast<uint64_t>(eng.trailLength);
                if (ImGui::SliderScalar(name.c_str(), ImGuiDataType_U64, &tmp, &trmin, &trmax)) {
                    eng.trailLength = static_cast<size_t>(tmp);
                    for (auto& p : eng.kerrPhotons)  p.trail.reserve(eng.trailLength);
                    for (auto& p : eng.testPhotons)  p.trail.reserve(eng.trailLength);
                }
            }},
            {"Zoom", [](Engine& eng, const std::string& name) {
                double zoommin = 0.01, zoommax = 2.0;
                ImGui::SliderScalar(name.c_str(), ImGuiDataType_Double, &eng.zoom, &zoommin, &zoommax, "%.3f");
            }},
            {"Switch Mode", [](Engine& eng, const std::string&) {
                if (ImGui::RadioButton("Kerr-Modus", reinterpret_cast<int*>(&eng.mode), static_cast<int>(Mode::Kerr))) {
                    eng.resetPhotons(eng.photonCount);
                }
                ImGui::SameLine();
                if (ImGui::RadioButton("Test-Modus", reinterpret_cast<int*>(&eng.mode), static_cast<int>(Mode::Test))) {
                    eng.resetPhotons(eng.photonCount);
                }
            }},
            {"Reset Photons", [](Engine& eng, const std::string& name) {
                if (ImGui::Button(name.c_str())) {
                    eng.resetPhotons(eng.photonCount);
                }
            }},
            {"End Controlpanel", [](Engine& eng, const std::string&) {
                ImGui::Text("Photons: %d  TrailLen: %zu", eng.photonCount, eng.trailLength);
                ImGui::Text("Spin: %.2f     Zocom: %.2f", eng.a_spin, eng.zoom);
            }},
            {"Show Horizons", [](Engine& eng, const std::string&) {
                ImGui::Checkbox("Show Horizons", &eng.showHorizons);
            }},
            {"Show Ergosphere", [](Engine& eng, const std::string&) {
                ImGui::Checkbox("Show Ergosphere", &eng.showErgosphere);
            }}
    };
}

    inline void Controlpanel(Engine& eng) {
        ImGui::Begin("Kontrolle");

        static PanelList panels = createPanels();
        for (auto& [name, fn] : panels) {
            fn(eng, name);
        }

        ImGui::End();
    }
}