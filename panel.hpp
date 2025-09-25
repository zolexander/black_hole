#include <string>
#include <map>
#include <functional>
#include "engine.hpp"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <kerrstate.hpp>
namespace BlackholeSim
{
    using PanelList = std::vector<std::pair<std::string, std::function<void(Engine &)>>>;

    inline PanelList createPanels()
    {
        return {
            {"Photonen", [](Engine &eng)
             {
                 if (ImGui::SliderInt("Photonen", &eng.photonCount, 10, 2000))
                 {
                     eng.resetPhotons(eng.photonCount);
                 }
             }},
            {"Spin a/M", [](Engine &eng)
             {
                 double amin = 0.0, amax = 1.0;
                 if (ImGui::SliderScalar("Spin a/M", ImGuiDataType_Double, &eng.a_spin, &amin, &amax))
                 {
                     // nothing else needed here; integrator reads a_spin each step
                 }
             }},
            {"Trail Length", [](Engine &eng)
             {
                 size_t trmin = 10, trmax = 5000;
                 uint64_t tmp = (uint64_t)eng.trailLength;
                 if (ImGui::SliderScalar("Trail", ImGuiDataType_U64, &tmp, &trmin, &trmax))
                 {
                     eng.trailLength = (size_t)tmp;
                     for (auto &p : eng.kerrPhotons)
                         p.trail.reserve(eng.trailLength);
                     for (auto &p : eng.testPhotons)
                         p.trail.reserve(eng.trailLength);
                 }
             }},
            {"Zoom", [](Engine &eng)

             {
                 double zoommin = 0.01, zoommax = 2.0;
                 ImGui::SliderScalar("Zoom", ImGuiDataType_Double, &eng.zoom, &zoommin, &zoommax, "%.3f");
             }},
            {"Switch Mode", [](Engine &eng)
             {
                 if (ImGui::RadioButton("Kerr-Modus", (int *)&eng.mode, (int)Mode::Kerr))
                 {
                     eng.resetPhotons(eng.photonCount);
                 }
                 ImGui::SameLine();
                 if (ImGui::RadioButton("Test-Modus", (int *)&eng.mode, (int)Mode::Test))
                 {
                     eng.resetPhotons(eng.photonCount);
                 }
             }},
            {"Reset Photons", [](Engine &eng)
             {
                 if (ImGui::Button("Reset Photonen"))
                 {
                     eng.resetPhotons(eng.photonCount);
                 }
             }},
            {"End Controlpanel", [](Engine &eng)
             {
                 ImGui::Text("Photons: %d  TrailLen: %zu", eng.photonCount, eng.trailLength);
                 ImGui::Text("Spin: %.2f     Zoom: %.2f", eng.a_spin, eng.zoom);
             }}};
    }

    // Panel-Aufruf
    inline void Controlpanel(Engine &eng)
    {
        ImGui::Begin("Kontrolle");

        static PanelList panels = createPanels();
        for (auto &[name, fn] : panels)
        {
            fn(eng);
        }

        ImGui::End();
    }
}