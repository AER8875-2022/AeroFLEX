#pragma once

#include "common_aeroflex.hpp"
#include "vlm/vlm.hpp"
#include "structure/structure.hpp"
#include "aerotools.h"
#include <vector>

namespace aero {

struct Settings {
    double tolerance = 1e-3;
};

class Aero {

public:
    Settings settings;

    std::vector<double> residuals;
    std::atomic<int> iters = 0;
    GUIHandler &gui;


    // Modules
    vlm::VLM &vlm;
    structure::Structure &structure;

    Aero(GUIHandler &gui, vlm::VLM &vlm, structure::Structure &structure);
    void input();
    void solve();

private:
    interpolation pos;
    interpolation_f force;
};

}
