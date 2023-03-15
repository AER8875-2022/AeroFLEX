
#ifndef __AEROELASTICITY__
#define __AEROELASTICITY__

#include "common_aeroflex.hpp"
#include "vlm/vlm.hpp"
#include "structure/structure.hpp"
#include <vector>

namespace aero {

struct Settings {
    // TODO:
};

class Aero {

public:
    // TODO: STORE SOLUTION
    Settings settings;

    std::vector<double> residuals;
    std::atomic<int> iters = 0;
    GUIHandler &gui;
    structure Structure;

    // Modules
    vlm::VLM &vlm;
    structure::Structure &structure;

public:
    Aero(GUIHandler &gui, vlm::VLM &vlm, structure::Structure &structure);
    void input();
    void solve();

};

}

#endif
