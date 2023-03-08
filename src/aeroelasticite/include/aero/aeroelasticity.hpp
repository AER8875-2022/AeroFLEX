
#ifndef __AEROELASTICITY__
#define __AEROELASTICITY__

#include "common_aeroflex.hpp"
#include <vector>

namespace aero {

struct Settings {
    // TODO:
};

class Aero {

public:
    // TODO: SOLUTION
    Settings settings;

    std::vector<double> residuals;
    std::atomic<int> iters = 0;
    GUIHandler &gui;

public:
    Aero(GUIHandler &gui);
    void input();
    void solve();

};

}

#endif
