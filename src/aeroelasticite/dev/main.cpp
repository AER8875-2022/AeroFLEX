//
// Created by bocan on 2023-02-16.
//
#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include "aero/aerotools.h"
#include "vlm/model.hpp"
#include "aero/aeroelasticity.hpp"




int main() {




    GUIHandler gui;
    vlm::VLM vlm(gui);
    structure::Structure structure(gui);

    aero::Aero elast(gui, vlm, structure);
    elast.input();
    elast.solve();

    return 0;
}
