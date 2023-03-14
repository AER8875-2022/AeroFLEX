
#include "common_aeroflex.hpp"
#include <aero/aeroelasticity.hpp>

using namespace structure;


Aero::Aero(GUIHandler &gui) : gui(gui) {};

void Aero::input() {
    // TODO
    structure.input();
    auto mapStructni= FEM.indexation_switch;
    auto mapStruct= FEM.Grid_MAP;




}

void Aero::solve() {
    // TODO
}
