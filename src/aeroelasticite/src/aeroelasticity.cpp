
#include "common_aeroflex.hpp"
#include <aero/aeroelasticity.hpp>

using namespace structure;


Aero::Aero(GUIHandler &gui) : gui(gui) {};

void Aero::input() {
    // TODO
    structure.input();
    interpolation pos;
    auto mapStructni= FEM.indexation_switch;
    auto mapStruct= FEM.Grid_MAP;
    dispInterpol(pos, model.wingStations, model.wings, model.vortexRings, model.nodes, mapStruct, mapStructni);







}

void Aero::solve() {
    // TODO
    structure.solve();
    computeVLMDispalecement(pos, );
}
