
#include "common_aeroflex.hpp"
#include <aero/aeroelasticity.hpp>

using namespace structure;
using namespace aero;


Aero::Aero(GUIHandler &gui, vlm::VLM &vlm, structure::Structure &structure)
    : gui(gui), vlm(vlm), structure(structure) {};

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
