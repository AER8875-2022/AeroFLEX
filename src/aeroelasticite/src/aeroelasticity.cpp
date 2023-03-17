
#include "common_aeroflex.hpp"
#include <aero/aeroelasticity.hpp>

using namespace structure;
using namespace aero;


Aero::Aero(GUIHandler &gui, vlm::VLM &vlm, structure::Structure &structure)
    : gui(gui), vlm(vlm), structure(structure) {};

void Aero::input() {
    // TODO
    structure.input();
    vlm.input();
    interpolation pos;
    interpolation_f force;
    auto mapStructni= FEM.indexation_switch;
    auto mapStruct= FEM.Grid_MAP;
    dispInterpol(pos, object.wingStations, object.wings, object.vortexRings, object.nodes, mapStruct, mapStructni);
    LoadInterpol(force, object.wingStations, object.wings, object.vortexRings, object.nodes, mapStruct, mapStructni);








}

void Aero::solve() {
    // TODO

    vlm.solve();

    structure.solve(computeVLMForce(force));

    object.updateGeometry(computeVLMDispalecement(pos,object.wingStations,object.wings,object.vortexRings,object.nodes ));
}
