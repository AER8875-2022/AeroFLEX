
#include "common_aeroflex.hpp"
#include <aero/aeroelasticity.hpp>
#include "aero/aerotools.h"
#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include "vlm/model.hpp"
#include "structure/structure.hpp"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "vlm/vlm.hpp"
#include <cmath>

using namespace structure;
using namespace aero;
using namespace vlm;


Aero::Aero(GUIHandler &gui, vlm::VLM &vlm, structure::Structure &structure)
    : gui(gui), vlm(vlm), structure(structure) {};

void Aero::input() {
    // TODO
    structure.input();
    vlm.initialize();
    aero::interpolation pos;
    aero::interpolation_f force;
    auto mapStructni= structure.FEM.indexation_switch;
    auto mapStruct= structure.FEM.Grid_MAP;
    aero::DispInterpol(pos, vlm.object.wingStations, vlm.object.wings, vlm.object.vortexRings, vlm.object.nodes, mapStruct, mapStructni);
    aero::LoadInterpol(force, vlm.object.wingStations, vlm.object.wings, vlm.object.vortexRings, vlm.object.nodes, mapStruct, mapStructni);








}

void Aero::solve() {
    // TODO

    vlm.solve();
    interpolation_f force;
    interpolation pos;
    structure.FEM.set_Load_Vector_From_Vector(ComputeStructureForces(force,vlm.object.wingStations));
    structure.solve();

    vlm.object.updateGeometry(computeVLMDispalecement(pos,vlm.object.wingStations,vlm.object.wings,vlm.object.vortexRings,vlm.object.nodes,structure.Solutions ));
}
