
#include "common_aeroflex.hpp"
#include <aero/aeroelasticity.hpp>
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


    auto& mapStructni= structure.FEM.indexation_switch;
    auto& mapStruct= structure.FEM.Grid_MAP;

    std::cout << "vlm and structure ini check"  << std::endl;
    aero::DispInterpol(pos, vlm.object.wingStations, vlm.object.wings, vlm.object.vortexRings, vlm.object.nodes, mapStruct, mapStructni);
    std::cout << "dispInt check"  << std::endl;
    aero::LoadInterpol(force, vlm.object.wingStations, vlm.object.wings, vlm.object.vortexRings, vlm.object.nodes, mapStruct, mapStructni);
    std::cout << "loadInt check"  << std::endl;
}

void Aero::solve() {
    std::cout << "begin first iter"  << std::endl;

    //first one to initilaize cl
    vlm.solve();
    std::cout << "vlm solved"  << std::endl;
    for(int i=0; i<vlm.object.wingStations.size(); i++){
        std::cout << "station: " << i << std::endl;
        for(int j=0; j<vlm.object.wingStations[i].get_forces().size(); j++) {
            std::cout << vlm.object.wingStations[i].get_forces()[j] << std::endl;
        }}
    Eigen::VectorXd forcestruct= ComputeStructureForces(force, vlm.object.wingStations);
    //forcestructsize
    std::cout << "forcestructsize: " << forcestruct.size() << std::endl;


    structure.FEM.set_Load_Vector_From_Vector(forcestruct);
    //print force
    //for (int i=0; i<forcestruct.size(); i++){
        //std::cout << "force: " << forcestruct[i] << std::endl;
   // }

    std::cout << "structure forces computed"  << std::endl;
    structure.solve();
    std::cout << "structure solved"  << std::endl;
    double old_cl= vlm.object.wings[0].get_cl();


    vlm.object.updateGeometry(computeVLMDispalecement(pos, vlm.object.wingStations, vlm.object.wings, vlm.object.vortexRings,vlm.object.nodes, structure.Solutions));
    std::cout << "end first iter " << std::endl;
    double tol = 0.01;
    do {

        while (gui.signal.pause) std::this_thread::sleep_for(std::chrono::milliseconds(100));

        vlm.solve();

        structure.FEM.set_Load_Vector_From_Vector(ComputeStructureForces(force, vlm.object.wingStations));
        structure.solve();

        vlm.object.updateGeometry(computeVLMDispalecement(pos, vlm.object.wingStations, vlm.object.wings, vlm.object.vortexRings,vlm.object.nodes, structure.Solutions));
        double new_cl = vlm.object.wings[0].get_cl();
        tol = std::abs(new_cl-old_cl);
        std::cout << "Aero tol: " << tol << std::endl;
        old_cl = new_cl;
    }while ( tol > settings.tolerance && !gui.signal.stop);
}
