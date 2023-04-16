
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
            std::cout << vlm.object.wingStations[i].get_forces()[j] <<"  " ;
        }
        std::cout<<std::endl;
        }
    Eigen::VectorXd forcestruct= ComputeStructureForces(force, vlm.object.wingStations);
    //forcestructsize
    std::cout << "forcestructsize: " << forcestruct.size() << std::endl;
    int s=forcestruct.size()/6;
    for (int i=0; i<s; i++)
    {
        std::cout<< i <<" :";

        std::cout<< forcestruct[6*i];
        std::cout<< forcestruct[6*i+1];
        std::cout<< forcestruct[6*i+2];
        std::cout<< forcestruct[6*i+3];
        std::cout<< forcestruct[6*i+4];
        std::cout<< forcestruct[6*i+5]<<std::endl;

    }

    structure.FEM.set_Load_Vector_From_Vector(forcestruct);
    //print force
   /* for (int i=0; i<forcestruct.size(); i++){
        std::cout << "force: " << forcestruct[i] << std::endl;
    }*/

    std::cout << "structure forces computed"  << std::endl;
    structure.solve();
    std::cout << "structure solved"  << std::endl;
    double old_cl= vlm.object.wings[0].get_cl();
    std::cout << "old_cl: " << old_cl << std::endl;


    vlm.object.updateGeometry(computeVLMDispalecement(pos, vlm.object.wingStations, vlm.object.wings, vlm.object.vortexRings,vlm.object.nodes, structure.Solutions));
    std::cout << "end first iter " << std::endl;
    //print solutions

    double tol = 0.01;
    int iter=1;
    do {

        while (gui.signal.pause) std::this_thread::sleep_for(std::chrono::milliseconds(100));

        vlm.solve();

        structure.FEM.set_Load_Vector_From_Vector(ComputeStructureForces(force, vlm.object.wingStations));
        structure.solve();


        vlm.object.updateGeometry(computeVLMDispalecement(pos, vlm.object.wingStations, vlm.object.wings, vlm.object.vortexRings,vlm.object.nodes, structure.Solutions));
        double new_cl = vlm.object.wings[0].get_cl();
        //std::cout << "cl " << new_cl << std::endl;

        tol = std::abs(new_cl-old_cl);
        //std::cout << "Aero tol: " << tol << std::endl;
        old_cl = new_cl;
        iter=iter+1;
    }while ( tol > settings.tolerance && !gui.signal.stop);

    std::cout<<"Nombre d'itérations: "<< iter << std::endl;


 for (int i=0; i<structure.Solutions.size(); i++) {
     std::cout<< "========================== " <<std::endl;


     std::cout<< "Itération : " << i <<std::endl;

     std::cout<< "========================== " <<std::endl;
     for (int j=9; j<structure.Solutions[i].size()/6;j++)
     {
         //double d=99999;
         if(i==0)
         {
         double x= structure.Solutions[i][6*j];
         double y= structure.Solutions[i][6*j+1];
         double z= structure.Solutions[i][6*j+2];
        double d = std::sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        std::cout << "solution: " << d << std::endl;
         }
         else{
             double x= structure.Solutions[i][6*j]-structure.Solutions[i-1][6*j];
         double y= structure.Solutions[i][6*j+1]-structure.Solutions[i-1][6*j+1];
         double z= structure.Solutions[i][6*j+2]-structure.Solutions[i-1][6*j+2];
         double d = std::sqrt(pow(x,2)+pow(y,2)+pow(z,2));
         std::cout << "solution: " << d << std::endl;

         }

        

        
     }

    }
}