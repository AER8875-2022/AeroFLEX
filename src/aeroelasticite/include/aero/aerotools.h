//
// Created by bocan on 2023-02-16.
#pragma once




// Path: aerotools.cpp
// Compare this snippet from main.cpp:
namespace aero{
    struct interpolation{std::vector<double> node;
        std::vector<double> weight;};

    void DispInterpol(interpolation &pos,std::vector<surface::wingStation> wingStations,
                      std::vector<surface::wing> wings,std::vector<element::vortexRing> vortexRings,std::vector<Vector3d> nodes,map<int,Vector3d> mapStruct,map<int, int> mapStructni);
    auto computeVLMDispalecement(interpolation pos,std::vector<surface::wingStation> wingStations,
                                 std::vector<surface::wing> wings,std::vector<element::vortexRing> vortexRings,std::vector<Vector3d> nodes);
    
    struct interpolation_f {std::vector<double> point_fa; std::vector<double> point_fs;
                        std::vector<double> poids;};
    
    //void LoadInterpol(interpolation_f &force, std::vector<surface::wingStation> wingStations,
     //             std::vector<surface::wing> wings,std::vector<element::vortexRing> vortexRings,std::vector<Vector3d> nodes,map<int,Vector3d> mapStruct,map<int, int> mapStructni)

    //std::vector<double> ComputeStructureForces(interpolation_f force,Matrix<double, 6, 1> forces_to_inertial_frame);

    std::vector<double> crossProduct(const std::vector<double>& i, const std::vector<double>& j);

    double norme(std::vector<double> v);


}