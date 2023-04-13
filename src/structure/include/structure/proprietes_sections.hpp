#ifndef STRUCTURE_PROPRIETES_SECTIONS_HPP
#define STRUCTURE_PROPRIETES_SECTIONS_HPP

#include "Eigen/Dense"
#include <vector>
#include <cmath>
#include <iostream>
#include "structure/materials.hpp"

namespace structure {

class PBAR{

public:
    double A ;             //Aire de section
    double E ;             //Module de Young
    double G ;             //Module de rigidité en torsion
    double Iz;             //Second moment de surface en z
    double Iy;             //Second moment de surface en z
    double J ;             //Second moment de surface en 


    PBAR(){};
    PBAR(double a,double iz,double iy,double j, MAT1 mat1){
        A = a;
        E  = mat1.E;
        G  = mat1.G;
        Iz = iz;
        Iy = iy;
        J  = j ;
    }

};

}

#endif
