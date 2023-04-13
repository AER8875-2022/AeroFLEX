#ifndef STRUCTURE_MATERIALS_HPP
#define STRUCTURE_MATERIALS_HPP

#include "Eigen/Dense"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <map>

namespace structure {

class MAT1
{
    public:
    double E ;             //Module de Young
    double G ;             //Module de rigidité en torsion
    
    MAT1(){};
    MAT1(double e,double g)
    {
        E=e;
        G=g;
    };

    
};

}

#endif
