#ifndef STRUCTURE_LOADS_HPP
#define STRUCTURE_LOADS_HPP

#include "Eigen/Dense"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <map>

namespace structure {

class FORCE{
    public:
    int               Node_ID;  //Node dans le code imposé à quelque chose         
    double               Norm;  //DDL imposé à zéro : 123456
    Eigen::Vector3d Direction;

    FORCE(){};

    FORCE(int node_id, double norm, Eigen::Vector3d direction){

        Node_ID   = node_id;
        Norm      = norm;
        Direction = direction/direction.norm();
    }

    Eigen::Vector3d get_xyz_force()
    {
        return Norm*Direction;
    }

};

class MOMENT{
    public:
    int               Node_ID;  //Node dans le code imposé à quelque chose         
    double               Norm;  //DDL imposé à zéro : 123456
    Eigen::Vector3d Direction;

    MOMENT(){};

    MOMENT(int node_id, double norm, Eigen::Vector3d direction){

        Node_ID   = node_id;
        Norm      = norm;
        Direction = direction/direction.norm();
    }

    Eigen::Vector3d get_xyz_moment()
    {
        return Norm*Direction;
    }



};

}

#endif
