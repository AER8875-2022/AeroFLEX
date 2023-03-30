#ifndef STRUCTURE_SPC_HPP
#define STRUCTURE_SPC_HPP

#include "Eigen/Dense"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <map>

namespace structure {

class SPC1{
    public:
    int          Node_ID;  //Node dans le code imposé à quelque chose         
    std::string     CODE;  //DDL imposé à zéro : 123456

    SPC1(){};

    SPC1(int node_id, std::string c){

        Node_ID = node_id;
        CODE = c;
    }

};

}

#endif
