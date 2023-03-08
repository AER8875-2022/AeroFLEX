//
// Created by bocan on 2023-02-16.
//

#include "aerotools.h"


#include <iostream>
#include<algorithm>
#include<array>
#include<vector>
#include"model.hpp"
#include<Eigen/Dense>
#include <Eigen/Geometry>


namespace aero{
    struct interpolation{
        std::vector<double> node;
        std::vector<double> weight;

    };


    auto DispInterpol( interpolation &pos ){
        for (k = 0; k < ; ++k)
        {
            for (i = 0; i < ; ++i)
            {

                auto wstation= wingStations[wings[0].get_wingStationsIDs()[k]];
                auto vring = vortexRings[wstation.get_vortexIDs()[i]];
                int numPointsVort =sizeof(vring.getnodesIDs())/sizeof(vring.getnodesIDs()[0]);


                int bestDistance[2]={0,0};
                int bestPoint[2]={0,0};
                for (p=0; p<numPointsVort; ++p)
                {
                    auto node=nodes[vring.get_nodeIDs()[p]];

                    //compute the distance between the node on vlm and nodes on structure

                    for (s=0 ; s<numPointsStruct; ++s) {
                        auto nodeStruct = nodes[struct.get_nodeIDs()[s]];
                        distance = sqrt(
                                pow(nodeStruct[0]-node[0], 2) + pow(nodeStruct[1]-node[1], 2) +
                                pow(nodeStruct[2]-node[2], 2));
                        if (distance < max(bestDistance[0], bestDistance[1])) {
                            if (bestDistance[0] < bestDistance[1]) {
                                bestDistance[0] = distance;
                                bestPoint[0] = s;
                            } else {
                                bestDistance[1] = distance;
                                bestPoint[1] = s;
                            }


                        }

                    }

                    //projection of the node on the structure
                    auto nodeStruct1 = nodes[struct.get_nodeIDs()[bestPoint[0]]];
                    auto nodeStruct2 = nodes[struct.get_nodeIDs()[bestPoint[1]]];

                    pos.node.push_back(struct.get_nodeIDs()[bestPoint[0]);
                    pos.node.push_back(struct.get_nodeIDs()[bestPoint[1]);

                    auto x1= nodeStruct1-nodeStruct2 ;
                    auto x2= nodeStruct1-node ;
                    double projection = (x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2])/pow(x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2],0.5);

                    pos.weight.push_back(projection);
                    pos.weight.push_back(1-projection);




                }
            }


        }
    }

    function computeVLMDispalecement() {
        for (k = 0; k<; ++k) {
            for (i = 0; i<; ++i) {

                auto wstation = wingStations[wings[0].get_wingStationsIDs()[k]];
                auto vring = vortexRings[wstation.get_vortexIDs()[i]];
                int numPointsVort = vring.size();
                for (p=0; p<numPointsVort; ++p)
                {
                    auto node = nodes[vring.get_nodeIDs()[p]];
                    auto nodeStruct1 = nodes[pos.node[2*p]];
                    auto nodeStruct2 = nodes[pos.node[2*p+1]];
                    auto weight1 = pos.weight[2*p];
                    auto weight2 = pos.weight[2*p+1];

                    Eigen::Vector3d disp;
                    disp[0]=weight1*nodeStruct1[0]+weight2*nodeStruct2[0];
                    disp[1]=weight1*nodeStruct1[1]+weight2*nodeStruct2[1];
                    disp[2]=weight1*nodeStruct1[2]+weight2*nodeStruct2[2];
                    //prise en compte de la rotation
                    Eigen::vector4d q1(1,2,3,4), q2(1,2,3,4), q3(1,2,3,4);

                    

                    node[0]=node[0]+disp[0];
                    node[1]=node[1]+disp[1];
                    node[2]=node[2]+disp[2];

                }
            }
        }
    }





}
