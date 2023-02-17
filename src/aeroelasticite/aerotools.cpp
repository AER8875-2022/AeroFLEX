//
// Created by bocan on 2023-02-16.
//

#include "aerotools.h"


#include <iostream>
#include<algorithm>
#include<array>
#include<vector>
#include"model.hpp"


namespace aero{
    auto DispInterpol(){
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
                        if (distance < min(bestDistance[0], bestDistance[1])) {
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

                    // A faire : resoudre le systeme de 2 eq(produit scalaire et colinearité) avec Eigen pour trouver la projection du noeud sur la structure

                }
            }


        }
    }






    return 0;
}
