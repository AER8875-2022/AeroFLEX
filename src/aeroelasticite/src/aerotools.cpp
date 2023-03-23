//
// Created by bocan on 2023-02-16.
//

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

using namespace Eigen;
using namespace vlm;
using namespace structure;
using namespace std;


namespace aero{
    struct interpolation;
    struct interpolation_f;


    void DispInterpol( interpolation &pos,std::vector<vlm::surface::wingStation> wingStations,
                       std::vector<vlm::surface::wing> wings,std::vector<element::vortexRing> vortexRings,
                      std::vector<vlm::surface::wing> wings,std::vector<vlm::element::vortexRing> vortexRings,std::vector<Vector3d> nodes,Map<int,Vector3d> mapStruct,map<int, int> mapStructni) {
        for (int k = 0; k <wings.size() ; ++k)
        {
            for (int i = 0; i <wingStations.size() ; ++i)
            {

                auto wstation= wingStations[wings[0].get_stationsIDs()[k]];
                auto vring = vortexRings[wstation.get_vortexIDs()[i]];




                for (int p=0; p<vring.size(); ++p)
                {
                    auto node=nodes[vring.get_nodeIDs()[p]];
                    double bestDistance[2]={10000,10000};
                    double bestPoint[2]={10000,10000};
                    //compute the distance between the node on vlm and nodes on structure

                    for (auto s : mapStructni) {
                        auto nodeStruct = mapStruct[s.first];
                        double distance = sqrt(
                                pow(nodeStruct[0]-node[0], 2) + pow(nodeStruct[1]-node[1], 2) +
                                pow(nodeStruct[2]-node[2], 2));
                        if (distance < std::max(bestDistance[0], bestDistance[1])) {
                            if (distance< bestDistance[0]) {
                                bestDistance[1] = bestDistance[0];
                                bestDistance[0] = distance;
                                bestPoint[1] = bestPoint[0];
                                bestPoint[0] = s.first;
                            } else {
                                bestDistance[1] = distance;
                                bestPoint[1] = s.first;
                            }


                        }

                    }

                    //projection of the node on the structure
                    auto nodeStruct1 = mapStruct[bestPoint[0]];
                    auto nodeStruct2 = mapStruct[bestPoint[1]];

                    pos.node.push_back(bestPoint[0]);
                    pos.node.push_back(bestPoint[1]);

                    double x1[3];
                    double x2[3];
                    for (int i=0; i<3; ++i) {
                        x1[i] = nodeStruct2[i]-nodeStruct1[i];
                        x2[i] = node[i]-nodeStruct1[i];
                    }
                    double projection = (x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2])/pow(x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2],0.5);

                    pos.weight.push_back(1-projection);
                    pos.weight.push_back(projection);




                }
            }


        }
    }

    auto computeVLMDispalecement(interpolation pos,std::vector<vlm::surface::wingStation> wingStations,
                                 std::vector<vlm::surface::wing> wings,std::vector<vlm::element::vortexRing> vortexRings,std::vector<Vector3d> nodes,Map<int,Vector3d> mapStruct,map<int, int> mapStructni) {
        for (int k = 0; k<wings.size(); ++k) {
            for (int i = 0; i<wingStations.size(); ++i) {

                auto wstation = wingStations[wings[0].get_wingStationsIDs()[k]];
                auto vring = vortexRings[wstation.get_vortexIDs()[i]];
                int numPointsVort = vring.size();
                for (int p=0; p<numPointsVort; ++p)
                {
                    auto node = nodes[vring.get_nodeIDs()[p]];
                    auto nodeStruct1 = mapStruct[pos.node[2*p]];
                    auto nodeStruct2 = mapStruct[pos.node[2*p+1]];
                    auto weight1 = pos.weight[2*p];
                    auto weight2 = pos.weight[2*p+1];

                    Vector3d disp;
                    disp[0]=weight1*nodeStruct1[0]+weight2*nodeStruct2[0];
                    disp[1]=weight1*nodeStruct1[1]+weight2*nodeStruct2[1];
                    disp[2]=weight1*nodeStruct1[2]+weight2*nodeStruct2[2];
                    //prise en compte de la rotation



                    

                    node[0]=node[0]+disp[0];
                    node[1]=node[1]+disp[1];
                    node[2]=node[2]+disp[2];

                }
            }
        }
        return node ;
    }
   
        
      

    vector<double> crossProduct(const vector<double>& i, const vector<double>& j)
    {
             if (i.size() !=3 || j.size() != 3){
                 throw std::runtime_error("Vectors must have dimension 3.");
             }
             vector<double> result(3);
             result[0] = i[1]*j[2] - i[2]*j[1];
             result[1] = i[2]*j[0] - i[0]*j[2];
             result[2] = i[0]*j[1] - i[1]*j[0];
             return result;
    }
         
         
    double norme(vector<double> v)
    {
        double sum = 0;
        for(int i=0; i<v.size(); i++)
        {
         sum += pow(v[i], 2);
        }
        return sqrt(sum);
    }





}
