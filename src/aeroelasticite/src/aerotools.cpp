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

using namespace Eigen;
using namespace vlm;
using namespace structure;


namespace aero{
    struct interpolation;
    struct interpolation_f;


    void DispInterpol( interpolation &pos,std::vector<surface::wingStation> wingStations,
                       std::vector<surface::wing> wings,std::vector<element::vortexRing> vortexRings,
                      std::vector<surface::wing> wings,std::vector<element::vortexRing> vortexRings,std::vector<Vector3d> nodes,map<int,Vector3d> mapStruct,map<int, int> mapStructni) {
        for (k = 0; k <wings.size() ; ++k)
        {
            for (i = 0; i <wingStations.size() ; ++i)
            {

                auto wstation= wingStations[wings[0].get_wingStationsIDs()[k]];
                auto vring = vortexRings[wstation.get_vortexIDs()[i]];




                for (p=0; p<vring.size(); ++p)
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

                    pos.node.push_back(mapStructni[bestPoint[0]);
                    pos.node.push_back(mapStructni[bestPoint[1]);

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

    auto computeVLMDispalecement(interpolation pos,std::vector<surface::wingStation> wingStations,
                                 std::vector<surface::wing> wings,std::vector<element::vortexRing> vortexRings,std::vector<Vector3d> nodes) {
        for (k = 0; k<wings.size(); ++k) {
            for (i = 0; i<wingStations.size(); ++i) {

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
    void LoadInterpol(interpolation_f &force, std::vector<surface::wingStation> wingStations,
                       std::vector<surface::wing> wings,std::vector<element::vortexRing> vortexRings,std::vector<Vector3d> nodes,map<int,Vector3d> mapStruct,map<int, int> mapStructni)
   {
        // Loading VLM forceacting Point
        
        int n=wings.size(); // nombre de wingstations=size (wings)
        int m=wingStations.size(); // nombre de vortexRings sur une wingstation
        int i=0 // numéro du premier vortexring d une wingstation
        
        for (int j=0; j<n; ++j)
        {
            auto wingstation = wingStations[wings[0].get_wingStationsIDs()[j]];
            
            auto vortexring= vortexRings[wstation.get_vortexIDs()[i]];
            force.point_fa.push.back(vortexRing[i].ForcesActingPoint(nodes, vortexRings));
            
            i= i + m;
        }
        
        //    Loading STRUCTURE connectivity

        for (int k=0, k<map.size(), ++k)
        {
            force.point_fs.push.back(vector<double>(3));

        }


        for (auto s : mapStructni) 
        {
             auto nodeStruct = mapStruct[s.first];
             force.point_fs[3*s]=nodeStruct[0];
             force.point_fs[3*s+1]=nodeStruct[1];
             force.point_fs[3*s+2]=nodeStruct[2];
        }
        
        
        
        // Projeter le point_fa sur la droite et calcul de la distance
        vector <vector<double>> point_fp;
        vector <double> dist(n);
        vector <double> v(3);
        
        vector <double> u(3); // vecteur directeur de la poutre
        u[0]= force.point_fs[0] - force.point_fs[3];
        u[1]= force.point_fs[1] - force.point_fs[4];
        u[2]= force.point_fs[2] - force.point_fs[5];
        
        for (int i=0; i<n; ++i){
            point_fp.push.back(vector<double> (3));
            double t;
            t= (u[0](force.point_fa[3*i]-force.point_fs[0]) + u[1](force.point_fa[3*i+1]-force.point_fs[1]) + 
                   u[2](force.point_fa[3*i+2]-force.point_fs[2]))/pow(norme(u),2);
            
            point_fp[i][0]=force.point_fs[0] + t*u[0];
            point_fp[i][1]=force.point_fs[1] + t*u[1];
            point_fp[i][2]=force.point_fs[2] + t*u[2];
            
            v[0]= force.point_fa[3*i]   -   point_fp[i][0];
            v[1]= force.point_fa[3*i+1] -   point_fp[i][1];
            v[2]= force.point_fa[3*i+2] -   point_fp[i][2];
            dist[i]=norme(v);
            
        }
        
        // Calcul des coefficients d'interpolation
        ///vector <vector> forces_s[m][6];
        double dist_0;
        double dist_1;
        double dist_2;
        
    
        
        for (int i=0; i<n; ++i){
            ///auto forces = model.forces_to_inertial(i);
            force.poids.pushback(vector<double>(3));
            force.poids[i][0]=i;
            
            
            for (int j=0; j<map.size(); ++j){
                v[0]= force.point_fs[3*j] - point_fs[3*(j+1)];
                v[1]= force.point_fs[3*j+1] - point_fp[3*(j+1)+1];
                v[2]= force.point_fs[3*j+2] - point_fp[3*(j+1)+2];
                dist_0=norme(v);
                
                v[0]= force.point_fs[3*j]   - point_fp[i][0];
                v[1]= force.point_fs[3*j+1] - point_fp[i][1];
                v[2]= force.point_fs[3*j+2] - point_fp[i][2];
                dist_1=norme(v);
                
                v[0]= force.point_fs[3*(j+1)]   - point_fp[i][0];
                v[1]= force.point_fs[3*(j+1)+1] - point_fp[i][1];
                v[2]= force.point_fs[3*(j+1)+2] - point_fp[i][2];
                dist_2=norme(v);
                
                if (dist_1<=dist_0 && j != map.size()-1)
                {
                    epsilon_l= dist_2/dist_0;
                    epsilon_r= dist_1/dist_0;
                    force.poids[3*i]=j;
                    force.poids[3*i+1]=epsilon_l;
                    force.poids[3*i+2]=epsilon_r;
                    
                }
                else if (j==map.size()-1){
                    force.poids[3*i]=j;
                    force.poids[3*i+1]=1;
                    force.poids[3*i+2]=0;
                    
                }
            }
            
        }
   }

   auto ComputeStructureForces(interpolation_f force,Matrix<double, 6, 1> forces_to_inertial_frame)

   {
        n=force.point_fa.size();
        m=force.point_fs.size();
        
        vector <double> forces_s(6*m,0);
        vector <double> M(3);
        vector <double> r(3);
        
        for (int i=0; i<n; ++i)
        {

            auto forces = object.forces_to_inertial_frame(i);
            j=force.poids[3*i];
            if (j!=m-1)
            {
            
                forces_s[6*j]   += force.poids[3*i+1] * forces[i][0];
                forces_s[6*j+1] += force.poids[3*i+1] * forces[i][1];
                forces_s[6*j+2] += force.poids[3*i+1] * forces[i][2];
                
                forces_s[6*(j+1)]   += force.poids[3*i+2] * forces[i][0];
                forces_s[6*(j+1)+1] += force.poids[3*i+2] * forces[i][1];
                forces_s[6*(j+1)+2] += force.poids[3*i+2] * forces[i][2];
                
                M[0]=forces[i][3];
                M[1]=forces[i][4];
                M[2]=forces[i][5];
                
                r[0]=force.point_fa[3*i] -   force.point_fs[3*j];
                r[1]=force.point_fa[3*i+1] - force.point_fs[3*j+1];
                r[2]=force.point_fa[3*i+2] - force.point_fs[3*j+2];
                
                M_s=crossProduct (r,M);
                
                forces_s[6*j+3]   += force.poids[3*i+1] * (M_s[0] + forces[i][3]);
                forces_s[6*j+4]   += force.poids[3*i+1] * (M_s[1] + forces[i][4]);
                forces_s[6*j+5]   += force.poids[3*i+1] * (M_s[2] + forces[i][5]);
                
                r[0]=force.point_fa[3*i] -   force.point_fs[3*(j+1)];
                r[1]=force.point_fa[3*i+1] - force.point_fs[3*(j+1)+1];
                r[2]=force.point_fa[3*i+2] - force.point_fs[3*(j+1)+2];
                
                M_s=crossProduct (r,M);
                
                forces_s[6*(j+1)+3]   += force.poids[3*i+2] * (M_s[0] + forces[i][3]);
                forces_s[6*(j+1)+4]   += force.poids[3*i+2] * (M_s[1] + forces[i][4]);
                forces_s[6*(j+1)+5]   += force.poids[3*i+2] * (M_s[2] + forces[i][5]);
            
            }

            else if (j==m-1)
            {
                forces_s[6*j]   += force.poids[i][1] * forces[i][0];
                forces_s[6*j+1] += force.poids[i][1] * forces[i][1];
                forces_s[6*j+2] += force.poids[i][1] * forces[i][2];

                M[0]=forces[i][3];
                M[1]=forces[i][4];
                M[2]=forces[i][5];
                
                r[0]=force.point_fa[3*i] -   force.point_fs[3*j];
                r[1]=force.point_fa[3*i+1] - force.point_fs[3*j+1];
                r[2]=force.point_fa[3*i+2] - force.point_fs[3*j+2];
                
                M_s=crossProduct (r,M);
                
                forces_s[6*j+3]   += force.poids[i][1] * (M_s[0] + forces[i][3]);
                forces_s[6*j+4]   += force.poids[i][1] * (M_s[1] + forces[i][4]);
                forces_s[6*j+5]   += force.poids[i][1] * (M_s[2] + forces[i][5]);
                
            }
            
            
        }
        return forces_s;
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
