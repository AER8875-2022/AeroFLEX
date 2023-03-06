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
    
auto LoadInterpol(){
    
    // Loading VLM forceacting Point
    vector<vector <double>> point_fa;
    int n=size.wings(); // nombre de wingstations=size (wings)
    int m=size.wingstation; // nombre de vortexRings sur une wingstation
    int i=0 // numéro du premier vortexring d une wingstation
    
    for (int j=0; j<=n; ++j){
        auto wingstation = wingStations[wings[0].get_wingStationsIDs()[j]];
        m=sizeowingstation);
        
        auto vortexring= vortexRings[wstation.get_vortexIDs()[i]];
        point_fa.push.back(vortexRing[i].ForcesActingPoint(nodes, vortexRings));
        
        i= i + m;
        }
    
    //    Loading STRUCTURE connectivity
    vector<vector <double>> point_fs;
    
    

    // Projeter le point_fa sur la droite et calcul de la distance
    vector<vector<double>> point_fp [n];
    vector <double> dist[n];
    vector <double> v[3]
    
    vector <double> u(3); // vecteur directeur de la poutre
    u[1]= point_fs[2][1] - point_fs[1][1];
    u[2]= point_fs[2][2] - point_fs[1][2];
    u[3]= point_fs[2][3] - point_fs[1][3];
    
    for (int i=0; i<=n; ++i){
        double t;
        t= (u[1](point_fa[i][1]-point_fs[1][1]) + u[2](point_fa[i][2]-point_fs[1][2]) + u[1](point_fa[i][1]-point_fs[1][1]))/norme(u);
        
        point_fp[i][1]=point_fs[1][1] + t*u[1];
        point_fp[i][2]=point_fs[1][2] + t*u[2];
        point_fp[i][3]=point_fs[1][3] + t*u[3];
        
        v[1]= point_fa[i][1] - point_fp[i][1];
        v[2]= point_fa[i][2] - point_fp[i][2];
        v[3]= point_fa[i][3] - point_fp[i][3];
        dist[i]=norme(v);
        
     }
    
    // Calcul des coefficients d'interpolation
    vector <vector> forces_s[m][6];
    double dist_0;
    double dist_1;
    double dist_2;
    
    for (int i=0; i<n; ++i){
        auto forces = model.forces_to_inertial(i);

        
        for (int j=0; j<m; ++j){
            v[1]= point_fs[j][1] - point_fs[j+1][1];
            v[2]= point_fs[j][2] - point_fp[j+1][2];
            v[3]= point_fs[j][3] - point_fp[j+1][3];
            dist_0=norme(v);
            
            v[1]= point_fs[j][1] - point_fp[i][1];
            v[2]= point_fs[j][2] - point_fp[i][2];
            v[3]= point_fs[j][3] - point_fp[i][3];
            dist_1=norme(v);
            
            v[1]= point_fs[j+1][1] - point_fp[i][1];
            v[2]= point_fs[j+1][2] - point_fp[i][2];
            v[3]= point_fs[j+1][3] - point_fp[i][3];
            dist_2=norme(v);
            
            if dist_1<=dist_0{
                epsilon_l= dist_2/dist_0;
                epsilon_r= dist_1/dist_0;
                
                forces_s[j][1]= forces_s[j][1] + epsilon_l*forces[i][1];
                forces_s[j][2]= forces_s[j][2] + epsilon_l*forces[i][2];
                forces_s[j][3]= forces_s[j][3] + epsilon_l*forces[i][3];
                forces_s[j][4]= forces_s[j][4] + epsilon_l*forces[i][4]*dist[i];
                forces_s[j][5]= forces_s[j][5] + epsilon_l*forces[i][5]*dist[i];
                forces_s[j][6]= forces_s[j][6] + epsilon_l*forces[i][6]*dist[i];
                
                forces_s[j+1][1]= forces_s[j+1][1] + epsilon_l*forces[i][1];
                forces_s[j+1][2]= forces_s[j+1][2] + epsilon_l*forces[i][2];
                forces_s[j+1][3]= forces_s[j+1][3] + epsilon_l*forces[i][3];
                forces_s[j+1][4]= forces_s[j+1][4] + epsilon_l*forces[i][4]*dist[i];
                forces_s[j+1][5]= forces_s[j+1][5] + epsilon_l*forces[i][5]*dist[i];
                forces_s[j+1][6]= forces_s[j+1][6] + epsilon_l*forces[i][6]*dist[i];
                
            }
            
        }
     
    }
      
}
        double dotProduct(const vector<double>& a, const vector<double>& b) {
            if (a.size() != b.size()) {
                throw std::runtime_error("Vectors must have the same dimension.");
            }
            double result = 0.0;
            for (size_t i = 0; i < a.size(); i++) {
                result += a[i] * b[i];
            }
            return result;
         }
         vector<double> crossProduct(const vector<double>& i, const vector<double>& j){
             if (i.size() !=3 || j.size() != 3){
                 throw std::runtime_error("Vectors must have dimension 3.");
             }
             vector<double> result(3);
             result[0] = i[1]*j[2] - i[2]*j[1];
             result[1] = i[2]*j[0] - i[0]*j[2];
             result[2] = i[0]*j[1] - i[1]*j[0];
             return result;
         }
         double norme(vector<double> v){
             double sum = 0;
             for(int i=0; i<v.size(); i++){
                 sum += pow(v[i], 2);
             }
             return sqrt(sum);
         }
