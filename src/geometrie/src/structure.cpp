#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "geometrie/structure.hpp"
#include "geometrie/geometry.hpp"

using namespace std;

//Calcul du point centre
//Centre placé au quart de corde
vector<vector<double>> get_center_point(vector<vector<double>> X, vector<vector<double>> Y, vector<vector<double>> Su, vector<vector<double>> Sl){
    vector<vector<double>> center(X.size(), vector<double>(3,0));
    
    double dim = X[0].size()-1;
    
    for (int i=0; i<X.size(); i++){
        center[i][0] = X[i][0] + 0.25*(X[i][dim]-X[i][0]);
        center[i][1] = Y[i][0] + 0.25*(Y[i][dim]-Y[i][0]);
        center[i][2] = Su[i][0] + 0.25*(Su[i][dim]-Su[i][0]);

        }
    return center;
    
}

vector<double> vecteur_normal(vector<double> p1, vector<double> p2){
    vector<double> p3 = p1;
    p3[0]= p3[0]+1;
    vector<double> normal(3,0), v1(3,0), v2(3,0), normal_unitaire(3,0);
    double norme;
    for (int i=0; i<= v1.size(); i++){
        v1[i]=p2[i]-p1[i];
        v2[i]=p3[i]-p1[i];
    }
    
    normal[0]=(v1[1]*v2[2])-(v2[1]*v1[2]);
    normal[1]=-((v1[0]*v2[2])-(v2[0]*v1[2]));
    normal[2]=-(v1[0]*v2[1])-(v2[0]*v1[1]);
    
    norme=sqrt(pow(normal[0],2)+pow(normal[1],2)+pow(normal[2],2));
    
    normal_unitaire[0]=normal[0]/norme;
    normal_unitaire[1]=normal[1]/norme;
    normal_unitaire[2]=normal[2]/norme;
    return normal_unitaire;
}


vector<vector<vector<vector<double>>>> get_geometry(Body wing){
    vector<vector<vector<vector<double>>>> surfaces = wing.get_paired_body_surfaces();
    return surfaces;
}

vector<tuple<int,vector<double>,vector<double>,vector<double>>> maillage_structure(Body wing){
    vector<tuple<int,vector<double>,vector<double>,vector<double>>> element;
    vector<double> normal(3,0),pt_normal(3,0);
    
    vector<vector<vector<vector<double>>>> surfaces = get_geometry(wing);
    
    vector<vector<double>> X_U, X_L, Y_U, Y_L, SU, SL;
    for (int m=0; m < surfaces[0][0].size()-1; m++){
        X_U.push_back(surfaces[0][0][m]);
        X_L.push_back(surfaces[0][1][m]);
        Y_U.push_back(surfaces[0][2][m]);
        X_L.push_back(surfaces[0][3][m]);
        SU.push_back(surfaces[0][4][m]);
        SL.push_back(surfaces[0][5][m]);

        
    }
    
    if (surfaces.size()>1){
        for (int k=0; k < surfaces.size()-1; k++){
            vector<vector<vector<double>>> surface = surfaces[k+1];
            for (int m=0; m < surface[0].size()-1; m++){
                X_U.push_back(surface[0][m+1]);
                X_L.push_back(surface[1][m+1]);
                Y_U.push_back(surface[2][m+1]);
                X_L.push_back(surface[3][m+1]);
                SU.push_back(surface[4][m+1]);
                SL.push_back(surface[5][m+1]);
                
            }
        }
    }
    
    

        
    vector<vector<double>> centre = get_center_point(X_U, Y_U, SU, SL);
    
    
    
    
    
    ofstream myfile;
    myfile.open ("Point_maillage_structure.txt");
    


        
        myfile<<"$##### NOEUDS #####"<<endl;
        myfile<<setprecision(5)<<fixed;
        for(int i=0; i<=(centre.size()-1); i++){
            myfile<<"GRID, "<<i+100<<", ,"<<setprecision(5)<<centre[i][0]<<", "<<setprecision(5)<<centre[i][0]<<", "<<setprecision(5)<<centre[i][2]<<", ,"<<"345"<<endl;
        }
        myfile <<""<<endl;
        myfile <<"$##### Connectivity #####"<<endl;
        for(int i=0; i<=(centre.size()-2); i++){
            normal=vecteur_normal(centre[i+1], centre[i]);
            pt_normal[0]=centre[i][0]+normal[0];
            pt_normal[1]=centre[i][1]+normal[1];
            pt_normal[2]=centre[i][2]+normal[2];
            myfile<<setprecision(5)<<fixed;
            myfile<<"CBAR, "<<i+1<<", "<<"10, "<<i+100<<", "<<i+101<<", "<<normal[0]<<", "<<normal[1]<<", "<<normal[2]<<endl;
            
            element.push_back(make_tuple(i,centre[i],centre[i+1],pt_normal));

        }
    
    myfile<<""<<endl;
    myfile<<"$##### Propriétés de sections #####"<<endl;
    myfile<<"PBAR, "<<"10, "<<"11, "<<"1.0, "<<"0.1, "<<"0.1, "<<"0.1"<<endl;
    myfile<<"MAT1, "<<"11, "<<"140.E9, "<<"70.E9,"<<endl;
    myfile <<""<<endl;
    myfile<<"$##### Contrainte #####"<<endl;
    myfile<<"SPC1, "<<"100, "<<"123456, "<<"100,"<<endl;
    myfile <<""<<endl;
    myfile<<"$##### Charges #####"<<endl;
    myfile<<"MOMENT, "<<"103, "<<"112, "<<", "<<"3.E4, "<<"0.0, "<<"0.0, "<<"1.0"<<endl;
    myfile <<""<<endl;
    myfile<<"$"<<endl;
    myfile<<"ENDDATA "<<endl;

    myfile.close();
    
    return element;
}



