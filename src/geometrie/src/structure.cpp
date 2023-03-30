#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <tuple>
#include <array>
#include "geometrie/structure.hpp"
#include "geometrie/geometry.hpp"


//Calcul du point centre
//Centre placé au quart de corde
std::vector<std::vector<double>> get_center_point(std::vector<std::vector<double>> X, std::vector<std::vector<double>> Y, std::vector<std::vector<double>> Su, std::vector<std::vector<double>> Sl){
    std::vector<std::vector<double>> center(X.size(), std::vector<double>(3,0));

    double dim = X[0].size()-1;

    for (int i=0; i<X.size(); i++){
        center[i][0] = X[i][0] + 0.25*(X[i][dim]-X[i][0]);
        center[i][1] = Y[i][0] + 0.25*(Y[i][dim]-Y[i][0]);
        center[i][2] = Su[i][0] + 0.25*(Su[i][dim]-Su[i][0]);

        }
    return center;

}

void vecteur_normal(std::array<double, 3> &normal, std::vector<double> p1, std::vector<double> p2){
    std::vector<double> v1(3,0), v2(3,0.);
    double norme;
    for (int i=0; i<= v1.size(); i++){
        v1[i]=p2[i]-p1[i];
    }
    v2[0]=1;

    normal[0]=(v1[1]*v2[2])-(v2[1]*v1[2]);
    normal[1]=-((v1[0]*v2[2])-(v2[0]*v1[2]));
    normal[2]=-(v1[0]*v2[1])-(v2[0]*v1[1]);

    norme=sqrt(pow(normal[0],2)+pow(normal[1],2)+pow(normal[2],2));

    normal[0] /= norme;
    normal[1] /= norme;
    normal[2] /= norme;
}

std::vector<std::vector<std::vector<std::vector<double>>>> get_geometry(Body wing){
    std::vector<std::vector<std::vector<std::vector<double>>>> surfaces = wing.get_paired_body_surfaces();
    return surfaces;
}

std::vector<std::tuple<int,std::vector<double>,std::vector<double>,std::vector<double>>> maillage_structure(Body wing){
    std::vector<std::tuple<int,std::vector<double>,std::vector<double>,std::vector<double>>> element;
    std::vector<double> pt_normal(3,0);

    std::vector<std::vector<std::vector<std::vector<double>>>> surfaces = get_geometry(wing);

    std::vector<std::vector<double>> X_U, X_L, Y_U, Y_L, SU, SL;
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
            std::vector<std::vector<std::vector<double>>> surface = surfaces[k+1];
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




    std::vector<std::vector<double>> centre = get_center_point(X_U, Y_U, SU, SL);





    std::ofstream myfile;
    myfile.open ("Point_maillage_structure.txt");




        myfile<<"$##### NOEUDS #####"<<std::endl;
        myfile<<std::setprecision(5)<<std::fixed;
        for(int i=0; i<=(centre.size()-1); i++){
            myfile<<"GRID, "<<i+100<<", ,"<<std::setprecision(5)<<centre[i][0]<<", "<<std::setprecision(5)<<centre[i][0]<<", "<<std::setprecision(5)<<centre[i][2]<<", ,"<<"345"<<std::endl;
        }
        myfile <<""<<std::endl;
        myfile <<"$##### Connectivity #####"<<std::endl;
        std::array<double, 3> normal = {0., 0., 0.};
        for(int i=0; i<=(centre.size()-2); i++){
            vecteur_normal(normal, centre[i+1], centre[i]);
            pt_normal[0]=centre[i][0]+normal[0];
            pt_normal[1]=centre[i][1]+normal[1];
            pt_normal[2]=centre[i][2]+normal[2];
            myfile<<std::setprecision(5)<<std::fixed;
            myfile<<"CBAR, "<<i+1<<", "<<"10, "<<i+100<<", "<<i+101<<", "<<normal[0]<<", "<<normal[1]<<", "<<normal[2]<<std::endl;

            element.push_back(std::make_tuple(i,centre[i],centre[i+1],pt_normal));
        }

    myfile<<""<<std::endl;
    myfile<<"$##### Propriétés de sections #####"<<std::endl;
    myfile<<"PBAR, "<<"10, "<<"11, "<<"1.0, "<<"0.1, "<<"0.1, "<<"0.1"<<std::endl;
    myfile<<"MAT1, "<<"11, "<<"140.E9, "<<"70.E9,"<<std::endl;
    myfile <<""<<std::endl;
    myfile<<"$##### Contrainte #####"<<std::endl;
    myfile<<"SPC1, "<<"100, "<<"123456, "<<"100,"<<std::endl;
    myfile <<""<<std::endl;
    myfile<<"$##### Charges #####"<<std::endl;
    myfile<<"MOMENT, "<<"103, "<<"112, "<<", "<<"3.E4, "<<"0.0, "<<"0.0, "<<"1.0"<<std::endl;
    myfile <<""<<std::endl;
    myfile<<"$"<<std::endl;
    myfile<<"ENDDATA "<<std::endl;

    myfile.close();

    return element;
}



