//
//  main.cpp
//  Module_geometrie
//
//  Created by Nader El shami on 2023-02-14.
//
#include <cstddef>
#include <fstream>
# include <iostream>
# include <vector>
# include <cmath>
# include <array>
#include "geometrie/euler.hpp"
#include "geometrie/geometry.hpp"
#include "geometrie/structure.hpp"

using namespace std;
int main(int argc, const char * argv[]) {

    std::string body_type = "General";
    Body WING_RIGHT(body_type);
    // WING_RIGHT.add_wing_surface_cst(
    //     100,                            // nx
    //     5,                             // ny
    //     10,                             // b, envergure totale
    //     0,                              // d
    //     {1.2, 0.7, 0.35, 0.3, 0.25},    // c vecteur pour la corde, 2
    //     {0.001},                        // z_te_half propre cst
    //     {0.05, 0.05, 0.01},             // r_le propre cst
    //     {0.4, 0.1},                     // beta propre cst, angle bord de fuite
    //     {0, 2},                         // x_le
    //     {0.05, 0.1, 0.4},               // z_n
    //     {0.05, 0.2}                     // delta_alpha_t, twist 
    // );
    WING_RIGHT.add_wing_surface_naca(
        100,                             // nx
        10,                            // ny
        10,                             // b
        0,                              // d
        {1.2, 0.7, 0.35, 0.3, 0.25},    // c
        {0},                         // m, 0 (1seul)
        {0},                         // p, 0 (1seul)
        {12},                       // t, 12 (1seul)
        {0, 2},                         // x_le
        {0.05, 0.1, 0.4},               // z_n
        {0.05, 0.2}                     // delta_alpha_t
    );

    Body WING_LEFT(body_type);
    WING_LEFT.add_wing_surface_cst(
        10,                            // nx
        5,                             // ny
        10,                             // b
        0,                              // d
        {1.2, 0.7, 0.35, 0.3, 0.25},    // c
        {0.001},                        // z_te_half
        {0.05, 0.05, 0.01},             // r_le
        {0.4, 0.1},                     // beta
        {0, 2},                         // x_le
        {0.05, 0.1, 0.4},               // z_n
        {0.05, 0.2}                     // delta_alpha_t
    );
    /*WING_LEFT.add_wing_surface_naca(
        70,                             // nx
        100,                            // ny
        10,                             // b
        0,                              // d
        {1.2, 0.7, 0.35, 0.3, 0.25},    // c
        {2, 4},                         // m
        {4, 2},                         // p
        {15, 23},                       // t
        {0, 2},                         // x_le
        {0.05, 0.1, 0.4},               // z_n
        {0.05, 0.2}                     // delta_alpha_t
    );*/
    WING_LEFT.mirror_body();

    WING_RIGHT.change_all_distributions("partial", "cartesian");
    std::vector<std::vector<std::vector<std::vector<double>>>> WR_surfaces = WING_RIGHT.get_paired_body_surfaces();    
    std::vector<std::vector<std::vector<std::vector<double>>>> WL_surfaces = WING_LEFT.get_paired_body_surfaces();    
    // std::vector<std::vector<std::vector<double>>> WR_surface_0 = WR_surfaces[0];
    // std::vector<std::vector<double>> X = WR_surface_0[0];
    // std::vector<std::vector<double>> Y = WR_surface_0[2];
    // std::vector<std::vector<double>> SU = WR_surface_0[4];
    // std::vector<std::vector<double>> SL = WR_surface_0[5];

    //  tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> A = WING_RIGHT.build_airfoil();
    // vector<vector<double>> X=get<0>(A);
    // vector<vector<double>> Y=get<1>(A);
    // vector<vector<double>> SU=get<2>(A);
    // vector<vector<double>> SL=get<3>(A);

    // generer(X, Y, SU, SL);  
    // vector for 1st element size 

    //MESH RANS

    std::vector<double> disc{100,150,200};
    cout<<disc[0]<<endl;
    std::vector<string> file_name{"Airfoil_coarse.msh","Airfoil_normal.msh","Airfoil_fine.msh"};
    cout<<file_name[0]<<endl;
    std::vector<std::vector<std::vector<std::vector<double>>>> surfaces = WR_surfaces;
    bool RANS = false; // True = solver RANS, False = solver Euler
    for (int i=0; i < surfaces.size(); i++){
        for (int j=0; j<3; j++){
            generer(surfaces[i][0], surfaces[i][1], surfaces[i][2], surfaces[i][4], surfaces[i][5], disc[j], file_name[j], RANS);
        }
    }

    vector<tuple<int,vector<double>,vector<double>,vector<double>>> element = maillage_structure(WING_RIGHT);


    
    //std::vector<std::vector<std::vector<double>>> xu = surfaces[i][0][0];
    //cout<<xu<<endl; 
    // int n = surfaces[0][0][0].size();
    // cout<<n<<endl; 
    // int c = surfaces[0][1][0].size();
    // cout<<c<<endl; 
    //int n = 10;
    // double maxi = *max_element(surfaces[0][0][0].begin(), surfaces[0][0][0].end());
    // cout<<maxi<<endl; 
    // double mini = *max_element(surfaces[0][1][0].begin(), surfaces[0][1][0].end());
    // cout<<mini<<endl; 
    // double min = *min_element(surfaces[0][1][0].begin(), surfaces[0][1][0].end());
    // cout<<min<<endl;
    // double mm=max(maxi,min);
    // cout<<mm<<endl; 
    // double chord = mm-min;
    // double c4_x=min+(chord/4);
    // cout<<c4_x<<endl; 

  

    return 0;
    }
