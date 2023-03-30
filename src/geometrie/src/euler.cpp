#include "geometrie/euler.hpp"
#include "geometrie/geometry.hpp"
#include <algorithm>

void generer(std::vector<std::vector<double>> xu, std::vector<std::vector<double>> xl, std::vector<std::vector<double>> y, std::vector<std::vector<double>> su, std::vector<std::vector<double>> sl, double disc, std::string file_name, bool RANS){

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal",0); //Mute GMSH on terminal

  gmsh::model::add("Euler");

    // initialization
    int indx=0;               // slice selection
    int n = xu[indx].size();  // number of points

    int indx_pt = 1;          // point index initialisation
    int indx_line = 1;        // line index initialisation
    int indx_phys=1;          // Curve index initialisation

    // Quarter Chord Info
    auto max_xu = std::max_element(xu[indx].begin(), xu[indx].end());
    auto max_xl = std::max_element(xu[indx].begin(), xu[indx].end());
    double max_x=std::max(*max_xu,*max_xl);
    auto min_x = std::min_element(xu[indx].begin(), xu[indx].end());
    double chord = max_x-*min_x;
    double c4_x=*min_x+(chord/4);
    double su_int=1;
    double sl_int=1;
    bool xu_int = true;
    bool xl_int = true;
    for (int i=0; i < n; i++){
        if (xu[indx][i] > c4_x && xu_int) {
            su_int=su[indx][i-1]+(c4_x-xu[indx][i-1])*((su[indx][i]-su[indx][i-1])/(xu[indx][i]-xu[indx][i-1]));
            xu_int = false;
            }
        if (xl[indx][i] > c4_x && xl_int) {
            sl_int=sl[indx][i-1]+(c4_x-xl[indx][i-1])*((sl[indx][i]-sl[indx][i-1])/(xl[indx][i]-xl[indx][i-1]));
            xl_int = false;
            }
    }
    double c4_z=(su_int+sl_int)/2;

    // First element size
    double lc = chord/disc;         // mesh element size

    // Point definition
  for (int i=0; i < n; i++){
    gmsh::model::geo::addPoint(xu[indx][i], su[indx][i], 0, lc, indx_pt);
    if (i != n-1) {
        gmsh::model::geo::addPoint(xl[indx][n-1-i], sl[indx][n-1-i], 0, lc, n+indx_pt);
    }
    indx_pt++;
    }
    indx_pt=2*n;

    // Airfoil Points Tag List
    std::vector<int> point_list;
    for (int i = 1; i < 2*n; i++)
        point_list.push_back(i);

    // Fan Points Tag List
    std::vector<double> fan_point_list;
    fan_point_list.push_back(n);
    fan_point_list.push_back(n+1);


    // Line definition
    for (int i = 1; i < 2*n; i++){
        if (i==2*n-1){
            gmsh::model::geo::addLine(i, 1, indx_line);
        }
        else {
            gmsh::model::geo::addLine(i, i+1, indx_line);
        }
        indx_line++;
    }
    //std::vector<int> new_list{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

    // Line tag list
    std::vector<int> lin_list;
    for (int i = 1; i < 2*n; i++)
        lin_list.push_back(i);

    std::vector<double> lin_list2;
    for (int i = 1; i < 2*n; i++)
        lin_list2.push_back(i);



    // center Point (quarter chord)
    gmsh::model::geo::addPoint(c4_x, c4_z, 0, lc, indx_pt);
    indx_pt++;

    // Create Circle Points
    gmsh::model::geo::addPoint(c4_x, c4_z+50*chord, 0,lc*500, indx_pt);
    indx_pt++;
    gmsh::model::geo::addPoint(c4_x, c4_z-50*chord, 0,lc*500, indx_pt);
    indx_pt++;

    // Create Circle curves
     gmsh::model::geo::addCircleArc(indx_pt-1, indx_pt-3, indx_pt-2, indx_line);
     indx_line++;
     gmsh::model::geo::addCircleArc(indx_pt-2, indx_pt-3, indx_pt-1, indx_line);
     indx_line++;

    //Curve loop definition
    gmsh::model::geo::addCurveLoop(lin_list, 101);                    // airfoil
    gmsh::model::geo::addCurveLoop({indx_line-1, indx_line-2}, 102);  //domain

    // Surface definition, airfoil
    gmsh::model::geo::addPlaneSurface({102, 101}, 1);

    // Physical curve
    gmsh::model::addPhysicalGroup(1, lin_list, indx_phys,"Airfoil");          //physical curve, airfoil
    indx_phys++;
    gmsh::model::addPhysicalGroup(1, {indx_line-1, indx_line-2}, indx_phys,"Farfield");          //physical curve, farfield
    indx_phys++;
    gmsh::model::addPhysicalGroup(2, {1}, indx_phys, "Internal");    //phyisical surface, airfoil
    indx_phys++;

    // Boundary Layer for RANS solver
    if (RANS) {
         gmsh::model::mesh::field::add("BoundaryLayer", 1);
        gmsh::model::mesh::field::setNumbers(1, "CurvesList", lin_list2);
        gmsh::model::mesh::field::setNumber(1, "Thickness", 0.01);
        gmsh::model::mesh::field::setNumber(1, "Size", 0.0001);
        gmsh::model::mesh::field::setNumber(1, "Ratio", 1.25);
        gmsh::model::mesh::field::setNumber(1, "Quads", 1);
        gmsh::model::mesh::field::setNumbers(1, "FanPointsList", fan_point_list);
        gmsh::model::mesh::field::setAsBoundaryLayer(1);
    }




    gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(2);

  gmsh::write(file_name);

  gmsh::finalize();

}
