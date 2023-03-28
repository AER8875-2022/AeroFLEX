
#ifndef __GEOM__
#define __GEOM__

#include <cstddef>
#include <fstream>
# include <iostream>
# include <vector>
# include <cmath>
# include <array>

#include "common_aeroflex.hpp"
//les mettre ici ?
#include "geometrie/euler.hpp"
#include "geometrie/geometry.hpp"
#include "geometrie/structure.hpp"
#include "tinyconfig.hpp"

namespace geom {

struct Settings {
    //General
    double cr = 1.0;          //Chord Root
    double ct = 1.0;          //Chords Tip
    double envergure = 30.0;   //total Span

    //Angles
    double twist = 0.0;       //twist angle
    double fleche = 0.0;      //sweep angle
    double dihedre = 0.0;     //diherdral angle

    // Wing position
    double P_beam = 0;      //Wing structural beam position
    double P_aile = 0;      //Wing position on fuselage

    //NACA 4-digit
    double m = 0;
    double p = 0;
    double t = 12;

    //CST
    double z_te = 0.001;        //Trailing edge distance betweem 2 last points
    double r_le = 0.05;        //Leading edge radius
    double Beta = 0.2;        //Leading edge angle

    //Winglet;
    bool Winglet = 0;       //1-> yes, 0->no

    //Airfoil Type
    bool S_type = 0;        //1-> CST, 0->NACA

    //import
    void import_config_file(tiny::config &io);
    void export_config_file(tiny::config &io);
};

class Geom {

public:
    GUIHandler &gui;
    Settings settings;

public:
    Geom(GUIHandler &gui);
    void Geom_gen();
};

}

#endif
