
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

namespace geom {

struct Settings {
    //General
    double cr;          //Chord Root
    double ct;          //Chords Tip
    double envergure;   //total Span

    //Angles
    double twist;       //twist angle
    double fleche;      //sweep angle
    double dihedre;     //diherdral angle

    // Wing position
    double P_beam;      //Wing structural beam position
    double P_aile;      //Wing position on fuselage

    //NACA 4-digit
    double m;
    double p;
    double t;

    //CST
    double z_te;        //Trailing edge distance betweem 2 last points
    double r_le;        //Leading edge radius
    double Beta;        //Leading edge angle

    //Winglet;
    bool Winglet;       //1-> yes, 0->no

    //Airfoil Type
    bool S_type;        //1-> CST, 0->NACA
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
