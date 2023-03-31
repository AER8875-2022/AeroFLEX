
#ifndef __GEOM__
#define __GEOM__

#include <vector>
#include <string>

#include "common_aeroflex.hpp"
//les mettre ici ?

#include "tinyconfig.hpp"

namespace geom {

struct Settings {

    //Airfoil Type
    int S_type = 1;         //1-> CST, 0->NACA
    std::vector<std::string> solver_options = {"NACA", "CST"};

    inline std::string solver_type() {return solver_options.at(S_type);};
    inline void set_solver_type(const std::string& type_) {
        if (type_ == "NACA") {
            S_type = 0;
        } else if (type_ == "CST") {
            S_type = 1;
        }
    };

    // Winglet
    int Winglet = 0;       //1-> yes, 0->no
    std::vector<std::string> Winglet_options = {"OUI", "NON"};

    inline std::string Winglet_type() {return Winglet_options.at(Winglet);};
    inline void set_Winglet_type(const std::string& Winglet_) {
        if (Winglet_ == "NON") {
            Winglet = 0;
        } else if (Winglet_ == "OUI") {
            Winglet = 1;
        }
    };

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

    //Propriétés du matériaux
    double E = 1000000000.0;
    double G = 1000000.0;

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