/*

        ____     ____  ___    _   _______
       /  _/    / __ \/   |  / | / / ___/
       / /_____/ /_/ / /| | /  |/ /\__ \
     _/ /_____/ _, _/ ___ |/ /|  /___/ /
    /___/    /_/ |_/_/  |_/_/ |_//____/

    Implicit RANS Solver
    Description: One file include for I-RANS functionality
    By: Alexis Angers

*/
#pragma once

#include <rans/multigrid.h>
#include <rans/post.h>
#include "common_aeroflex.hpp"
#include "database/database.hpp"

#include <string>

namespace rans {

class Rans {
    public:
    // This is a mess because of no separation of IO and computation...
    std::vector<double> residuals = {1.0};
    std::atomic<int> iters = 0;

    CpProfile profile;

    Settings settings;

    std::vector<mesh> ms;
    std::vector<double> alphas;

    bool mesh_loaded = false;

    GUIHandler &gui;

    Rans(GUIHandler &gui) : gui(gui) {};

    void input();
    void solve();
    void solve_airfoil(const std::string& airfoil, database::airfoil& db);
    void compute_alphas();

    template<class T> void run();
    template<class T> void run_airfoil(const std::string& airfoil, database::airfoil& db);
};

void Rans::compute_alphas() {
    alphas.clear();
    for (double i = settings.alpha_start; i <= settings.alpha_end; i += settings.alpha_step) {
        alphas.push_back(i);
    }
}

void Rans::input() {
    if (!mesh_loaded) {
        for (const auto& mesh_name : settings.meshes) {
            ms.push_back(mesh(mesh_name));
        }
        mesh_loaded = true;
    }
};

template<class T>
void Rans::run() {
    multigrid<T> multi(ms, settings, gui, residuals, iters, profile);

    rans::solver& s = multi.run(true);

    save(settings.outfilename, s);
    std::cout << "Saved results to file " << settings.outfilename << "\n" << std::endl;
}

template<class T>
void Rans::run_airfoil(const std::string& airfoil, database::airfoil& db) {
    gui.msg.push("[RANS] Solving airfoil: " + airfoil);
    ms.clear();
    // TODO: check with geom for the naming convention
    // For the moment we will only load 1 mesh
    std::cout<<"Nom du profil RANS :";
    std::cout<<airfoil<<std::endl;
    
    ms.push_back(mesh("./"+airfoil+"_coarse.msh"));
    ms.push_back(mesh("./"+airfoil+"_normal.msh"));
    ms.push_back(mesh("./"+airfoil+"_fine.msh"));

    settings.bcs["farfield"].vars_far.angle = db.alpha[0] * 0.01745;
    multigrid<T> multi(ms, settings, gui, residuals, iters, profile);
    multi.solvers[0].init();

    for (auto& alpha: db.alpha) {
        gui.msg.push("[RANS] Solving for alpha = " + std::to_string(alpha) + " deg.");
        settings.bcs["farfield"].vars_far.angle = alpha * 0.01745; // deg to rad
        for (auto& s: multi.solvers) {
            s.set_bcs(settings.bcs);
        };
        rans::solver& s = multi.run(false);
        wallProfile wp = get_wall_profile(s, "wall");
        db.cl.push_back(wp.cl);
        db.cd.push_back(wp.cd);
        db.cmy.push_back(wp.cm);
        save(airfoil + "_" + std::to_string(alpha) + ".vtu", s);
        if (gui.signal.stop) break;
    }
}

void Rans::solve() {
    if (settings.solver_type() == "implicit") {
        run<implicitSolver>();
    } else if (settings.solver_type() == "explicit") {
        run<explicitSolver>();
    }
};

void Rans::solve_airfoil(const std::string& airfoil, database::airfoil& db) {
    iters = 0;
    if (settings.solver_type() == "implicit") {
        run_airfoil<implicitSolver>(airfoil, db);
    } else if (settings.solver_type() == "explicit") {
        run_airfoil<explicitSolver>(airfoil, db);
    }
};

}
