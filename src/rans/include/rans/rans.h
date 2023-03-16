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
    std::vector<double> residuals;
    //wallProfile wp;
    std::atomic<int> iters = 0;
    Settings settings;

    std::vector<mesh> ms;
    bool mesh_loaded = false;

    GUIHandler &gui;

    Rans(GUIHandler &gui) : gui(gui) {};

    void input();
    void solve();
    void solve_airfoil(const std::string& airfoil, database::airfoil& db);

    template<class T> void run();
    template<class T> void run_airfoil(const std::string& airfoil, database::airfoil& db);
};

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
    multigrid<T> multi(ms, settings, gui, residuals, iters);
    gui.event.rans_preprocess = true;

    rans::solver& s = multi.run(true);
    gui.event.rans_solve = true;

    save(settings.outfilename, s);
    gui.event.rans_postprocess = true;
    std::cout << "Saved results to file " << settings.outfilename << "\n" << std::endl;
}

template<class T>
void Rans::run_airfoil(const std::string& airfoil, database::airfoil& db) {
    gui.msg.push("[RANS] Solving airfoil: " + airfoil);
    ms.clear();
    // TODO: check with geom for the naming convention
    // For the moment we will only load 1 mesh
    ms.push_back(mesh("../../../../examples/rans/" + airfoil + ".msh"));
    multigrid<T> multi(ms, settings, gui, residuals, iters);
    multi.solvers[0].init();

    for (auto& alpha: db.alpha) {
        gui.msg.push("[RANS] Solving for alpha = " + std::to_string(alpha) + "°");
        settings.bcs["Farfield"].vars_far.angle = alpha;
        for (auto& s: multi.solvers) {
            s.set_bcs(settings.bcs);
        };
        rans::solver& s = multi.run(false);
        wallProfile wp = get_wall_profile(s, "Airfoil");
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