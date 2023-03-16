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

namespace rans {

class Rans {
    public:
    // This is a mess because of no separation of IO and computation...
    std::vector<double> residuals;
    std::atomic<int> iters = 0;
    Settings settings;

    std::vector<mesh> ms;
    bool mesh_loaded = false;

    GUIHandler &gui;

    Rans(GUIHandler &gui) : gui(gui) {};

    void input();
    void solve();
    template<class T> void run();
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

void Rans::solve() {
    if (settings.solver_type() == "implicit") {
        run<implicitSolver>();
    } else if (settings.solver_type() == "explicit") {
        run<explicitSolver>();
    }
};

}