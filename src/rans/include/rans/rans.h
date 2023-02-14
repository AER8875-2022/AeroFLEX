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

namespace rans {

class Rans {
    public:
    // This is a mess because of no separation of IO and computation...
    std::vector<double> residuals;
    std::atomic<int> iters = 0;
    Settings data;

    std::vector<mesh> ms;
    bool mesh_loaded = false;

    SignalHandler &signal_gui;

    Rans(SignalHandler &signal_gui) : signal_gui(signal_gui) {};

    void input();
    void solve();
    template<class T> void run();
};

void Rans::input() {
    if (!mesh_loaded) {
        for (const auto& mesh_name : data.meshes) {
            ms.push_back(mesh(mesh_name));
        }
        mesh_loaded = true;
    }
};

template<class T>
void Rans::run() {
    multigrid<T> multi(ms, data.g, data.second_order, signal_gui, residuals, iters);

    rans::solver& s = multi.run(data.vars_far, true, data.relaxation);

    save(data.outfilename, s);
    std::cout << "Saved results to file " << data.outfilename << "\n" << std::endl;
}

void Rans::solve() {
    if (data.solver_type == "implicit") {
        run<implicitSolver<inviscid>>();
    } else if (data.solver_type == "explicit") {
        run<explicitSolver<inviscid>>();
    }
};

}