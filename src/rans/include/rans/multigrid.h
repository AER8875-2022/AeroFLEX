/*

        ____     ____  ___    _   _______
       /  _/    / __ \/   |  / | / / ___/
       / /_____/ /_/ / /| | /  |/ /\__ \
     _/ /_____/ _, _/ ___ |/ /|  /___/ /
    /___/    /_/ |_/_/  |_/_/ |_//____/

    Implicit RANS Solver
    Description: Multigrid class
    By: Alexis Angers

*/
#pragma once
#include <rans/solver.h>
#include <atomic>
#include <algorithm>
#include <thread>

#include <common_aeroflex.hpp>

namespace rans {

// Implementation of Full Multigrid (FMG) initialization

template<class solverType>
class multigrid {
    SignalHandler &signal_gui;
    std::vector<double> &residuals;
    std::atomic<int> &iters;

    std::vector<solverType> solvers;

    std::vector<Eigen::SparseMatrix<double>> mappers;

    Eigen::SparseMatrix<double> gen_mapper(const uint i);

    int run_solver(solverType& s, double start_cfl=40, const uint max_iters=300, double tolerance=1e-6, const double relaxation=1);

public:

    multigrid(std::vector<mesh> ms, std::map<std::string, boundary_condition> bcs, gas g, bool second_order, SignalHandler &signal_gui, std::vector<double> &residuals, std::atomic<int> &iters);

    solverType& run(const bool reinit=true, const double relaxation=1);

};

template<class solverType>
multigrid<solverType>::multigrid(std::vector<mesh> ms, std::map<std::string, boundary_condition> bcs, gas g, bool second_order, SignalHandler &signal_gui, std::vector<double> &residuals, std::atomic<int> &iters) : signal_gui(signal_gui), residuals(residuals), iters(iters) {

    solvers.reserve(ms.size());
    for (auto& mi : ms) {
        solvers.push_back(solverType(mi, g, "laminar"));
        solvers[solvers.size()-1].set_bcs(bcs);
        solvers[solvers.size()-1].set_second_order(second_order);
    }

    #ifdef RANS_DEBUG
        #ifdef _OPENMP
            std::cout << "\nMultigrid : Eigen using " << Eigen::nbThreads() << " threads" << std::endl;
        #endif
    #endif

    std::cout << "\nMultigrid : Precompute " << ms.size()-1 << " matrices" << std::endl;

    if (ms.size() > 1) {
        mappers.resize(ms.size()-1);

        for (uint i=0; i<ms.size()-1; ++i) {
            std::cout << " - " << i+1 << " ";
            mappers[i] = gen_mapper(i);
            std::cout << " done\n" << std::flush;
        }
    }
}


template<class solverType>
Eigen::SparseMatrix<double> multigrid<solverType>::gen_mapper(const uint i) {

    auto& coarse = solvers[i].get_mesh();
    auto& fine = solvers[i+1].get_mesh();

    const uint m = fine.cellsAreas.size();
    const uint n = coarse.cellsAreas.size();

    std::vector<double> scales(m);

    const uint print_interval = (m-1)/9;

    for (uint i=0; i<m; ++i) {
        // Loop over rows, meaning fine mesh elements

        const double& xi = fine.cellsCentersX[i];
        const double& yi = fine.cellsCentersY[i];

        scales[i] = 0;

        for (uint j=0; j<n; ++j) {
            // Loop over columns, meaning coarse mesh elements
            const double& xj = coarse.cellsCentersX[j];
            const double& yj = coarse.cellsCentersY[j];

            // Compute distance between these two elements
            const double d2 = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj);

            const double& r2 = coarse.cellsAreas[j];

            if (d2 < 2*r2) {
                // Include this cell
                const double si = 1/std::max(0.1*sqrt(r2), sqrt(d2));
                scales[i] += si;
            }
        }

        if (i % print_interval == 0) std::cout << "." << std::flush;
    }

    Eigen::SparseMatrix<double> mapper(4*m, 4*n);

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(6*4*n);


    for (uint i=0; i<m; ++i) {
        // Loop over rows, meaning fine mesh elements

        const double& xi = fine.cellsCentersX[i];
        const double& yi = fine.cellsCentersY[i];

        for (uint j=0; j<n; ++j) {
            // Loop over columns, meaning coarse mesh elements
            const double& xj = coarse.cellsCentersX[j];
            const double& yj = coarse.cellsCentersY[j];

            // Compute distance between these two elements
            const double d2 = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj);

            const double& r2 = coarse.cellsAreas[j];

            if (d2 < 2*r2) {
                // Include this cell
                const double si = 1/std::max(0.1*sqrt(r2), sqrt(d2));

                for (int k=0; k<4; ++k) {
                    tripletList.push_back(Eigen::Triplet<double>(4*i+k, 4*j+k, si/scales[i]));
                }
            }
        }

        if (i % print_interval == 0) std::cout << "." << std::flush;
    }

    mapper.setFromTriplets(tripletList.begin(), tripletList.end());

    return mapper;
}



template<>
int multigrid<explicitSolver>::run_solver(
    explicitSolver& s,
    const double start_cfl,
    const uint max_iters,
    double tolerance,
    const double relaxation
) {

    double cfl = start_cfl;

    // Iterations loop

    int i = 0;
    double err = 0;
    double err_last = 0;

    double err_0 = 0;
    do {
        while (signal_gui.pause) std::this_thread::sleep_for(std::chrono::milliseconds(100));

        if (i % s.get_print_interval() == 0) std::cout << "Iteration " << i << " " << std::flush;
        s.set_cfl(cfl);

        s.fill();
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        int ok = s.compute();
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        if (ok == 0) {
            err_last = err;
            err = s.solve(relaxation);
        } else
            err = -1;
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        if (i == 0) err_0 = err;

        if (err > 0) {
            err /= err_0;
        }

        if (i % s.get_print_interval() == 0) std::cout << " Residual = " << err << std::endl;


        if (err < 0) return 1;
        // Store every 10 residual
        if (i % 10 == 0) {
            residuals[iters] = err;
            iters++;
        }

        i++;
    } while ((err > tolerance) && (i < max_iters) && !signal_gui.stop);

    return 0;
}

template<>
int multigrid<implicitSolver>::run_solver(
    implicitSolver& s,
    const double start_cfl,
    const uint max_iters,
    double tolerance,
    const double relaxation
) {

    double cfl = start_cfl;

    // Iterations loop

    int i = 0;
    double err = 0;
    double err_last = 0;

    double err_0 = s.get_uniform_residual();
    do {
        while (signal_gui.pause) std::this_thread::sleep_for(std::chrono::milliseconds(100));

        if (i % s.get_print_interval() == 0) std::cout << "Iteration " << i << " " << std::flush;
        s.set_cfl(cfl);

        s.fill();
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        int ok = s.compute();
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        if (ok == 0) {
            err_last = err;
            err = s.solve(relaxation, err_0 * tolerance);
        } else
            err = -1;
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        err /= err_0;

        if (err >= 0) {

            if (i > 0) {
                cfl = cfl * std::max(0.1, std::min(2., (err_last)/(err)));
            }
        }

        if (i % s.get_print_interval() == 0) std::cout << " Residual = " << err << std::endl;
        i++;

        if (err < 0) return 1;
        residuals[iters] = err;
        iters++;
        // do-while c'est cursed
        // Do smt if max iters is reached
    } while ((err > tolerance) && (i < max_iters) && !signal_gui.stop);

    return 0;
}

template<>
explicitSolver& multigrid<explicitSolver>::run(const bool reinit, const double relaxation) {
    const uint max_iters = 30000;
    residuals.resize(max_iters / 10);
    if (reinit) solvers[0].init();

    for (uint i=0; i<solvers.size(); ++i) {
        iters = 0;
        std::fill(residuals.begin(), residuals.end(), 0.0);

        // Multigrid loop

        std::cout << "\nMultigrid : Stage " << i+1 << "/" << solvers.size() << "\n" << std::endl;

        if (i > 0) {
            // Map last solution to current grid
            solvers[i-1].bcs_from_internal();
            solvers[i].get_q() = mappers[i-1] * solvers[i-1].get_q();
            solvers[i].refill_bcs();
        }

        double start_cfl = 1.25;
        double tolerance = (i == 0)|(i == (solvers.size()-1)) ? 1e-4 : 1e-2;

        int state = run_solver(solvers[i], start_cfl, max_iters, tolerance, relaxation);

        if (state) return solvers[i];
    }

    std::cout << std::endl;

    return solvers[solvers.size()-1];
}

template<>
implicitSolver& multigrid<implicitSolver>::run(const bool reinit, const double relaxation) {
    const uint max_iters = 300;
    residuals.resize(max_iters);

    if (reinit) solvers[0].init();
    solvers[0].refill_bcs();

    for (uint i=0; i<solvers.size(); ++i) {
        iters = 0;
        std::fill(residuals.begin(), residuals.end(), 0.0);

        // Multigrid loop

        std::cout << "\nMultigrid : Stage " << i+1 << "/" << solvers.size() << "\n" << std::endl;

        if (i > 0) {
            // Map last solution to current grid
            solvers[i-1].bcs_from_internal();
            solvers[i].get_q() = mappers[i-1] * solvers[i-1].get_q();
            solvers[i].refill_bcs();
        }

        double start_cfl = i == 0 ? 10 : 10;
        double tolerance = (i == 0)|(i == (solvers.size()-1)) ? 1e-4 : 1e-4;

        int state = run_solver(solvers[i], start_cfl, max_iters, tolerance, relaxation);

        if (state || signal_gui.stop) return solvers[i];
    }

    std::cout << std::endl;

    return solvers[solvers.size()-1];
}

}