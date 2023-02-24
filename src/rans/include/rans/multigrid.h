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
    GUIHandler &gui;
    Settings &settings;
    std::vector<double> &residuals;
    std::atomic<int> &iters;

    std::vector<solverType> solvers;

    std::vector<Eigen::SparseMatrix<double>> mappers;

    Eigen::SparseMatrix<double> gen_mapper(const uint i);

    int run_solver(solverType& s);

public:

    multigrid(std::vector<mesh> ms, Settings &settings, GUIHandler &gui, std::vector<double> &residuals, std::atomic<int> &iters);

    solverType& run(const bool reinit=true);

};

template<class solverType>
multigrid<solverType>::multigrid(std::vector<mesh> ms, Settings &settings, GUIHandler &gui, std::vector<double> &residuals, std::atomic<int> &iters)
: gui(gui), settings(settings), residuals(residuals), iters(iters)
{

    solvers.reserve(ms.size());
    for (auto& mi : ms) {
        solvers.push_back(solverType(mi, settings.g, settings.viscosity_model));
        solvers[solvers.size()-1].set_bcs(settings.bcs);
        solvers[solvers.size()-1].set_second_order(settings.second_order);
        solvers[solvers.size()-1].set_gradient_scheme(settings.gradient_scheme);
        solvers[solvers.size()-1].set_limiter_k(settings.limiter_k);
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
    explicitSolver& s
) {

    double cfl = settings.start_cfl;

    // Iterations loop

    int i = 0;
    double err = 0;
    double err_last = 0;

    double err_0 = s.get_uniform_residual();
    do {
        while (gui.signal.pause) std::this_thread::sleep_for(std::chrono::milliseconds(100));

        if (i % s.get_print_interval() == 0) std::cout << "Iteration " << i << " " << std::flush;
        s.set_cfl(cfl);

        s.fill();
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        int ok = s.compute();
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        if (ok == 0) {
            err_last = err;
            err = s.solve(settings.relaxation);
        } else
            err = -1;
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

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
    } while ((err > settings.tolerance) && (i < settings.max_iterations) && !gui.signal.stop);

    return 0;
}

template<>
int multigrid<implicitSolver>::run_solver(
    implicitSolver& s
) {

    double cfl = settings.start_cfl;

    // Iterations loop

    int i = 0;
    double err = 0;
    double err_last = 0;

    double err_0 = s.get_uniform_residual();
    do {
        while (gui.signal.pause) std::this_thread::sleep_for(std::chrono::milliseconds(100));

        if (i % s.get_print_interval() == 0) std::cout << "Iteration " << i << " " << std::flush;
        s.set_cfl(cfl);

        s.fill();
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        int ok = s.compute();
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        if (ok == 0) {
            err_last = err;
            err = s.solve(settings.relaxation, err_0 * settings.tolerance, settings.rhs_iterations);
        } else
            err = -1;
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        err /= err_0;


        //cfl = std::max(settings.start_cfl, std::min(settings.max_cfl, cfl * err_last / err));
        cfl = std::min(settings.start_cfl + (i+1)*settings.slope_cfl, settings.max_cfl);

        if (i % s.get_print_interval() == 0) std::cout << " Residual = " << err << std::endl;
        i++;

        if (err < 0) return 1;
        residuals[iters] = err;
        iters++;
        // do-while c'est cursed
        // Do smt if max iters is reached
    } while ((err > settings.tolerance) && (i < settings.max_iterations) && !gui.signal.stop);

    return 0;
}

template<>
explicitSolver& multigrid<explicitSolver>::run(const bool reinit) {
    const uint max_iters = settings.max_iterations;
    residuals.resize(solvers.size() * (max_iters / 10));

    if (reinit) solvers[0].init();
    solvers[0].refill_bcs();

    for (uint i=0; i<solvers.size(); ++i) {

        // Multigrid loop

        std::cout << "\nMultigrid : Stage " << i+1 << "/" << solvers.size() << "\n" << std::endl;

        if (i > 0) {
            // Map last solution to current grid
            solvers[i-1].bcs_from_internal();
            solvers[i].get_q() = mappers[i-1] * solvers[i-1].get_q();
            solvers[i].refill_bcs();
        }

        int state = run_solver(solvers[i]);

        if (state) return solvers[i];
    }

    std::cout << std::endl;

    return solvers[solvers.size()-1];
}

template<>
implicitSolver& multigrid<implicitSolver>::run(const bool reinit) {
    const uint max_iters = settings.max_iterations;
    residuals.resize(solvers.size() * max_iters);

    if (reinit) solvers[0].init();
    solvers[0].refill_bcs();

    for (uint i=0; i<solvers.size(); ++i) {
        // Multigrid loop

        std::cout << "\nMultigrid : Stage " << i+1 << "/" << solvers.size() << "\n" << std::endl;

        if (i > 0) {
            // Map last solution to current grid
            solvers[i-1].bcs_from_internal();
            solvers[i].get_q() = mappers[i-1] * solvers[i-1].get_q();
            solvers[i].refill_bcs();
        }

        int state = run_solver(solvers[i]);

        if (state || gui.signal.stop) return solvers[i];
    }

    std::cout << std::endl;

    return solvers[solvers.size()-1];
}

}