/*

        ____     ____  ___    _   _______
       /  _/    / __ \/   |  / | / / ___/
       / /_____/ /_/ / /| | /  |/ /\__ \
     _/ /_____/ _, _/ ___ |/ /|  /___/ /
    /___/    /_/ |_/_/  |_/_/ |_//____/

    Implicit RANS Solver
    Description: Rans solver class
    By: Alexis Angers

*/
#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <memory>

#include <rans/mesh.h>
#include <rans/core.h>
#include <rans/physics.h>

#ifdef RANS_MATRIX_FREE
    #include <rans/ransMatrix.h>
#endif




namespace rans {




class solver {

#ifdef RANS_DEBUG_PUBLIC
public:
#else
protected:
#endif

    // Gas & mesh
    gas g;
    mesh m;

    // q = rho, rhou, rhov, rhoe
    // nu is the turbulence model variable
    Eigen::VectorXd q;
    Eigen::VectorXd nu;

    // w -> updates computed by the linear solver step
    Eigen::VectorXd qW;
    Eigen::VectorXd nuW;

    // The time step
    Eigen::VectorXd dt;

    // Gradients and limiters (q only)
    Eigen::VectorXd gx;
    Eigen::VectorXd gy;
    Eigen::VectorXd limiters;
    Eigen::VectorXd q_min;
    Eigen::VectorXd q_max;

    bool second_order;

    std::vector<std::unique_ptr<flux>> edges_flux_functions;


    void set_mesh_and_gas(const mesh& m_in, const gas& g_in);

    void set_walls_from_internal(Eigen::VectorXd&);

    void calc_dt();
    void calc_gradients(const Eigen::VectorXd&);
    void calc_limiters(const Eigen::VectorXd&);

    double cfl = 1;

    uint print_interval;

public:

    solver() {}

    solver(const mesh& m_in, const gas& g_in);

    solver(const solver& s) : solver(s.get_cmesh(), s.get_gas()) {}

    void set_cfl(const double& cfl_in);

    void init(farfield_conditions);


    double get_uniform_residual(farfield_conditions);

    void farfield(farfield_conditions);

    Eigen::VectorXd& get_q();
    gas& get_gas() {return g;}
    const gas& get_gas() const {return g;}

    void zero_bcs();
    void bcs_from_internal();

    solution get_solution();

    void set_second_order(const bool x=true);

    mesh& get_mesh() {return m;}
    const mesh& get_cmesh() const {return m;}
    uint get_print_interval() {return print_interval;}


    virtual void fill() {}
    virtual int compute() {return -1;}
    virtual double solve(const double relaxation=1, const double tol=0) {return -1;}

};


void solver::set_mesh_and_gas(const mesh& m_in, const gas& g_in) {
    second_order = true;

    m = m_in;
    g = g_in;

    int N = m.cellsAreas.size();

    // Fill edges flux functions
    std::vector<int> edges_flux_to_add(m.edgesCells.cols());
    std::fill(edges_flux_to_add.begin(), edges_flux_to_add.end(), 0);

    for (int b=0; b<m.boundaryEdges.size(); ++b) {
        const uint& e = m.boundaryEdges[b];
        if (m.boundaryEdgesPhysicals[b] == "farfield") {
            edges_flux_to_add[e] = 1;
        } else if (m.boundaryEdgesPhysicals[b] == "wall") {
            edges_flux_to_add[e] = 2;
        }
    }

    edges_flux_functions.reserve(m.edgesCells.cols());
    for (int e=0; e<m.edgesCells.cols(); ++e) {
        double& nx = m.edgesNormalsX[e];
        double& ny = m.edgesNormalsY[e];
        if (edges_flux_to_add[e] == 0) {
            edges_flux_functions.emplace_back(std::make_unique<internal_flux<inviscid>>(g, nx, ny));
        } else if (edges_flux_to_add[e] == 1) {
            edges_flux_functions.emplace_back(std::make_unique<farfield_flux<inviscid>>(g, nx, ny));
        } else if (edges_flux_to_add[e] == 2) {
            edges_flux_functions.emplace_back(std::make_unique<slip_wall_flux<inviscid>>(g, nx, ny));
        }
    }

    q.resize(4*N);
    qW.resize(4*N);

    nu.resize(N);
    nuW.resize(N);

    dt.resize(N);

    gx.resize(4*N);
    gy.resize(4*N);
    limiters.resize(4*N);
    q_min.resize(4*N);
    q_max.resize(4*N);
}

solver::solver(const mesh& m_in, const gas& g_in) {
    set_mesh_and_gas(m_in, g_in);
}


void solver::set_cfl(const double& cfl_in) {
    cfl = cfl_in;
}

void solver::set_second_order(const bool x) {
    second_order = x;
}


Eigen::VectorXd& solver::get_q() {
    return q;
}

void solver::zero_bcs() {
    for (uint i=m.nRealCells; i<m.cellsAreas.size(); ++i) {
        for (int k=0; k<4; ++k) {
            q(4*i + k) = 0;
        }
    }
}

void solver::bcs_from_internal() {
    for (uint b=0; b<m.boundaryEdges.size(); ++b) {
        const uint e = m.boundaryEdges[b];

        int cell0 = m.edgesCells(e, 0);
        int cell1 = m.edgesCells(e, 1);

        for (int k=0; k<4; ++k) {
            q(4*cell1 + k) = q(4*cell0 + k);
        }
    }
}

void solver::set_walls_from_internal(Eigen::VectorXd& q_) {
    for (uint b=0; b<m.boundaryEdges.size(); ++b) {
        if (m.boundaryEdgesPhysicals[b] == "wall") {
            const uint e = m.boundaryEdges[b];

            int cell0 = m.edgesCells(e, 0);
            int cell1 = m.edgesCells(e, 1);

            for (int k=0; k<4; ++k) {
                q_(4*cell1 + k) = q_(4*cell0 + k);
            }
        }
    }
}


void solver::calc_dt() {
    #pragma omp parallel for
    for (int i=0; i<dt.size(); ++i) dt(i) = 0;

    // Reset rho matrix to zero
    for (int e=0; e<m.edgesCells.cols(); ++e) {

        // Get this edges cells index
        int cell0 = m.edgesCells(e, 0);
        int cell1 = m.edgesCells(e, 1);

        int k0 = 4*cell0;
        int k1 = 4*cell1;

        Eigen::VectorXd n(2);
        n(0) = m.edgesNormalsX[e];
        n(1) = m.edgesNormalsY[e];

        double eig;

        if (edges_flux_functions[e]->two_sided) {
            const double V_L = (q(k0+1)*n(0) + q(k0+2)*n(1))/q(k0);
            const double V_R = (q(k1+1)*n(0) + q(k1+2)*n(1))/q(k1);

            const double p_L = (g.gamma-1)*(q(k0+3) - 0.5/q(k0)*(q(k0+1)*q(k0+1)+q(k0+2)*q(k0+2)));
            const double p_R = (g.gamma-1)*(q(k1+3) - 0.5/q(k1)*(q(k1+1)*q(k1+1)+q(k1+2)*q(k1+2)));

            const double eig_L = sqrt(p_L*g.gamma / q(k0)) + abs(V_L);
            const double eig_R = sqrt(p_R*g.gamma / q(k1)) + abs(V_R);
            eig = std::max(eig_L, eig_R);
        } else {
            // Use the internal field
            const double V_L = (q(k0+1)*n(0) + q(k0+2)*n(1))/q(k0);

            const double p_L = (g.gamma-1)*(q(k0+3) - 0.5/q(k0)*(q(k0+1)*q(k0+1)+q(k0+2)*q(k0+2)));

            eig = sqrt(p_L*g.gamma / q(k0)) + abs(V_L);
        }

        dt(cell0) += eig * m.edgesLengths[e];
        if (edges_flux_functions[e]->two_sided)
            dt(cell1) += eig * m.edgesLengths[e];
    }

    #pragma omp parallel for
    for (int i=0; i<dt.size(); ++i) {
        dt(i) = cfl * m.cellsAreas[i] / dt(i);
    }
}


void solver::calc_gradients(const Eigen::VectorXd& q_) {

    #pragma omp parallel for
    for (int i=0; i<gx.size(); ++i) {gx(i) = 0; gy(i) = 0;}

    for (uint e=0; e<m.edgesNodes.cols(); ++e) {
        const uint& i = m.edgesCells(e, 0);
        const uint& j = m.edgesCells(e, 1);

        const double dxif = m.edgesCentersX[e] - m.cellsCentersX[i];
        const double dyif = m.edgesCentersY[e] - m.cellsCentersY[i];
        const double dif = sqrt(dxif*dxif + dyif*dyif);

        const double dxij = m.cellsCentersX[i] - m.cellsCentersX[j];
        const double dyij = m.cellsCentersY[i] - m.cellsCentersY[j];
        const double dij = sqrt(dxij*dxij + dyij*dyij);

        const double geom_factor = dif / dij;

        if (i != j) {
            for (uint k=0; k<4; ++k) {
                const double fk = (q_(4*i+k)*(1.0 - geom_factor) + q_(4*j+k) * geom_factor) * m.edgesLengths[e];

                gx(4*i+k) += fk * m.edgesNormalsX[e];
                gy(4*i+k) += fk * m.edgesNormalsY[e];

                gx(4*j+k) -= fk * m.edgesNormalsX[e];
                gy(4*j+k) -= fk * m.edgesNormalsY[e];
            }
        } else {
            for (uint k=0; k<4; ++k) {
                const double fk = q(4*i+k) * m.edgesLengths[e];

                gx(4*i+k) += fk * m.edgesNormalsX[e];
                gy(4*i+k) += fk * m.edgesNormalsY[e];
            }
        }
    }

    #pragma omp parallel for
    for (int i=0; i<m.nRealCells; ++i) {
        gx.segment(4*i, 4) /= m.cellsAreas[i];
        gy.segment(4*i, 4) /= m.cellsAreas[i];
    }
    #pragma omp parallel for
    for (int i=m.nRealCells; i<m.cellsAreas.size(); ++i) {
        gx.segment(4*i, 4) *= 0;
        gy.segment(4*i, 4) *= 0;
    }
}



void solver::calc_limiters(const Eigen::VectorXd& q_) {
    #pragma omp parallel for
    for (int i=0; i<limiters.size(); ++i) {limiters(i) = 1; q_min(i) = q_(i); q_max(i) = q_(i);}

    // Compute qmin and qmax
    for (uint e=0; e<m.edgesNodes.cols(); ++e) {
        const uint& i = m.edgesCells(e, 0);
        const uint& j = m.edgesCells(e, 1);

        for (uint k=0; k<4; ++k) {
            q_min(4*i+k) = std::min(q_min(4*i+k), q_(4*j+k));
            q_min(4*j+k) = std::min(q_min(4*j+k), q_(4*i+k));
            q_max(4*i+k) = std::max(q_max(4*i+k), q_(4*j+k));
            q_max(4*j+k) = std::max(q_max(4*j+k), q_(4*i+k));
        }
    }
    // Compute limiters
    std::array<uint, 2> ids;
    for (uint e=0; e<m.edgesNodes.cols(); ++e) {
        ids[0] = m.edgesCells(e, 0);
        ids[1] = m.edgesCells(e, 1);

        for (auto& id : ids) {
            if (id < m.nRealCells) {
                const double dx = m.edgesCentersX[e] - m.cellsCentersX[id];
                const double dy = m.edgesCentersY[e] - m.cellsCentersY[id];
                const double sqrt_area = sqrt(m.cellsAreas[id]);

                for (uint k=0; k<4; ++k) {
                    double dqg = gx(4*id+k)*dx + gy(4*id+k)*dy;

                    double delta_max = q_max(4*id+k) - q_(4*id+k);
                    double delta_min = q_min(4*id+k) - q_(4*id+k);

                    const double Ka = 5. * sqrt_area;
                    const double K3a = Ka * Ka * Ka;
                    const double dMaxMin2 = (delta_max - delta_min)*(delta_max - delta_min);

                    double lim = 1.0;

                    #ifdef RANS_MICHALAK_LIMITER
                        double sig;
                        if (dMaxMin2 <= K3a) {
                            sig = 1.;
                        } else if (dMaxMin2 <= 2*K3a) {
                            const double y = (dMaxMin2/K3a - 1.0);
                            sig = 2.0*y*y*y - 3.0*y*y + 1.0;
                        } else {
                            sig = 0.;
                        }
                        if (sig < 1.0) {
                            if (dqg > 1e-14) {
                                lim = michalak_limiter(delta_max/dqg);
                            } else if (dqg < -1e-14) {
                                lim = michalak_limiter(delta_min/dqg);
                            } else {
                                lim = 1.0;
                            }
                        }
                        lim = sig + (1.0 - sig)*lim;
                    #else
                        if (dqg > 1e-16) {
                            lim = 1/dqg * ((delta_max*delta_max + K3a)*dqg + 2*dqg*dqg*delta_max) / (delta_max*delta_max + 2*dqg*dqg + delta_max*dqg + K3a);
                        } else if (dqg < -1e-16) {
                            lim = 1/dqg * ((delta_min*delta_min + K3a)*dqg + 2*dqg*dqg*delta_min) / (delta_min*delta_min + 2*dqg*dqg + delta_min*dqg + K3a);
                        } else {
                            lim = 1.0;
                        }
                    #endif

                    limiters(4*id+k) = std::min(limiters(4*id+k), lim);
                }
            }
        }
    }

}



void solver::init(farfield_conditions var_init) {
    int N = m.nRealCells;

    auto [rho, rhou, rhov, rhoe] = var_init.get_conservative(g);

    #pragma omp parallel for
    for (int i=0; i<N; ++i) {
        q(4*i) = rho;
        q(4*i+1) = rhou;
        q(4*i+2) = rhov;
        q(4*i+3) = rhoe;
    }
}



void solver::farfield(farfield_conditions var_far) {

    auto [rho, rhou, rhov, rhoe] = var_far.get_conservative(g);

    #pragma omp parallel for
    for (int b=0; b<m.boundaryEdges.size(); ++b) {
        if (m.boundaryEdgesPhysicals[b] == "farfield") {
            int cb = m.edgesCells(m.boundaryEdges[b], 1);
            int k = cb*4;

            q(k) = rho;
            q(k+1) = rhou;
            q(k+2) = rhov;
            q(k+3) = rhoe;
        }
    }
}



double solver::get_uniform_residual(farfield_conditions var_far) {
    // Get residual assuming a uniform flow field corresponding to vars_far

    // Reset qW to zero
    #pragma omp parallel for
    for (int i=m.nRealCells; i<m.cellsAreas.size(); ++i) {
        // Dirichlet for boundary cell
        for (int j=0; j<4; ++j) {
            qW.coeffRef(4*i+j) = 0;
        }
    }

    auto [rho, rhou, rhov, rhoe] = var_far.get_conservative(g);
    Eigen::VectorXd q_far(4);
    q_far(0) = rho;
    q_far(1) = rhou;
    q_far(2) = rhov;
    q_far(3) = rhoe;

    // Compute rho vector of fluxes
    for (int e=0; e<m.edgesCells.cols(); ++e) {

        // Get this edges cells index
        int cell0 = m.edgesCells(e, 0);
        int cell1 = m.edgesCells(e, 1);

        int k0 = 4*cell0;
        int k1 = 4*cell1;


        Eigen::VectorXd n(2);
        n(0) = m.edgesNormalsX[e];
        n(1) = m.edgesNormalsY[e];

        // Set flux
        {
            const Eigen::VectorXd this_flux = (*edges_flux_functions[e])(q_far, q_far) * m.edgesLengths[e];
            for (int i=0; i<4; ++i) {
                qW(k0+i) -= this_flux(i);
                if (edges_flux_functions[e]->two_sided)
                    qW(k1+i) += this_flux(i);
            }
        }
    }

    // Return norm of rho vector
    return qW.norm();
}


solution solver::get_solution() {
    int N = m.cellsAreas.size();

    solution s;
    s.gamma = g.gamma;
    s.R = g.R;

    s.rho.resize(N);
    s.rhou.resize(N);
    s.rhov.resize(N);
    s.rhoe.resize(N);

    #pragma omp parallel for
    for (int i=0; i<N; ++i) {
        s.rho(i) = q(4*i);
        s.rhou(i) = q(4*i+1);
        s.rhov(i) = q(4*i+2);
        s.rhoe(i) = q(4*i+3);
    }

    return s;
}
























class explicitSolver : public solver {

    std::vector<double> alpha = {0.25, 0.5, 1.};
    Eigen::VectorXd qk;

    void calc_residual(const Eigen::VectorXd&);

public:

    explicitSolver(const mesh& m_in, const gas& g_in) {
        set_mesh_and_gas(m_in, g_in);
        print_interval = 100;
        qk.resize(q.size());
    }

    explicitSolver(const explicitSolver& s) : explicitSolver(s.get_cmesh(), s.get_gas()) {}

    void fill() {}
    int compute() {return 0;}
    double solve(const double relaxation=1, const double tol=0);

};


void explicitSolver::calc_residual(const Eigen::VectorXd& q_) {
    #pragma omp parallel for
    for (int i=0; i<qW.size(); ++i) {
        qW(i) = 0;
    }

    for (int e=0; e<m.edgesCells.cols(); ++e) {

        // Get this edges cells index
        int cell0 = m.edgesCells(e, 0);
        int cell1 = m.edgesCells(e, 1);

        int k0 = 4*cell0;
        int k1 = 4*cell1;

        Eigen::VectorXd q_L = q_.segment(k0, 4);
        Eigen::VectorXd q_R = q_.segment(k1, 4);

        Eigen::VectorXd n(2);
        n(0) = m.edgesNormalsX[e];
        n(1) = m.edgesNormalsY[e];

        // Use linear interpolation (second order)
        Eigen::VectorXd this_flux;
        if (second_order) {
            const double d0x = m.edgesCentersX[e] - m.cellsCentersX[cell0];
            const double d0y = m.edgesCentersY[e] - m.cellsCentersY[cell0];

            const double d1x = m.edgesCentersX[e] - m.cellsCentersX[cell1];
            const double d1y = m.edgesCentersY[e] - m.cellsCentersY[cell1];

            const Eigen::VectorXd q_L_o2 = q_L + ( gx.segment(k0,4)*d0x + gy.segment(k0,4)*d0y ).cwiseProduct(limiters.segment(k0,4));
            const Eigen::VectorXd q_R_o2 = q_R + ( gx.segment(k1,4)*d1x + gy.segment(k1,4)*d1y ).cwiseProduct(limiters.segment(k1,4));

            this_flux = (*edges_flux_functions[e])(q_L_o2, q_R_o2) * m.edgesLengths[e];
        } else {
            this_flux = (*edges_flux_functions[e])(q_L, q_R) * m.edgesLengths[e];
        }
        for (int i=0; i<4; ++i) {
            qW(k0+i) -= this_flux(i);
            if (edges_flux_functions[e]->two_sided)
                qW(k1+i) += this_flux(i);
        }
    }

    #pragma omp parallel for
    for (int i=0; i<m.nRealCells; ++i) {
        for (int j=0; j<4; ++j)
            qW(4*i+j) /= m.cellsAreas[i];
    }
}


double explicitSolver::solve(const double relaxation, const double tol) {
    calc_dt();

    #pragma omp parallel for
    for (int i=0; i<qk.size(); ++i) qk(i) = q(i);

    for (const double& ai : alpha) {
        if (second_order) {
            set_walls_from_internal(qk);
            calc_gradients(qk);
            calc_limiters(qk);
        }
        calc_residual(qk);
        for (uint i=0; i<m.nRealCells; ++i) {
            for (uint j=0; j<4; ++j)
                qk(4*i+j) = q(4*i+j) + qW(4*i+j) * dt(i) * ai * relaxation;
        }
    }

    #pragma omp parallel for
    for (int i=0; i<qk.size(); ++i) q(i) = qk(i);

    return qW.norm();
}
























class implicitSolver : public solver {

    // The right hand side vectors
    Eigen::VectorXd RhoVector;
    Eigen::VectorXd TurbVector;


    void fillRhoLHS();
    void fillRhoRHS();


    #ifdef RANS_MATRIX_FREE

        RansMatrix RhoMatrix;
        Eigen::SparseMatrix<double, Eigen::RowMajor> TurbMatrix;

        Eigen::GMRES<RansMatrix, Eigen::IdentityPreconditioner> rho_solver;

    #else

        #ifdef RANS_SPARSELU
            Eigen::SparseMatrix<double> RhoMatrix;
            Eigen::SparseMatrix<double> TurbMatrix;
        #else
            Eigen::SparseMatrix<double, Eigen::RowMajor> RhoMatrix;
            Eigen::SparseMatrix<double, Eigen::RowMajor> TurbMatrix;
        #endif

        #ifdef RANS_SPARSELU
            Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> rho_solver;
        #else
            #ifdef RANS_BICGSTAB
                Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>> rho_solver;
            #else
                Eigen::GMRES<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>> rho_solver;
            #endif
        #endif

    #endif


public:

    implicitSolver(const mesh& m_in, const gas& g_in) {
        set_mesh_and_gas(m_in, g_in);
        print_interval = 1;

        // Linear solver parameters
        #ifndef RANS_BICGSTAB
            rho_solver.setTolerance(1e-2);
            rho_solver.preconditioner().setFillfactor(3);
            rho_solver.preconditioner().setDroptol(1e-4);

            rho_solver.setMaxIterations(500);
        #endif

        int N = m.cellsAreas.size();

        #ifdef RANS_MATRIX_FREE
            RhoMatrix.set_mesh(&m);
            RhoMatrix.set_gas(&g);
            RhoMatrix.set_q(&q);
            RhoMatrix.set_dt(&dt);
            RhoMatrix.set_flux_functions(&edges_flux_functions);
        #else
            RhoMatrix.resize(4*N, 4*N);
        #endif

        TurbMatrix.resize(N, N);

        RhoVector.resize(4*N);


        #ifndef RANS_MATRIX_FREE
            std::vector<Eigen::Triplet<double>> tripletList;
            tripletList.reserve(16*4*N);

            // Internal field
            for (int e=0; e<m.edgesCells.cols(); ++e) {

                // Get this edges cells index
                int cell0 = m.edgesCells(e, 0);
                int cell1 = m.edgesCells(e, 1);

                int k0 = 4*cell0;
                int k1 = 4*cell1;

                // Set Rho matrix elements to zero for this cell pair
                for (int i=0; i<4; ++i) {
                    for (int j=0; j<4; ++j) {
                        tripletList.push_back(Eigen::Triplet<double>(k0+i, k0+j, 0));
                        tripletList.push_back(Eigen::Triplet<double>(k0+i, k1+j, 0));
                        tripletList.push_back(Eigen::Triplet<double>(k1+i, k0+j, 0));
                        tripletList.push_back(Eigen::Triplet<double>(k1+i, k1+j, 0));
                    }
                }
            }
            RhoMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
        #endif
    }

    implicitSolver(const implicitSolver& s) : implicitSolver(s.get_cmesh(), s.get_gas()) {}


    void fill();
    int compute();
    double solve(const double relaxation=1, const double tol=0);

};



void implicitSolver::fill() {
    // Fill the jacobian matrix
    fillRhoLHS();
}



void implicitSolver::fillRhoLHS() {

    calc_dt();

    // Reset rho matrix to zero
    #ifndef RANS_MATRIX_FREE
        #pragma omp parallel for
        for (int e=0; e<m.edgesCells.cols(); ++e) {

            // Get this edges cells index
            int cell0 = m.edgesCells(e, 0);
            int cell1 = m.edgesCells(e, 1);

            int k0 = 4*cell0;
            int k1 = 4*cell1;

            // Reset Rho matrix elements to zero for this cell pair
            for (int i=0; i<4; ++i) {
                for (int j=0; j<4; ++j) {
                    RhoMatrix.coeffRef(k0+i, k0+j) = 0;
                    RhoMatrix.coeffRef(k0+i, k1+j) = 0;
                    RhoMatrix.coeffRef(k1+i, k0+j) = 0;
                    RhoMatrix.coeffRef(k1+i, k1+j) = 0;
                }
            }
        }

        // Fill internal matrix elements from time discretization
        #pragma omp parallel for
        for (int i=0; i<m.cellsAreas.size(); ++i) {
            int k = 4*i;
            for (int j=0; j<4; ++j) {
                RhoMatrix.coeffRef(k+j, k+j) = m.cellsAreas[i]/dt(i);
            }
        }
    #endif

    // Fill internal matrix elements from flux jacobian
    for (int e=0; e<m.edgesCells.cols(); ++e) {

        // Get this edges cells index
        int cell0 = m.edgesCells(e, 0);
        int cell1 = m.edgesCells(e, 1);

        int k0 = 4*cell0;
        int k1 = 4*cell1;

        Eigen::VectorXd q_L = q.segment(k0, 4);
        Eigen::VectorXd q_R = q.segment(k1, 4);

        Eigen::VectorXd n(2);
        n(0) = m.edgesNormalsX[e];
        n(1) = m.edgesNormalsY[e];


        #ifndef RANS_MATRIX_FREE
        {
            const Eigen::MatrixXd Jacobian = calc_convective_jacobian(*edges_flux_functions[e], q_L, q_R) * m.edgesLengths[e];

            // Set jacobian elements
            for (int i=0; i<4; ++i) {
                for (int j=0; j<4; ++j) {
                    RhoMatrix.coeffRef(k0+i, k0+j) += Jacobian(i, j);
                    RhoMatrix.coeffRef(k0+i, k1+j) += Jacobian(i, j+4);

                    if (edges_flux_functions[e]->two_sided) {
                        RhoMatrix.coeffRef(k1+i, k0+j) += Jacobian(i+4, j);
                        RhoMatrix.coeffRef(k1+i, k1+j) += Jacobian(i+4, j+4);
                    }
                }
            }
        }
        #endif
    }

    // Also fill bcs
    #pragma omp parallel for
    for (int ci=m.nRealCells; ci<m.cellsAreas.size(); ++ci) {
        int k = 4*ci;

        #ifndef RANS_MATRIX_FREE
        for (int i=0; i<4; ++i) {
            for (int j=0; j<4; ++j) {
                RhoMatrix.coeffRef(k+i, k+j) = 0;
            }
            RhoMatrix.coeffRef(k+i, k+i) = 1;
        }
        #endif
    }


}




void implicitSolver::fillRhoRHS() {

    calc_dt();

    if (second_order) {
        set_walls_from_internal(q);
        calc_gradients(q);
        calc_limiters(q);
    }

    // Reset rho vector
    #pragma omp parallel for
    for (int i=0; i<RhoVector.size(); ++i) {
        RhoVector(i) = 0;
    }


    // Fill internal matrix elements from flux jacobian
    for (int e=0; e<m.edgesCells.cols(); ++e) {

        // Get this edges cells index
        int cell0 = m.edgesCells(e, 0);
        int cell1 = m.edgesCells(e, 1);

        int k0 = 4*cell0;
        int k1 = 4*cell1;

        Eigen::VectorXd q_L = q.segment(k0, 4);
        Eigen::VectorXd q_R = q.segment(k1, 4);

        Eigen::VectorXd n(2);
        n(0) = m.edgesNormalsX[e];
        n(1) = m.edgesNormalsY[e];

        // Set flux
        {
            // Use linear interpolation (second order)
            Eigen::VectorXd this_flux;
            if (second_order) {
                const double d0x = m.edgesCentersX[e] - m.cellsCentersX[cell0];
                const double d0y = m.edgesCentersY[e] - m.cellsCentersY[cell0];

                const double d1x = m.edgesCentersX[e] - m.cellsCentersX[cell1];
                const double d1y = m.edgesCentersY[e] - m.cellsCentersY[cell1];

                const Eigen::VectorXd q_L_o2 = q_L + ( gx.segment(k0,4)*d0x + gy.segment(k0,4)*d0y ).cwiseProduct(limiters.segment(k0,4));
                const Eigen::VectorXd q_R_o2 = q_R + ( gx.segment(k1,4)*d1x + gy.segment(k1,4)*d1y ).cwiseProduct(limiters.segment(k1,4));

                this_flux = (*edges_flux_functions[e])(q_L_o2, q_R_o2) * m.edgesLengths[e];
            } else {
                this_flux = (*edges_flux_functions[e])(q_L, q_R) * m.edgesLengths[e];
            }
            for (int i=0; i<4; ++i) {
                RhoVector(k0+i) -= this_flux(i);
                if (edges_flux_functions[e]->two_sided)
                    RhoVector(k1+i) += this_flux(i);
            }
        }
    }


    #pragma omp parallel for
    for (int ci=m.nRealCells; ci<m.cellsAreas.size(); ++ci) {
        int k = 4*ci;
        // Dirichlet for boundary cell
        for (int j=0; j<4; ++j) {
            RhoVector.coeffRef(k+j) = 0;
        }
    }
}








int implicitSolver::compute() {
    rho_solver.compute(RhoMatrix);

    if (rho_solver.info() != 0) {
        return -1;
    }
    return 0;
}


double implicitSolver::solve(
    const double relaxation,
    const double tol
) {
    // Solve rho, rho_u, rho_v, rho_e
    fillRhoRHS();

    double err = RhoVector.norm();
    if (err < tol) return err;
    double err_ini = err;

    qW = rho_solver.solve(RhoVector);

    if (rho_solver.info() != 0) {
        return -1;
    } else {
        q += qW*relaxation;
    }

    for (uint i=0; i<10; ++i) {
        fillRhoRHS();

        err = RhoVector.norm();
        if (err < tol) return err;
        if (err > 10*err_ini) return err;

        qW = rho_solver.solve(RhoVector);
        if (rho_solver.info() != 0) {
            return -1;
        } else {
            q += qW*relaxation;
        }
    }

    #ifndef RANS_SPARSELU
        #ifdef RANS_DEBUG
            std::cout << "(iters=" << rho_solver.iterations() << ")";
        #endif
    #endif

    fillRhoRHS();
    return RhoVector.norm();
}








}

