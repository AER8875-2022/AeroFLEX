
#pragma once
#include <tuple>

#include <Eigen/Dense>
#include <rans/mesh.h>

namespace rans {



class cell_handler {

    mesh& m;
    Eigen::VectorXd& q;
    std::vector<std::unique_ptr<flux>>& flux_functions;
    std::vector<Eigen::Matrix2d>& leastSquaresMatrices;

    const bool second_order;
    const uint id;

public:

    cell_handler(
        mesh& m, 
        Eigen::VectorXd& q, 
        std::vector<std::unique_ptr<flux>>& flux_functions, 
        std::vector<Eigen::Matrix2d>& leastSquaresMatrices, 
        bool second_order,
        uint& id
    ) : m(m), q(q), flux_functions(flux_functions), leastSquaresMatrices(leastSquaresMatrices), second_order(second_order), id(id) 
    {}

    uint size() {return m.cellsIsTriangle[id] ? 3 : 4;}
    uint get_id() const {return id;}

    Eigen::MatrixXd gradients();

    Eigen::VectorXd limiters(const Eigen::MatrixXd& grad);

};


Eigen::MatrixXd cell_handler::gradients() {
    // For now, assume least squares gradient scheme
    Eigen::MatrixXd grads(4, 2);

    if (id >= m.nRealCells) return grads;

    Eigen::MatrixXd d(size(), 2);
    Eigen::Vector2d grad_i;

    for (uint j=0; j<size(); ++j) {

        const uint e = m.cellsEdges(id, j);

        const uint cell_p = id;
        const uint cell_n = m.edgesCells(e, 0) == id ? m.edgesCells(e, 1) : m.edgesCells(e, 0);

        d(j, 0) = m.cellsCentersX[cell_n] - m.cellsCentersX[cell_p];
        d(j, 1) = m.cellsCentersY[cell_n] - m.cellsCentersY[cell_p];
    }

    Eigen::MatrixXd deltas_k(size(), 4);

    for (uint j=0; j<size(); ++j) {
        const uint e = m.cellsEdges(id, j);

        const uint cell_p = id;
        const uint cell_n = m.edgesCells(e, 0) == id ? m.edgesCells(e, 1) : m.edgesCells(e, 0);

        Eigen::VectorXd q_L = q.segment(4*cell_p, 4);
        Eigen::VectorXd q_R = flux_functions[e]->vars(q_L, q.segment(4*cell_n, 4), Eigen::VectorXd::Zero(4), Eigen::VectorXd::Zero(4));
        
        for (uint k=0; k<4; ++k) {
            deltas_k(j, k) = q_R(k) - q_L(k);
        }
    }

    for (uint k=0; k<4; ++k) {
        grad_i = leastSquaresMatrices[id] * d.transpose() * deltas_k.col(k);
        grads(k, 0) = grad_i(0);
        grads(k, 1) = grad_i(1);
    }

    return grads;
}



Eigen::VectorXd cell_handler::limiters(
    const Eigen::MatrixXd& grad
) {
    Eigen::VectorXd limiters(4);

    if (id >= m.nRealCells) return limiters;

    for (int i=0; i<4; ++i) limiters(i) = 1.;

    // Compute qmin and qmax
    Eigen::VectorXd q_min = q.segment(4*id, 4);
    Eigen::VectorXd q_max = q.segment(4*id, 4);

    for (uint j=0; j<size(); ++j) {
        const uint e = m.cellsEdges(id, j);

        const uint cell_n = m.edgesCells(e, 0) == id ? m.edgesCells(e, 1) : m.edgesCells(e, 0);

        for (uint k=0; k<4; ++k) {
            q_min(k) = std::min(q_min(k), q(4*cell_n+k));
            q_max(k) = std::max(q_max(k), q(4*cell_n+k));
        }
    }

    // Compute limiters
    for (uint j=0; j<size(); ++j) {
        const uint e = m.cellsEdges(id, j);

        const double dx = m.edgesCentersX[e] - m.cellsCentersX[id];
        const double dy = m.edgesCentersY[e] - m.cellsCentersY[id];
        const double sqrt_area = sqrt(m.cellsAreas[id]);

        for (uint k=0; k<4; ++k) {
            double dqg = grad(k,0)*dx + grad(k,1)*dy;

            double delta_max = q_max(k) - q(4*id+k);
            double delta_min = q_min(k) - q(4*id+k);

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

            limiters(k) = std::min(limiters(k), lim);
        }
    }

    return limiters;
}







class edge_handler {

    mesh& m;
    Eigen::VectorXd& q;
    std::vector<std::unique_ptr<flux>>& flux_functions;
    std::vector<Eigen::Matrix2d>& leastSquaresMatrices;

    const bool second_order;

    const uint id;

public:

    edge_handler(
        mesh& m, 
        Eigen::VectorXd& q, 
        std::vector<std::unique_ptr<flux>>& flux_functions, 
        std::vector<Eigen::Matrix2d>& leastSquaresMatrices, 
        bool& second_order,
        uint& id
    ) : m(m), q(q), flux_functions(flux_functions), leastSquaresMatrices(leastSquaresMatrices), second_order(second_order), id(id) 
    {}

    uint get_id() const {return id;}

    Eigen::MatrixXd average_gradients(const Eigen::MatrixXd& grad0, const Eigen::MatrixXd& grad1);

    Eigen::VectorXd get_flux();

    Eigen::MatrixXd jacobian();

};

Eigen::MatrixXd edge_handler::average_gradients(const Eigen::MatrixXd& grad0, const Eigen::MatrixXd& grad1) {
    Eigen::MatrixXd grad(4, 2);

    uint cell0 = m.edgesCells(id, 0);
    uint cell1 = m.edgesCells(id, 1);

    if (cell0 == cell1) {
        for (uint i=0; i<4; ++i) {
            grad(i,0) = grad0(i, 0);
            grad(i,1) = grad0(i, 1);
        }
        return grad;
    }
    
    double tij[2];
    tij[0] = m.cellsCentersX[cell1] - m.cellsCentersX[cell0];
    tij[1] = m.cellsCentersY[cell1] - m.cellsCentersY[cell0];

    const double lij = sqrt(tij[0]*tij[0] + tij[1]*tij[1]);

    tij[0] /= lij;
    tij[1] /= lij;

    // Directional derivative
    double grad_dir[4];
    for (uint i=0; i<4; ++i) {
        grad_dir[i] = (q(4*cell1+i) - q(4*cell0+i))/lij;
    }

    // Arithmetic average gradient
    double gradx_bar[4];
    double grady_bar[4];
    for (uint i=0; i<4; ++i) {
        gradx_bar[i] = (grad0(i,0) + grad1(i,0))*0.5;
        grady_bar[i] = (grad0(i,1) + grad1(i,1))*0.5;
    }

    // Corrected gradient
    for (uint i=0; i<4; ++i) {
        const double grad_dot = gradx_bar[i]*tij[0] + grady_bar[i]*tij[1];
        grad(i,0) = gradx_bar[i] - (grad_dot - grad_dir[i])*tij[0];
        grad(i,1) = grady_bar[i] - (grad_dot - grad_dir[i])*tij[1];
    }

    return grad;
}


Eigen::VectorXd edge_handler::get_flux() {
    uint cell0_id = m.edgesCells(id, 0);
    uint cell1_id = m.edgesCells(id, 1);

    /**/
    cell_handler cell0(m, q, flux_functions, leastSquaresMatrices, second_order, cell0_id);
    cell_handler cell1(m, q, flux_functions, leastSquaresMatrices, second_order, cell1_id);

    Eigen::MatrixXd grad0 = cell0.gradients();
    Eigen::MatrixXd grad1 = cell1.gradients();

    Eigen::MatrixXd grads = average_gradients(grad0, grad1);
    /**/

    Eigen::VectorXd q0 = q.segment(4*cell0_id, 4);
    Eigen::VectorXd q1 = q.segment(4*cell1_id, 4);

    // Linear interpolate with limiters
    /*
    if (second_order) {
        Eigen::VectorXd lim0 = cell0.limiters(grad0);
        Eigen::VectorXd lim1 = cell1.limiters(grad1);
        for (int i=0; i<4; ++i) {
            {
                const double dx = m.edgesCentersX[id] - m.cellsCentersX[cell0.get_id()];
                const double dy = m.edgesCentersY[id] - m.cellsCentersY[cell0.get_id()];
                q0(i) += lim0(i)*(grad0(i,0)*dx + grad0(i,1)*dy)*0.1;
            }
            {
                const double dx = m.edgesCentersX[id] - m.cellsCentersX[cell1.get_id()];
                const double dy = m.edgesCentersY[id] - m.cellsCentersY[cell1.get_id()];
                q1(i) += lim1(i)*(grad1(i,0)*dx + grad1(i,1)*dy)*0.1;
            }
        }
    }
    /**/

    return flux_functions[id]->operator()(
        q0, q1,
        grads.block(0,0,4,1), grads.block(0,1,4,1)
    );
}

Eigen::MatrixXd edge_handler::jacobian() {
    Eigen::MatrixXd J(8, 8);

    uint cell0_id = m.edgesCells(id, 0);
    uint cell1_id = m.edgesCells(id, 1);

    const Eigen::VectorXd f = get_flux();

    for (int i=0; i<4; ++i) {
        // Part of cell 0
        {
            const double update = std::max(1e-6, abs(q(4*cell0_id+i))*1e-6);
            q(4*cell0_id+i) += update;
            const Eigen::VectorXd fp = get_flux();
            q(4*cell0_id+i) -= update;

            J.block(0, i, 4, 1) = (fp  - f)/update;
            J.block(4, i, 4, 1) = -J.block(0, i, 4, 1);
        }
        // Part of cell 1
        {
            const double update = std::max(1e-6, abs(q(4*cell1_id+i))*1e-6);
            q(4*cell1_id+i) += update;
            const Eigen::VectorXd fp = get_flux();
            q(4*cell1_id+i) -= update;

            J.block(0, i+4, 4, 1) = (fp  - f)/update;
            J.block(4, i+4, 4, 1) = -J.block(0, i+4, 4, 1);
        }
    }

    return J;
}






}

