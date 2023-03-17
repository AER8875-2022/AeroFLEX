
#pragma once
#include <tuple>

#include <Eigen/Dense>
#include <rans/mesh.h>



class edge_handler {

    mesh& m;
    Eigen::VectorXd& q;
    std::vector<std::unique_ptr<flux>>& flux_functions;

    const uint id;

public:

    edge_handler(mesh& m, Eigen::VectorXd& q, std::vector<std::unique_ptr<flux>>& flux_functions, uint& id)
    : m(m), q(q), flux_functions(flux_functions), id(id) 
    {}

    Eigen::MatrixXd average_gradients(const Eigen::MatrixXd& grad0, const Eigen::MatrixXd& grad1);

    Eigen::VectorXd flux();

    Eigen::MatrixXd jacobian();

};

Eigen::MatrixXd edge_handler::average_gradients(const Eigen::MatrixXd& grad0, const Eigen::MatrixXd& grad1) {
    Eigen::MatrixXd grad(4,2);

    return grad;
}


Eigen::VectorXd edge_handler::flux() {
    uint cell0_id = m.edgesCells(id, 0);
    uint cell1_id = m.edgesCells(id, 1);

    cell_handler cell0(m, q, flux_functions, cell0_id);
    cell_handler cell1(m, q, flux_functions, cell1_id);

    Eigen::MatrixXd grad0 = cell0.gradients();
    Eigen::MatrixXd grad1 = cell1.gradients();

    grads = average_gradients(grad0, grad1);

    return flux_functions[id]->operator()(
        q.segment(4*cell0_id, 4), q.segment(4*cell1_id, 4),
        grads.block(0,0,4,1), grads.block(0,1,4,1)
    );
}

Eigen::MatrixXd edge_handler::jacobian() {
    Eigen::MatrixXd J(8,8);

    uint cell0_id = m.edgeCells(id, 0);
    uint cell1_id = m.edgeCells(id, 1);

    const Eigen::VectorXd f = flux();
    for (int i=0; i<4; ++i) {
        // Part of cell 0
        {
            const double update = std::max(1e-6, abs(q(4*cell0_id+i))*1e-6);
            q(4*cell0_id+i) += update;
            const Eigen::VectorXd fp = flux();
            q(4*cell0_id+i) -= update;

            J.block(0, i, 4, 1) = (fp  - f)/update;
            J.block(4, i, 4, 1) = -J.block(0, i, 4, 1);
        }
        // Part of cell 1
        {
            const double update = std::max(1e-6, abs(q(4*cell1_id+i))*1e-6);
            q(4*cell1_id+i) += update;
            const Eigen::VectorXd fp = flux();
            q(4*cell1_id+i) -= update;

            J.block(0, i+4, 4, 1) = (fp  - f)/update;
            J.block(4, i+4, 4, 1) = -J.block(0, i+4, 4, 1);
        }
    }

    return J;
}


class cell_handler {

    mesh& m;
    Eigen::VectorXd& q;
    std::vector<std::unique_ptr<flux>>& flux_functions;

    const uint id;

public:

    cell_handler(mesh& m, Eigen::VectorXd& q, std::vector<std::unique_ptr<flux>>& flux_functions, uint& id)
    : m(m), q(q), flux_functions(flux_functions), id(id) 
    {}

    uint size() {return m.cellsIsTriangle[id] ? 3 : 4;}

    Eigen::MatrixXd gradients();

};


Eigen::MatrixXd cell_handler::gradients() {
    Eigen::MatrixXd grads(4,2);


    for (uint edge_id=0; edge_id<size(); ++edge_id) {
        
    }

    return grads;
}