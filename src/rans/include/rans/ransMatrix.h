/*

        ____     ____  ___    _   _______
       /  _/    / __ \/   |  / | / / ___/
       / /_____/ /_/ / /| | /  |/ /\__ \
     _/ /_____/ _, _/ ___ |/ /|  /___/ /
    /___/    /_/ |_/_/  |_/_/ |_//____/

    Implicit RANS Solver
    Description: RANS matrix for matrix-free solvers
    By: Alexis Angers

*/
#pragma once


#include <rans/core.h>
#include <rans/physics.h>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>




namespace rans {

struct edgeInfo {
    uint cell0;
    uint cell1;
    flux *flux_function;
    Eigen::VectorXd n;
    double length;

    edgeInfo(uint& a, uint& b, flux* c, Eigen::VectorXd& d, double& e) {
        cell0 = a;
        cell1 = b;
        flux_function = c;
        n = d;
        length = e;
    }
};

}



class RansMatrix;
//using Eigen::SparseMatrix;

namespace Eigen {
    namespace internal {
    // RansMatrix looks-like a SparseMatrix, so let's inherits its traits:
    template<>
    struct traits<RansMatrix> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
    {};
    }
}

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class RansMatrix : public Eigen::EigenBase<RansMatrix> {
    public:
    rans::mesh *m;
    std::vector<std::unique_ptr<rans::flux>> *edges_flux_functions;

    rans::gas* g;
    Eigen::VectorXd *q;
    Eigen::VectorXd *dt;

    public:
    // Required typedefs, constants, and method:
    typedef double Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;
    enum {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    Index rows() const { return m->cellsAreas.size()*4; }
    Index cols() const { return m->cellsAreas.size()*4; }

    template<typename Rhs>
    Eigen::Product<RansMatrix,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
        return Eigen::Product<RansMatrix,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
    }

    // Custom API:
    RansMatrix() : mp_mat(0) {}

    void attachMyMatrix(const Eigen::SparseMatrix<double> &mat) {
        mp_mat = &mat;
    }

    const Eigen::SparseMatrix<double> my_matrix() const { return *mp_mat; }

    uint get_n_cells() const {return m->cellsAreas.size();}
    uint get_n_edges() const {return m->edgesLengths.size();}

    rans::edgeInfo get_edge(const uint e) const {
        Eigen::VectorXd n(2);
        n(0) = m->edgesNormalsX[e];
        n(1) = m->edgesNormalsY[e];
        return rans::edgeInfo(m->edgesCells(e, 0), m->edgesCells(e, 1), &(*edges_flux_functions->operator[](e)), n, m->edgesLengths[e]);
    }

    void set_flux_functions(std::vector<std::unique_ptr<rans::flux>>* ffx) {edges_flux_functions = ffx;}

    void set_mesh(rans::mesh *mi) {m = mi;}

    void set_gas(rans::gas *gi) {g = gi;}

    void set_q(Eigen::VectorXd *qi) {q = qi;}

    void set_dt(Eigen::VectorXd *dti) {dt = dti;}

private:
    const Eigen::SparseMatrix<double> *mp_mat;
};


// Implementation of RansMatrix * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
    namespace internal {

    template<typename Rhs>
    struct generic_product_impl<RansMatrix, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
    : generic_product_impl_base<RansMatrix,Rhs,generic_product_impl<RansMatrix,Rhs> >
    {
        typedef typename Product<RansMatrix,Rhs>::Scalar Scalar;

        template<typename Dest>
        static void scaleAndAddTo(Dest& dst, const RansMatrix& lhs, const Rhs& rhs, const Scalar& alpha)
        {
            // This method should implement "dst += alpha * lhs * rhs" inplace,
            // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
            #ifndef NDEBUG
                assert(alpha==Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);
            #endif

            // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
            // but let's do something fancier (and less efficient):

            for (uint i=0; i<lhs.m->nRealCells; ++i) {
                dst.segment(4*i, 4) = rhs.segment(4*i, 4) * lhs.m->cellsAreas[i]/(*lhs.dt)(i);
            }
            for (uint i=lhs.m->nRealCells; i<lhs.m->cellsAreas.size(); ++i) {
                dst.segment(4*i, 4) = rhs.segment(4*i, 4);
            }

            const double d = rhs.dot(rhs);
            const double h = 1e-6;

            for (uint e=0; e<lhs.get_n_edges(); ++e) {

                auto [i, j, ffx, n, l] = lhs.get_edge(e);

                Eigen::VectorXd q_L = lhs.q->segment(4*i, 4);
                Eigen::VectorXd q_R = lhs.q->segment(4*j, 4);

                Eigen::VectorXd W_L = rhs.segment(4*i, 4);
                Eigen::VectorXd W_R = rhs.segment(4*j, 4);

                Eigen::VectorXd FB = (*ffx)(q_L, q_R);

                dst.segment(4*i, 4) += ((*ffx)(q_L + h*W_L, q_R) - FB)/h;

                if (ffx->two_sided) {

                    dst.segment(4*j, 4) -= ((*ffx)(q_L, q_R + h*W_R) - FB)/h;
                }
            }
        }
    };

    }
}


