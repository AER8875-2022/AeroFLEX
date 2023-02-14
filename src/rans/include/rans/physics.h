/*

        ____     ____  ___    _   _______
       /  _/    / __ \/   |  / | / / ___/
       / /_____/ /_/ / /| | /  |/ /\__ \
     _/ /_____/ _, _/ ___ |/ /|  /___/ /
    /___/    /_/ |_/_/  |_/_/ |_//____/

    Implicit RANS Solver
    Description: Flux and its jacobian matrix
    By: Alexis Angers

*/
#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdexcept>

#include <rans/core.h>



namespace rans {


// Crude evaluation using first-order scalar dissipation upwind scheme
class flux {
public:
    gas& g;
    double& nx;
    double& ny;

    bool two_sided;

    flux(gas& g, double& nx, double& ny) : g(g), nx(nx), ny(ny) {}

    virtual Eigen::VectorXd operator()(
        const Eigen::VectorXd& q_L,
        const Eigen::VectorXd& q_R,
        const double& nu_L = 0,
        const double& nu_R = 0
    ) {
        throw std::logic_error( "default flux called" );
        return Eigen::VectorXd(4);
    }

};

// Classes for template specializations
class inviscid {
};
class viscous {
};


template <class kind>
class internal_flux : public flux {
private:
    double entropy_correction(const double& l, const double& d) const {
        return l > d ? l : (l*l + d*d)/(2*d);
    }
public:
    internal_flux(gas& g, double& nx, double& ny) : flux(g, nx, ny) {
        two_sided = true;
    }
    inline Eigen::VectorXd operator()(
        const Eigen::VectorXd& q_L,
        const Eigen::VectorXd& q_R,
        const double& nu_L = 0,
        const double& nu_R = 0
    ) {
        throw std::logic_error( "default flux called" );
        return Eigen::VectorXd(4);
    }
};


template<>
inline Eigen::VectorXd internal_flux<inviscid>::operator()(
    const Eigen::VectorXd& q_L,
    const Eigen::VectorXd& q_R,
    const double& nu_L,
    const double& nu_R
) {
    Eigen::VectorXd f(4);

    const double& rho_L = q_L(0);
    const double& rho_u_L = q_L(1);
    const double& rho_v_L = q_L(2);
    const double& rho_e_L = q_L(3);

    const double& rho_R = q_R(0);
    const double& rho_u_R = q_R(1);
    const double& rho_v_R = q_R(2);
    const double& rho_e_R = q_R(3);

    const double V_L = (rho_u_L*nx + rho_v_L*ny)/rho_L;
    const double V_R = (rho_u_R*nx + rho_v_R*ny)/rho_R;

    const double p_L = (g.gamma - 1)*(rho_e_L - 0.5/rho_L*(rho_u_L*rho_u_L + rho_v_L*rho_v_L));
    const double p_R = (g.gamma - 1)*(rho_e_R - 0.5/rho_R*(rho_u_R*rho_u_R + rho_v_R*rho_v_R));
    
    f(0) = (V_L*rho_L + V_R*rho_R)*0.5;
    f(1) = (V_L*rho_u_L + p_L*nx + V_R*rho_u_R + p_R*nx)*0.5;
    f(2) = (V_L*rho_v_L + p_L*ny + V_R*rho_v_R + p_R*ny)*0.5;
    f(3) = (V_L*(rho_e_L + p_L) + V_R*(rho_e_R + p_R))*0.5;

    // Upwind flux
    const double pL = p_L;
    const double pR = p_R;

    // Roe variables
    const double uL = q_L(1)/q_L(0);
    const double uR = q_R(1)/q_R(0);
    const double vL = q_L(2)/q_L(0);
    const double vR = q_R(2)/q_R(0);

    const double srhoL = sqrt(q_L(0));
    const double srhoR = sqrt(q_R(0));
    const double rho = srhoR*srhoL;
    const double u = (uL*srhoL + uR*srhoR)/(srhoL + srhoR);
    const double v = (vL*srhoL + vR*srhoR)/(srhoL + srhoR);
    const double h = ((q_L(3) + pL)/q_L(0)*srhoL + (q_R(3) + pR)/q_R(0)*srhoR)/(srhoL + srhoR);
    const double q2 = u*u + v*v;
    const double c = sqrt( (g.gamma - 1.) * (h - 0.5*q2) );
    const double V = u*nx + v*ny;
    const double VR = uR*nx + vR*ny;
    const double VL = uL*nx + vL*ny;

    // Entropy correction
    const double lambda_cm = entropy_correction(abs(V-c), 0.05*c);
    const double lambda_c  = entropy_correction(abs(V), 0.05*c);
    const double lambda_cp = entropy_correction(abs(V+c), 0.05*c);

    const double kF1 = lambda_cm*((pR-pL) - rho*c*(VR-VL))/(2.*c*c);
    const double kF234_0 = lambda_c*((q_R(0) - q_L(0)) - (pR-pL)/(c*c));
    const double kF234_1 = lambda_c*rho;
    const double kF5 = lambda_cp*((pR-pL) + rho*c*(VR-VL))/(2*c*c);

    f(0) -= 0.5*(kF1            + kF234_0                                                       + kF5);
    f(1) -= 0.5*(kF1*(u-c*nx) + kF234_0*u      + kF234_1*(uR - uL - (VR-VL)*nx)             + kF5*(u+c*nx));
    f(2) -= 0.5*(kF1*(v-c*ny) + kF234_0*v      + kF234_1*(vR - vL - (VR-VL)*ny)             + kF5*(v+c*ny));
    f(3) -= 0.5*(kF1*(h-c*V)    + kF234_0*q2*0.5 + kF234_1*(u*(uR-uL) + v*(vR-vL) - V*(VR-VL))  + kF5*(h+c*V)); 

    return f;
}


template <class kind>
class slip_wall_flux : public flux {
protected:
    internal_flux<kind> invf;
public:
    slip_wall_flux(gas& g, double& nx, double& ny) : flux(g, nx, ny), invf(g, nx, ny) {
        two_sided = false;
    }

    inline Eigen::VectorXd operator()(
        const Eigen::VectorXd& q_L,
        const Eigen::VectorXd& _q_bc,
        const double& nu_L = 0,
        const double& nu_R = 0
    ) {

        Eigen::VectorXd q_R(4);

        // q_L is the internal field

        const double& rho = q_L(0);
        const double& rho_u = q_L(1);
        const double& rho_v = q_L(2);
        const double& rho_e = q_L(3);

        const double rhoV = rho_u*nx + rho_v*ny;

        const double rho_u_rev = rho_u - 2.*rhoV*nx;
        const double rho_v_rev = rho_v - 2.*rhoV*ny;

        q_R(0) = q_L(0);
        q_R(1) = rho_u_rev;
        q_R(2) = rho_v_rev;
        q_R(3) = q_L(3);

        return invf(q_L, q_R);
    }
};


template <class kind>
class farfield_flux : public flux {
protected:
    internal_flux<kind> invf;
public:
    farfield_flux(gas& g, double& nx, double& ny) : flux(g, nx, ny), invf(g, nx, ny) {
        two_sided = false;
    }

    inline Eigen::VectorXd operator()(
        const Eigen::VectorXd& q_L,
        const Eigen::VectorXd& q_bc,
        const double& nu_L = 0,
        const double& nu_R = 0
    ) {

        Eigen::VectorXd q_R(4);

        // q_R contains bc info
        // q_L is the internal field

        const double& rho = q_L(0);
        const double& rho_u = q_L(1);
        const double& rho_v = q_L(2);
        const double& rho_e = q_L(3);

        const double& bc_rho = q_bc(0);
        const double& bc_rhou = q_bc(1);
        const double& bc_rhov = q_bc(2);
        const double& bc_rhoe = q_bc(3);

        const double bc_u = bc_rhou / bc_rho;
        const double bc_v = bc_rhov / bc_rho;
        const double bc_p = (g.gamma - 1)*(bc_rhoe - 0.5/bc_rho*(bc_rhou*bc_rhou + bc_rhov*bc_rhov));

        const double p = (g.gamma - 1)*(rho_e - 0.5/rho*(rho_u*rho_u + rho_v*rho_v));

        const double c = sqrt(g.gamma * p / rho);
        const double mach = sqrt(rho_u*rho_u + rho_v*rho_v)/(rho*c);
        const double inlet_outlet = rho_u*nx + rho_v*ny;

        if (mach > 1) {
            // Supersonic
            if (inlet_outlet < 0) {
                // Inlet
                q_R(0) = bc_rho;
                q_R(1) = bc_rho * bc_u;
                q_R(2) = bc_rho * bc_v;
                q_R(3) = bc_p/(g.gamma - 1) + 0.5*bc_rho*(bc_u*bc_u + bc_v*bc_v);
            } else {
                // Outlet
                q_R(0) = rho;
                q_R(1) = rho_u;
                q_R(2) = rho_v;
                q_R(3) = rho_e;
            }
        } else {
            // Subsonic
            // a variables -> outside
            double pa = bc_p;
            double rhoa = bc_rho;
            double ua = bc_u;
            double va = bc_v;
            
            // d variables -> inside
            double pd = p;
            double rhod = rho;
            double ud = rho_u / rho;
            double vd = rho_v / rho;

            double rho0 = rho;
            double c0 = c;

            if (inlet_outlet < 0) {
                // Inlet
                double pb = 0.5*(pa + pd - rho0*c0*(nx*(ua - ud) + ny*(va - vd)));
                q_R(0) = rhoa + (pb - pa)/(c0*c0);
                q_R(1) = q_R(0) * (ua - nx*(pa - pb)/(rho0*c0));
                q_R(2) = q_R(0) * (va - ny*(pa - pb)/(rho0*c0));
                q_R(3) = pb / (g.gamma - 1) + 0.5/q_R(0)*(q_R(1)*q_R(1) + q_R(2)*q_R(2));
            } else {
                // Outlet
                double pb = pa;
                q_R(0) = rhod + (pb - pd)/(c0*c0);
                q_R(1) = q_R(0) * (ud + nx*(pd - pb)/(rho0*c0));
                q_R(2) = q_R(0) * (va + ny*(pd - pb)/(rho0*c0));
                q_R(3) = pb / (g.gamma - 1) + 0.5/q_R(0)*(q_R(1)*q_R(1) + q_R(2)*q_R(2));
            }
        }

        return invf(q_L, q_R);
    }
};



template <typename _flux>
inline Eigen::MatrixXd calc_convective_jacobian(
    _flux& flux_function,
    Eigen::VectorXd& q_L,
    Eigen::VectorXd& q_R,
    const double& nu_L = 0,
    const double& nu_R = 0
) {
    Eigen::MatrixXd J(8, 8);

    const Eigen::VectorXd f = flux_function(q_L, q_R, nu_L, nu_R);

    for (int i=0; i<4; ++i) {
        // Part of q_L
        {
            const double update = std::max(1e-6, abs(q_L(i))*1e-6);

            q_L(i) += update;
            const Eigen::VectorXd fp = flux_function(q_L, q_R, nu_L, nu_R);
            q_L(i) -= update;

            J.block(0, i, 4, 1) = (fp  - f)/update;
            J.block(4, i, 4, 1) = -J.block(0, i, 4, 1);
        }

        // Part of q_R
        {
            const double update = std::max(1e-6, abs(q_R(i))*1e-6);

            q_R(i) += update;
            const Eigen::VectorXd fp = flux_function(q_L, q_R, nu_L, nu_R);
            q_R(i) -= update;

            J.block(0, i+4, 4, 1) = (fp  - f)/update;
            J.block(4, i+4, 4, 1) = -J.block(0, i+4, 4, 1);
        }
    }

    return J;
}


#define RANS_YT 2.0

inline double michalak_limiter(const double& y) {

    if (y >= RANS_YT) {
        return 1.0;
    } else {
        constexpr double a = 1.0/(RANS_YT*RANS_YT) - 2.0/(RANS_YT*RANS_YT*RANS_YT);
        constexpr double b = -3.0/2.0*a*RANS_YT - 0.5/RANS_YT;
        return a*y*y*y + b*y*y + y;
    }
}

inline double michalak_limiter_derivative(const double& y) {

    if (y >= RANS_YT) {
        return 0;
    } else {
        constexpr double a = 1.0/(RANS_YT*RANS_YT) - 2.0/(RANS_YT*RANS_YT*RANS_YT);
        constexpr double b = -3.0/2.0*a*RANS_YT - 0.5/RANS_YT;
        return 3*a*y*y + 2*b*y + 1;
    }
}

inline double venkatakrishnan_limiter(const double& y) {
    return (y*y + 2*y)/(y*y + y + 2);
}



}
