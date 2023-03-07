/*

        ____     ____  ___    _   _______
       /  _/    / __ \/   |  / | / / ___/
       / /_____/ /_/ / /| | /  |/ /\__ \
     _/ /_____/ _, _/ ___ |/ /|  /___/ /
    /___/    /_/ |_/_/  |_/_/ |_//____/

    Implicit RANS Solver
    Description: Core objects
    By: Alexis Angers

*/
#pragma once
#include <cmath>
#include <string>
#include <map>

using uint = unsigned int;

namespace rans {


struct gas {
    double R=0.71428571428;
    double mu_L=1e-5;
    double Pr_L=0.72;
    double Pr_T=0.9;
    double gamma=1.4;

    gas(double R_=0.71428571428, double mu_L_=1e-5, double Pr_L_=0.72, double Pr_T_=0.9, double gamma_=1.4) {
        R = R_;
        mu_L = mu_L_;
        Pr_L = Pr_L_;
        Pr_T = Pr_T_;
        gamma = gamma_;
    }
};


struct conservative_variables {
    double rho;
    double rhou;
    double rhov;
    double rhoe;

    conservative_variables(double rho_, double rhou_, double rhov_, double rhoe_) {
        rho = rho_; rhou = rhou_; rhov = rhov_; rhoe = rhoe_;
    }
};


struct boundary_variables {
    double mach=0;
    double angle=0;
    double T=1;
    double p=1.;

    boundary_variables() {}

    boundary_variables(double machi, double anglei=0, double Ti=1, double pi=1.) {
        mach = machi; angle = anglei; T = Ti; p = pi;
    }

    conservative_variables get_conservative(const gas& g) const {
        double c = sqrt(g.gamma * g.R * T);

        double u = mach * c * std::cos(angle);
        double v = mach * c * std::sin(angle);

        double rho = p / (g.R * T);

        double rhoE = p / (g.gamma - 1) + 0.5*rho*(u*u + v*v);
        return conservative_variables(rho, rho*u, rho*v, rhoE);
    }
};


struct adim_scale {

    double corde;

    double rho_scale;
    double rhou_scale;
    double rhov_scale;
    double rhoe_scale;

    adim_scale() {}
    adim_scale(const boundary_variables& var_far, const double& corde, const gas& g) {
        auto [rho, rhou, rhov, rhoe] = var_far.get_conservative(g);

        double ur = sqrt(rhou*rhou + rhov*rhov)/rho;

        rho_scale = corde / (rho * ur);
        rhou_scale = corde / (rho * ur * ur);
        rhov_scale = corde / (rho * ur * ur);
        rhoe_scale = corde / (rho * ur * ur * ur);
    }

};



struct boundary_condition {
    std::string bc_type;
    boundary_variables vars_far;
};



struct solution {
    double gamma;
    double R;

    Eigen::VectorXd rho;
    Eigen::VectorXd rhou;
    Eigen::VectorXd rhov;
    Eigen::VectorXd rhoe;

    double p(const int& i) const {
        return (gamma - 1)*(rhoe(i) - 0.5/rho(i)*(rhou(i)*rhou(i) + rhov(i)*rhov(i)));
    }

    double u(const int& i) const {
        return rhou(i)/rho(i);
    }
    double v(const int& i) const {
        return rhov(i)/rho(i);
    }
    double T(const int& i) const {
        return p(i)/(R*rho(i));
    }
    double c(const int& i) const {
        return sqrt(gamma * p(i)/rho(i));
    }
    double mach(const int& i) const {
        return sqrt(rhou(i)*rhou(i) + rhov(i)*rhov(i))/(rho(i) * c(i));
    }
};

struct Settings {
    gas g;

    std::vector<std::string> meshes;

    std::map<std::string, boundary_condition> bcs;

    std::string solver_type = "implicit";
    bool second_order = true;
    double relaxation = 0.8;

    std::string outfilename;

    int read_failure = 1;
};

}