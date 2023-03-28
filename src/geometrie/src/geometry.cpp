#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "geometrie/geometry.hpp"

#ifndef M_PI
#define M_PI 3.141592653589793115997963468544185161590576171875
#endif

// Initialize the methods above to avoid issues when calling
std::vector<double> bernstein_coefficients(double z_te, double r_le, double beta, double c, int order, double scale);
std::vector<double> bernstein_function(std::vector<double> eps, std::vector<double> coefficients);
std::vector<double> c_function(std::vector<double> eps, std::vector<double> ns);
std::vector<double> distribute_points(double start, double end, int nb, std::string mode);
std::vector<double> linspace(double start, double end, int nb);



CSTAirfoil::CSTAirfoil() {
    this->is_used = false;
    this->type = "CST";
    this->C = {1.0};
    this->Z_TE = {0.0};
    this->R_LE = {0.0};
    this->BETA = {0.0};
    this->COEFFICIENTS = {1.0, 1.0};

    // Class function coefficients
    this->N1 = {0.5};
    this->N2 = {1.0};
    this->NC1 = 0.0;
    this->NC2 = 0.0;
};

CSTAirfoil::CSTAirfoil(std::vector<double> c, std::vector<double> z_te, std::vector<double> r_le, 
std::vector<double> beta, std::vector<double> coefficients, 
std::vector<double> n1, std::vector<double> n2, double nc1, double nc2) {
    this->is_used = true;
    this->type = "CST";
    this->C = c;
    this->Z_TE = z_te;
    this->R_LE = r_le;
    this->BETA = beta;
    this->COEFFICIENTS = coefficients;

    // Class function coefficients
    this->N1 = n1;
    this->N2 = n2;
    this->NC1 = nc1;
    this->NC2 = nc2;
};

void CSTAirfoil::flip() {
    std::reverse(C.begin(),C.end());
    std::reverse(Z_TE.begin(),Z_TE.end());
    std::reverse(R_LE.begin(),R_LE.end());
    std::reverse(BETA.begin(),BETA.end());
    std::reverse(N1.begin(),N1.end());
    std::reverse(N2.begin(),N2.end());
    double nc1 = this->NC2;
    double nc2 = this->NC1;
    this->NC1 = nc1;
    this->NC2 = nc2;
};

// Getters
bool CSTAirfoil::get_is_used() {return this->is_used;};
std::vector<double> CSTAirfoil::get_C() {return this->C;};
std::vector<double> CSTAirfoil::get_Z_TE() {return this->Z_TE;};
std::vector<double> CSTAirfoil::get_R_LE() {return this->R_LE;};
std::vector<double> CSTAirfoil::get_BETA() {return this->BETA;};
std::vector<double> CSTAirfoil::get_COEFFICIENTS() {return this->COEFFICIENTS;};
std::vector<double> CSTAirfoil::get_N1() {return this->N1;};
std::vector<double> CSTAirfoil::get_N2() {return this->N2;};
double CSTAirfoil::get_NC1() {return this->NC1;};
double CSTAirfoil::get_NC2() {return this->NC2;};

std::vector<double> CSTAirfoil::get_airfoil_profile(std::vector<double> eps_list, double zeta_t, std::vector<double> shape_coefficients, std::vector<double> class_coefficients, double eta_class) {
    std::vector<double> y;
    std::vector<double> shape_function;
    std::vector<double> class_function;
    double p;

    shape_function = bernstein_function(eps_list, shape_coefficients);
    class_function = c_function(eps_list, class_coefficients);
    for (int k=0; k < eps_list.size(); k++) {
        p = (shape_function[k] * class_function[k] * eta_class) + eps_list[k] * zeta_t;
        y.push_back(p);
    };
    return y;
};

std::vector<std::vector<std::vector<double>>> CSTAirfoil::build_surface(std::vector<std::vector<double>> eps, std::vector<std::vector<double>> eta, std::string surface_id) {
    std::vector<std::vector<std::vector<double>>> res;

    if (eps.size() == eta.size()) {
        std::vector<std::vector<double>> x;
        std::vector<std::vector<double>> s;

        for (int i=0; i < eps.size(); i++) {
            std::vector<double> eps_i;
            std::vector<double> eta_i;
            std::vector<double> c;
            std::vector<double> z_te;
            std::vector<double> r_le;
            std::vector<double> beta;
            std::vector<double> n1;
            std::vector<double> n2;
            std::vector<double> class_eta;

            eps_i = eps[i];
            eta_i = eta[i];
            c = bernstein_function(eta_i, this->C);
            z_te = bernstein_function(eta_i, this->Z_TE);
            r_le = bernstein_function(eta_i, this->R_LE);
            beta = bernstein_function(eta_i, this->BETA);
            n1 = bernstein_function(eta_i, this->N1);
            n2 = bernstein_function(eta_i, this->N2);
            class_eta = c_function(eta_i, std::vector<double> {this->NC1, this->NC2});

            if (eps[i].size() == eta[i].size()) {
                std::vector<double> x_i;
                std::vector<double> s_i;

                for (int j=0; j < eps[i].size(); j++) {
                    double eps_ij;
                    double x_ij;
                    double s_ij;
                    eps_ij = eps_i[j];

                    // If COEFFICIENTS attribute is empty, use r_le and beta
                    std::vector<double> a;
                    if (this->COEFFICIENTS.size() == 0) {a = {sqrt(2 * r_le[j] / c[j]), tan(beta[j]) + (z_te[j] / c[j])};}
                    else {a = this->COEFFICIENTS;};

                    x_ij = c[j] * eps_ij;
                    // s_ij = c[j] * get_airfoil_profile(std::vector<double> {eps_ij}, z_te[j] / c[j], a, std::vector<double> {n1[j], n2[j]}, class_eta[j])[0];
                    s_ij = c[j] * get_airfoil_profile(std::vector<double> {eps_ij}, 0, a, std::vector<double> {n1[j], n2[j]}, class_eta[j])[0];
                    if (surface_id == "Lower") {s_ij *= -1.0;};
                    s_ij += eps_ij * z_te[j];

                    x_i.push_back(x_ij);
                    s_i.push_back(s_ij);
                };
                x.push_back(x_i);
                s.push_back(s_i);
            };
        };
        res = {x, s};
    };
    return res;
};


// Default constructor, set to NACA0012 profile
NACAAirfoil::NACAAirfoil() {
    this->is_used = false;
    this->type = "NACA";
    this->C = {1.0};
    this->M = {0.0};
    this->P = {0.0};
    this->T = {12.0};
};

NACAAirfoil::NACAAirfoil(std::vector<double> c, std::vector<double> m, std::vector<double> p, std::vector<double> t) {
    this->is_used = true;
    this->type = "NACA";
    this->C = c;
    this->M = m;
    this->P = p;
    this->T = t;
};

// Getters
bool NACAAirfoil::get_is_used() {return this->is_used;};
std::vector<double> NACAAirfoil::get_C() {return this->C;};
std::vector<double> NACAAirfoil::get_M() {return this->M;};
std::vector<double> NACAAirfoil::get_P() {return this->P;};
std::vector<double> NACAAirfoil::get_T() {return this->T;};

void NACAAirfoil::flip() {
    std::reverse(C.begin(),C.end());
    std::reverse(M.begin(),M.end());
    std::reverse(P.begin(),P.end());
    std::reverse(T.begin(),T.end());
};
        
double NACAAirfoil::get_thickness(double eps_val, double t_val) {
    double yt_val;
    yt_val = 5 * t_val * (0.2969 * sqrt(eps_val) - 0.1260 * eps_val - 0.3516 * pow(eps_val, 2) +
    0.2843 * pow(eps_val, 3) - 0.1015 * pow(eps_val, 4));
    return yt_val;
};

double NACAAirfoil::get_camber(double eps_val, double m_val, double p_val) {
    double yc_val;
    if (0 <= eps_val < p_val) {
        yc_val = (m_val / pow(p_val, 2)) * (2 * p_val * eps_val - pow(eps_val, 2));
    } else if (p_val < eps_val <= 1) {
        yc_val = (m_val / pow(1 - p_val, 2)) * ((1 - 2 * p_val) + 2 * p_val * eps_val - pow(eps_val, 2));
    } else {
        yc_val = m_val;
    };
    return yc_val;
};

double NACAAirfoil::get_theta(double eps_val, double m_val, double p_val) {
    double dyc_dx;
    double theta_val;
    if (0 <= eps_val < p_val) {
        dyc_dx = (2 * m_val / pow(p_val, 2)) * (p_val - eps_val);
    } else if (p_val < eps_val <= 1) {
        dyc_dx = (2 * m_val / pow(1 - p_val, 2)) * (p_val - eps_val);
    } else {
        dyc_dx = 0;
    };
    theta_val = atan(dyc_dx);
    return theta_val;
};

std::vector<std::vector<std::vector<double>>> NACAAirfoil::build_surface(std::vector<std::vector<double>> eps, std::vector<std::vector<double>> eta, std::string surface_id) {
    std::vector<std::vector<std::vector<double>>> res;
    if (eps.size() == eta.size()) {
        std::vector<std::vector<double>> x;
        std::vector<std::vector<double>> s;

        for (int i=0; i < eps.size(); i++) {
            std::vector<double> c = bernstein_function(eta[i], this->C);
            std::vector<double> m = bernstein_function(eta[i], this->M);
            std::vector<double> p = bernstein_function(eta[i], this->P);
            std::vector<double> t = bernstein_function(eta[i], this->T);
            for (int k=0; k < m.size(); k++) {m[k] *= 0.01;};
            for (int k=0; k < p.size(); k++) {p[k] *= 0.1;};
            for (int k=0; k < t.size(); k++) {t[k] *= 0.01;};

            if (eps[i].size() == eta[i].size()) {
                std::vector<double> x_i;
                std::vector<double> s_i;

                for (int j=0; j < eta[i].size(); j++) {
                    double yc_ij = get_camber(eps[i][j], m[j], p[j]);
                    double theta_ij = get_theta(eps[i][j], m[j], p[j]);
                    double yt_ij = get_thickness(eps[i][j], t[j]);

                    // Flip yt if surface is lower
                    if (surface_id == "Lower") {yt_ij *= -1.0;};

                    x_i.push_back(c[j] * (eps[i][j] - yt_ij * sin(theta_ij)));
                    s_i.push_back(c[j] * (yc_ij + yt_ij * cos(theta_ij)));
                };
                x.push_back(x_i);
                s.push_back(s_i);
            };
        };
        res = {x, s};
    };
    return res;
};



Surface::Surface(int nx, int ny, double b, double d, std::vector<double> c, 
std::vector<double> x_le, std::vector<double> z_n, std::vector<double> delta_alpha_t, 
double nc1, double nc2, double ny1, double ny2, std::string identification, 
CSTAirfoil cst_airfoil, NACAAirfoil naca_airfoil) {
    this->Nx = nx;
    this->Ny = ny;
    this->CST_AIRFOIL = cst_airfoil;
    this->NACA_AIRFOIL = naca_airfoil;

    this->B = b;
    this->D = d;
    this->C = c;
    this->X_LE = x_le;
    this->Z_N = z_n;
    this->DELTA_ALPHA_T = delta_alpha_t;
    this->NC1 = nc1;
    this->NC2 = nc2;
    this->NY1 = ny1;
    this->NY2 = ny2;
    this->ID = identification;

    // Get adequate distributions
    this->x_distribution = "cartesian";
    this->y_distribution = "cartesian";
    this->is_flipped = false;

    // Set eps and eta distributions
    set_distributions();
};

double Surface::get_B() {return this->B;};
double Surface::get_D() {return this->D;};
std::vector<double> Surface::get_C() {return this->C;};

void Surface::flip_on_xz_plane() {
    std::reverse(C.begin(),C.end());
    std::reverse(X_LE.begin(),X_LE.end());
    std::reverse(Z_N.begin(),Z_N.end());
    std::reverse(DELTA_ALPHA_T.begin(),DELTA_ALPHA_T.end());
            
    this->D = - (this->B + this->D);
    this->CST_AIRFOIL.flip();
    this->NACA_AIRFOIL.flip();
    double nc1;
    double nc2;
    nc1 = this->NC2;
    nc2 = this->NC1;
    this->NC1 = nc1;
    this->NC2 = nc2;
};

void Surface::change_distributions(std::string x_dist, std::string y_dist) {
    this->x_distribution = x_dist;
    this->y_distribution = y_dist;
    set_distributions();
};

std::vector<std::vector<std::vector<double>>> Surface::build_complete_surface(std::string surface_type) {
    std::vector<std::vector<std::vector<double>>> new_airfoil;
    std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> r;
    // Pick correct Airfoil type to use. CST has priority over NACA.
    if (CST_AIRFOIL.get_is_used()) {
        new_airfoil = this->CST_AIRFOIL.build_surface(this->eps, this->eta, this->ID);
    } else if (NACA_AIRFOIL.get_is_used()) {
        new_airfoil = this->NACA_AIRFOIL.build_surface(this->eps, this->eta, this->ID);
    } else {
        std::cout << "No airfoil has been specified" << std::endl;
    };
    u = new_airfoil[0];
    r = new_airfoil[1];

    std::vector<std::vector<double>> x;
    std::vector<std::vector<double>> y;
    std::vector<std::vector<double>> s;
    for (int i=0; i < this->eps.size(); i++) {
        std::vector<double> x_i;
        std::vector<double> y_i;
        std::vector<double> s_i;

        std::vector<double> c = bernstein_function(eta[i], this->C);
        std::vector<double> x_le = bernstein_function(eta[i], this->X_LE);
        std::vector<double> z_n = bernstein_function(eta[i], this->Z_N);
        std::vector<double> delta_alpha_t = bernstein_function(eta[i], this->DELTA_ALPHA_T);
        std::vector<double> class_eta_y = c_function(eta[i], std::vector<double> {this->NY1, this->NY2});
        // std::transform(class_eta_y.begin(), class_eta_y.end(), class_eta_y.begin(), 
        // [](double element) { return element *= c_function_max(std::vector<double> {this->NY1, this->NY2}); });

        for (int j=0; j < this->eps[i].size(); j++) {
            double x_ij;
            double r_ij;
            double y_ij;
            double s_ij;

            // Rotate with alpha_delta_t and z_n
            x_ij = + (cos(M_PI * delta_alpha_t[j]) * u[i][j]) + (sin(M_PI * delta_alpha_t[j]) * r[i][j]);
            r_ij = - (sin(M_PI * delta_alpha_t[j]) * u[i][j]) + (cos(M_PI * delta_alpha_t[j]) * r[i][j]) + z_n[j];

            if (surface_type == "Revolution") {
                x_ij = x_ij + x_le[j];
                y_ij = (this->D + this->B * cos(2 * M_PI * eta[i][j]) * r_ij) / 2;
                s_ij = (this->D + this->B * sin(2 * M_PI * eta[i][j]) * r_ij) / 2;
            } else {
                x_ij = (class_eta_y[j] * x_ij) + (x_le[j] + (c[j] / 2) * (1 - class_eta_y[j]));
                y_ij = (this->D + this->B * eta[i][j]) / 2;
                s_ij = r_ij;
            };

            x_i.push_back(x_ij);
            y_i.push_back(y_ij);
            s_i.push_back(s_ij);
        };
        x.push_back(x_i);
        y.push_back(y_i);
        s.push_back(s_i);
    };
    return {x, y, s};
};

void Surface::set_distributions() {
    std::string mode = x_distribution;
    std::string shape = y_distribution;

    std::vector<double> eps_base;
    std::vector<double> eta_base;
                
    if (mode == "distribute") {
        eps_base = distribute_points(0., 1., Nx, mode);
        eta_base = distribute_points(0., 1., Ny, mode);
    } else if (mode == "partial") {
        eps_base = distribute_points(0., 1., Nx, mode);
        eta_base = linspace(0., 1., Ny);
    } else if (mode == "crushed") {
        eps_base = distribute_points(0., 1., Nx, mode);
        eta_base = linspace(0. + 0.005, 1. - 0.005, Ny);
    } else if (mode == "sphere") {
        eps_base = distribute_points(0., 1., Nx, mode);
        eta_base = distribute_points(0., 1., Ny, mode);
    } else {
        eps_base = linspace(0., 1., Nx);
        eta_base = linspace(0., 1., Ny);
    };

    std::vector<std::vector<double>> eps;
    std::vector<std::vector<double>> eta;

    for (int i = 0; i < eta_base.size(); i++) {
        std::vector<double> eps_i;
        eps_i = eps_base;
        std::vector<double> eta_i(eps_base.size(), eta_base[i]);

        eps.push_back(eps_i);
        eta.push_back(eta_i);
    };

    if (shape == "rounded") {
        for (int i = 0; i < eps.size(); i++) {
            for (int j = 0; j < eps[i].size(); j++) {
                double x_ij;
                double y_ij;
                double eps_ij;
                double eta_ij;
                y_ij = eta[i][j] - 0.5;
                x_ij = eps[i][j] - 0.5;

                // Correct positions for equal distance on circumference
                if (abs(x_ij) < abs(y_ij)) {
                    if (y_ij != 0) {
                        double correction;
                        correction = abs(tan(y_ij));
                        if (y_ij - x_ij * correction != 0) {
                            double delta_x;
                            delta_x = correction * (y_ij * y_ij - x_ij * x_ij) / (y_ij - x_ij * correction);
                            delta_x = (x_ij / y_ij) * delta_x;
                            x_ij = x_ij - delta_x;
                        };
                    };
                } else {
                    if (x_ij != 0) {
                        double correction;
                        correction = abs(tan(x_ij));
                        if (x_ij - y_ij * correction != 0) {
                            double delta_y;
                            delta_y = correction * (x_ij * x_ij - y_ij * y_ij) / (x_ij - y_ij * correction);
                            delta_y = (y_ij / x_ij) * delta_y;
                            y_ij = y_ij - delta_y;
                        };
                    };
                };

                double x_ij_new;
                double y_ij_new;
                if (x_ij * x_ij + y_ij * y_ij == 0) {
                    x_ij_new = 0;
                    y_ij_new = 0;
                } else {
                    y_ij_new = y_ij * std::max(abs(x_ij), abs(y_ij)) / sqrt(x_ij * x_ij + y_ij * y_ij);
                    if (abs(y_ij_new) == 0.5) {
                        x_ij_new = 0;
                    } else {
                        x_ij_new = x_ij * std::max(abs(x_ij), abs(y_ij)) / sqrt((x_ij * x_ij + y_ij * y_ij) * (1 - 4 * (y_ij_new * y_ij_new)));
                    };
                };

                eta_ij = y_ij_new + 0.5;
                eps_ij = x_ij_new + 0.5;

                eta[i][j] = std::max(0., std::min(1., eta_ij));
                eps[i][j] = std::max(0., std::min(1., eps_ij));
            };
        };
    };
    this->eps = eps;
    this->eta = eta;
};



Body::Body(std::string body_type) {
    this->type = body_type;
    this->is_flipped = false;
};

std::vector<Surface> Body::get_surfaces() {return this->surfaces;};

void Body::set_to_general_solid() {this->type = "General";};
void Body::set_to_revolution_solid() {this->type = "Revolution";};

void Body::add_surface_cst(int nx, int ny, double b, double d, std::vector<double> c, std::vector<double> z_te, 
std::vector<double> r_le, std::vector<double> beta, std::vector<double> coefficients, 
std::vector<double> x_le, std::vector<double> z_n, std::vector<double> delta_alpha_t,
std::vector<double> n1, std::vector<double> n2, double nc1, double nc2, double ny1, double ny2, std::string identification) {
    CSTAirfoil cst_airfoil(c, z_te, r_le, beta, coefficients, n1, n2, nc1, nc2);
    NACAAirfoil naca_airfoil;
    Surface new_surface(nx, ny, b, d, c, x_le, z_n, delta_alpha_t, nc1, nc2, ny1, ny2, identification, cst_airfoil, naca_airfoil);
    this->surfaces.push_back(new_surface);
};

void Body::add_surface_naca(std::vector<double> m, std::vector<double> p, std::vector<double> t, int nx, int ny, 
double b, double d, std::vector<double> c, std::vector<double> x_le, std::vector<double> z_n, 
std::vector<double> delta_alpha_t, double nc1, double nc2, double ny1, double ny2, std::string identification) {
    CSTAirfoil cst_airfoil;
    NACAAirfoil naca_airfoil(c, m, p, t);
    Surface new_surface(nx, ny, b, d, c, x_le, z_n, delta_alpha_t, nc1, nc2, ny1, ny2, identification, cst_airfoil, naca_airfoil);
    this->surfaces.push_back(new_surface);
};

void Body::add_wing_surface_cst(int nx, int ny, double b, double d, std::vector<double> c, 
std::vector<double> z_te_half, std::vector<double> r_le, std::vector<double> beta, 
std::vector<double> x_le, std::vector<double> z_n, std::vector<double> delta_alpha_t) {
    // Get current number of surfaces in Body
    // New surfaces will be stored at indices nb (upper) and nb+1 (lower)
    int nb = this->surfaces.size();
    // Set associativity of surfaces
    this->associativity.push_back(std::vector<int> {nb, nb+1});

    // Append surfaces to list
    std::vector<double> empty;
    add_surface_cst(nx, ny, b, d, c, z_te_half, r_le, beta, empty, x_le, z_n, delta_alpha_t, {0.5}, {1.0}, 0, 0, 0, 0, "Upper");
    
    std::transform(z_te_half.begin(), z_te_half.end(), z_te_half.begin(), [](double element) { return element *= -1; });
    std::transform(beta.begin(), beta.end(), beta.begin(), [](double element) { return element *= -1; });

    for (int i=0; i < z_te_half.size(); i++) {std::cout << z_te_half[i] << " ";};
    std::cout << std::endl;

    add_surface_cst(nx, ny, b, d, c, z_te_half, r_le, beta, empty, x_le, z_n, delta_alpha_t, {0.5}, {1.0}, 0, 0, 0, 0, "Lower");
};

void Body::add_general_surface_cst(int nx, int ny, double b, double d, std::vector<double> c, 
std::vector<double> z_te_u, std::vector<double> z_te_l, std::vector<double> coefficients_u, std::vector<double> coefficients_l, 
std::vector<double> x_le, std::vector<double> z_n, std::vector<double> delta_alpha_t, 
std::vector<double> n1, std::vector<double> n2, double nc1, double nc2, double ny1, double ny2) {
    // Get current number of surfaces in Body
    // New surfaces will be stored at indices nb (upper) and nb+1 (lower)
    int nb = this->surfaces.size();
    // Set associativity of surfaces
    this->associativity.push_back(std::vector<int> {nb, nb+1});

    // Append surfaces to list
    std::vector<double> empty;
    add_surface_cst(nx, ny, b, d, c, z_te_u, empty, empty, coefficients_u, x_le, z_n, delta_alpha_t, n1, n2, nc1, nc2, ny1, ny2, "Upper");
    add_surface_cst(nx, ny, b, d, c, z_te_l, empty, empty, coefficients_l, x_le, z_n, delta_alpha_t, n1, n2, nc1, nc2, ny1, ny2, "Lower");
};

void Body::add_wing_surface_naca(int nx, int ny, double b, double d, std::vector<double> c, 
std::vector<double> m, std::vector<double> p, std::vector<double> t,
std::vector<double> x_le, std::vector<double> z_n, std::vector<double> delta_alpha_t) {
    // Get current number of surfaces in Body
    // New surfaces will be stored at indices nb (upper) and nb+1 (lower)
    int nb = this->surfaces.size();
    // Set associativity of surfaces
    this->associativity.push_back(std::vector<int> {nb, nb+1});

    // Append surfaces to list
    add_surface_naca(m, p, t, nx, ny, b, d, c, x_le, z_n, delta_alpha_t, 0, 0, 0, 0, "Upper");
    add_surface_naca(m, p, t, nx, ny, b, d, c, x_le, z_n, delta_alpha_t, 0, 0, 0, 0, "Lower");
};
    
void Body::mirror_body() {
    // Update is_flipped state
    if (this->is_flipped) {this->is_flipped = false;}
    else {this->is_flipped = true;};
    // Flip surfaces
    for (int i=0; i < this->surfaces.size(); i++) {this->surfaces[i].flip_on_xz_plane();};
};

void Body::change_all_distributions(std::string x_dist, std::string y_dist) {
    // Change distributions of all surfaces in body
    for (int i=0; i < this->surfaces.size(); i++) {this->surfaces[i].change_distributions(x_dist, y_dist);};
};

std::vector<std::vector<std::vector<std::vector<double>>>> Body::get_paired_body_surfaces() {
    std::vector<std::vector<std::vector<std::vector<double>>>> paired_surfaces;
    for (int i=0; i < this->associativity.size(); i++) {
        int su_id = this->associativity[i][0];
        int sl_id = this->associativity[i][1];
        std::string surface_type = this->type;
        std::vector<std::vector<std::vector<double>>> surfaces_u = this->surfaces[su_id].build_complete_surface(surface_type);
        std::vector<std::vector<std::vector<double>>> surfaces_l = this->surfaces[sl_id].build_complete_surface(surface_type);
        std::vector<std::vector<std::vector<double>>> paired_surface{surfaces_u[0], surfaces_l[0], surfaces_u[1], surfaces_l[1], surfaces_u[2], surfaces_l[2]};
        paired_surfaces.push_back(paired_surface);
    };
    return paired_surfaces;
};

std::vector<std::vector<std::vector<std::vector<double>>>> Body::get_body_surfaces() {
    std::vector<std::vector<std::vector<std::vector<double>>>> all_surfaces;
    for (int i=0; i < this->surfaces.size(); i++) {
        std::vector<std::vector<std::vector<double>>> new_surface = this->surfaces[i].build_complete_surface(this->type);
        all_surfaces.push_back(new_surface);
    };
    return all_surfaces;
};



std::vector<double> linspace(double start, double end, int nb) {
    std::vector<double> res(nb, 0);
    //for (int i=0; i < nb; i++) {res[i] = i * (end - start) / (nb - 1);};
    for (int i=0; i < nb; i++) {res[i] = start + i * (end - start) / (nb - 1);};
    return res;
};


std::vector<double> distribute_points(double start, double end, int nb, std::string mode) {
    std::vector<double> points_init;
    std::vector<double> points_final;
    double point_new;
    points_init = linspace(-1, 1, nb);

    for (int i=0; i < points_init.size(); i++) {
        if (mode == "sphere") {
            point_new = start + (end - start) * (0.5 * (1 + sin(M_PI * points_init[i] / 2)));
        }
        else {
            point_new = points_init[i] / (1 + pow(points_init[i], 2));
            point_new = start + (end - start) * (0.5 + point_new);
        };
        // Make sure new points are bounded
        point_new = std::max(start, std::min(end, point_new));
        points_final.push_back(point_new);
    };
    return points_final;
};




std::vector<double> c_function(std::vector<double> eps, std::vector<double> ns) {
    std::vector<double> y;
    for (int i=0; i < eps.size(); i++) {
        if (ns[0] == 0 && eps[i] == 0) {y.push_back(1);}
        else if (ns[1] == 0 && eps[i] == 1) {y.push_back(1);}
        else {y.push_back(pow(eps[i], ns[0]) * pow((1 - eps[i]), ns[1]));};
    };
    return y;
};


double factorial(int n) {
    double res;
    res = 1.;
    for (int i=0; i < n; i++) {res *= (i+1)*1.0;};
    return res;
};


std::vector<double> bernstein_function(std::vector<double> eps, std::vector<double> coefficients) {
    std::vector<double> res(eps.size(), 0);
    int order;
    double bern_i;
    order = coefficients.size();
    for (int i=0; i < order; i++) {
        bern_i = factorial(order - 1) / (factorial(i) * factorial(order - 1 - i));
        for (int j=0; j < eps.size(); j++) {
            res[j] += bern_i * coefficients[i] * pow(eps[j], i) * pow(1 - eps[j], order - 1 - i);
        };
    };
    return res;
};


std::vector<double> bernstein_coefficients(double z_te, double r_le, double beta, double c, int order, double scale) {
    std::vector<double> a(order, scale);
    a[0] = (sqrt(2 * r_le / c));
    a[order - 1] = tan(beta) + (z_te / c);
    return a;
};

