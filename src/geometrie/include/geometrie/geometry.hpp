#ifndef geometry_HPP
#define geometry_HPP

#include <string>
#include <vector>
//
class CSTAirfoil
{
private:
    bool is_used;
    std::string type;
    std::vector<double> C;
    std::vector<double> Z_TE;
    std::vector<double> R_LE;
    std::vector<double> BETA;
    std::vector<double> COEFFICIENTS;
    std::vector<double> N1;
    std::vector<double> N2;
    double NC1;
    double NC2;
public:
    CSTAirfoil();
    CSTAirfoil(std::vector<double> c, std::vector<double> z_te, std::vector<double> r_le, 
        std::vector<double> beta, std::vector<double> coefficients, 
        std::vector<double> n1, std::vector<double> n2, double nc1, double nc2);  
    bool get_is_used();
    std::vector<double> get_C();
    std::vector<double> get_Z_TE();
    std::vector<double> get_R_LE();
    std::vector<double> get_BETA();
    std::vector<double> get_COEFFICIENTS();
    std::vector<double> get_N1();
    std::vector<double> get_N2();
    double get_NC1();
    double get_NC2();

    void flip();
    std::vector<double> get_airfoil_profile(std::vector<double> eps_list, double zeta_t, std::vector<double> shape_coefficients, std::vector<double> class_coefficients, double eta_class);
    std::vector<std::vector<std::vector<double>>> build_surface(std::vector<std::vector<double>> eps, std::vector<std::vector<double>> eta, std::string surface_id);
};


class NACAAirfoil 
{
private:
    bool is_used;
    std::string type;
    std::vector<double> C;
    std::vector<double> M;
    std::vector<double> P;
    std::vector<double> T;

public:
    NACAAirfoil();
    NACAAirfoil(std::vector<double> c, std::vector<double> m, std::vector<double> p, std::vector<double> t);

    bool get_is_used();
    std::vector<double> get_C();
    std::vector<double> get_M();
    std::vector<double> get_P();
    std::vector<double> get_T();

    void flip();
    double get_thickness(double eps_val, double t_val);
    double get_camber(double eps_val, double m_val, double p_val);
    double get_theta(double eps_val, double m_val, double p_val);
    std::vector<std::vector<std::vector<double>>> build_surface(std::vector<std::vector<double>> eps, std::vector<std::vector<double>> eta, std::string surface_id);
};


class Surface 
{
private:
    int Nx;
    int Ny;

    // airfoildata
    CSTAirfoil CST_AIRFOIL;
    NACAAirfoil NACA_AIRFOIL;

    double B;
    double D;
    std::vector<double> C;
    std::vector<double> X_LE;
    std::vector<double> Z_N;
    std::vector<double> DELTA_ALPHA_T;
    double NC1;
    double NC2;
    double NY1;
    double NY2;
    std::string ID;

    std::string x_distribution;
    std::string y_distribution;
    bool is_flipped;
    std::vector<std::vector<double>> eps;
    std::vector<std::vector<double>> eta;

public:
    Surface(int nx, int ny, double b, double d, std::vector<double> c, 
    std::vector<double> x_le, std::vector<double> z_n, std::vector<double> delta_alpha_t, 
    double nc1, double nc2, double ny1, double ny2, std::string identification, 
    CSTAirfoil cst_airfoil = CSTAirfoil{}, NACAAirfoil naca_airfoil = NACAAirfoil{});

    double get_B();
    double get_D();
    std::vector<double> get_C();
    std::string get_id();

    void flip_on_xz_plane();
    void change_distributions(std::string x_dist, std::string y_dist);
    std::vector<std::vector<std::vector<double>>> build_complete_surface(std::string surface_type);
    void set_distributions();
};


class Body 
{
private:
    std::string type;
    std::vector<Surface> surfaces;
    std::vector<std::vector<int>> associativity;
    bool is_flipped;

public:
    Body(std::string body_type);
    void clear();
    std::vector<Surface> get_surfaces();
    void set_to_general_solid();
    void set_to_revolution_solid();

    void add_surface_cst(int nx, int ny, double b, double d, std::vector<double> c, std::vector<double> z_te, 
    std::vector<double> r_le, std::vector<double> beta, std::vector<double> coefficients, 
    std::vector<double> x_le, std::vector<double> z_n, std::vector<double> delta_alpha_t,
    std::vector<double> n1, std::vector<double> n2, double nc1, double nc2, double ny1, double ny2, std::string identification);

    void add_surface_naca(std::vector<double> m, std::vector<double> p, std::vector<double> t, int nx, int ny, 
    double b, double d, std::vector<double> c, std::vector<double> x_le, std::vector<double> z_n, 
    std::vector<double> delta_alpha_t, double nc1, double nc2, double ny1, double ny2, std::string identification);

    void add_wing_surface_cst(int nx, int ny, double b, double d, std::vector<double> c, 
    std::vector<double> z_te_half, std::vector<double> r_le, std::vector<double> beta, 
    std::vector<double> x_le, std::vector<double> z_n, std::vector<double> delta_alpha_t);

    void add_general_surface_cst(int nx, int ny, double b, double d, std::vector<double> c, 
    std::vector<double> z_te_u, std::vector<double> z_te_l, std::vector<double> coefficients_u, std::vector<double> coefficients_l, 
    std::vector<double> x_le, std::vector<double> z_n, std::vector<double> delta_alpha_t, 
    std::vector<double> n1, std::vector<double> n2, double nc1, double nc2, double ny1, double ny2);

    void add_wing_surface_naca(int nx, int ny, double b, double d, std::vector<double> c, 
    std::vector<double> m, std::vector<double> p, std::vector<double> t,
    std::vector<double> x_le, std::vector<double> z_n, std::vector<double> delta_alpha_t);
    
    void mirror_body();
    void change_all_distributions(std::string x_dist, std::string y_dist);

    std::vector<std::vector<std::vector<std::vector<double>>>> get_paired_body_surfaces();
    std::vector<std::vector<std::vector<std::vector<double>>>> get_body_surfaces();
    std::vector<std::vector<int>> get_associativity();
    bool get_is_flipped();
};


//
#endif
