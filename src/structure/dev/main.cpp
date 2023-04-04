#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <structure/proprietes_sections.hpp>
#include <structure/element.hpp>
#include <structure/model.hpp>
#include <structure/structure.hpp>
#include <common_aeroflex.hpp>

using namespace structure;

int main(int argc, char **argv){

    // Reading config file
    tiny::config config;
    config.read(argv[1]);

    // Importing settings
    Settings settings;
    settings.import_config_file(config);
    
    // Dummy GUI variables
    GUIHandler gui;
    std::atomic<int> iters = 0;
    std::vector<double> residuals;
    std::cout<<settings.Force_follower+0<<std::endl;

    MODEL M1(settings.Mesh_file_path, gui, iters, residuals, settings.Force_follower);
    Eigen::VectorXd dep;

    if(!settings.get_solve_type().compare("NONLINEAR"))
    {   
        int Load_step = settings.N_step;
        double       tol = settings.Tolerance;
        double      amor = settings.Damping;

        dep = M1.get_NonLin_Solve(Load_step, tol, amor);
    }
    else
    {
        dep = M1.get_Lin_Solve();

    }

    std::string arg1 = settings.Mesh_file_path;
    std::cout<<"Test"<<std::endl;

    int L_string = arg1.length();
    
    std::string end = arg1.substr (L_string-9,L_string);
    std::string File_name = "Result"+end;

    std::cout<<"Résultats écrit dans le fichier : "<< File_name<<std::endl; 
    
    double L_elem = 300./M1.Nbr_Element;
    std::ofstream out;
    out.open(File_name);
    out << "x y z rx ry rz\n";
    
    for (int i=0; i< (dep.size()/6); ++i) 
    {   
        double x_ini = 0.; 
        double y_ini = 0.; 
        double z_ini = 0.;

        for (auto const& [id_user, id_code] : M1.indexation_switch) {
            if (id_code == i) {
                Eigen::Vector3d Coord = M1.Grid_MAP[id_user];
                x_ini = Coord(0);
                y_ini = Coord(1);
                z_ini = Coord(2);
                break;
            }
        }
        
        double x     =  x_ini + dep(i*6);
        double y     =  y_ini + dep(i*6 + 1);
        double z     =  z_ini + dep(i*6 + 2);
        double rx    =  dep(i*6 + 3);
        double ry    =  dep(i*6 + 4);
        double rz    =  dep(i*6 + 5);

        out << std::setprecision(16)<<x << " " << y <<" "<< z << " " << rx <<" "<< ry <<" "<< rz <<"\n";
    }
    out.close();



    return 0;
}

