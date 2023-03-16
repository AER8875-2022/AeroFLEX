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

    MODEL M1(settings.Mesh_file_path, gui, iters, residuals);
    Eigen::VectorXd dep;

    if(!settings.Solve_type.compare("NLS"))
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
    std::cout<<L_elem<<std::endl;
    std::ofstream out;
    out.open(File_name);
    out << "x y z rx ry rz\n";
    
    for (int i=0; i< (dep.size()/6); ++i) 
    {
        double x     =  L_elem*i+ dep(i*6);
        double y     =  dep(i*6 + 1);
        double z     =  dep(i*6 + 2);
        double rx    =  dep(i*6 + 3);
        double ry    =  dep(i*6 + 4);
        double rz    =  dep(i*6 + 5);

        out << std::setprecision(16)<<x << " " << y <<" "<< z << " " << rx <<" "<< ry <<" "<< rz <<"\n";
    }
    out.close();



    return 0;
}

