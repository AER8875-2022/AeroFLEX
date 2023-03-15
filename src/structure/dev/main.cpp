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
    
    std::cout << "You have entered " << argc << " arguments:" << "\n";

    // Dummy GUI variables
    GUIHandler gui;
    std::atomic<int> iters = 0;
    std::vector<double> residuals;

    
    MODEL M1(argv[1], gui, iters, residuals);
    Eigen::VectorXd dep;
    std::string S ("NLS");
 
    if(!S.compare(argv[2]))    //     /home/olivier/Structure-Dev/examples/Moment.txt NLS 10 1E-10 0.5
    {   
        int Load_step = std::stoi(argv[3]);
        double       tol = std::stod(argv[4]);
        double      amor = std::stod(argv[5]);

        dep = M1.get_NonLin_Solve(Load_step, tol, amor);
    }
    else
    {
        dep = M1.get_Lin_Solve();

    }

    std::string arg1 = argv[1];
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

