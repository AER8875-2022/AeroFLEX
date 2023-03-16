
#ifndef __STRUCTURE__
#define __STRUCTURE__

#include "common_aeroflex.hpp"
#include <iostream>
#include <string>
#include <exception>
#include <vector>
#include <atomic>
#include "structure/model.hpp"
#include "tinyconfig.hpp"

namespace structure {

struct Settings {
  std::string Mesh_file_path;
  std::string Solve_type;
  double Tolerance;
  int N_step;
  double Damping;

  void import_config_file(tiny::config &config) {

    try {

      Mesh_file_path = config.get<std::string>("structure-io", "mesh_file");
      Solve_type = config.get<std::string>("structure-solver", "type");
      Tolerance = config.get<double>("structure-solver", "tolerance");
      N_step = config.get<int>("structure-solver", "n_steps");
      Damping = config.get<double>("structure-solver", "damping");

    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
    }

  }

};

class Structure {

public:
  Settings settings;
  std::vector<Eigen::VectorXd> Solutions;
  std::vector<double> residuals;
  std::atomic<int> iters = 0;
  GUIHandler &gui;
  MODEL FEM;

public:
  Structure(GUIHandler &gui): gui(gui), FEM(gui, iters, residuals) {}
  void input()
  {
    FEM.read_data_file(settings.Mesh_file_path);
    FEM.set_K_global();
    FEM.set_Load_Vector_From_Load_Objects();      
    FEM.set_K_Final_sparse();
  }
   
  void solve(Eigen::VectorXd New_F){
    
    FEM.set_Load_Vector_From_Vector(New_F);

    if (settings.Solve_type == "Linear")
    {
      Solutions.push_back(FEM.get_Lin_Solve());
    }
    else if(settings.Solve_type == "Non-Linear")
    {
      Solutions.push_back(FEM.get_NonLin_Solve(settings.N_step, settings.Tolerance, settings.Damping));
    }


  }
};

}

#endif
