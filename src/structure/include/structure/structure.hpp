
#ifndef __STRUCTURE__
#define __STRUCTURE__

#include "common_aeroflex.hpp"
#include <ios>
#include <iostream>
#include <string>
#include <sstream>
#include <exception>
#include <vector>
#include <atomic>
#include "structure/model.hpp"
#include "tinyconfig.hpp"

namespace structure {

struct Settings {

  std::string Mesh_file_path = "../../../../examples/structure/Moment.txt";
  double Tolerance = 1e-6;
  int N_step = 10;
  double Damping = 1.0;

  std::vector<std::string> Solve_type_options = {"NONLINEAR", "LINEAR"};
  int Solve_type = 0;

  void set_solve_type(const std::string &solve_type) {
    if (!solve_type.compare("NONLINEAR"))
      this->Solve_type = 0;
    else if (!solve_type.compare("LINEAR"))
      this->Solve_type = 1;
  }

  std::string get_solve_type() {
    return Solve_type_options[Solve_type];
  }

  void import_config_file(tiny::config &config) {
      Mesh_file_path = config.get<std::string>("structure-io", "mesh_file");
      Tolerance = config.get<double>("structure-solver", "tolerance");
      N_step = config.get<int>("structure-solver", "n_steps");
      Damping = config.get<double>("structure-solver", "damping");
      set_solve_type(config.get<std::string>("structure-solver", "type"));
  }

  void export_config_file(tiny::config &config) {
    // Appending sections to config file
    config.sections.push_back("structure-io");
    config.sections.push_back("structure-solver");

    std::stringstream tolerance, damping;
    tolerance << std::scientific << this->Tolerance;
    damping << std::scientific << this->Damping;

    // Adding to config map
    // [structure-io]
    config.config["structure-io"]["mesh_file"] = Mesh_file_path;

    // [structure-solver]
    config.config["structure-solver"]["tolerance"] = tolerance.str();
    config.config["structure-solver"]["n_steps"] = N_step;
    config.config["structure-solver"]["damping"] = damping.str();
    config.config["structure-solver"]["type"] = get_solve_type();
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
   
  void solve() {

    // Resetting previous case
    reset();

    // Allocating memory for residuals
    // TODO: Setting a max_iter number as imput param!
    residuals.reserve(100);
    
    // FEM.set_Load_Vector_From_Vector(New_F);

    if (!settings.get_solve_type().compare("LINEAR"))
    {
      Solutions.push_back(FEM.get_Lin_Solve());
    }
    else if(!settings.get_solve_type().compare("NONLINEAR"))
    {
      Solutions.push_back(FEM.get_NonLin_Solve(settings.N_step, settings.Tolerance, settings.Damping));
    }

  }

private:
    void reset() {
      iters = 0;
      residuals.clear();
    }
};

}

#endif
