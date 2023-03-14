
#ifndef __STRUCTURE__
#define __STRUCTURE__

#include "common_aeroflex.hpp"
#include <vector>
#include <atomic>
#include "structure/model.hpp"

namespace structure {

struct Settings {
  // TODO  Tol, n_step, amor --> Done
  std::string Mesh_file_path;
  std::string Solve_type;
  double Tolerance;
  double N_step;
  double Damping; 
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
  void input() // TODO --> read_mesh
  {
    FEM.read_data_file(settings.Mesh_file_path);
    FEM.set_K_global();
    FEM.set_load_vector();        
    FEM.set_K_Final_sparse();
  }
   
  void solve(Eigen::VectorXd New_F){
    
    //FEM.set_FullLoadVector(New_F);

    if (settings.Solve_type == "Linear")
    {
      Solutions.push_back(FEM.get_LinSolve());
    }
    else if(settings.Solve_type == "Non-Linear")
    {
      Solutions.push_back(FEM.get_NonLinSolve(settings.N_step, settings.Tolerance, settings.Damping));
    }


  } // TODO --> set_loads , Lin_solve , NonLin_solve
};

}

#endif
