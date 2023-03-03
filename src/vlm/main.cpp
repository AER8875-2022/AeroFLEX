
#include "vlm/input.hpp"
#include "vlm/model.hpp"
#include "vlm/solver.hpp"
#include "vlm/utils.hpp"
#include <atomic>
#include <iostream>

int main(int argc, char **argv) {

  vlm::utils::printArtwork();

  // Verifying and handling input arguments
  if (argc < 2) {
    std::cerr << "\033[1;31m==>ERROR: expected config file as argument \033[0m"
              << std::endl;
    std::cerr << "\033[1;36m==>USAGE>> vlm.out path/to/config/file \033[0m"
              << std::endl;
    return 1;
  } else if (argc > 2) {
    std::cerr << "\033[1;33m==>WARNING: Extra arguments are ignored \033[0m"
              << std::endl;
    std::cerr << "==>\033[1;36mUSAGE>> vlm.out path/to/config/file\033[0m"
              << std::endl;
  }

  // #############################
  // LOADING SIMULATION PARAMETERS
  // #############################
  std::cout << "==>Loading simulation parameters...";
  auto [io, sim, solvP] = vlm::input::importConfigFile(argv[1]);
  std::cout << "\033[1;36mDone\033[0m" << std::endl;

  vlm::utils::printCaseInfo(sim, io, solvP);

  // #############################
  // IMPORTING MESH FILE
  // #############################
  std::cout << "==>Loading mesh file...";
  auto mesh = vlm::input::importMeshFile(io);
  std::cout << "\033[1;36mDone\033[0m" << std::endl;

  vlm::utils::printMeshInfo(mesh);

  // #############################
  // INITIALIZING VLM MODEL
  // #############################
  std::cout << "==>Initializing vortex lattice model...";
  vlm::model object;
  object.initialize(mesh, sim, io);
  std::cout << "\033[1;36mDone\033[0m" << std::endl;

  // #############################
  // INITIALIZING VISCOUS DATABASE
  // #############################
  database::table database;

  if (!solvP.type.compare("NONLINEAR")) {
    std::cout << "==>Initializing viscous database...";
    database.importFromFile(io.databaseFile, solvP);
  }

  std::cout << "\033[1;36mDone\033[0m" << std::endl;

  // #############################
  // CALLING SOLVER
  // #############################
  std::cout << "==>Solving...\n";

  std::atomic<int> iters = 0;
  std::vector<double> residuals;

  vlm::solver::base *solver;
  vlm::solver::linear::steady linear(iters, residuals);
  vlm::solver::nonlinear::steady nonlinear(iters, residuals);

  if (!solvP.type.compare("LINEAR")) {
    linear.initialize(solvP, object, database::table());
    solver = &linear;
  } else if (!solvP.type.compare("NONLINEAR")) {
    nonlinear.initialize(solvP, object, database);
    solver = &nonlinear;
  } else {
    std::cerr << "\033[1;31m==>ERROR: Unknown VLM solver \033[0m" << std::endl;
    return 1;
  }

  solver->solve(object);
  std::cout << "...\033[1;36mDone\033[0m" << std::endl;
  std::cout << "\n\033[1;32mSimulation done! \033[0m" << std::endl;

  return 0;
}
