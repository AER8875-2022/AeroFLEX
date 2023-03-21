
#include "vlm/vlm.hpp"
#include <atomic>
#include <iostream>

int main(int argc, char **argv) {

  vlm::info::printArtwork();

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
  vlm::Settings settings;
  tiny::config config;
  config.read(argv[1]);
  settings.import_config_file(config);
  std::cout << "\033[1;36mDone\033[0m" << std::endl;

  vlm::info::printCaseInfo(settings);

  // #############################
  // IMPORTING MESH FILE
  // #############################
  std::cout << "==>Loading mesh file...";
  auto mesh = vlm::input::importMeshFile(settings.io);
  std::cout << "\033[1;36mDone\033[0m" << std::endl;

  vlm::info::printMeshInfo(mesh);

  // #############################
  // INITIALIZING VLM MODEL
  // #############################
  std::cout << "==>Initializing vortex lattice model...";
  vlm::model object;
  object.initialize(mesh, settings);
  std::cout << "\033[1;36mDone\033[0m" << std::endl;

  // #############################
  // INITIALIZING VISCOUS DATABASE
  // #############################
  database::table database;

  if (!settings.solver.get_type().compare("NONLINEAR")) {
    std::cout << "\n";
    std::cout << "==>Initializing viscous database...";
    database.importAirfoils(settings.io.databaseFile);
    database.importLocations(settings.io.locationFile);
    if (!database.check()) {
      std::cout << "\033[1;31m==>ERROR: One or more airfoils were not found! - "
                   "Aborting \033[0m"
                << std::endl;
    }
    std::cout << "\033[1;36mDone\033[0m" << std::endl;
  }

  std::cout << std::endl;

  // #############################
  // CALLING SOLVER
  // #############################
  std::cout << "==>Solving...\n";

  std::atomic<int> iters = 0;
  std::vector<double> residuals;
  GUIHandler gui; // Empty GUI Handler

  vlm::solver::base *solver;
  vlm::solver::linear::steady linear(settings.solver, iters, residuals, gui);
  vlm::solver::nonlinear::steady nonlinear(settings.solver, iters, residuals,
                                           gui);

  if (!settings.solver.get_type().compare("LINEAR")) {
    linear.initialize(object, database::table());
    solver = &linear;
  } else if (!settings.solver.get_type().compare("NONLINEAR")) {
    nonlinear.initialize(object, database);
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
