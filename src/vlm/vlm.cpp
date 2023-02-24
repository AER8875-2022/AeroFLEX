
#include "vlm.hpp"

/*
//Calling solver from CLASS VLM
VLM::VLM(SignalHandler signal_gui) : signal_gui(signal_gui) {}

void VLM::input(){

  //importing the model
  std::cout << "==>Loading simulation parameters...";
  [io, sim, solvP] = vlm::input::importConfigFile(argv[1]);
  std::cout << "\033[1;36mDone\033[0m" << std::endl;


  std::cout << "==>Loading mesh file...";
  auto mesh = vlm::input::importMeshFile(io);
  std::cout << "\033[1;36mDone\033[0m" << std::endl;


  //creating the model
  std::cout << "==>Initializing vortex lattice model...";
  object(mesh, sim, io);
  std::cout << "\033[1;36mDone\033[0m" << std::endl;

}

VLM::solve(solvP){

  std::cout << "==>Solving...\n";
  auto solver = vlm::solver::initializeSolver(vlm::input::solvP, object, 1.0);
  solver->solve(object);
  std::cout << "...\033[1;36mDone\033[0m" << std::endl;

  }

  */
