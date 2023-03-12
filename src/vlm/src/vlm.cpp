
#include "vlm/vlm.hpp"
#include "common_aeroflex.hpp"

using namespace vlm;

VLM::VLM(GUIHandler &gui) : gui(gui) {};

void VLM::input() {
  if (!is_initialized) {
    auto mesh = vlm::input::importMeshFile(settings.io);
    model object;
    object.initialize(mesh, settings.sim, settings.io);
    is_initialized = true;
  }
};

void VLM::solve() {
  solver::base *solver;
  solver::linear::steady linear(iters, residuals, gui);
  solver::nonlinear::steady nonlinear(iters, residuals, gui);

  // Clear previous solution
  reinitialize();

  if (!settings.solver.type.compare("LINEAR")) {
    linear.initialize(settings.solver, object, database::table());
    solver = &linear;
  } else if (!settings.solver.type.compare("NONLINEAR")) {
    nonlinear.initialize(settings.solver, object, database);
    solver = &nonlinear;
  } else {
    return;
  }

  solver->solve(object);
};

void VLM::reinitialize() {
  // Reset model
  object.resetWake();
  object.clear();
  // Reset iterations
  residuals.clear();
  iters = 0;
}
