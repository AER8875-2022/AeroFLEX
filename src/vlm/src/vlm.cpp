
#include "vlm/vlm.hpp"
#include "common_aeroflex.hpp"
#include "vlm/info.hpp"
#include "vlm/input.hpp"

using namespace vlm;

VLM::VLM(GUIHandler &gui) : gui(gui) {};

void VLM::input() {
  if (!is_initialized) {
    auto mesh = vlm::input::importMeshFile(settings.io);
    object.initialize(mesh, settings);
    is_initialized = true;
  }
};

void VLM::solve() {

  gui.msg.push("[VLM] Solving VLM with solver: " + settings.solver.get_type());

  solver::base *solver;
  solver::linear::steady linear(settings.solver, iters, residuals, gui);
  solver::nonlinear::steady nonlinear(settings.solver, iters, residuals, gui);

  // Clear previous solution
  // reinitialize();

  if (!settings.solver.get_type().compare("LINEAR")) {
    linear.initialize(object, database::table());
    solver = &linear;
  } else if (!settings.solver.get_type().compare("NONLINEAR")) {
    nonlinear.initialize(object, database);
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
