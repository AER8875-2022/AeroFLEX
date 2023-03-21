
#include "vlm/vlm.hpp"
#include "common_aeroflex.hpp"
#include "vlm/info.hpp"
#include "vlm/input.hpp"

using namespace vlm;

VLM::VLM(GUIHandler &gui) : gui(gui) {};

void VLM::initialize() {
  if (!is_initialized) {
    auto mesh = vlm::input::importMeshFile(settings.io);
    object.initialize(mesh, settings);
    is_initialized = true;
  }
};

void VLM::solve() {

  gui.msg.push("[VLM] Solving VLM case " + settings.io.baseName + " with solver " + settings.solver.get_type());

  // Check for validity of database
  if (!settings.solver.get_type().compare("NONLINEAR") && !database.check()) {
    gui.msg.push("[VLM] ERROR: One or more airfoils have not been found - Aborting");
    return;
  }

  // Clear previous solution
  reinitialize();

  // Allocating memory for residuals
  residuals.reserve(settings.solver.max_iter);

  solver::base *solver;
  solver::linear::steady linear(settings.solver, iters, residuals, gui);
  solver::nonlinear::steady nonlinear(settings.solver, iters, residuals, gui);

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

  gui.msg.push("[VLM] Simulation completed in " + std::to_string(iters) + " iterations");
  gui.msg.push("[VLM] Exported results to " + settings.io.outDir);
};

void VLM::reinitialize() {
  // Reset model
  object.resetWake();
  object.clear();
  object.updateGeometry(object.nodes);
  // Reset iterations
  residuals.clear();
  iters = 0;
}
