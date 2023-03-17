
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

  gui.msg.push("[VLM] Solving VLM with solver " + settings.solver.get_type());
  gui.msg.push("[VLM] with base name " + settings.io.baseName);

  solver::base *solver;
  solver::linear::steady linear(settings.solver, iters, residuals, gui);
  solver::nonlinear::steady nonlinear(settings.solver, iters, residuals, gui);

  // Clear previous solution
  reinitialize();

  // Allocating memory for residuals
  residuals.reserve(settings.solver.max_iter);

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
