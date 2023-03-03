
#include "vlm.hpp"
#include "common_aeroflex.hpp"

using namespace vlm;

VLM::VLM(GUIHandler &gui) : gui(gui){};

void VLM::input() {
  if (!is_initialized) {
    auto mesh = vlm::input::importMeshFile(data.io);
    model object;
    object.initialize(mesh, data.sim, data.io);
    is_initialized = true;
  }
};

void VLM::solve() {
  solver::base *solver;
  solver::linear::steady linear(iters, residuals);
  solver::nonlinear::steady nonlinear(iters, residuals);

  if (!data.solver.type.compare("LINEAR")) {
    linear.initialize(data.solver, object, database::table());
    solver = &linear;
  } else if (!data.solver.type.compare("NONLINEAR")) {
    nonlinear.initialize(data.solver, object, database);
    solver = &nonlinear;
  } else {
    return;
  }

  solver->solve(object);
};
