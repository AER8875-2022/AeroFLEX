
#include "vlm/solver.hpp"
#include "database/database.hpp"
#include "vlm/output.hpp"
#include <any>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <thread>

using namespace vlm;
using namespace solver;
using namespace Eigen;

linear::steady::steady(std::atomic<int> &iter, std::vector<double> &residuals,
                       GUIHandler &gui)
    : iter(iter), residuals(residuals), gui(gui) {}

void linear::steady::initialize(const input::solverParam &solvP,
                                const model &object, const database::table &) {
  this->solvP = solvP;
  this->lhs =
      MatrixXd::Zero(object.vortexRings.size() + object.doubletPanels.size(),
                     object.vortexRings.size() + object.doubletPanels.size());
  this->rhs =
      VectorXd::Zero(object.vortexRings.size() + object.doubletPanels.size());
  this->DoubletMatrixINF =
      MatrixXd::Zero(object.doubletPanels.size(), object.doubletPanels.size());
}

void linear::steady::solve(model &object) {
  // Building wake
  auto freeStream = object.sim.freeStream();
  freeStream.normalize();
  object.initializeWake(100.0 * object.sim.cref);

  // Initializing linear solver
  BiCGSTAB<MatrixXd> system;

  // Building left hand side
  buildLHS(object);
  system.compute(lhs);

  // Building linear system
  buildRHS(object);
  // Solving
  VectorXd gamma = system.solve(rhs);

  // Saving solution to model for further uses
  saveSolution(object, gamma);

  // Computing forces
  computeForces(object);
  computeInducedDrag(object);

  // Exporting solution
  exportSolution(object);

  std::cout << "\t Residual = " << system.error() << std::endl;
  residuals.push_back(system.error());
  iter++;
}

void linear::steady::saveSolution(model &object, const VectorXd &gamma) {
#pragma omp parallel
  {
    // Saving gamma to vortex rings
#pragma omp for
    for (int id = 0; id < object.vortexRings.size(); id++) {
      auto &vortex = object.vortexRings[id];

      vortex.updateGamma(gamma[vortex.get_globalIndex()]);
    }
// Saving gamma to wake panels
#pragma omp for
    for (int id = 0; id < object.wakePanels.size(); id++) {
      auto &wake = object.wakePanels[id];

      wake.updateGamma(gamma[wake.get_globalIndex()]);
    }
  }
}

void linear::steady::computeInducedDrag(model &object) {
  double drag = 0.0;
  // Looping over wake panels
#pragma omp parallel for reduction(- : drag)
  for (int influencedID = 0; influencedID < object.wakePanels.size();
       influencedID++) {
    auto &influencedWake = object.wakePanels[influencedID];

    Vector3d v = Vector3d::Zero();
    for (int influencingID = 0; influencingID < object.wakePanels.size();
         influencingID++) {
      auto &influencingWake = object.wakePanels[influencingID];
      v += influencingWake.influence(influencedWake.get_collocationPoint());
    }
    Vector3d dl = influencedWake.leadingEdgeDl();
    drag -= 0.5 * object.sim.rho * influencedWake.get_gamma() *
            v.dot(influencedWake.get_normal() * dl.norm());
  }
  object.cd += drag / (object.sim.dynamicPressure() * object.sim.sref);
}

void linear::steady::computeForces(model &object) {
  object.cl = 0.0;
  object.cd = 0.0;
  object.cm = Vector3d::Zero();
  // Computing forces at each wing station
#pragma omp parallel for
  for (int stationID = 0; stationID < object.wingStations.size(); stationID++) {
    auto &station = object.wingStations[stationID];
    station.computeForces(object.sim);
  }
  // Updating global forces
  for (int wingID = 0; wingID < object.wings.size(); wingID++) {
    auto &wing = object.wings[wingID];

    wing.computeForces(object.sim);
    object.cl += wing.get_cl() * wing.get_area() / object.sim.sref;
    object.cm += wing.get_cm() * wing.get_area() / object.sim.sref;
  }
}

void linear::steady::exportSolution(const model &object) {
  output::exportForces(object, 0);
  output::exportSurfacesVTU(object, 0);
  output::exportWakeVTU(object, 0);
}

void linear::steady::buildLHS(const model &object) {
  // Panel scan loop
  const int size_Vortex = object.vortexRings.size();
  // influence de tous les panneaux sur les VortexRings
#pragma omp parallel
  {
#pragma omp for
    for (int influencedVortexID = 0;
         influencedVortexID < object.vortexRings.size(); influencedVortexID++) {
      auto &influencedVortex = object.vortexRings[influencedVortexID];

      // influence Vortex -> Vortex
      for (int influencingVortexID = 0;
           influencingVortexID < object.vortexRings.size();
           influencingVortexID++) {
        auto &influencingVortex = object.vortexRings[influencingVortexID];

        auto v1 = influencingVortex.influence(
            influencedVortex.get_collocationPoint());
        lhs(influencedVortex.get_globalIndex(),
            influencingVortex.get_globalIndex()) =
            v1.dot(influencedVortex.get_normal());
      }

      // influence Doublets -> Vortex
      for (int influencingDoubletsID = 0;
           influencingDoubletsID < object.doubletPanels.size();
           influencingDoubletsID++) {
        auto &influencingDoublets = object.doubletPanels[influencingDoubletsID];

        auto v2 = influencingDoublets.influence(
            influencedVortex
                .get_collocationPoint()); // modifier la fonction
                                          // doublets.influence dans panel.cpp
        lhs(influencedVortex.get_globalIndex(),
            size_Vortex + influencingDoublets.get_globalIndex()) =
            v2.dot(influencedVortex.get_normal());
      }
    }
    // influence de tous les panneaux sur les doublets
#pragma omp for
    for (int influencedDoubletsID = 0;
         influencedDoubletsID < object.doubletPanels.size();
         influencedDoubletsID++) {
      auto &influencedDoublets = object.doubletPanels[influencedDoubletsID];

      // influence Vortex -> Doublets
      for (int influencingVortexID = 0;
           influencingVortexID < object.vortexRings.size();
           influencingVortexID++) {
        auto &influencingVortex = object.vortexRings[influencingVortexID];

        auto v3 = influencingVortex.influence(
            influencedDoublets
                .get_center()); // modifier la fonction doublets.influence dans
                                // panel.cpp

        lhs(size_Vortex + influencedDoublets.get_globalIndex(),
            influencingVortex.get_globalIndex()) =
            v3.dot(influencedDoublets.get_normal());
      }

      // influence Doublets -> Doublets
      for (int influencingDoubletsID = 0;
           influencingDoubletsID < object.doubletPanels.size();
           influencingDoubletsID++) {
        auto &influencingDoublets = object.doubletPanels[influencingDoubletsID];

        auto v4 = influencingDoublets.influence(
            influencedDoublets
                .get_center()); // modifier la fonction doublets.influence dans
                                // panel.cpp
        lhs(size_Vortex + influencedDoublets.get_globalIndex(),
            size_Vortex + influencingDoublets.get_globalIndex()) =
            v4.dot(influencedDoublets.get_normal());
        DoubletMatrixINF(influencedDoublets.get_globalIndex(),
                         influencingDoublets.get_globalIndex()) =
            v4.dot(influencedDoublets.get_normal());
      }
    }

// Adding contribution of wake
#pragma omp for
    for (int influencedVortexID = 0;
         influencedVortexID < object.vortexRings.size(); influencedVortexID++) {
      auto &influencedVortex = object.vortexRings[influencedVortexID];

      for (int influencingWakeID = 0;
           influencingWakeID < object.wakePanels.size(); influencingWakeID++) {
        auto &influencingWake = object.wakePanels[influencingWakeID];

        auto v =
            influencingWake.influence(influencedVortex.get_collocationPoint());
        lhs(influencedVortex.get_globalIndex(),
            influencingWake.get_globalIndex()) +=
            v.dot(influencedVortex.get_normal());
      }
    }
#pragma omp for
    for (int influencedDoubletsID = 0;
         influencedDoubletsID < object.doubletPanels.size();
         influencedDoubletsID++) {
      auto &influencedDoublets = object.doubletPanels[influencedDoubletsID];

      for (int influencingWakeID = 0;
           influencingWakeID < object.wakePanels.size(); influencingWakeID++) {
        auto &influencingWake = object.wakePanels[influencingWakeID];

        auto v = influencingWake.influence(influencedDoublets.get_center());
        lhs(size_Vortex + influencedDoublets.get_globalIndex(),
            influencingWake.get_globalIndex()) +=
            v.dot(influencedDoublets.get_normal());
      }
    }
  }
}

void linear::steady::buildRHS(const model &object) {
  VectorXd rhs_VLM = VectorXd::Zero(object.vortexRings.size());
  VectorXd sources = VectorXd::Zero(object.doubletPanels.size());
#pragma omp parallel
  {
#pragma omp for
    for (int vortexID = 0; vortexID < object.vortexRings.size(); vortexID++) {
      auto &vortex = object.vortexRings[vortexID];

      rhs_VLM(vortex.get_globalIndex()) =
          -object.sim.freeStream().dot(vortex.get_normal());
    }

    // Building sources vector
#pragma omp for
    for (int doubletID = 0.0; doubletID < object.doubletPanels.size();
         doubletID++) {
      auto &doubs = object.doubletPanels[doubletID];

      sources(doubs.get_globalIndex()) =
          -object.sim.freeStream().dot(doubs.get_normal());
    }
  }

  VectorXd rhs_PanMethod = DoubletMatrixINF * sources;
  // concatenate rhs vectors
  rhs << rhs_VLM, rhs_PanMethod;
}
// -------------------------------------------

// -------------------------------------------

nonlinear::steady::steady(std::atomic<int> &iter,
                          std::vector<double> &residuals, GUIHandler &gui)
    : linear::steady(iter, residuals, gui) {}

void nonlinear::steady::initialize(const input::solverParam &solvP,
                                   const model &object,
                                   const database::table &database) {
  this->solvP = solvP;
  this->lhs =
      MatrixXd::Zero(object.vortexRings.size() + object.doubletPanels.size(),
                     object.vortexRings.size() + object.doubletPanels.size());
  this->rhs =
      VectorXd::Zero(object.vortexRings.size() + object.doubletPanels.size());
  this->DoubletMatrixINF =
      MatrixXd::Zero(object.doubletPanels.size(), object.doubletPanels.size());
  this->database = database;
}

void nonlinear::steady::buildRHS(const model &object) {
  VectorXd rhs_VLM = VectorXd::Zero(object.vortexRings.size());
  VectorXd sources = VectorXd::Zero(object.doubletPanels.size());
#pragma omp parallel
  {
#pragma omp for
    for (int vortexID = 0; vortexID < object.vortexRings.size(); vortexID++) {
      auto &vortex = object.vortexRings[vortexID];

      rhs_VLM(vortex.get_globalIndex()) =
          -object.sim.freeStream(vortex.local_aoa).dot(vortex.get_normal());
    }
    // Building sources vector
#pragma omp for
    for (int doubletID = 0; doubletID < object.doubletPanels.size();
         doubletID++) {
      auto &doubs = object.doubletPanels[doubletID];

      sources(doubs.get_globalIndex()) =
          -object.sim.freeStream().dot(doubs.get_normal());
    }
  }

  VectorXd rhs_PanMethod = DoubletMatrixINF * sources;
  // concatenate rhs vectors
  rhs << rhs_VLM, rhs_PanMethod;
}

void nonlinear::steady::solve(model &object) {
  // Building wake
  auto freeStream = object.sim.freeStream();
  freeStream.normalize();
  object.initializeWake(100.0 * object.sim.cref);

  // Initializing linear solver
  BiCGSTAB<MatrixXd> system;

  // Building LHS
  buildLHS(object);
  system.compute(lhs);

  // Initializing iteration
  double residual;

  // Main solving loop
  do {
    while (gui.signal.pause)
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    // Step 1 : Solving VLM
    buildRHS(object);
    VectorXd gamma = system.solve(rhs);
    // Saving solution to model for further uses
    saveSolution(object, gamma);
    // Computing 3D lift coefficient
    iterateLift(object);
    // Step 2: One iteration of aoa correction
    residual = iterate(object);

    // gui.msg.push("Iteration " + std::to_string(iter));
    // gui.msg.push("\t Residual = " + std::to_string(residual));

    std::cout << "Iteration " << iter << std::endl;
    std::cout << "\t Residual = " << residual << std::endl;

    residuals.push_back(residual);
    iter++;
  } while ((residual > solvP.tolerance) && (iter <= solvP.max_iter) &&
           (!gui.signal.stop));
  // Compute viscous forces
  computeForces(object);
  // Computing 3D effects on drag
  computeInducedDrag(object);
  // Exporting solution
  exportSolution(object);
}

double nonlinear::steady::iterate(model &object) {
  double residual = 0.0;

  for (int wingID = 0; wingID < object.wings.size(); wingID++) {
    auto &wing = object.wings[wingID];

#pragma omp parallel for reduction(+ : residual)
    for (int i = 0; i < wing.get_stationIDs().size(); i++) {
      auto stationID = wing.get_stationIDs()[i];

      double root =
          object.wingStations[wing.get_stationIDs().front()].get_spanLoc();

      auto &station = object.wingStations[stationID];
      double spanLoc = (station.spanLoc - root) / wing.get_span();

      // Step 3 : Effective angle of attack
      double cl_inv = station.cl;
      double aoa_eff =
          cl_inv / VLM_CL_ALPHA_DEG - station.local_aoa + object.sim.aoa;

      // Step 4 : Viscous lift interpolation
      double cl_visc = database.cl(aoa_eff, wing.get_globalIndex(), spanLoc);

      // Print error if extrapolation is detected
      if (std::abs(cl_visc) < 1e-15) {
        std::cerr << "\033[1;31mERROR: Extrapolation detected\033[0m"
                  << std::endl;
        exit(1);
      };

      // Step 5 : Applying aoa correction
      double dalpha = solvP.relaxation * (cl_visc - cl_inv) / VLM_CL_ALPHA_DEG;
      station.updateLocalAoa(dalpha);

      residual += (cl_visc - cl_inv) / cl_visc * (cl_visc - cl_inv) / cl_visc;
    }
  }

  return (std::sqrt(residual));
}

void nonlinear::steady::iterateLift(model &object) {
  object.cl = 0.0;
#pragma omp parallel for
  for (int stationID = 0; stationID < object.wingStations.size(); stationID++) {
    auto &station = object.wingStations[stationID];

    station.computeForces(object.sim);
  }
}

void nonlinear::steady::computeForces(model &object) {
  object.cl = 0.0;
  object.cd = 0.0;
  object.cm = Vector3d::Zero();
  // Interpolating viscous forces at each wing station
#pragma omp parallel for
  for (int wingID = 0; wingID < object.wings.size(); wingID++) {
    auto &wing = object.wings[wingID];

    double root =
        object.wingStations[wing.get_stationIDs().front()].get_spanLoc();

    for (int i = 0; i < wing.get_stationIDs().size(); i++) {
      auto stationID = wing.get_stationIDs()[i];
      auto &station = object.wingStations[stationID];

      double spanLoc = (station.spanLoc - root) / wing.get_span();
      // Computing effective aoa to interpolate
      double aoa_eff =
          station.cl / VLM_CL_ALPHA_DEG - station.local_aoa + object.sim.aoa;

      // Interpolating all coefficients at each station
      auto [cl, cd, cmy] =
          database.coefficients(aoa_eff, wing.get_globalIndex(), spanLoc);
      // Lever used to transfer to 2D moment to 3D moment at specified origin
      double lever = object.sim.origin()(0) - station.forceActingPoint()(0);
      // Updating station's force coefficients
      station.cl = cl;
      station.cd = cd;
      station.cm(1) = cmy + lever / station.chord * cl;
    }
    // Updating global forces
    wing.computeForces(object.sim);
    object.cl += wing.get_cl() * wing.get_area() / object.sim.sref;
    object.cd += wing.get_cd() * wing.get_area() / object.sim.sref;
    object.cm += wing.get_cm() * wing.get_area() / object.sim.sref;
  }
}
