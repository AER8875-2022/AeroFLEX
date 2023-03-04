
#include "vlm/model.hpp"
#include <iostream>
#include <omp.h>
#include <unordered_set>

using namespace vlm;
using namespace Eigen;

void model::initialize(const input::meshData &mesh, const input::simParam &sim,
                       const input::ioParam &io) {
  // Initializing inputs
  this->sim = sim;
  this->io = io;
  // Allocating memory for each objects
  nodes.reserve(mesh.nodes.size());
  vortexRings.reserve(mesh.vortexIDs.size());
  doubletPanels.reserve(mesh.doubletIDs.size());
  wingStations.reserve(mesh.stationIDs.size());
  wings.reserve(mesh.wingIDs.size());
  patches.reserve(mesh.patchIDs.size());
  wakeNodes.reserve(2 * mesh.stationIDs.size());
  wakePanels.reserve(mesh.stationIDs.size());
  // Creating model from mesh object
  initializeMesh(mesh);
}

void model::initializeWake(const double wakeLength) {
  for (auto &station : wingStations) {
    station.generateWake(wakeLength, sim, wakeNodes, wakePanels);
  }
}

void model::initializeMesh(const input::meshData &mesh) {
  build(mesh);
  for (auto &wing : wings) {
    wing.initialize(sim);
  }
  for (auto &patch : patches) {
    patch.initialize(sim);
  }
}

void model::updateGeometry(const std::vector<Vector3d> &nodes) {
  this->nodes = nodes;
  for (auto &wing : wings) {
    wing.updateGeometry();
  }
  for (auto &patch : patches) {
    patch.updateGeometry();
  }
}

void model::resetWake() {
  wakePanels.clear();
  wakeNodes.clear();
}

void model::clear() {
  for (auto &vortex : vortexRings) {
    vortex.updateGamma(0.0);
  }
  // TODO: Update Doublet solution
  for (auto &wingStation : wingStations) {
    wingStation.resetLocalAoa(sim);
  }
}

void model::build(const input::meshData &mesh) {
  // Building nodes
  for (auto &[key, node] : mesh.nodes) {
    nodes.push_back(node);
  }
  // Building vortex rings
  for (auto &[key, ids] : mesh.vortexIDs) {
    vortexRings.push_back(element::vortexRing(key, ids, nodes, 1.0));
  }
  // Building doublet panels
  for (auto &[key, ids] : mesh.doubletIDs) {
    doubletPanels.push_back(element::doubletPanel(key, ids, nodes, 1.0));
  }
  // Building stations
  for (auto &[key, ids] : mesh.stationIDs) {
    wingStations.push_back(surface::wingStation(key, ids, vortexRings));
  }
  // Building wings
  for (auto &[key, ids] : mesh.wingIDs) {
    wings.push_back(surface::wing(key, ids, wingStations));
  }
  // Building patches
  for (auto &[key, ids] : mesh.patchIDs) {
    patches.push_back(surface::patch(key, ids, doubletPanels));
  }
}

double model::get_cl() const { return cl; }

double model::get_cd() const { return cd; }

Vector3d model::get_cm() const { return cm; }

Matrix<double, 6, 1>
model::forces_to_inertial_frame(const int stationID) const {
  auto &station = wingStations[stationID];
  Vector3d lift = station.get_cl() * sim.liftAxis();
  Vector3d drag = station.get_cd() * sim.streamAxis();
  Vector3d force = (lift + drag) * sim.dynamicPressure() * station.get_area();
  Vector3d moment = station.get_cm() * sim.dynamicPressure() *
                    station.get_area() * station.get_chord();
  return {force(0), force(1), force(2), moment(0), moment(1), moment(2)};
}
