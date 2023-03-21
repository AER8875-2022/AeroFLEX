
#include "vlm/model.hpp"
#include <iostream>
#include <unordered_set>

using namespace vlm;
using namespace Eigen;

void model::initialize(const input::meshData &mesh, const Settings &settings) {
  // Initializing inputs
  this->sim = settings.sim;
  this->io = settings.io;
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
    vortex.updateGamma(1.0);
  }
  // TODO: Update Doublet solution
  for (auto &wingStation : wingStations) {
    wingStation.resetLocalStream(sim);
  }
}

void model::build(const input::meshData &mesh) {

  // Clearing previous mesh
  nodes.clear();
  vortexRings.clear();
  doubletPanels.clear();
  wingStations.clear();
  wings.clear();
  patches.clear();

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

double model::get_cy() const { return cy; }

Vector3d model::get_cm() const { return cm; }
