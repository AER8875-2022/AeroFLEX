
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
    station.generateWake(wakeLength, nodes, vortexRings, sim, wakeNodes,
                         wakePanels);
  }
}

void model::initializeMesh(const input::meshData &mesh) {
  build(mesh);
  for (auto &wing : wings) {
    wing.initialize(nodes, wingStations, vortexRings, sim);
  }
  for (auto &patch : patches) {
    patch.initialize(nodes, doubletPanels, sim);
    patch.ScanNeighbor(doubletPanels);
  }
}

void model::updateGeometry(const std::vector<Vector3d> &nodes) {
  this->nodes = nodes;
  for (auto &wing : wings) {
    wing.updateGeometry(nodes, wingStations, vortexRings);
  }
  for (auto &patch : patches) {
    patch.updateGeometry(nodes, doubletPanels);
  }
}

void model::resetWake() {
  wakePanels.clear();
  wakeNodes.clear();
}

void model::build(const input::meshData &mesh) {
  // Building nodes
  for (auto &[key, node] : mesh.nodes) {
    nodes.push_back(node);
  }
  // Building wings
  for (auto &[key, ids] : mesh.wingIDs) {
    wings.push_back(surface::wing(key, ids));
  }
  // Building patches
  for (auto &[key, ids] : mesh.patchIDs) {
    patches.push_back(surface::patch(key, ids));
  }
  // Building stations
  for (auto &[key, ids] : mesh.stationIDs) {
    wingStations.push_back(surface::wingStation(key, ids));
  }
  // Building vortex rings
  for (auto &[key, ids] : mesh.vortexIDs) {
    vortexRings.push_back(element::vortexRing(key, ids, 1.0));
  }
  // Building doublet panels
  for (auto &[key, ids] : mesh.doubletIDs) {
    doubletPanels.push_back(element::doubletPanel(key, ids, 1.0));
  }
}

double model::get_cl() const { return cl; }

double model::get_cd() const { return cd; }

Vector3d model::get_cm() const { return cm; }
