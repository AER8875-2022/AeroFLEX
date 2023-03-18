
#include "vlm/surface.hpp"
#include "vlm/input.hpp"

using namespace vlm;
using namespace surface;
using namespace Eigen;

wingStation::wingStation(const int globalIndex,
                         const std::vector<int> &vortexIDs,
                         std::vector<element::vortexRing> &vortices)
    : globalIndex(globalIndex), vortexIDs(vortexIDs), vortices(vortices),
      area(0.0), cl(0.0), cm(Vector3d::Zero()) {}

void wingStation::initialize(const input::simParam &sim) {
  for (auto &vortexID : vortexIDs) {
    vortices[vortexID].initialize(sim);
  }
  computeArea();
  computeChordLength();
  // Setting the local aoa to the geometric aoa
  local_aoa = sim.aoa;
}

void wingStation::updateGeometry() {
  for (auto &vortexID : vortexIDs) {
    vortices[vortexID].updateGeometry();
  }
  computeArea();
  computeChordLength();
}

void wingStation::generateWake(const double wakeLength,
                               const input::simParam &sim,
                               std::vector<Vector3d> &wakeNodes,
                               std::vector<element::vortexRing> &wakePanels) {
  auto &trailPanel = vortices[vortexIDs.back()];
  auto &nodes = trailPanel.panel.nodes;
  // Adding trailing edge nodes
  auto p1 = nodes[trailPanel.get_nodeIDs()[3]];
  auto p2 = nodes[trailPanel.get_nodeIDs()[2]];
  wakeNodes.push_back(p1);
  wakeNodes.push_back(p2);
  // Generating new wake nodes
  auto &wakeOrientation = sim.freeStream().normalized();
  wakeNodes.push_back(p2 + wakeLength * wakeOrientation);
  wakeNodes.push_back(p1 + wakeLength * wakeOrientation);
  // Generating new wake panel
  int globalIndex = trailPanel.get_globalIndex();
  int id1 = wakeNodes.size() - 4;
  int id2 = wakeNodes.size() - 3;
  int id3 = wakeNodes.size() - 2;
  int id4 = wakeNodes.size() - 1;
  auto wakePanel =
      element::vortexRing(globalIndex, {id1, id2, id3, id4}, wakeNodes, 1.0);
  wakePanel.initialize(sim);
  wakePanels.push_back(wakePanel);
}

Vector3d wingStation::forceActingPoint() const {
  return (vortices[vortexIDs.front()].forceActingPoint());
}

void wingStation::updateLocalAoa(const double dalpha) {
  local_aoa += dalpha;
  for (auto &vortexID : vortexIDs) {
    vortices[vortexID].local_aoa += dalpha;
    vortices[vortexID].computeCollocationPoint();
  }
}

void wingStation::resetLocalAoa(const input::simParam &sim) {
  local_aoa = sim.aoa;
  for (auto &vortexID : vortexIDs) {
    vortices[vortexID].local_aoa = sim.aoa;
    vortices[vortexID].computeCollocationPoint();
  }
}

double wingStation::get_globalIndex() const { return globalIndex; }

double wingStation::get_area() const { return area; }

double wingStation::get_chord() const { return chord; }

double wingStation::get_spanLoc() const { return spanLoc; }

std::vector<int> wingStation::get_vortexIDs() const { return vortexIDs; }

double wingStation::get_cl() const { return cl; }

double wingStation::get_cd() const { return cd; }

Vector3d wingStation::get_cm() const { return cm; }

void wingStation::computeArea() {
  area = 0.0;
  for (auto &vortexID : vortexIDs) {
    area += vortices[vortexID].get_area();
  }
}

void wingStation::computeChordLength() {
  auto &leadingEdge = vortices[vortexIDs.front()];
  auto &trailingEdge = vortices[vortexIDs.back()];
  auto &nodes = leadingEdge.panel.nodes;

  double d1 = nodes[trailingEdge.get_nodeIDs()[3]](0) -
              nodes[leadingEdge.get_nodeIDs()[0]](0);
  double d2 = nodes[trailingEdge.get_nodeIDs()[2]](0) -
              nodes[leadingEdge.get_nodeIDs()[1]](0);
  chord = 0.5 * (d1 + d2);
}

void wingStation::to_local(Vector3d &vector) {

  // Computing local referential
  Vector3d leadingEdge = vortices[vortexIDs.front()].leadingEdgeDl();
  Vector3d x_local = Vector3d::UnitX();
  Vector3d z_local = x_local.cross(leadingEdge).normalized();
  Vector3d y_local = z_local.cross(x_local).normalized();

  // Rotating freeStream
  Matrix3d rotationMatrix{{x_local(0), x_local(1), x_local(2)},
                          {y_local(0), y_local(1), y_local(2)},
                          {z_local(0), z_local(1), z_local(2)}};

  vector = rotationMatrix * vector;
}

void wingStation::to_global(Vector3d &vector) {

  // Computing local referential
  Vector3d leadingEdge = vortices[vortexIDs.front()].leadingEdgeDl();
  Vector3d x_local = Vector3d::UnitX();
  Vector3d z_local = x_local.cross(leadingEdge).normalized();
  Vector3d y_local = z_local.cross(x_local).normalized();

  // Rotating freeStream
  Matrix3d rotationMatrix{{x_local(0), x_local(1), x_local(2)},
                          {y_local(0), y_local(1), y_local(2)},
                          {z_local(0), z_local(1), z_local(2)}};

  vector = rotationMatrix.inverse() * vector;
}

Vector3d wingStation::liftAxis(const input::simParam &sim) {
  Vector3d leadingEdge = vortices[vortexIDs.front()].leadingEdgeDl();
  return (sim.freeStream(local_aoa).cross(leadingEdge).normalized());
}

void wingStation::computeForces(const input::simParam &sim) {

  Vector3d stream = sim.freeStream(local_aoa);

  to_local(stream);

  double previousGamma = 0.0;
  force = Vector3d::Zero();
  moment = Vector3d::Zero();

  for (auto &vortexID : vortexIDs) {
    auto &vortex = vortices[vortexID];

    // Computing distances
    Vector3d dl = vortex.leadingEdgeDl();
    to_local(dl);
    Vector3d lever = sim.origin() - vortex.forceActingPoint();

    // Local forces
    Vector3d local_force = stream.cross(dl);
    local_force *= sim.rho * (vortex.get_gamma() - previousGamma);
    to_global(local_force);

    Vector3d local_moment = local_force.cross(lever);

    // Incrementing total force
    force += local_force;
    moment += local_moment;

    // Panel coefficient
    vortex.cf = local_force / (sim.dynamicPressure() * sim.sref);
    vortex.cm = local_moment / (sim.dynamicPressure() * sim.sref * sim.cref);

    // Setting gamma reference for next panel
    previousGamma = vortex.get_gamma();
  }

  // Oriented in section referential
  cl_local = force.dot(liftAxis(sim)) / (sim.dynamicPressure() * area);
  // Oriented in inertial referential
  cl = force.dot(sim.liftAxis()) / (sim.dynamicPressure() * area);
  cm = moment / (sim.dynamicPressure() * area * sim.cref);
}

// ----------------------------

wing::wing(const int globalIndex, const std::vector<int> &stationIDs,
           std::vector<wingStation> &stations)
    : globalIndex(globalIndex), stationIDs(stationIDs), stations(stations),
      area(0.0), cl(0.0), cm(Vector3d::Zero()) {}

void wing::initialize(const input::simParam &sim) {
  for (auto &stationID : stationIDs) {
    stations[stationID].initialize(sim);
  }
  computeArea();
  computeSpan();
}

void wing::updateGeometry() {
  for (auto &stationID : stationIDs) {
    stations[stationID].updateGeometry();
  }
  computeArea();
  computeSpan();
}

double wing::get_globalIndex() const { return globalIndex; }

double wing::get_area() const { return area; }

double wing::get_span() const { return span; }

std::vector<int> wing::get_stationIDs() const { return stationIDs; }

double wing::get_cl() const { return cl; }

double wing::get_cd() const { return cd; }

Vector3d wing::get_cm() const { return cm; }

void wing::computeArea() {
  area = 0.0;
  for (auto &stationID : stationIDs) {
    area += stations[stationID].get_area();
  }
}

void wing::computeSpan() {
  span = 0.0;
  for (auto &stationID : stationIDs) {
    auto &station = stations[stationID];
    Vector3d leadingEdge = station.vortices[station.vortexIDs.front()].leadingEdgeDl();
    span += leadingEdge.norm();
    station.spanLoc = span - 0.5*leadingEdge.norm();
  }
}

void wing::computeForces(const input::simParam &sim) {
  cl = 0.0;
  cd = 0.0;
  cm = Vector3d::Zero();

  for (auto &stationID : stationIDs) {
    auto &station = stations[stationID];
    cl += station.get_cl() * station.get_area() / area;
    cd += station.get_cd() * station.get_area() / area;
    cm += station.get_cm() * station.get_area() / area * station.get_chord() /
          sim.cref;
  }
}

// ----------------------------

patch::patch(const int globalIndex, const std::vector<int> &doubletIDs,
             std::vector<element::doubletPanel> &doubletPanels)
    : globalIndex(globalIndex), doubletIDs(doubletIDs),
      doubletPanels(doubletPanels), area(0.0) {}

void patch::initialize(const input::simParam &sim) {
  for (auto &doubletID : doubletIDs) {
    doubletPanels[doubletID].initialize(sim);
  }
  computeArea();
}

void patch::updateGeometry() {
  for (auto &doubletID : doubletIDs) {
    doubletPanels[doubletID].updateGeometry();
  }
  computeArea();
}

double patch::get_globalIndex() const { return globalIndex; }

double patch::get_area() const { return area; }

std::vector<int> patch::get_doubletIDs() const { return doubletIDs; }

void patch::computeArea() {
  area = 0.0;
  for (auto &doubletID : doubletIDs) {
    area += doubletPanels[doubletID].get_area();
  }
}
