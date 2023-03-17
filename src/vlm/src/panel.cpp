
#include "vlm/panel.hpp"

using namespace vlm;
using namespace geom;
using namespace element;
using namespace Eigen;

// -------------------

panel::panel(const std::vector<int> &nodeIDs, std::vector<Vector3d> &nodes)
    : nodeIDs(nodeIDs), nodes(nodes), area(0.0), normal(Vector3d::Zero()),
      center(Vector3d::Zero()) {
  this->edges.reserve(nodeIDs.size());
}

void panel::updateGeometry() {
  for (auto &edge : edges) {
    edge.updateGeometry(nodes);
  }
  computeCenter();
  computeArea();
  computeNormal();
}

void panel::initialize() {
  if (edges.empty()) {
    for (size_t i = 0; i != nodeIDs.size() - 1; i++) {
      edges.push_back(geom::edgeLine(nodeIDs[i], nodeIDs[i + 1]));
    }
    edges.push_back(geom::edgeLine(nodeIDs.back(), nodeIDs.front()));
  }
  updateGeometry();
}

std::vector<int> panel::get_nodeIDs() const { return nodeIDs; }

double panel::get_area() const { return area; }

Vector3d panel::get_normal() const { return normal; }

Vector3d panel::get_center() const { return center; }

std::vector<geom::edgeLine> panel::get_edges() const { return edges; }

void panel::computeCenter() {
  center = Vector3d::Zero();
  for (auto &nodeID : nodeIDs) {
    auto node = nodes[nodeID];
    center += node;
  }
  center /= (double)nodeIDs.size();
}

void panel::computeNormal() {
  normal = Vector3d::Zero();
  for (size_t i = 0; i != nodeIDs.size() - 1; i++) {
    normal -= edges[i].get_dl().cross(edges[i + 1].get_dl());
  }
  normal -= edges.back().get_dl().cross(edges.front().get_dl());
  normal.normalize();
}

void panel::computeArea() {
  area = 0.0;
  Vector3d cross, centerEdge = Vector3d::Zero();
  for (size_t i = 0; i != nodeIDs.size(); i++) {
    auto node = nodes[nodeIDs[i]];
    centerEdge = center - node;
    cross = edges[i].get_dl().cross(centerEdge);
    area += 0.5 * cross.norm();
  }
}

// -------------------

vortexRing::vortexRing(const int globalIndex, const std::vector<int> &nodeIDs,
                       std::vector<Vector3d> &nodes, const double gamma)
    : globalIndex(globalIndex), gamma(gamma), cl(0.0), cm(Vector3d::Zero()),
      panel(nodeIDs, nodes) {
  vortices.reserve(4);
}

void vortexRing::initialize(const input::simParam &sim) {
  panel.initialize();
  for (size_t i = 0; i != panel.nodeIDs.size(); i++) {
    vortices.push_back(fil::vortexLine(1.0, sim.coreRadius));
  }
  local_aoa = sim.aoa;
  computeCollocationPoint();
}

Vector3d vortexRing::forceActingPoint() const {
  return (0.5 *
          (panel.nodes[panel.nodeIDs[0]] + panel.nodes[panel.nodeIDs[1]]));
}

Vector3d vortexRing::leadingEdgeDl() {
  return (panel.nodes[panel.nodeIDs[1]] - panel.nodes[panel.nodeIDs[0]]);
}

void vortexRing::computeCollocationPoint() {
  // High aoa correction term according to x axis
  double k = local_aoa * M_PI / (180.0 * std::sin(local_aoa * M_PI / 180.0));
  collocationPoint = k * panel.center + (1 - k) * forceActingPoint();
}

Vector3d vortexRing::influence(const Vector3d &collocationPoint) const {
  Vector3d v = Vector3d::Zero();
  for (size_t i = 0; i != vortices.size(); i++) {
    v += vortices[i].influence(collocationPoint, panel.nodes, panel.edges[i]);
  }
  return v;
}

Vector3d vortexRing::streamInfluence(const Vector3d &collocationPoint) const {
  Vector3d v = Vector3d::Zero();
  for (size_t i = 1; i != vortices.size(); i += 2) {
    v += vortices[i].influence(collocationPoint, panel.nodes, panel.edges[i]);
  }
  return v;
}

void vortexRing::updateGeometry() {
  panel.updateGeometry();
  // computeCollocationPoint();
}

void vortexRing::updateGamma(const double gamma) {
  this->gamma = gamma;
  for (auto &vortex : vortices) {
    vortex.gamma = gamma;
  }
}

double vortexRing::get_area() const { return panel.area; }

Vector3d vortexRing::get_normal() const { return panel.normal; }

Vector3d vortexRing::get_center() const { return panel.center; }

std::vector<int> vortexRing::get_nodeIDs() const { return panel.nodeIDs; }

int vortexRing::get_globalIndex() const { return globalIndex; }

double vortexRing::get_gamma() const { return gamma; }

Matrix<double, 6, 1> vortexRing::get_forces() const { return forces; }

double vortexRing::get_cl() const { return cl; }

double vortexRing::get_cd() const { return cd; }

Vector3d vortexRing::get_cm() const { return cm; }

Vector3d vortexRing::get_collocationPoint() const { return collocationPoint; }

// -------------------

doubletPanel::doubletPanel(const int globalIndex,
                           const std::vector<int> &nodeIDs,
                           std::vector<Vector3d> &nodes, const double sigma)
    : globalIndex(globalIndex), sigma(sigma), panel(nodeIDs, nodes) {

  Doublets_vortices.reserve(4);
}

void doubletPanel::initialize(const input::simParam &sim) {
  panel.initialize();
  // Localreference.reserve(3); // potentiellement temporaire
  LocalCoordinate();
  for (size_t i = 0; i != panel.nodeIDs.size(); i++) {
    Doublets_vortices.push_back(fil::vortexLine(1.0, sim.coreRadius));
  }
  /* //PAS UTILE POUR L'INSTANT
  //half median length of the panel (could be place in a new function)
  SMP=
  ((nodes[panel.nodeIDs[3]]+nodes[panel.nodeIDs[2]])/2-panel.center).norm();
  SMQ=((nodes[panel.nodeIDs[0]]+nodes[panel.nodeIDs[3]])/2-panel.center).norm();
  */
}

void doubletPanel::updateGeometry() { panel.updateGeometry(); }

void doubletPanel::LocalCoordinate() {

  Eigen::Vector3d m =
      ((panel.nodes[panel.nodeIDs[1]] + panel.nodes[panel.nodeIDs[2]]) / 2 -
       panel.center)
          .normalized(); // could be simplefied with panel centerpoint of the
                         // edge
  Eigen::Vector3d l = panel.normal.cross(m);
  Localreference = {l, m, panel.normal}; // erreur potentiel dans la notation
  // std::cout << Localreference <<std::endl;
}

void doubletPanel::updateNodes(const std::vector<Vector3d> &nodes) {
  for (size_t n = 0; n < panel.nodeIDs.size(); n++) {
    Eigen::Vector3d tempo = nodes[panel.nodeIDs[n]] - panel.center;
    ProjectedNodes.push_back(
        {tempo.dot(Localreference[0]), tempo.dot(Localreference[1]),
         tempo.dot(Localreference[2])}); // placing it into a different labels
                                         // (temporary)
  }
}
// Solving the influence of the doublets
Vector3d doubletPanel::influence(const Vector3d &collocationPoint) const {
  Vector3d d = Vector3d::Zero();
  // la partie en dessous n'est pas encore adapté
  for (size_t i = 0; i != Doublets_vortices.size(); i++) {
    d += Doublets_vortices[i].influence(collocationPoint, panel.nodes,
                                        panel.edges[i]);
  }
  return d;
}

double doubletPanel::get_area() const { return panel.area; }

Vector3d doubletPanel::get_normal() const { return panel.normal; }

Vector3d doubletPanel::get_center() const { return panel.center; }

std::vector<int> doubletPanel::get_nodeIDs() const { return panel.nodeIDs; }

int doubletPanel::get_globalIndex() const { return globalIndex; }

std::array<Vector3d, 3> doubletPanel::get_LocalCoordinate() const {
  return Localreference;
}

double doubletPanel::get_sigma() const { return sigma; }
