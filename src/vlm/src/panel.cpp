
#include "vlm/panel.hpp"
#include "iostream"

using namespace vlm;
using namespace geom;
using namespace element;
using namespace Eigen;

// -------------------

panel::panel(const std::vector<int> &nodeIDs)
    : nodeIDs(nodeIDs), area(0.0), normal(Vector3d::Zero()),
      center(Vector3d::Zero()) {
  this->edges.reserve(nodeIDs.size());
}

void panel::updateGeometry(const std::vector<Vector3d> &nodes) {
  for (auto &edge : edges) {
    edge.updateGeometry(nodes);
  }
  computeCenter(nodes);
  computeArea(nodes);
  computeNormal(nodes);
}

void panel::initialize(const std::vector<Vector3d> &nodes) {
  for (size_t i = 0; i != nodeIDs.size() - 1; i++) {
    edges.push_back(geom::edgeLine(nodeIDs[i], nodeIDs[i + 1]));
    edge_center.push_back(nodes[nodeIDs[i + 1]] - ((nodes[nodeIDs[i + 1]] - nodes[nodeIDs[i]])/2));
  }
  edges.push_back(geom::edgeLine(nodeIDs.back(), nodeIDs.front()));
  edge_center.push_back(nodes[nodeIDs.front()] - ((nodes[nodeIDs.front()] - nodes[nodeIDs.back()])/2));
  updateGeometry(nodes);
}

std::vector<int> panel::get_nodeIDs() const { return nodeIDs; }

double panel::get_area() const { return area; }

Vector3d panel::get_normal() const { return normal; }

Vector3d panel::get_center() const { return center; }

std::vector<geom::edgeLine> panel::get_edges() const { return edges; }

void panel::computeCenter(const std::vector<Vector3d> &nodes) {
  center = Vector3d::Zero();
  for (auto &nodeID : nodeIDs) {
    auto node = nodes[nodeID];
    center += node;
  }
  center /= (double)nodeIDs.size();
}

void panel::computeNormal(const std::vector<Vector3d> &nodes) {
  normal = Vector3d::Zero();
  for (size_t i = 0; i != nodeIDs.size() - 1; i++) {
    normal -= edges[i].get_dl().cross(edges[i + 1].get_dl());
  }
  normal -= edges.back().get_dl().cross(edges.front().get_dl());
  normal.normalize();
}

void panel::computeArea(const std::vector<Vector3d> &nodes) {
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
                       const double gamma)
    : globalIndex(globalIndex), gamma(gamma), cl(0.0), cm(Vector3d::Zero()),
      panel(nodeIDs) {
  vortices.reserve(4);
}

void vortexRing::initialize(const std::vector<Vector3d> &nodes,
                            const input::simParam &sim) {
  panel.initialize(nodes);
  for (size_t i = 0; i != panel.nodeIDs.size(); i++) {
    vortices.push_back(fil::vortexLine(1.0, sim.coreRadius));
  }
  local_aoa = sim.aoa;
  computeCollocationPoint(nodes);
  LocalCoordinate(nodes);
}

Vector3d
vortexRing::forceActingPoint(const std::vector<Vector3d> &nodes) const {
  return (0.5 * (nodes[panel.nodeIDs[0]] + nodes[panel.nodeIDs[1]]));
}

Vector3d vortexRing::leadingEdgeDl(const std::vector<Vector3d> &nodes) {
  return (nodes[panel.nodeIDs[1]] - nodes[panel.nodeIDs[0]]);
}

void vortexRing::computeCollocationPoint(const std::vector<Vector3d> &nodes) {
  // High aoa correction term according to x axis
  double k = local_aoa * M_PI / (180.0 * std::sin(local_aoa * M_PI / 180.0));
  collocationPoint = k * panel.center + (1 - k) * forceActingPoint(nodes);
}

void vortexRing::LocalCoordinate(const std::vector<Vector3d> &nodes) {

  Eigen::Vector3d m =
      ((nodes[panel.nodeIDs[1]] + nodes[panel.nodeIDs[2]]) / 2 - panel.center)
          .normalized(); // could be simplefied with panel centerpoint of the
                         // edge
  Eigen::Vector3d l = panel.normal.cross(m);
  Localreference = {l, m, panel.normal};
}

//influence des panneaux sur l'aile (Condition de Neumann)
Vector3d vortexRing::influence_wing(const Vector3d &collocationPoint,
                               const std::vector<Vector3d> &nodes) const {
  Vector3d v = Vector3d::Zero();
  for (size_t i = 0; i != vortices.size(); i++) {
    v += vortices[i].influence(collocationPoint, nodes, panel.edges[i]);
  }
  return v;
}
//influence des panneaux sur la surface non-portante/fuselage (Condition de Dirichlet)
double vortexRing::influence_patch(const Vector3d &collocationPoint,
                               const std::vector<Vector3d> &nodes) const {
  //Vector3d v = Vector3d::Zero();
  double v;
  for (size_t i = 0; i != vortices.size(); i++) {
    v += vortices[i].influence_patch(collocationPoint, nodes, panel.edges[i], Localreference, panel.center);
  }
  return v;
}
Vector3d vortexRing::streamInfluence(const Vector3d &collocationPoint,
                                     const std::vector<Vector3d> &nodes) const {
  Vector3d v = Vector3d::Zero();
  for (size_t i = 1; i != vortices.size(); i += 2) {
    v += vortices[i].influence(collocationPoint, nodes, panel.edges[i]);
  }
  return v;
}

void vortexRing::updateGeometry(const std::vector<Vector3d> &nodes) {
  panel.updateGeometry(nodes);
  computeCollocationPoint(nodes);
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

double vortexRing::get_cl() const { return cl; }

double vortexRing::get_cd() const { return cd; }

Vector3d vortexRing::get_cm() const { return cm; }

Vector3d vortexRing::get_collocationPoint() const { return collocationPoint; }

// -------------------

doubletPanel::doubletPanel(const int globalIndex,
                           const std::vector<int> &nodeIDs, const double sigma)
    : globalIndex(globalIndex), sigma(sigma), panel(nodeIDs) {

  //initializing Neighbor vector size and initiale value (-1)
  NeighborPanel_IDs.reserve(panel.nodeIDs.size());
  //std::fill(NeighborPanel_IDs.begin(), NeighborPanel_IDs.end(), -1);
  for (int w=0 ; w<panel.nodeIDs.size(); w++){
    NeighborPanel_IDs.push_back(-1);
  }
 //not used to be removed during clean up
  nondirectNeighbor_IDs.reserve(panel.nodeIDs.size());
  Doublets_vortices.reserve(panel.nodeIDs.size());
}

void doubletPanel::initialize(const std::vector<Vector3d> &nodes,
                              const input::simParam &sim) {
  panel.initialize(nodes);
  // Localreference.reserve(3); // potentiellement temporaire
  LocalCoordinate(nodes);
  for (size_t i = 0; i != panel.nodeIDs.size(); i++) {
    Doublets_vortices.push_back(fil::vortexLine(1.0, sim.coreRadius));
  }
  //half median length of the panel (not used removed during the clean up)
  /*
  T = (nodes[panel.nodeIDs[3]]+nodes[panel.nodeIDs[2]])/2-panel.center;
  SMP=((nodes[panel.nodeIDs[3]]+nodes[panel.nodeIDs[2]])/2-panel.center).norm();
  SMQ=((nodes[panel.nodeIDs[0]]+nodes[panel.nodeIDs[3]])/2-panel.center).norm();
  */
}

void doubletPanel::updateGeometry(const std::vector<Vector3d> &nodes) {
  panel.updateGeometry(nodes);
}

void doubletPanel::updateMu(const double mu) {
  this->mu = mu;
  for (auto &vortice : Doublets_vortices) {
    vortice.mu = mu;
  }
  //std::cout<< "panel number " << globalIndex << std::endl;
  //std::cout << "value of mu :" << mu << std::endl;
}

void doubletPanel::LocalCoordinate(const std::vector<Vector3d> &nodes) {

  Eigen::Vector3d m =
      ((nodes[panel.nodeIDs[1]] + nodes[panel.nodeIDs[2]]) / 2 - panel.center)
          .normalized(); // could be simplefied with panel centerpoint of the
                         // edge
  Eigen::Vector3d l = panel.normal.cross(m);
  Localreference = {l, m, panel.normal};
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
// Solving the influence of the wing (Neumann) 
// à modifier pour influence_wing
Vector3d doubletPanel::influence_wing(const Vector3d &collocationPoint,
                                 const std::vector<Vector3d> &nodes) const {
  Vector3d d = Vector3d::Zero();
  // la partie en dessous n'est pas encore adapté
  for (size_t i = 0; i != Doublets_vortices.size(); i++) {
    d +=
        Doublets_vortices[i].influence(collocationPoint, nodes, panel.edges[i]);
  }
  return d;
}

// Solving the influence of the doublets
double doubletPanel::influence_patch(const Vector3d &collocationPoint,
                                 const std::vector<Vector3d> &nodes) const {
  //Vector3d d = Vector3d::Zero();
  double d;
  // la partie en dessous n'est pas encore adapté
  for (size_t i = 0; i != Doublets_vortices.size(); i++) {
    d +=
        Doublets_vortices[i].influence_patch(collocationPoint, nodes, panel.edges[i], Localreference, panel.center);
  }
  return d;
}

Vector3d doubletPanel::ProjectingCollocation(const Vector3d &collocationPoint) const{
  return {collocationPoint.dot(Localreference[0]), collocationPoint.dot(Localreference[1]),
         collocationPoint.dot(Localreference[2])};
}
void doubletPanel::storing_velocity(Vector3d velocity) {
  local_velocity = velocity;
}

double doubletPanel::get_area() const { return panel.area; }

Vector3d doubletPanel::get_normal() const { return panel.normal; }

double doubletPanel::get_mu() const { return mu; }

double doubletPanel::get_cp() const { return cp; }

std::vector<int> doubletPanel::get_neighbor() const { return NeighborPanel_IDs; }

Vector3d doubletPanel::get_center() const { return panel.center; }

std::vector<geom::edgeLine> doubletPanel::get_edges() const { return panel.edges; }

std::vector<Vector3d> doubletPanel::get_edge_center() const { return panel.edge_center; }

std::vector<int> doubletPanel::get_nodeIDs() const { return panel.nodeIDs; }

int doubletPanel::get_globalIndex() const { return globalIndex; }

std::array<Vector3d, 3> doubletPanel::get_LocalCoordinate() const {
  return Localreference;
}


double doubletPanel::get_sigma() const { return sigma; } //inutile
