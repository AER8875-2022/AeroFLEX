
#include "vlm/panel.hpp"
#include "iostream"
#include <tuple>

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
  LocalCoordinate(nodes);
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

void panel::LocalCoordinate(const std::vector<Vector3d> &nodes) {
  Eigen::Vector3d m = ((nodes[nodeIDs[1]] + nodes[nodeIDs[2]]) / 2 - center).normalized();
  Eigen::Vector3d l = m.cross(normal);
  Localreference = {l, m, normal};
}

std::vector<int> panel::get_nodeIDs() const { return nodeIDs; }

double panel::get_area() const { return area; }

Vector3d panel::get_normal() const { return normal; }

std::array<Vector3d, 3> panel::get_LocalCoordinate() const { return Localreference; };

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
   std::tuple<double, bool, int> v_local_and_inside_sign;
  double v = 0.0;
  bool inside;
  double v_local;
  int sign;
  const double rlim = 1e-12; 
  for (size_t i = 0; i != vortices.size(); i++) {
     //globalindex added for debugging remove when done
    v_local_and_inside_sign = vortices[i].influence_patch(collocationPoint, nodes, panel.edges[i], panel.Localreference, panel.center, rlim, globalIndex);
    std::tie(v_local, inside, sign) = v_local_and_inside_sign;
    v += v_local;

    //pas fini
  }
  return v/(-4 * M_PI);
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

std::array<Vector3d, 3> vortexRing::get_LocalCoordinate() const { return panel.Localreference; };

// -------------------

doubletPanel::doubletPanel(const int globalIndex,
                           const std::vector<int> &nodeIDs)
    : globalIndex(globalIndex), panel(nodeIDs) {

  NeighborPanel_IDs.reserve(panel.nodeIDs.size());
  for (size_t w=0 ; w<panel.nodeIDs.size(); w++){
    NeighborPanel_IDs.push_back(-1);
  }
  nondirectNeighbor_IDs.reserve(panel.nodeIDs.size()); //not used to be removed during clean up
  Doublets_vortices.reserve(panel.nodeIDs.size());
}

void doubletPanel::initialize(const std::vector<Vector3d> &nodes,
                              const input::simParam &sim) {
  panel.initialize(nodes);
  for (size_t i = 0; i != panel.nodeIDs.size(); i++) {
    Doublets_vortices.push_back(fil::vortexLine(1.0, sim.coreRadius));
  }
}

void doubletPanel::updateGeometry(const std::vector<Vector3d> &nodes) {
  panel.updateGeometry(nodes);
}

void doubletPanel::updateMu(const double mu) {
  this->mu = mu;
  for (auto &vortice : Doublets_vortices) {
    vortice.mu = mu;
  }
    std::cout<< "panel number " << globalIndex << std::endl;
    std::cout << "value of mu :" << mu << std::endl;
}
//not used yet
void doubletPanel::updateNodes(const std::vector<Vector3d> &nodes) {
  for (size_t n = 0; n < panel.nodeIDs.size(); n++) {
    Eigen::Vector3d tempo = nodes[panel.nodeIDs[n]] - panel.center;
    ProjectedNodes.push_back(
        {tempo.dot(panel.Localreference[0]), tempo.dot(panel.Localreference[1]),
         tempo.dot(panel.Localreference[2])});
  }
}

// Solving the influence of the wing (Neumann) 
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
  std::tuple<double, bool, int> C_local_and_bool_sign;
  double C = 0.0;
  bool inside;
  double C_local;
  const double rlim = 1e-12;
  int sign;
  
  for (size_t i = 0; i != Doublets_vortices.size(); i++) {
    //globalindex added for debugging remove when done
    C_local_and_bool_sign=
        Doublets_vortices[i].influence_patch(collocationPoint, nodes, panel.edges[i], panel.Localreference, panel.center, rlim, globalIndex);
        std::tie(C_local, inside, sign) = C_local_and_bool_sign;
    C+= C_local;
  }
  //std::cout << C/(-4 * M_PI) << std::endl;
  Vector3d local_center_point = {panel.center.dot(panel.Localreference[0]), panel.center.dot(panel.Localreference[1]), panel.center.dot(panel.Localreference[2])};
  if (std::abs((collocationPoint - local_center_point)[2]) < rlim){ //.dot(panel.Localreference[2])
    if (inside){
      C = sign * 2 * M_PI; //necessite l'ajout du signe
    }
  }
  return C/(-4 * M_PI);
}

// Solving the influence of the sources
double doubletPanel::influence_sources(const Vector3d &collocationPoint,
                                 const std::vector<Vector3d> &nodes) const {
  double B;
  const double rlim = 1e-12; 
  for (size_t i = 0; i != Doublets_vortices.size(); i++) {
    B +=
        Doublets_vortices[i].influence_sources(collocationPoint, nodes, panel.edges[i], panel.Localreference, panel.center, rlim);
  }

  //std::cout << "panel num" << globalIndex << std::endl; 
  //std::cout << "influence source :" << B << std::endl; 
  return B/(-4 * M_PI);
}

Vector3d doubletPanel::ProjectingToLocal(const Vector3d &vector) const{
  return {vector.dot(panel.Localreference[0]), vector.dot(panel.Localreference[1]),
         vector.dot(panel.Localreference[2])};
}

void doubletPanel::storing_velocity(Vector3d velocity, double normalvelocity) {
  local_velocity = velocity;
  velocity_div_vinf = normalvelocity;
  
}
void doubletPanel::storing_sigma(double sigma1){
  sigma = sigma1;
}

double doubletPanel::get_area() const { return panel.area; }

Vector3d doubletPanel::get_normal() const { return panel.normal; }

double doubletPanel::get_mu() const { return mu; }

double doubletPanel::get_cp() const { return cp; }

std::vector<int> doubletPanel::get_neighbor() const { return NeighborPanel_IDs; }

Vector3d doubletPanel::get_center() const { return panel.center; }

Vector3d doubletPanel::get_local_velocity() const { return local_velocity; }

std::vector<geom::edgeLine> doubletPanel::get_edges() const { return panel.edges; }

std::vector<Vector3d> doubletPanel::get_edge_center() const { return panel.edge_center; }

std::vector<int> doubletPanel::get_nodeIDs() const { return panel.nodeIDs; }

int doubletPanel::get_globalIndex() const { return globalIndex; }

std::array<Vector3d, 3> doubletPanel::get_LocalCoordinate() const { return panel.Localreference; }

Vector3d doubletPanel::get_segment_normal() const { return segment_normal; } //for troubleshooting

Vector3d doubletPanel::get_localstream() const { return localstream; } //for troubleshooting

double doubletPanel::get_sigma() const { return sigma; } //inutile

double doubletPanel::get_velocity_div_vinf() const { return velocity_div_vinf; } //for troubleshooting
