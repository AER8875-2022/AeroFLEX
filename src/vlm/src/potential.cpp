
#include "vlm/potential.hpp"
#include <cmath>
#include <iostream>

using namespace vlm;
using namespace geom;
using namespace fil;
using namespace Eigen;

// -------------------

edgeLine::edgeLine(const int n1, const int n2)
    : n1(n1), n2(n2), dl(Vector3d::Zero()) {}

void edgeLine::updateGeometry(const std::vector<Vector3d> &nodes) {
  auto p1 = nodes[n1];
  auto p2 = nodes[n2];
  computeDl(p1, p2);
}

int edgeLine::get_n1() const { return n1; }

int edgeLine::get_n2() const { return n2; }

Vector3d edgeLine::get_dl() const { return dl; }

void edgeLine::computeDl(const Vector3d &p1, const Vector3d &p2) {
  dl = p2 - p1;
}

// -------------------

vortexLine::vortexLine(const double gamma, const double coreRadius)
    : gamma(gamma), coreRadius(coreRadius) {}

Vector3d vortexLine::influence(const Vector3d &collocationPoint,
                               const std::vector<Vector3d> &nodes,
                               const edgeLine &edge) const {
  auto [r1, r2] =
      distToPoint(collocationPoint, nodes[edge.get_n1()], nodes[edge.get_n2()]);
  // Geometric computations
  auto r1r2 = r1.cross(r2);
  auto norm = r1r2.norm();
  auto dlr1 = edge.get_dl().dot(r1);
  auto dlr2 = edge.get_dl().dot(r2);

  // Computing minimal distance from segment to point
  double h = norm / (dlr1 / r1.norm() - dlr2 / r2.norm());

  // Computing direction of induced velocity
  Vector3d direction = inducedVelOrientation(edge.get_dl(), r1);

  // Computing Lamb-Oseen viscous correction
  double lo = 1.0 - std::exp(-(h / coreRadius) * (h / coreRadius));

  return (gamma / (4.0 * M_PI * h) * lo * direction);
}

double vortexLine::influence_patch(const Vector3d &collocationPoint,
                     const std::vector<Vector3d> &nodes,
                     const geom::edgeLine &edge, 
                     const std::array<Vector3d, 3> &Localreference,
                     const Vector3d &center_point) const {
  
  if (nodes[edge.get_n1()]==nodes[edge.get_n2()]){
    //Pour les points coincidents, l'influence est nul
    //std::cout<< "same" << std::endl;
    return 0;
  } else {
  //il serait mieux de changer le système de référence en matrice plutôt que vecteur
  //Projecting the nodes in the local reference system
  Vector3d nodes1 = {nodes[edge.get_n1()].dot(Localreference[0]), nodes[edge.get_n1()].dot(Localreference[1]) , nodes[edge.get_n1()].dot(Localreference[2])};
  Vector3d nodes2 = {nodes[edge.get_n2()].dot(Localreference[0]), nodes[edge.get_n2()].dot(Localreference[1]) , nodes[edge.get_n2()].dot(Localreference[2])};
  auto [a, b] = distToPoint(collocationPoint, nodes1, nodes2);
  // (l,m,normal)
  auto SM = edge.get_dl().dot(Localreference[1]);
  auto SL = edge.get_dl().dot(Localreference[0]);
  auto A = a.norm(); 
  auto B = b.norm();
  auto AL = a.dot(Localreference[0]);
  auto AM = a.dot(Localreference[1]);
  auto BM = b.dot(Localreference[1]);
  auto Al= AM*SL - AL*SM;

//std::cout<< "new value" << std::endl;
//std::cout<< Al << std::endl;
//std::cout<< Al2 << std::endl; //l'un est le négatif de l'autre
  //auto Rck = center_point; //Erreur possible (Centre local ou **centre global**)
  Vector3d local_center_point = {center_point.dot(Localreference[0]), center_point.dot(Localreference[1]), center_point.dot(Localreference[2])};
  auto Pjk = collocationPoint - local_center_point;
  auto PN = Pjk.dot(Localreference[2]);

  auto PA= (PN*PN) * SL + Al * AM;
  auto PB= (PN*PN) * SL + Al * BM;

  auto RNUM = SM * PN * (B * PA - A * PB);
  auto DNOM = PA * PB + (PN*PN) * A * B * (SM*SM);

  //std::cout << "vector new" << std::endl;
  //std::cout << DNOM << std::endl;

  //atan is in radians
  return {atan2(RNUM, DNOM)};
  }
}

double vortexLine::get_gamma() const { return gamma; }

std::tuple<Vector3d, Vector3d> vortexLine::distToPoint(const Vector3d &collocationPoint, const Vector3d &p1,
                        const Vector3d &p2) const {
  Vector3d r1 = collocationPoint - p1;
  Vector3d r2 = collocationPoint - p2;
  return {r1, r2};
}

Vector3d vortexLine::inducedVelOrientation(const Vector3d &dl,
                                           const Vector3d &r) const {
  auto direction = dl.cross(r);
  return (direction.normalized());
}
