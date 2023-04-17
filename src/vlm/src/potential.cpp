
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

std::tuple<double, bool, int> vortexLine::influence_patch(const Vector3d &collocationPoint,
                         const std::vector<Vector3d> &nodes,
                         const geom::edgeLine &edge, 
                         const std::array<Vector3d, 3> &Localreference,
                         const Vector3d &center_point, const double &rlim, const int &globalIndex) const {
  
  if ((nodes[edge.get_n1()] - nodes[edge.get_n2()]).norm() < rlim){
    return {0.0, false, -1};
  } else  {
  bool inside = true;
  int sign = 1;
  Vector3d nodes1 = {nodes[edge.get_n1()].dot(Localreference[0]), nodes[edge.get_n1()].dot(Localreference[1]), nodes[edge.get_n1()].dot(Localreference[2])};
  Vector3d nodes2 = {nodes[edge.get_n2()].dot(Localreference[0]), nodes[edge.get_n2()].dot(Localreference[1]), nodes[edge.get_n2()].dot(Localreference[2])};
  Vector3d local_center_point = {center_point.dot(Localreference[0]), center_point.dot(Localreference[1]), center_point.dot(Localreference[2])};
  auto [a, b] = distToPoint(collocationPoint - local_center_point, nodes1 - local_center_point , nodes2 - local_center_point);

  auto SM = (edge.get_dl()).dot(Localreference[1]);
  auto SL = (edge.get_dl()).dot(Localreference[0]);
  auto A = a.norm(); 
  auto B = b.norm();
  auto AL = a[0]; //.dot(Localreference[0]);
  auto AM = a[1]; //.dot(Localreference[1]);
  auto BM = b[1]; //.dot(Localreference[1]);
  auto Al= AM*SL - AL*SM;
  
  auto Pjk = collocationPoint - local_center_point;
  auto PN = Pjk[2]; //.dot(Localreference[2]);
  auto PA = (PN*PN) * SL + Al * AM;
  auto PB = (PN*PN) * SL + Al * BM;

  //std::cout<< "new value Panel " << globalIndex << " : " << std::endl;
  //std::cout<< PA << std::endl;

  if (PN >= 0.0) {sign = 1.0;} else {sign = -1.0;}

  auto RNUM = SM * PN * (B * PA - A * PB);
  auto DNOM = PA * PB + (PN*PN) * A * B * (SM*SM);

  auto Cx = edge.get_dl().dot(Localreference[0])/(edge.get_dl().norm());
  auto Sy = edge.get_dl().dot(Localreference[1])/(edge.get_dl().norm());
  auto R12 = a[0] * Sy - a[1] * Cx;

  if (R12<0) {
    inside = false;
  }

  if (std::abs(PN) < rlim){
    return {0, inside, sign};
  }
  else {
    return {atan2(RNUM, DNOM), inside, sign}; //atan is in radians
  }  
  }
}

double vortexLine::influence_sources(const Vector3d &collocationPoint,
                     const std::vector<Vector3d> &nodes,
                     const geom::edgeLine &edge, 
                     const std::array<Vector3d, 3> &Localreference,
                     const Vector3d &center_point, const double &rlim) const {
  
  if ((nodes[edge.get_n1()] - nodes[edge.get_n2()]).norm() < rlim){
    return 0;
  } 
  else {
    //int sign;
    Vector3d nodes1 = {nodes[edge.get_n1()].dot(Localreference[0]), nodes[edge.get_n1()].dot(Localreference[1]) , nodes[edge.get_n1()].dot(Localreference[2])};
    Vector3d nodes2 = {nodes[edge.get_n2()].dot(Localreference[0]), nodes[edge.get_n2()].dot(Localreference[1]) , nodes[edge.get_n2()].dot(Localreference[2])};
    Vector3d local_center_point = {center_point.dot(Localreference[0]), center_point.dot(Localreference[1]), center_point.dot(Localreference[2])};
    // local reference (l,m,normal)
    auto [a, b] = distToPoint(collocationPoint - local_center_point, nodes1 - local_center_point, nodes2 - local_center_point);
    
    auto SM = (edge.get_dl()).dot(Localreference[1]);
    auto SL = (edge.get_dl()).dot(Localreference[0]);
    auto A = a.norm(); 
    auto B = b.norm();
    auto AL = a[0];//.dot(Localreference[0]);
    auto AM = a[1];//.dot(Localreference[1]);
    auto BM = b[1];//.dot(Localreference[1]);
    auto Al= AM*SL - AL*SM;

    auto Pjk = collocationPoint - local_center_point;
    auto PN = Pjk[2]; //.dot(Localreference[2]);
    auto PA = (PN*PN) * SL + Al * AM;
    auto PB = (PN*PN) * SL + Al * BM;

    //if (PN >= 0.0) {sign = 1.0;} else {sign = -1.0;}

    auto RNUM = SM * PN * (B * PA - A * PB);
    auto DNOM = PA * PB + (PN*PN) * A * B * (SM*SM);

    double Cjk;
    if (std::abs(PN) < rlim){
      Cjk = 0;
    }
    else {
      Cjk = atan2(RNUM, DNOM); //atan is in radians
    }

    double GL;
    if (std::abs(A + B - edge.get_dl().norm()) < rlim){
      GL = 0;
    } else {
      GL = (1/edge.get_dl().norm()) * std::log(std::abs((A + B + edge.get_dl().norm())/(A + B - edge.get_dl().norm())));
    }

  //std::cout<< "new value" << std::endl;
  //std::cout<< GL << std::endl;
  return (Al * GL - PN * Cjk);
  }
}

double vortexLine::get_gamma() const { return gamma; }

std::tuple<Vector3d, Vector3d> vortexLine::distToPoint(const Vector3d &collocationPoint, const Vector3d &p1,
                        const Vector3d &p2) const {
  Vector3d r1 = collocationPoint - p1; //for the panel method P1 is the second point (since it anti-clockwise)
  Vector3d r2 = collocationPoint - p2; //for the panel method P2 is the first point (since it anti-clockwise)
  return {r1, r2};
}

Vector3d vortexLine::inducedVelOrientation(const Vector3d &dl,
                                           const Vector3d &r) const {
  auto direction = dl.cross(r);
  return (direction.normalized());
}
