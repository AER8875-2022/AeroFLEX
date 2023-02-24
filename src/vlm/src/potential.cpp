
#include "vlm/potential.hpp"
#include <cmath>

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

double vortexLine::get_gamma() const { return gamma; }

std::tuple<Vector3d, Vector3d>
vortexLine::distToPoint(const Vector3d &collocationPoint, const Vector3d &p1,
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
