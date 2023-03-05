
#ifndef VLM_POTENTIAL_HPP
#define VLM_POTENTIAL_HPP

#include "Eigen/Dense"
#include <tuple>
#include <vector>

#ifndef M_PI
#define M_PI 3.141592653589793115997963468544185161590576171875
#endif

/** @file potential.hpp */

using namespace Eigen;

// Forward declaration of vortexRing class
namespace vlm::element {
class vortexRing;
}

namespace vlm {

namespace geom {

/** @brief Object representing a geometric edge of a surface element */
class edgeLine {

  /** @brief ID of node 1 */
  int n1;

  /** @brief ID of node 2 */
  int n2;

  /** @brief Vector going from node 1 to node 2 */
  Vector3d dl;

public:
  /** @param n1 ID of node 1
   *  @param n2 ID of node 2 */
  edgeLine(const int n1, const int n2);

  /** @brief Method that updates the geometry from new moved nodes */
  void updateGeometry(const std::vector<Vector3d> &nodes);

  /** @brief Getter method for node 1 */
  int get_n1() const;

  /** @brief Getter method for node 2 */
  int get_n2() const;

  /** @brief Getter method for dl */
  Vector3d get_dl() const;

private:
  /** @brief Method that computes the dl vector
   *  @param p1 Coordinates of node 1
   *  @param p2 Coordinates of node 2 */
  void computeDl(const Vector3d &p1, const Vector3d &p2);
};

} // namespace geom

namespace fil {

/** @brief Object that represents a single vortex filament between two points */
class vortexLine {

  /** @brief Strength of the vortex filament */
  double gamma;

  /** @brief Core radius associated with current vortex filament */
  double coreRadius;

public:
  /** @param gamma Strength of the vortex filament
   *  @param coreRadius Core radius associated with current vortex filament */
  vortexLine(const double gamma, const double coreRadius);

  /** @brief Method computing the influenced of the current vortex ring element
   * on a given collocation point
   *  @param collocationPoint Point in a tridimensional space
   *  @param nodes Vector containing all the nodes coordinates
   *  @param edge Edge geometry associated with this vortex filament
   *  @return Induced velocity at collocation point */
  Vector3d influence(const Vector3d &collocationPoint,
                     const std::vector<Vector3d> &nodes,
                     const geom::edgeLine &edge) const;

  /** @brief Getter method for gamma */
  double get_gamma() const;

private:
  /** @brief Method computing the distance between the collocation point and
   * both of the filament's nodes
   *  @param collocationPoint Point in tridimensional space
   *  @param p1 Coordinates of node 1
   *  @param p2 Coordinates of node 2
   *  @return Tuple holding the distances to p1 and p2 */
  std::tuple<Vector3d, Vector3d> distToPoint(const Vector3d &collocationPoint,
                                             const Vector3d &p1,
                                             const Vector3d &p2) const;

  /** @brief Method computing the orientation of the induced velocity at
   * collocation point
   *  @param dl Length of the vortex filament
   *  @param r vector linking p1 or p2 and collocation point */
  Vector3d inducedVelOrientation(const Vector3d &dl, const Vector3d &r) const;

  friend element::vortexRing;
};

} // namespace fil
} // namespace vlm

#endif
