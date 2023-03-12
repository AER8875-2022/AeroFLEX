
#ifndef VLM_PANEL_HPP
#define VLM_PANEL_HPP

#include "Eigen/Dense"
#include "vlm/input.hpp"
#include "vlm/potential.hpp"

#ifndef M_PI
#define M_PI 3.141592653589793115997963468544185161590576171875
#endif

/** @file panel.hpp */

using namespace Eigen;

#ifndef VLM_SOLVER_CLASS
#define VLM_SOLVER_CLASS
namespace vlm::solver::linear {
class steady;
class unsteady;
} // namespace vlm::solver::linear
namespace vlm::solver::nonlinear {
class steady;
class unsteady;
} // namespace vlm::solver::nonlinear
#endif

#ifndef VLM_ELEMENT_CLASS
#define VLM_ELEMENT_CLASS
namespace vlm::element {
class vortexRing;
class doubletPanel;
} // namespace vlm::element
#endif

#ifndef VLM_SURFACE_CLASS
#define VLM_SURFACE_CLASS
namespace vlm::surface {
class wingStation;
class wing;
class patch;
} // namespace vlm::surface
#endif

namespace vlm {

namespace geom {

/** @brief Object that represents the geometry of a panel */
class panel {

  /** @brief Vector holding the IDs of the panel's nodes */
  std::vector<int> nodeIDs;

  /** @brief Reference to nodes */
  std::vector<Vector3d> &nodes;

  /** @brief Area of the panel */
  double area;

  /** @brief Vector normal to the panel */
  Vector3d normal;

  /** @brief Coordinates of the panel's center */
  Vector3d center;

  /** @brief Vector holding the geometry of the panel's edges */
  std::vector<geom::edgeLine> edges;

public:
  /** @param nodeIDs IDs of the nodes defining the panel
   *  @param nodes Nodes of the mesh */
  panel(const std::vector<int> &nodeIDs, std::vector<Vector3d> &nodes);

  /** @brief Method that updates the geometry of the panel from moved nodes */
  void updateGeometry();

  /** @brief Method that initializes the geometry of the panel */
  void initialize();

  /** @brief Getter method for nodeIDs */
  std::vector<int> get_nodeIDs() const;

  /** @brief Getter method for area */
  double get_area() const;

  /** @brief Getter method for normal */
  Vector3d get_normal() const;

  /** @brief Getter method for center */
  Vector3d get_center() const;

  /** @brief Getter method for edges */
  std::vector<geom::edgeLine> get_edges() const;

private:
  /** @brief Method that computes the center of the panel */
  void computeCenter();

  /** @brief Method that computes the vector normal to the panel */
  void computeNormal();

  /** @brief Method that computes the area of the panel */
  void computeArea();

  friend class element::vortexRing;
  friend class element::doubletPanel;
  friend class surface::wingStation;
};

} // namespace geom

namespace element {

/** @brief Object representing a vortex ring element */
class vortexRing {

  /** @brief Unique global index of the current element */
  int globalIndex;

  /** @brief Strength of the vortices of the current element */
  double gamma;

  /** @brief Local lift coefficient */
  double cl;

  /** @brief Local drag coefficient */
  double cd;

  /** @brief Local moment coefficients */
  Vector3d cm;

  /** @brief Wrapper for geometric computations */
  geom::panel panel;

  /** @brief Corrected geoemtric collocation point */
  Vector3d collocationPoint;

  /** @brief Vortex filaments bounding the vortex ring */
  std::vector<fil::vortexLine> vortices;

  /** @brief local angle of attack correction */
  double local_aoa;

public:
  /** @param globalIndex Unique global index of the current element
   *  @param nodeIDs Nodes defining the element
   *  @param nodes Nodes of the mesh
   *  @param gamma Strength of the vortices bouding the ring */
  vortexRing(const int globalIndex, const std::vector<int> &nodeIDs,
             std::vector<Vector3d> &nodes, const double gamma);

  /** @brief Method initializing the element's geometry
   *  @param sim Simulation parameters */
  void initialize(const input::simParam &sim);

  /** @brief Method computing the coordinates where to forces are acting
   *  @return Coordinates of the force acting point */
  Vector3d forceActingPoint() const;

  /** @brief Method computing the length of the trailing edge of the ring
   *  @return Vector defining the trailing edge of the ring */
  Vector3d leadingEdgeDl();

  /** @brief Method computing the influence of the current vortex ring on a
   * collocation point
   *  @param collocationPoint Point in tridimensional space
   *  @return Induced velocity vector at collocation point */
  Vector3d influence(const Vector3d &collocationPoint) const;

  /** @brief Method computing the streamwise influence (only) of the current
   * vortex ring on a collocation point
   *  @param collocationPoint Point in tridimensional space
   *  @return Induced velocity vector at collocation point */
  Vector3d streamInfluence(const Vector3d &collocationPoint) const;

  /** @brief Method updating the geometry of the current element based on new
   * moved nodes */
  void updateGeometry();

  /** @brief Saves and updates the strength of the current vortex ring
   *  @param gamma Strength to be saved to the current element */
  void updateGamma(const double gamma);

  /** @brief Method computing the collocation point corrected for high anlges of
   * attack */
  void computeCollocationPoint();

  /** @brief Getter method for area */
  double get_area() const;

  /** @brief Getter method for normal */
  Vector3d get_normal() const;

  /** @brief Getter method for center */
  Vector3d get_center() const;

  /** @brief Getter method for nodeIDs */
  std::vector<int> get_nodeIDs() const;

  /** @brief Getter method for globalIndex */
  int get_globalIndex() const;

  /** @brief Getter method for gamma */
  double get_gamma() const;

  /** @brief Getter method for cl */
  double get_cl() const;

  /** @brief Getter method for cd */
  double get_cd() const;

  /** @brief Getter method for cm */
  Vector3d get_cm() const;

  /** @brief Getter method for collocationPoint */
  Vector3d get_collocationPoint() const;

  friend class surface::wingStation;
  friend class solver::linear::steady;
  friend class solver::linear::unsteady;
  friend class solver::nonlinear::steady;
  friend class solver::nonlinear::unsteady;
};

class doubletPanel {
  int globalIndex;
  double sigma;
  std::array<Vector3d, 3> Localreference; // (l,m,normal)
  std::vector<Vector3d> ProjectedNodes;   // temporary label (to be removed)
  double SMP;                             // potentiellement inutile
  double SMQ;                             // potentiellement inutile
  std::vector<fil::vortexLine> Doublets_vortices;
  geom::panel panel;

public:
  doubletPanel(const int globalIndex, const std::vector<int> &nodeIDs,
               std::vector<Vector3d> &nodes, const double sigma);
  void initialize(const input::simParam &sim);
  Vector3d influence(const Vector3d &collocationPoint) const;
  void updateGeometry();
  void LocalCoordinate();
  void
  updateNodes(const std::vector<Vector3d>
                  &nodes); // to create later updating/projecting the position
                           // of the node onto the mean of the panel (co-planor
                           // panel in the localcoordinates)
  double get_area() const;
  Vector3d get_normal() const;
  Vector3d get_center() const;
  std::vector<int> get_nodeIDs() const;
  int get_globalIndex() const;
  std::array<Vector3d, 3> get_LocalCoordinate() const; // À ajouter dans le .cpp
  double get_sigma() const;

  friend class surface::patch;
  friend class solver::linear::steady;
  friend class solver::linear::unsteady;
  friend class solver::nonlinear::steady;
  friend class solver::nonlinear::unsteady;
};

} // namespace element
} // namespace vlm

#endif
