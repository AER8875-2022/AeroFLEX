
#ifndef PANEL_HPP
#define PANEL_HPP

#include "Eigen/Dense"
#include "vlm/input.hpp"
#include "vlm/potential.hpp"

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

  /** @brief Area of the panel */
  double area;

  /** @brief Vector normal to the panel */
  Vector3d normal;

  /** @brief Coordinates of the panel's center */
  Vector3d center;

  /** @brief Vector holding the local reference system of the panel (l,m,normal) */
  std::array<Vector3d, 3> Localreference;

  /** @brief Vector holding the global coordonates of each edge's center */
  std::vector<Vector3d> edge_center;

  /** @brief Vector holding the geometry of the panel's edges */
  std::vector<geom::edgeLine> edges;

public:
  /** @param nodeIDs IDs of the nodes defining the panel */
  panel(const std::vector<int> &nodeIDs);

  /** @brief Method that updates the geometry of the panel from moved nodes
   *  @param nodes Nodes of the mesh */
  void updateGeometry(const std::vector<Vector3d> &nodes);

  /** @brief Method that initializes the geometry of the panel
   *  @param nodes Nodes of the mesh */
  void initialize(const std::vector<Vector3d> &nodes);

  /** @brief Method that compute the local reference system of the panel (l,m,normal) */
  void LocalCoordinate(const std::vector<Vector3d> &nodes); 

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
  /** @brief Method that computes the center of the panel
   *  @param nodes Nodes of the mesh */
  void computeCenter(const std::vector<Vector3d> &nodes);

  /** @brief Method that computes the vector normal to the panel
   *  @param nodes Nodes of the mesh */
  void computeNormal(const std::vector<Vector3d> &nodes);

  /** @brief Method that computes the area of the panel
   *  @param nodes Nodes of the mesh */
  void computeArea(const std::vector<Vector3d> &nodes);

  /** @brief Getter method for the local coordinates */
  std::array<Vector3d, 3> get_LocalCoordinate() const;

  friend class element::vortexRing;
  friend class element::doubletPanel;
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
   *  @param gamma Strength of the vortices bouding the ring */
  vortexRing(const int globalIndex, const std::vector<int> &nodeIDs,
             const double gamma);

  /** @brief Method initializing the element's geometry
   *  @param nodes Nodes of the mesh
   *  @param sim Simulation parameters */
  void initialize(const std::vector<Vector3d> &nodes,
                  const input::simParam &sim);

  /** @brief Method computing the coordinates where to forces are acting
   *  @param nodes Nodes of the mesh
   *  @return Coordinates of the force acting point */
  Vector3d forceActingPoint(const std::vector<Vector3d> &nodes) const;

  /** @brief Method computing the length of the trailing edge of the ring
   *  @param nodes Nodes of the mesh
   *  @return Vector defining the trailing edge of the ring */
  Vector3d leadingEdgeDl(const std::vector<Vector3d> &nodes);

  /** @brief Method computing the influence of the current vortex ring on a
   * collocation point
   *  @param collocationPoint Point in tridimensional space
   *  @param nodes Nodes of the mesh
   *  @return Induced velocity vector at collocation point */
  Vector3d influence_wing(const Vector3d &collocationPoint,
                     const std::vector<Vector3d> &nodes) const;


  double influence_patch(const Vector3d &collocationPoint,
                     const std::vector<Vector3d> &nodes) const;

  /** @brief Method computing the streamwise influence (only) of the current
   * vortex ring on a collocation point
   *  @param collocationPoint Point in tridimensional space
   *  @param nodes Nodes of the mesh
   *  @return Induced velocity vector at collocation point */
  Vector3d streamInfluence(const Vector3d &collocationPoint,
                           const std::vector<Vector3d> &nodes) const;

  /** @brief Method updating the geometry of the current element based on new
   * moved nodes
   *  @param nodes Nodes of the mesh */
  void updateGeometry(const std::vector<Vector3d> &nodes);

  /** @brief Saves and updates the strength of the current vortex ring
   *  @param gamma Strength to be saved to the current element */
  void updateGamma(const double gamma);

  /** @brief Method computing the collocation point corrected for high angles of
   * attack
   *  @param nodes Nodes of the mesh */
  void computeCollocationPoint(const std::vector<Vector3d> &nodes);

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

  /** @brief Getter method for the local coordinates */
  std::array<Vector3d, 3> get_LocalCoordinate() const;

  friend class surface::wingStation;
  friend class solver::linear::steady;
  friend class solver::linear::unsteady;
  friend class solver::nonlinear::steady;
  friend class solver::nonlinear::unsteady;
};

/** @brief Object representing a doublet panel element */
class doubletPanel {
  
  /** @brief Unique global index of the current element */
  int globalIndex;

  /** @brief Strength of the doublet of the current element */
  double mu;

  /** @brief Vector containing the ID of each neighbor the current element */
  std::vector<int> NeighborPanel_IDs;

  /** @brief local velocity of the current element */
  Vector3d local_velocity;

  /** @brief Local pressure coefficient */
  double cp;

  double local_area;

  std::vector<fil::vortexLine> Doublets_vortices;

  /** @brief Wrapper for geometric computations */
  geom::panel panel;

  std::vector<Vector3d> ProjectedNodes;   // NOT USED (to be removed)
  std::vector<int> nondirectNeighbor_IDs; //not used  (to be removed)
  Vector3d segment_normal; // only for toubleshooting
  Vector3d localstream;  // only for toubleshooting
  double velocity_div_vinf; //only for troobleshooting
  Vector3d T; //to be removed
  double SMP; //to be removed                          
  double SMQ; //to be removed

public:
  /** @param globalIndex Unique global index of the current element
   *  @param nodeIDs Nodes defining the element */
  doubletPanel(const int globalIndex, const std::vector<int> &nodeIDs);

  /** @brief Method initializing the element's geometry
   *  @param nodes Nodes of the mesh
   *  @param sim Simulation parameters */
  void initialize(const std::vector<Vector3d> &nodes,
                  const input::simParam &sim);

  /** @brief Method computing the influence of the current vortex ring on a
   * collocation point
   *  @param collocationPoint Point in tridimensional space
   *  @param nodes Nodes of the mesh
   *  @return Induced velocity vector at collocation point */
  Vector3d influence_wing(const Vector3d &collocationPoint,
                     const std::vector<Vector3d> &nodes) const;

  /** @brief Method computing the influence of the potentiel at a
   * collocation point
   *  @param collocationPoint Point in tridimensional space
   *  @param nodes Nodes of the mesh
   *  @return influence of the potentiel at a collocation point */
  double influence_patch(const Vector3d &collocationPoint,
                     const std::vector<Vector3d> &nodes) const;

  
  double influence_sources(const Vector3d &collocationPoint,
                                 const std::vector<Vector3d> &nodes) const;

  /** @brief Method updating the geometry of the current element based on new
   * moved nodes (not used yet to be added for future release)
   *  @param nodes Nodes of the mesh */              
  void updateGeometry(const std::vector<Vector3d> &nodes);


  void updateNodes(const std::vector<Vector3d>
                  &nodes); //not done

  /** @brief Project the vector to the local reference system
   *  @param vector projected to the current element */
  Vector3d ProjectingToLocal(const Vector3d &vector) const;

  void storing_velocity(Vector3d velocity, double normalvelocity); // for troobleshooting

  /** @brief Saves and updates the strength of the current doublets
   *  @param mu Strength to be saved to the current element */
  void updateMu(const double mu);

  /** @brief Getter method for edges */
  std::vector<geom::edgeLine> get_edges() const;
  
  /** @brief Getter method for the area */
  double get_area() const;

  /** @brief Getter method for the normal vector */
  Vector3d get_normal() const;

  /** @brief Getter method for the pressure coefficient */
  double get_cp() const;

  /** @brief Getter method for the normal */
  Vector3d get_center() const;

  /** @brief Getter method for the panel neighbor ID */
  std::vector<int> get_neighbor() const;

  /** @brief Getter method for the nodeIDs */
  std::vector<int> get_nodeIDs() const;

  /** @brief Getter method for each edge center point */
  std::vector<Vector3d> get_edge_center() const;

  /** @brief Getter method for globalIndex */
  int get_globalIndex() const;

  /** @brief Getter method for the local coordinates */
  std::array<Vector3d, 3> get_LocalCoordinate() const;

  /** @brief Getter method for the local velocity */
  Vector3d get_local_velocity() const;
  
  /** @brief Getter method for mu */
  double get_mu() const;

  /** @brief Getter method for local_area */
  double get_local_area() const;

Vector3d get_segment_normal() const; //only for troubleshooting
Vector3d get_localstream() const; //only for troubleshooting
double get_sigma() const; //only for troubleshooting
double get_velocity_div_vinf() const; //only for troubleshooting

  friend class surface::patch;
  friend class solver::linear::steady;
  friend class solver::linear::unsteady;
  friend class solver::nonlinear::steady;
  friend class solver::nonlinear::unsteady;
};

} // namespace element
} // namespace vlm

#endif
