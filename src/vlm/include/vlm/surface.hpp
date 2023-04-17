
#ifndef VLM_SURFACE_HPP
#define VLM_SURFACE_HPP

#include "Eigen/Dense"
#include "vlm/input.hpp"
#include "vlm/panel.hpp"

/** @file surface.hpp */

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

namespace vlm {

namespace surface {

/** @brief Object that holds the characteristics of a wing station/section */
class wingStation {

  /** @brief Unique global index of the current element */
  int globalIndex;

  /** @brief IDs of the vortices forming to the current element */
  std::vector<int> vortexIDs;

  /** @brief Area of the current element */
  double area;

  /** @brief Local chord of the current wing station */
  double chord;

  /** @brief Local lift coefficient */
  double cl = 0.0;

  /** @brief Local drag coefficient */
  double cd = 0.0;

  /** @brief Local moment coefficients */
  Vector3d cm = Vector3d::Zero();

  /** @brief local angle of attack correction */
  double local_aoa = 0.0;

  /** @brief Spanwise location of the current wing station */
  double spanLoc;

public:
  /** @param globalIndex Unique global index of the current element
   *  @param vortexIDs IDs of the vortices forming the current element */
  wingStation(const int globalIndex, const std::vector<int> &vortexIDs);

  /** @brief Method initializing the current element
   *  @param nodes Nodes of the mesh
   *  @param vortices Vortex ring elements
   *  @param sim Simulation parameters */
  void initialize(const std::vector<Vector3d> &nodes,
                  std::vector<element::vortexRing> &vortices,
                  const input::simParam &sim);

  /** @brief Method generating a wake panel for the current wing station
   *  @param wakeLength Length of the generated wake panel
   *  @param nodes Nodes of the mesh
   *  @param vortices Vortex ring elements
   *  @param sim Simulation parameters
   *  @param wakeNodes To be generated wake nodes
   *  @param wakePanels To be generated wake panels */
  void generateWake(const double wakeLength, const std::vector<Vector3d> &nodes,
                    const std::vector<element::vortexRing> &vortices,
                    const input::simParam &sim,
                    std::vector<Vector3d> &wakeNodes,
                    std::vector<element::vortexRing> &wakePanels);

  /** @brief Method updating the geometry of the current element base on new
   * moved nodes
   *  @param nodes Nodes of the mesh
   *  @param vortices Vortex ring elements */
  void updateGeometry(const std::vector<Vector3d> &nodes,
                      std::vector<element::vortexRing> &vortices);

  /** @brief Method computing the coordinates where to forces are acting
   *  @param nodes Nodes of the mesh
   *  @param vortices Vortex ring elements
   *  @return Coordinates of the force acting point */
  Vector3d
  forceActingPoint(const std::vector<Vector3d> &nodes,
                   const std::vector<element::vortexRing> &vortices) const;

  /** @brief Method updating the local angle of attack for the current wing
   * station
   *  @param dalpha Variation of the local angle of attack
   *  @param vortices Vortex ring elements */
  void updateLocalAoa(const double dalpha, const std::vector<Vector3d> &nodes,
                      std::vector<element::vortexRing> &vortices);

  /** @brief Method computing the inviscid local forces on the current element */
  void computeForces(const input::simParam &sim,
                     const std::vector<Vector3d> &nodes,
                     std::vector<element::vortexRing> &vortices);

  /** @brief Getter method for globalIndex */
  double get_globalIndex() const;

  /** @brief Getter method for area */
  double get_area() const;

  /** @brief Getter method for chord */
  double get_chord() const;

  /** @brief Getter method for spanLoc */
  double get_spanLoc() const;

  /** @brief Getter method for vortexIDS */
  std::vector<int> get_vortexIDs() const;

  /** @brief Getter method for cl */
  double get_cl() const;

  /** @brief Getter method for cd */
  double get_cd() const;

  /** @brief Getter method for cm */
  Vector3d get_cm() const;

private:
  /** @brief Method computing the area of the current element
   *  @param vortices Vortex ring elements */
  void computeArea(const std::vector<element::vortexRing> &vortices);

  /** @brief Method computing the local chord of the current element
   *  @param nodes Nodes of the mesh
   *  @param vortices Vortex ring elements */
  void computeChordLength(const std::vector<Vector3d> &nodes,
                          const std::vector<element::vortexRing> &vortices);

  friend class solver::nonlinear::steady;
  friend class solver::nonlinear::unsteady;
};

/** @brief Object holding the characteristics of a wing/lifting surface */
class wing {

  /** @brief Unique global index for the current element */
  int globalIndex;

  /** @brief IDs of the wing stations forming the current element */
  std::vector<int> stationIDs;

  /** @brief Area of the current element */
  double area;

  /** @brief Span of the wing */
  double span;

  /** @brief Lift coefficient of the current surface */
  double cl = 0.0;

  /** @brief Drag coefficient of the current surface */
  double cd = 0.0;

  /** @brief Moment coefficients of the current surface */
  Vector3d cm = Vector3d::Zero();

public:
  /** @param globalIndex Unique global index for the current element
   *  @param stationIDs IDs of the wing stations forming the current element */
  wing(const int globalIndex, const std::vector<int> &stationIDs);

  /** @brief Method initializing the geometry of the current element
   *  @param nodes Nodes of the mesh
   *  @param stations Wing station elements
   *  @param vortices Vortex ring elements
   *  @param sim Simulation parameters */
  void initialize(const std::vector<Vector3d> &nodes,
                  std::vector<wingStation> &stations,
                  std::vector<element::vortexRing> &vortices,
                  const input::simParam &sim);

  /** @brief Method updating the geometry of the current element from new moved
   * nodes
   *  @param nodes Nodes of the mesh
   *  @param stations Wing station elements
   *  @param vortices Vortex ring elements */
  void updateGeometry(const std::vector<Vector3d> &nodes,
                      std::vector<wingStation> &stations,
                      std::vector<element::vortexRing> &vortices);

  /** @brief Method computing the forces on the current surface
   *  @param stations Wing station elements
   *  @param sim Simulation parameters */
  void computeForces(const input::simParam &sim,
                     const std::vector<surface::wingStation> &stations);

  /** @brief Getter method for globalindex */
  double get_globalIndex() const;

  /** @brief Getter method for area */
  double get_area() const;

  /** @brief Getter method for span */
  double get_span() const;

  /** @brief Getter method for stationIDs */
  std::vector<int> get_stationIDs() const;

  /** @brief Getter method for cl */
  double get_cl() const;

  /** @brief Getter method for cd */
  double get_cd() const;

  /** @brief Getter method for cm */
  Vector3d get_cm() const;

private:
  void computeArea(const std::vector<wingStation> &stationIDs);
};

/** @brief Object that holds the characteristics of a patch/non-lifting surface */
class patch {
  
  /** @brief Unique global index for the current element */
  int globalIndex;

  /** @brief IDs of the doublet's panel forming to the current element */
  std::vector<int> doubletIDs;

  /** @brief Area of the current element */
  double area;

public:
  /** @param globalIndex Unique global index for the current element
   *  @param doubletIDs IDs of the wing stations forming the current element */
  patch(const int globalIndex, const std::vector<int> &doubletIDs);

  /** @brief Method initializing the geometry of the current element
   *  @param nodes Nodes of the mesh
   *  @param doublets doubletPanel of the element
   *  @param sim Simulation parameters */
  void initialize(const std::vector<Vector3d> &nodes,
                  std::vector<element::doubletPanel> &doublets,
                  const input::simParam &sim);

  /** @brief Method updating the geometry of the current element base on new
   * moved nodes (not used for future release)
   *  @param doublets doubletPanel of the element */
  void updateGeometry(const std::vector<Vector3d> &nodes,
                      std::vector<element::doubletPanel> &doublets);
  
  /** @brief Method searching the each panel neighbor 
   *  @param doublets doubletPanel of the element */
  void ScanNeighbor(std::vector<element::doubletPanel> &doublets);
  
  void Storing_nondirectPanel(std::vector<element::doubletPanel> &doublets); //not used to be removed during clean up

  /** @brief Method computing the local pressure on the current non lifting surface
   * @param doublets doubletPanel of the element
   * @param sim Simulation parameters */
  void computePressure(const input::simParam &sim,
              std::vector<element::doubletPanel> &doublets, std::vector<Vector3d> nodes);
  
  /** @brief Getter method for the globalIndex */
  double get_globalIndex() const;

  /** @brief Getter method for the area */
  double get_area() const;

  /** @brief Getter method for the doubletsIDs */
  std::vector<int> get_doubletIDs() const;

private:
  void computeArea(const std::vector<element::doubletPanel> &doubletIDs);
};

} // namespace surface
} // namespace vlm

#endif
