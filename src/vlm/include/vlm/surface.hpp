
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

#ifndef VLM_SURFACE_CLASS
#define VLM_SURFACE_CLASS
namespace vlm::surface {
class wingStation;
class wing;
class patch;
} // namespace vlm::surface
#endif

namespace vlm {

namespace surface {

/** @brief Object that holds the characteristics of a wing station/section */
class wingStation {

  /** @brief Unique global index of the current element */
  int globalIndex;

  /** @brief IDs of the vortices forming to the current element */
  std::vector<int> vortexIDs;

  /** @brief Reference to vortex elements */
  std::vector<element::vortexRing> &vortices;

  /** @brief Area of the current element */
  double area;

  /** @brief Local chord of the current wing station */
  double chord;

  /** @brief Force vector */
  Vector3d force;

  /** @rief Moment vector */
  Vector3d moment;

  /** @brief Lift coefficient */
  double cl = 0.0;

  /** @brief Drag coefficient */
  double cd = 0.0;

  /** @brief Side force coefficient */
  double cy = 0.0;

  /** @brief Moment coefficients */
  Vector3d cm = Vector3d::Zero();

  /** @brief Local lift coefficient */
  double cl_local = 0.0;

  /** @brief local angle of attack correction */
  double local_aoa = 0.0;

  /** @brief local flow stream */
  Vector3d local_stream;

  /** @brief Spanwise location of the current wing station */
  double spanLoc;

public:
  /** @param globalIndex Unique global index of the current element
   *  @param vortexIDs IDs of the vortices forming the current element
   *  @param vortices VortexRings of the mesh */
  wingStation(const int globalIndex, const std::vector<int> &vortexIDs,
              std::vector<element::vortexRing> &vortices);

  /** @brief Method initializing the current element
   *  @param sim Simulation parameters */
  void initialize(const input::simParam &sim);

  /** @brief Method to update the geometric characteristics of the surface */
  void updateGeometry();

  /** @brief Method generating a wake panel for the current wing station
   *  @param wakeLength Length of the generated wake panel
   *  @param sim Simulation parameters
   *  @param wakeNodes To be generated wake nodes
   *  @param wakePanels To be generated wake panels */
  void generateWake(const double wakeLength, const input::simParam &sim,
                    std::vector<Vector3d> &wakeNodes,
                    std::vector<element::vortexRing> &wakePanels);

  /** @brief Method computing the coordinates where to forces are acting
   *  @return Coordinates of the force acting point */
  Vector3d forceActingPoint() const;

  /** @brief Method updating the local angle of attack for the current wing
   * station
   *  @param dalpha Variation of the local angle of attack */
  void updateLocalStream(const double dalpha, const input::simParam &sim);

  /** @brief Method to reset the local angle of attack to the geometric one
   *  @param sim Simulation parameters */
  void resetLocalStream(const input::simParam &sim);

  /** @brief Function that returns the local lift axis of the section */
  Vector3d liftAxis();

  /** @brief Method computing the inviscid local forces on the current element
   */
  void computeForces(const input::simParam &sim);

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

  /** @brief Getter method for cy */
  double get_cy() const;

  /** @brief Getter method for cm */
  Vector3d get_cm() const;

private:
  /** @brief Method computing the area of the current element */
  void computeArea();

  /** @brief Method computing the local chord of the current element */
  void computeChordLength();

  /** @brief Method to obtain a vector in the local referential of the section
   */
  inline void to_local(Vector3d &vector);

  friend class wing;
  friend class solver::nonlinear::steady;
  friend class solver::nonlinear::unsteady;
};

/** @brief Object holding the characteristics of a wing/lifting surface */
class wing {

  /** @brief Unique global index for the current element */
  int globalIndex;

  /** @brief IDs of the wing stations forming the current element */
  std::vector<int> stationIDs;

  /** @brief Reference to wing stations */
  std::vector<wingStation> &stations;

  /** @brief Area of the current element */
  double area;

  /** @brief Span of the wing */
  double span;

  /** @brief Lift coefficient of the current surface */
  double cl = 0.0;

  /** @brief Drag coefficient of the current surface */
  double cd = 0.0;

  /** @brief Side force coefficient of the current surface */
  double cy = 0.0;

  /** @brief Moment coefficients of the current surface */
  Vector3d cm = Vector3d::Zero();

public:
  /** @param globalIndex Unique global index for the current element
   *  @param stationIDs IDs of the wing stations forming the current element
   *  @param stations Wing station of the mesh */
  wing(const int globalIndex, const std::vector<int> &stationIDs,
       std::vector<wingStation> &stations);

  /** @brief Method initializing the geometry of the current element
   *  @param sim Simulation parameters */
  void initialize(const input::simParam &sim);

  /** @brief Method to update the geometric characteristics of the surface */
  void updateGeometry();

  /** @brief Method computing the forces on the current surface
   *  @param sim Simulation parameters */
  void computeForces(const input::simParam &sim);

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

  /** @brief Getter method for cy */
  double get_cy() const;

  /** @brief Getter method for cm */
  Vector3d get_cm() const;

private:
  /** @brief Method that computes the area of the wing surface */
  void computeArea();

  /** @brief Method that computes the wing's span as well as the span location
   * of every of its wing stations */
  void computeSpan();
};

class patch {
  int globalIndex;
  std::vector<int> doubletIDs;
  std::vector<element::doubletPanel> &doubletPanels;
  double area;

public:
  patch(const int globalIndex, const std::vector<int> &doubletIDs,
        std::vector<element::doubletPanel> &doubletPanels);
  void initialize(const input::simParam &sim);
  void updateGeometry();
  double get_globalIndex() const;
  double get_area() const;
  std::vector<int> get_doubletIDs() const;

private:
  void computeArea();
};

} // namespace surface
} // namespace vlm

#endif
