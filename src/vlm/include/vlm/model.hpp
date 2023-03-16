
#ifndef VLM_MODEL_HPP
#define VLM_MODEL_HPP

#include "Eigen/Dense"
#include "vlm/input.hpp"
#include "vlm/panel.hpp"
#include "vlm/surface.hpp"

/** @file model.hpp */

using namespace Eigen;

namespace vlm {

/** @brief Object that holds all mesh elements and presents an interface to
 * manipulating the numerical model throughtout the simulation */
class model {
public:
  /** @brief Object holding simulation parameters */
  input::simParam sim;

  /** @brief Object holding input/output parameters */
  input::ioParam io;

  /** @brief Vector holding the nodes */
  std::vector<Vector3d> nodes;

  /** @brief Vector holding the vortex ring elements */
  std::vector<element::vortexRing> vortexRings;

  /** @brief Vector holding the doublet elements */
  std::vector<element::doubletPanel> doubletPanels;

  /** @brief Vector holding wing station characteristics */
  std::vector<surface::wingStation> wingStations;

  /** @brief Vector holding wing characteristics */
  std::vector<surface::wing> wings;

  /** @brief Vector holding patch characteristics */
  std::vector<surface::patch> patches;

  /** @brief Vector holding the automatically generated wake nodes */
  std::vector<Vector3d> wakeNodes;

  /** @brief Vector holding the automatically generated wake panels */
  std::vector<element::vortexRing> wakePanels;

private:
  /** @brief Global lift coefficient */
  double cl = 0.0;

  /** @brief Global drag coefficient */
  double cd = 0.0;

  /** @brief Global moment coefficients vector */
  Vector3d cm = Vector3d::Zero();

public:
  /** @brief Main method to to initialize the VLM model
   *  @param mesh Object holding information on the mesh
   *  @param sim Object holding simulation parameters
   *  @param io Object holding input/output parameters */
  void initialize(const input::meshData &mesh, const Settings &settings);

  /** @brief Wrapper method initializing wake elements for each lifting surface
   *  @param wakeLength Length of the generated wake in the x direction
   * (chordwise) */
  void initializeWake(const double wakeLength);

  /** @brief Method allowing the deformation of the mesh
   *  @param nodes Vector holding the nodes of the new deformed mesh */
  void updateGeometry(const std::vector<Vector3d> &nodes);

  /** @brief Method reinitializing the wake to a non existing state */
  void resetWake();

  /** @brief Method reinitializing the previous solution to zero */
  void clear();

  /** @brief Getter method for the lift coefficient */
  double get_cl() const;

  /** @brief Getter method for the drag coefficient */
  double get_cd() const;

  /** @brief Getter method for the moment coefficients vector */
  Vector3d get_cm() const;

private:
  /** @brief Method building the objects from data acquired from the mesh
   *  @param mesh Information on the mesh */
  void build(const input::meshData &mesh);

  /** @brief Method that computes the geometric properties and metrics of all
   * elements
   *  @param mesh Information on the mesh */
  void initializeMesh(const input::meshData &mesh);

  friend class solver::linear::steady;
  friend class solver::linear::unsteady;
  friend class solver::nonlinear::steady;
  friend class solver::nonlinear::unsteady;
};

} // namespace vlm

#endif
