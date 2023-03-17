
#ifndef VLM_INPUT_HPP
#define VLM_INPUT_HPP

#include "Eigen/Dense"
#include <map>
#include <string>
#include <tuple>
#include "tinyconfig.hpp"

#ifndef M_PI
#define M_PI 3.141592653589793115997963468544185161590576171875
#endif

/** @file input.hpp */

using namespace Eigen;

namespace vlm {

namespace input {

/** @brief Object holding case and physics-oriented parameters */
struct simParam {

  /** @brief Geometric angle of attack in degrees */
  double aoa = 5.0;

  /** @brief Geometric angle of side slip in degrees */
  double sideslip = 0.0;

  /** @brief Free stream magnitude velocity */
  double vinf = 1.0;

  /** @brief Density of the fluid */
  double rho = 1.0;

  /** @brief Reference chord length */
  double cref = 1.0;

  /** @brief Reference surface area */
  double sref = 32.0;

  /** @brief X component of origin to which the x and z moment are computed */
  double x0 = 0.0;

  /** @brief Y component of origin to which the x and z moment are computed */
  double y0 = 0.0;

  /** @brief Z component of origin to which the x and z moment are computed */
  double z0 = 0.0;

  /** @brief Viscous relaxation value applied on the vortex filament kernel */
  double coreRadius = 0.0;

  /** @brief Input options format of the database */
  std::vector<std::string> databaseFormat_options = {"NONE", "FILE"};

  /** @brief Input format of the database */
  int databaseFormat = 0;

  /** @brief Setter for databaseFormat */
  void set_databaseFormat(const std::string &item);

  /** @brief Getter for databaseFormat */
  std::string get_databaseFormat() const;

  /** @brief Method computing the free stream vector in the plane's referential
   *  @return Free stream vector */
  Vector3d freeStream() const;

  /** @brief Method computing the free stream vector in the plane's referential
   * and according to a local angle of attack
   *  @param alpha Local angle of attack
   *  @return Local free stream vector */
  Vector3d freeStream(const double alpha) const;

  /** @brief Method computing the stream flow orientation
   *  @return Stream axis vector */
  Vector3d streamAxis() const;

  /** @brief Method computing the stream flow orientation
   *  @param alpha Local angle of attack
   *  @return Stream axis vector */
  Vector3d streamAxis(const double alpha) const;

  /** @brief Method computing the lift axis
   *  @return Lift axis vector */
  Vector3d liftAxis() const;

  /** @brief Method computing the lift axis
   *  @param alpha Local angle of attack
   *  @return Lift axis vector */
  Vector3d liftAxis(const double alpha) const;

  /** @brief Method computing the dynamic pressure
   *  @return Dynamic pressure value */
  double dynamicPressure() const;

  /** @brief Method that returns the origin to which the moments are calculated
   *  @return Origin vector */
  Vector3d origin() const;

};

/** @brief Object holding input/output information */
struct ioParam {

  /** @brief Prefix that will be applied to all output files */
  std::string baseName = "vlm";

  /** @brief Path of the directory to which the solution will be outputted */
  std::string outDir = "vlm_out/";

  /** @brief Path to the mesh file */
  std::string meshFile = "../../../../examples/vlm/smith_wing.dat";

  /** @brief Path to the database file (if required) */
  std::string databaseFile = "../../../../examples/vlm/database_stall.dat";

  /** @brief Path to the location file for database evalutation */
  std::string locationFile = "../../../../examples/vlm/locations.dat";
};

/** @brief Object holding solver parameters */
struct solverParam {

  /** @brief Time domain options of the simulation */
  std::vector<std::string> timeDomain_options = {"STEADY", "UNSTEADY"};

  /** @brief Time domain of the simulation (steady or time-dependent) */
  int timeDomain = 0;

  /** @brief Setter for timeDomain */
  void set_timeDomain(const std::string &item);

  /** @brief Getter for timeDomain */
  std::string get_timeDomain() const;

  /** @brief Type of solver options utilised (linear or non linear) */
  std::vector<std::string> type_options = {"LINEAR", "NONLINEAR"};

  /** @brief Type of solver utilised (linear or non linear) */
  int type = 1;

  /** @brief Setter for type */
  void set_type(const std::string &item);

  /** @brief Getter for type */
  std::string get_type() const;

  /** @brief Tolerance of the non linear iterative solver */
  double tolerance = 1e-14;

  /** @brief Linear solver options type */
  std::vector<std::string> linearSolver_options = {"BICGSTAB", "DIRECT"};

  /** @brief Linear solver type */
  int linearSolver = 0;

  /** @brief Setter for linearSolver */
  void set_linearSolver(const std::string &item);

  /** @brief Getter for linearSolver */
  std::string get_linearSolver() const;

  /** @brief Relaxation for the iterative scheme */
  double relaxation = 1.0;

  /** @brief Max number of iterations for the iterative scheme */
  int max_iter = 200;
};

/** @brief Object holding the connectivity obtained from the mesh file */
struct meshData {

  /** @brief Map holding the nodes */
  std::map<unsigned int, Vector3d> nodes;

  /** @brief Map holding the IDs of the nodes belonging to a vortex ring*/
  std::map<unsigned int, std::vector<int>> vortexIDs;

  /** @brief Map holding the IDs of the nodes belonging to a doublet panel */
  std::map<unsigned int, std::vector<int>> doubletIDs;

  /** @brief Map holding the IDs of the vortex rings belonging to a wing station
   */
  std::map<unsigned int, std::vector<int>> stationIDs;

  /** @brief Map holding the IDs of the wing stations belonging to a wing
   * surface */
  std::map<unsigned int, std::vector<int>> wingIDs;

  /** @brief Map holding the IDs of the doublet panels belonging to a patch
   * surface */
  std::map<unsigned int, std::vector<int>> patchIDs;
};

/** @brief Function wrapper that import connectivity data specified in the mesh
 * file
 *  @param names input/output parameters object
 *  @return meshData object holding all the information on the current mesh */
meshData importMeshFile(const ioParam &names);

} // namespace input

/** @brief Main struct holding parameters on the current case */
struct Settings {

  /** @brief Object holding case and physics-oriented parameters */
  input::simParam sim;

  /** @brief Object holding input/output information */
  input::ioParam io;

  /** @brief Object holding solver parameters */
  input::solverParam solver;

  /** @brief Function to import parameters from config file */
  void import_config_file(tiny::config &config);

  /** @brief Function to export parameters to config file */
  void export_config_file(tiny::config &config);
};

} // namespace vlm

#endif
