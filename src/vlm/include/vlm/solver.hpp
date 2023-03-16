
#ifndef VLM_SOLVER_HPP
#define VLM_SOLVER_HPP

#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include "database/database.hpp"
#include "vlm/input.hpp"
#include "vlm/model.hpp"
#include <atomic>
#include <unordered_map>

#ifndef VLM_CL_ALPHA_DEG
#define VLM_CL_ALPHA_DEG                                                       \
  0.109662271123215082635482531259185634553432464599609375
#endif

using namespace Eigen;

/** @file solver.hpp */

namespace vlm {

namespace solver {

/** @brief Base solver class that defines a base for steady, unsteady and
 * linear/non-linear solvers */
class base {

protected:
  /** @brief Solver parameters */
  input::solverParam solvP;

  /** @brief Left hand side of the linear system */
  MatrixXd lhs;

  /** @brief Right hand side of the linear system */
  VectorXd rhs;

  /** @brief Source influence matrix */
  MatrixXd DoubletMatrixINF;

public:
  /** @brief Method to initialize the solver object
   *  @param solvP Struct holding the solver parameters
   *  @param object Model holding the VLM elements */
  virtual void initialize(const input::solverParam &solvP, const model &object,
                          const database::table &) = 0;

  /** @brief Main method to solve non linear VLM
   *  @param object Model holding the VLM elements */
  virtual void solve(model &object) = 0;

protected:
  /** @brief Pure virtual function to build the left hand side
   *  @param object Model holding the VLM elements */
  virtual void buildLHS(const model &object) = 0;

  /** @brief Pure virtual function to build the right hand side
   *  @param object Model holding the VLM elements */
  virtual void buildRHS(const model &object) = 0;

  /** @brief Pure virtual function to save solution into the model
   *  @param object Model holding the VLM elements
   *  @param X Vector holding the strengths of each vortex ring(gamma) and strengths of each doublets*/
  virtual void saveSolution(model &object, const VectorXd &X) = 0;

  /** @brief Pure virtual function to compute global induced drag
   *  @param object Model holding the VLM elements */
  virtual void computeInducedDrag(model &object) = 0;

  /** @brief Pure virtual function to compute inviscid forces
   *  @param object Model holding the VLM elements */
  virtual void computeForces(model &object) = 0;

  virtual void computePressure(model &object) = 0;

  /** @brief Pure virtual function to export the solution to output files
   *  @param object Model holding the VLM elements */
  virtual void exportSolution(const model &object) = 0;
};

namespace linear {

/** @brief Linear steady solver class */
class steady : public base {

public:
  /** @brief Iteration counter */
  std::atomic<int> &iter;

  /** @brief Residual vector */
  std::vector<double> &residuals;

public:
  steady(std::atomic<int> &iter, std::vector<double> &residuals);

  /** @brief Method to initialize the solver object
   *  @param solvP Struct holding the solver parameters
   *  @param object Model holding the VLM elements */
  virtual void initialize(const input::solverParam &solvP, const model &object,
                          const database::table &) override;

  /** @brief Main method to solve non linear VLM
   *  @param object Model holding the VLM elements */
  virtual void solve(model &object) override;

protected:
  /** @brief Method to build the left hand side of the linear system
   *  @param object Model holding the VLM elements */
  virtual void buildLHS(const model &object) override;

  /** @brief Method to build the right hand side of the linear system
   *  @param object Model holding the VLM elements */
  virtual void buildRHS(const model &object) override;

  /** @brief Method to save solution into the model
   *  @param object Model holding the VLM elements
   *  @param X Vector holding the strengths of each vortex ring(gamma) and strengths of each doublets*/
  virtual void saveSolution(model &object, const VectorXd &X) override;

  /** @brief Method to compute global induced drag
   *  @param object Model holding the VLM elements */
  virtual void computeInducedDrag(model &object) override;

  /** @brief Method to compute inviscid forces
   *  @param object Model holding the VLM elements */
  virtual void computeForces(model &object) override;

  virtual void computePressure(model &object) override;

  /** @brief Method to export the solution to output files
   *  @param object Model holding the VLM elements */
  virtual void exportSolution(const model &object) override;
};

/** @brief Linear unsteady solver */
class unsteady : public linear::steady {};

} // namespace linear

namespace nonlinear {

/** @brief Non linear steady solver */
class steady : public linear::steady {

  /** @brief Database/table that serves as an interpolation support for viscous
   * force correction */
  database::table database;

public:
  steady(std::atomic<int> &iter, std::vector<double> &residuals);

  /** @brief Method to initialize the solver object
   *  @param solvP Struct holding the solver parameters
   *  @param object Model holding the VLM elements */
  virtual void initialize(const input::solverParam &solvP, const model &object,
                          const database::table &database) override;

  /** @brief Main method to solve non linear VLM
   *  @param object Model holding the VLM elements */
  virtual void solve(model &object) override;

protected:
  /** @brief Method to build the right hand side of the linear system
   *  @param object Model holding the VLM elements */
  virtual void buildRHS(const model &object) override;

  /** @brief Method executing one iteration of angle of attack correction
   *  @param object Model holding the VLM elements
   *  @return Residual computed on the inviscid and viscous lift coefficients */
  virtual double iterate(model &object);

  /** @brief Method computing inviscid lift only
   *  @param object Model holding the VLM elements */
  virtual void iterateLift(model &object);

  /** @brief Method computing viscous forces
   *  @param object Model holding the VLM elements */
  virtual void computeForces(model &object) override;
};

/** @brief Non linear unsteady solver */
class unsteady : public linear::unsteady {};

} // namespace nonlinear
} // namespace solver
} // namespace vlm

#endif
