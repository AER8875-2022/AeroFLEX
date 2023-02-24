
#ifndef VLM_SOLVER_HPP
#define VLM_SOLVER_HPP

#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include "database/database.hpp"
#include "vlm/input.hpp"
#include "vlm/model.hpp"
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
  MatrixXd DoubletMatrixINF;

public:
  /** @brief Base constructor */
  base();

  /** @brief Constructor to initialize parameters and matrices
   *  @param solvP Solver parameters
   *  @param object Model holding the VLM elements */
  base(const input::solverParam &solvP, const model &object);

  /** @brief Pure virtual wrapper method to solve VLM algorithm
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
   *  @param gamma Vector holding the strengths of each vortex ring */
  virtual void saveSolution(model &object, const VectorXd &gamma) = 0;

  /** @brief Pure virtual function to compute global induced drag
   *  @param object Model holding the VLM elements */
  virtual void computeInducedDrag(model &object) = 0;

  /** @brief Pure virtual function to compute inviscid forces
   *  @param object Model holding the VLM elements */
  virtual void computeForces(model &object) = 0;

  /** @brief Pure virtual function to export the solution to output files
   *  @param object Model holding the VLM elements */
  virtual void exportSolution(const model &object) = 0;
};

namespace linear {

/** @brief Linear steady solver class */
class steady : public base {

public:
  /**
   *  @param solvP Solver parameters
   *  @param object Model holding the VLM elements */
  steady(const input::solverParam &solvP, const model &object);

  /** @brief Main method to solve the linear VLM problem */
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
   *  @param gamma Vector holding the strengths of each vortex ring */
  virtual void saveSolution(model &object, const VectorXd &gamma) override;

  /** @brief Method to compute global induced drag
   *  @param object Model holding the VLM elements */
  virtual void computeInducedDrag(model &object) override;

  /** @brief Method to compute inviscid forces
   *  @param object Model holding the VLM elements */
  virtual void computeForces(model &object) override;

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
  /**
   *  @param solvP solvP Solver parameters
   *  @param object Model holding the VLM elements
   *  @param database Interpolation table */
  steady(const input::solverParam &solvP, const model &object,
         const database::table database);

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

/** @brief Generates the appropriate solver based ont the input parameters
 *  @param solvP Solver parameters
 *  @param object Model holding the VLM elements
 *  @param eta Viscous slope parameter
 *  @return Pointer to appropriate solver class */
solver::base *initializeSolver(const input::solverParam &solvP,
                               const model &object,
                               const database::table &database);

} // namespace solver
} // namespace vlm

#endif
