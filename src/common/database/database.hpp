
#ifndef VISCOUS_DATABASE_HPP
#define VISCOUS_DATABASE_HPP

#include "libInterpolate/AnyInterpolator.hpp"
#include "libInterpolate/Interpolate.hpp"
#include "vlm/input.hpp"
#include "vlm/model.hpp"
#include <array>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#ifndef VLM_CL_ALPHA_DEG
/** @brief Value of the slope of the CL vs angle of attack function in 1/deg */
#define VLM_CL_ALPHA_DEG                                                       \
  0.109662271123215082635482531259185634553432464599609375
#endif

/** @file database.hpp */

namespace database {

/** @brief Bidimensional table that allows interpolation of viscous force
 * coefficients of an specific airfoil */
struct airfoil {

  /** @brief Array of 3 interpolators; One for each coefficient (0: CL, 1: CD,
   * 2: CMy) */
  std::array<_1D::AnyInterpolator<double>, 3> interpolator;

  /**
   *  @param alpha Vector holding the values of angle of attack at which the
   * viscous forces are evaluated
   *  @param cl Vector holding the values of the viscous lift coefficient
   *  @param cd Vector holding the values of the viscous drag coefficient
   *  @param cmy Vector holding the values of the viscous moment coefficient
   * according to the Y axis (spanwise)
   *  @param interpolationMethod Interpolation method based on the values of
   * interpolationType
   * */
  airfoil(std::vector<double> &alpha, std::vector<double> &cl,
          std::vector<double> &cd, std::vector<double> &cmy,
          const std::string interpolationMethod);

  /** @brief Interpolates the 3 viscous coefficients at specified angle of
   * attack
   *  @param alpha Angle of attack at which the forces are interpolated
   *  @return Tuple holding the values of all 3 interpolated coefficients (0:
   * CL, 1: CD, 2: CMy) */
  std::tuple<double, double, double> coefficients(const double alpha);

  /** @brief Interpolates the value of the lift coefficient at specified angle
   * of attack
   *  @param alpha Angle of attack at which the forces are interpolated
   *  @return Interpolated lift coefficient */
  double cl(const double alpha);
};

/** @brief Tridimensional table that allows the interpolaton of viscous force
 * coefficients on a specified surface and at a specifed angle of attack */
struct table {

  /** @brief Map holding the characteristics (bidimensional table) of a specific
   * airfoil. The key of the map corresponds to the ID of a wing */
  std::unordered_map<std::string, database::airfoil> airfoils;

  /** @brief Map holding the name of the airfoil at which the viscous forces are
   * evaluated with the RANS/Euler model. The key of the map corresponds to the
   * ID of a wing */
  std::unordered_map<int, std::vector<std::string>> sectionAirfoils;

  /** @brief Map holding the spanwise locations at which the viscous forces are
   * evaluated with the RANS/Euler model. The key of the map corresponds to the
   * ID of a wing */
  std::unordered_map<int, std::vector<double>> sectionSpanLocs;

  /** @brief Import database from a file
   *  @param path Path to the file containing the database
   *  @param solvP Struct holding parameters for the solver */
  void importFromFile(const std::string &path,
                      const vlm::input::solverParam &solvP);

  /** @brief Generate a database from a lift polar
   *  @param polar Lift polar equation as a string
   *  @param object Object holding all elements of the mesh
   *  @param solvP Struct holding parameters for the solver */
  void generateFromPolar(const std::string &polar, const vlm::model &object,
                         const vlm::input::solverParam &solvP);

  /** @brief Interpolates the 3 viscous coefficients at specified angle of
   * attack, surface and spanwise location
   *  @param alpha Angle of attack at which the forces are interpolated
   *  @param surfaceID ID of the wing surface on which the forces are evaluated
   *  @param spanLoc spanwise location where the forces are interpolated
   *  @return Tuple holding the values of all 3 interpolated coefficients (0:
   * CL, 1: CD, 2: CMy) */
  std::tuple<double, double, double> coefficients(const double alpha,
                                                  const double surfaceID,
                                                  const double spanLoc);

  /** @brief Interpolates the value of the lift coefficient at specified angle
   * of attack, surface and spanwise location
   *  @param alpha Angle of attack at which the forces are interpolated
   *  @param surfaceID ID of the wing surface on which the forces are evaluated
   *  @param spanLoc spanwise location where the forces are interpolated
   *  @return Interpolated lift coefficient */
  double cl(const double alpha, const double surfaceID, const double spanLoc);
};

} // namespace database

#endif
