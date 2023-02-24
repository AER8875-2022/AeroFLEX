
#include "database/database.hpp"
#include "exprtk.hpp"
#include <cmath>
#include <fstream>
#include <string>

database::airfoil::airfoil(std::vector<double> &alpha, std::vector<double> &cl,
                           std::vector<double> &cd, std::vector<double> &cmy,
                           const std::string interpolationMethod = "LAGRANGE") {
  // Selecting appropriate interpolation method
  if (!interpolationMethod.compare("LAGRANGE")) {
    this->interpolator[0] = _1D::LinearInterpolator<double>();
    this->interpolator[1] = _1D::LinearInterpolator<double>();
    this->interpolator[2] = _1D::LinearInterpolator<double>();
  } else if (!interpolationMethod.compare("SPLINE")) {
    this->interpolator[0] = _1D::CubicSplineInterpolator<double>();
    this->interpolator[1] = _1D::CubicSplineInterpolator<double>();
    this->interpolator[2] = _1D::CubicSplineInterpolator<double>();
  } else if (!interpolationMethod.compare("MONOTONIC")) {
    this->interpolator[0] = _1D::MonotonicInterpolator<double>();
    this->interpolator[1] = _1D::MonotonicInterpolator<double>();
    this->interpolator[2] = _1D::MonotonicInterpolator<double>();
  } else {
    this->interpolator[0] = _1D::LinearInterpolator<double>();
    this->interpolator[1] = _1D::LinearInterpolator<double>();
    this->interpolator[2] = _1D::LinearInterpolator<double>();
  }
  // Building interpolator
  this->interpolator[0].setData(alpha.size(), alpha.data(), cl.data());
  this->interpolator[1].setData(alpha.size(), alpha.data(), cd.data());
  this->interpolator[2].setData(alpha.size(), alpha.data(), cmy.data());
}

std::tuple<double, double, double>
database::airfoil::coefficients(const double alpha) {
  return {interpolator[0](alpha), interpolator[1](alpha),
          interpolator[2](alpha)};
}

double database::airfoil::cl(const double alpha) {
  return (interpolator[0](alpha));
}

// ---------------------------------

void database::table::importFromFile(const std::string &path,
                                     const vlm::input::solverParam &solvP) {

  std::ifstream file(path);
  if (!file.is_open()) {
    std::cerr << "\n\033[1;31m ->VLM ERROR: database file \"" << path
              << "\" not found! \033[0m" << std::endl;
    exit(1);
  }

  // Declaring output
  std::vector<std::string> line;
  const char delimiter = ' ';

  // Flags
  bool isAirfoil = false;
  bool isLocation = false;
  std::string airfoilName;

  int surfaceID;

  // coefficients & locations
  std::vector<double> alpha;
  std::vector<double> cl;
  std::vector<double> cd;
  std::vector<double> cmy;

  std::vector<std::string> sections;
  std::vector<double> spanLocs;

  std::string fileLine;
  while (std::getline(file, fileLine)) {

    // Ignoring empty lines
    if (!fileLine.size()) {
      continue;
    }

    line.clear();

    std::stringstream lineStream(fileLine);
    std::string word;

    // Splitting into an array
    while (std::getline(lineStream, word, delimiter)) {
      line.push_back(word);
    }

    // Ignoring comments
    if (!line[0].compare("#")) {
      continue;
    }

    // Airfoil declaration
    if (!line[0].compare("END") && !line[1].compare("AIRFOIL")) {
      isAirfoil = false;

      airfoils.insert(
          {airfoilName, airfoil(alpha, cl, cd, cmy, solvP.interpolation)});

      alpha.clear();
      cl.clear();
      cd.clear();
      cmy.clear();
      airfoilName.clear();
      continue;
    }
    if (isAirfoil) {
      alpha.push_back(std::stod(line[0]));
      cl.push_back(std::stod(line[1]));
      cd.push_back(std::stod(line[2]));
      cmy.push_back(std::stod(line[3]));
    }
    if (!line[0].compare("AIRFOIL")) {
      isAirfoil = true;
      airfoilName = line[1];
      continue;
    }

    // Location declaration
    if (!line[0].compare("END") && !line[1].compare("SURFACE")) {
      isLocation = false;

      sectionAirfoils.insert({surfaceID, sections});
      sectionSpanLocs.insert({surfaceID, spanLocs});

      sections.clear();
      spanLocs.clear();
      continue;
    }
    if (isLocation) {
      sections.push_back(line[0]);
      spanLocs.push_back(std::stod(line[1]));
    }
    if (!line[0].compare("SURFACE")) {
      isLocation = true;
      surfaceID = std::stoi(line[1]);
      continue;
    }
  }

  file.close();
}

void database::table::generateFromPolar(const std::string &polar,
                                        const vlm::model &object,
                                        const vlm::input::solverParam &solvP) {
  // Database data
  std::vector<double> alpha;
  std::vector<double> cl;
  std::vector<double> cd;
  std::vector<double> cmy;

  // Defining symbols
  exprtk::symbol_table<double> symbols;
  double x = 0.0;
  symbols.add_variable("x", x);
  symbols.add_constants();

  // Defining expression
  exprtk::expression<double> expression;
  expression.register_symbol_table(symbols);

  // Parsing expression
  exprtk::parser<double> parser;
  parser.compile(polar, expression);

  // Generating table with multiple aoa
  for (x = -15.; x <= 50.; x += 0.1) {
    const double cli = expression.value();
    alpha.push_back(x);
    cl.push_back(cli);
    cd.push_back(0);
    cmy.push_back(0);
  }
  airfoils.insert({"EQ", airfoil(alpha, cl, cd, cmy, solvP.interpolation)});

  // Generating evaluation points
  std::vector<std::string> sections = {"EQ", "EQ"};
  for (auto &wing : object.wings) {
    sectionAirfoils.insert({wing.get_globalIndex(), sections});
    sectionSpanLocs.insert({wing.get_globalIndex(), {0., 1.}});
  }
}

std::tuple<double, double, double>
database::table::coefficients(const double alpha, const double surfaceID,
                              const double spanLoc) {
  // Interpolating coefficients for a given alpha
  auto &sections = sectionAirfoils.at(surfaceID);
  auto &spanLocs = sectionSpanLocs.at(surfaceID);
  std::array<std::vector<double>, 3> coeff;
  coeff[0].reserve(spanLocs.size());
  coeff[1].reserve(spanLocs.size());
  coeff[2].reserve(spanLocs.size());
  for (auto &section : sections) {
    auto &airfoil = airfoils.at(section);
    auto [cl, cd, cmy] = airfoil.coefficients(alpha);
    coeff[0].push_back(cl);
    coeff[1].push_back(cd);
    coeff[2].push_back(cmy);
  }
  // Setting up spanwise interpolator
  std::array<_1D::LinearInterpolator<double>, 3> interpolator;
  interpolator[0].setData(spanLocs, coeff[0]);
  interpolator[1].setData(spanLocs, coeff[1]);
  interpolator[2].setData(spanLocs, coeff[2]);

  return {interpolator[0](spanLoc), interpolator[1](spanLoc),
          interpolator[2](spanLoc)};
}

double database::table::cl(const double alpha, const double surfaceID,
                           const double spanLoc) {
  // Interpolating coefficient for a given alpha
  auto &sections = sectionAirfoils.at(surfaceID);
  auto &spanLocs = sectionSpanLocs.at(surfaceID);
  std::vector<double> coeff;
  coeff.reserve(spanLocs.size());
  for (auto &section : sections) {
    auto &airfoil = airfoils.at(section);
    auto cl = airfoil.cl(alpha);
    coeff.push_back(cl);
  }
  // Setting up spanwise interpolator
  _1D::LinearInterpolator<double> interpolator;
  interpolator.setData(spanLocs, coeff);

  return (interpolator(spanLoc));
}
