
#include "database/database.hpp"
#include "mlinterp/mlinterp.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

database::airfoil::airfoil(std::vector<double> &alpha, std::vector<double> &cl,
                           std::vector<double> &cd, std::vector<double> &cmy)
    : alpha(alpha), cl(cl), cd(cd), cmy(cmy) {}

std::tuple<double, double, double>
database::airfoil::interpolate_coeff(const double alpha) {
  const int n_d[] = {int(this->alpha.size())};
  const double alpha_i[1] = {alpha};

  double cl_i[1], cd_i[1], cmy_i[1];

  mlinterp::interp(n_d, 1, &cl[0], cl_i, &this->alpha[0], alpha_i);
  mlinterp::interp(n_d, 1, &cd[0], cd_i, &this->alpha[0], alpha_i);
  mlinterp::interp(n_d, 1, &cmy[0], cmy_i, &this->alpha[0], alpha_i);

  return {cl_i[0], cd_i[0], cmy_i[0]};
}

double database::airfoil::interpolate_cl(const double alpha) {
  const int n_d[] = {int(this->alpha.size())};
  const double alpha_i[1] = {alpha};

  double cl_i[1];

  mlinterp::interp(n_d, 1, &cl[0], cl_i, &this->alpha[0], alpha_i);
  return (cl_i[0]);
}

// ---------------------------------

void database::table::importFromFile(const std::string &path,
                                     const vlm::input::solverParam &solvP) {

  std::ifstream file(path);
  if (!file.is_open()) {
    std::cerr << "\n\033[1;31m ->VLM ERROR: database file \"" << path
              << "\" not found! \033[0m" << std::endl;
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

      airfoils.insert({airfoilName, airfoil(alpha, cl, cd, cmy)});

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

std::tuple<double, double, double>
database::table::coefficients(const double alpha, const int surfaceID,
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
    auto [cl, cd, cmy] = airfoil.interpolate_coeff(alpha);
    coeff[0].push_back(cl);
    coeff[1].push_back(cd);
    coeff[2].push_back(cmy);
  }

  const int n_d[] = {int(spanLocs.size())};
  const double spanLoc_i[1] = {spanLoc};

  double cl_i[1], cd_i[1], cmy_i[1];

  mlinterp::interp(n_d, 1, &coeff[0][0], cl_i, &spanLocs[0], spanLoc_i);
  mlinterp::interp(n_d, 1, &coeff[1][0], cd_i, &spanLocs[0], spanLoc_i);
  mlinterp::interp(n_d, 1, &coeff[2][0], cmy_i, &spanLocs[0], spanLoc_i);

  return {cl_i[0], cd_i[0], cmy_i[0]};
}

double database::table::cl(const double alpha, const int surfaceID,
                           const double spanLoc) {
  // Interpolating coefficient for a given alpha
  auto &sections = sectionAirfoils.at(surfaceID);
  auto &spanLocs = sectionSpanLocs.at(surfaceID);
  std::vector<double> coeff;
  coeff.reserve(spanLocs.size());
  for (auto &section : sections) {
    auto &airfoil = airfoils.at(section);
    auto cl = airfoil.interpolate_cl(alpha);
    coeff.push_back(cl);
  }

  const int n_d[] = {int(spanLocs.size())};
  const double spanLoc_i[1] = {spanLoc};

  double cl_i[1];

  mlinterp::interp(n_d, 1, &coeff[0], cl_i, &spanLocs[0], spanLoc_i);

  return (cl_i[0]);
}
