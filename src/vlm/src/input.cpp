
#include "vlm/input.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <exception>

using namespace vlm;
using namespace input;
using namespace Eigen;

Vector3d simParam::freeStream() const {
  return Vector3d(
      {vinf * std::cos(aoa * M_PI / 180.0) * std::cos(sideslip * M_PI / 180.0),
       vinf * std::cos(aoa * M_PI / 180.0) * std::sin(sideslip * M_PI / 180.0),
       vinf * std::sin(aoa * M_PI / 180.0)});
}

Vector3d simParam::freeStream(const double alpha) const {
  return Vector3d({vinf * std::cos(alpha * M_PI / 180.0) *
                       std::cos(sideslip * M_PI / 180.0),
                   vinf * std::cos(alpha * M_PI / 180.0) *
                       std::sin(sideslip * M_PI / 180.0),
                   vinf * std::sin(alpha * M_PI / 180.0)});
}

Vector3d simParam::streamAxis() const { return (freeStream().normalized()); }

Vector3d simParam::streamAxis(const double alpha) const {
  return (freeStream(alpha).normalized());
}

Vector3d simParam::liftAxis() const {
  return (freeStream().cross(Vector3d::UnitY()).normalized());
}

Vector3d simParam::liftAxis(const double alpha) const {
  return (freeStream(alpha).cross(Vector3d::UnitY()).normalized());
}

double simParam::dynamicPressure() const { return (0.5 * rho * vinf * vinf); }

// -------------------------------------

meshData input::importMeshFile(const ioParam &names) {
  // Opening file and verifying if already exists
  std::ifstream file(names.meshFile);
  if (!file.is_open()) {
    std::cerr << "\n\033[1;31m ->VLM ERROR: mesh file \"" << names.meshFile
              << "\" not found! \033[0m" << std::endl;
  }

  // Mesh object
  meshData mesh;

  // Declaring output
  std::vector<std::string> line;
  const char delimiter = ' ';

  // Starting file reading
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

    // BEGINING LINE PARSING
    if (!line[0].compare("NODE")) {
      int nodeID = std::stoi(line[1]);
      Vector3d coord = {std::stod(line[2]), std::stod(line[3]),
                        std::stod(line[4])};
      // Adding to corresponding map
      mesh.nodes[nodeID] = coord;

    } else if (!line[0].compare("VORTEX")) {
      int vortexID = std::stoi(line[1]);
      std::vector<int> nodeIDs;
      for (size_t i = 2; i != line.size(); i++) {
        nodeIDs.push_back(std::stoi(line[i]));
      }
      // Adding to corresponding map
      mesh.vortexIDs[vortexID] = nodeIDs;

    } else if (!line[0].compare("DOUBLET")) {
      int doubletID = std::stoi(line[1]);
      std::vector<int> nodeIDs;
      for (size_t i = 2; i != line.size(); i++) {
        nodeIDs.push_back(std::stoi(line[i]));
      }
      // Adding to corresponding map
      mesh.doubletIDs[doubletID] = nodeIDs;

    } else if (!line[0].compare("WINGSTATION")) {
      int stationID = std::stoi(line[1]);
      std::vector<int> vortexIDs;
      for (size_t i = 2; i != line.size(); i++) {
        vortexIDs.push_back(std::stoi(line[i]));
      }
      // Adding to corresponding map
      mesh.stationIDs[stationID] = vortexIDs;

    } else if (!line[0].compare("WING") || !line[0].compare("LIFTINGSURFACE")) {
      int wingID = std::stoi(line[1]);
      std::vector<int> stationIDs;
      for (size_t i = 2; i != line.size(); i++) {
        stationIDs.push_back(std::stoi(line[i]));
      }
      // Adding to corresponding map
      mesh.wingIDs[wingID] = stationIDs;

    } else if (!line[0].compare("PATCH") ||
               !line[0].compare("NONLIFTINGSURFACE")) {
      int patchID = std::stoi(line[1]);
      std::vector<int> doubletIDs;
      for (size_t i = 2; i != line.size(); i++) {
        doubletIDs.push_back(std::stoi(line[i]));
      }
      // Adding to corresponding map
      mesh.patchIDs[patchID] = doubletIDs;
    }
  }
  file.close();
  // Looking for errors in mesh file extract
  return (mesh);
}

// ---------------------------

void Settings::import_config_file(tiny::config &config) {

  try {

  // [SIMULATION]
  sim.aoa = config.get<double>("vlm-simulation", "aoa", 5.0);
  sim.sideslip = config.get<double>("vlm-simulation", "sideslip", 0.0);
  sim.vinf = config.get<double>("vlm-simulation", "v_inf", 1.0);
  sim.rho = config.get<double>("vlm-simulation", "density", 1.0);
  sim.cref = config.get<double>("vlm-simulation", "c_ref", 1.0);
  sim.sref = config.get<double>("vlm-simulation", "s_ref", 1.0);
  auto xref = config.get<double>("vlm-simulation", "x_ref", 0.0);
  auto yref = config.get<double>("vlm-simulation", "y_ref", 0.0);
  auto zref = config.get<double>("vlm-simulation", "z_ref", 0.0);
  sim.origin = {xref, yref, zref};
  sim.coreRadius = config.get<double>("vlm-simulation", "lamb-oseen_radius", 0.0);
  sim.databaseFormat = config.get<std::string>("vlm-simulation", "database_format");

  // [IO]
  io.baseName = config.get<std::string>("vlm-io", "basename");
  io.outDir = config.get<std::string>("vlm-io", "output_dir");
  io.meshFile =
      config.get<std::string>("vlm-io", "mesh_file");
  io.databaseFile =
      config.get<std::string>("vlm-io", "database_file");
  io.locationFile =
      config.get<std::string>("vlm-io", "location_file");

  // [SOLVER]
  solver.timeDomain =
      config.get<std::string>("vlm-solver", "time_domain");
  solver.type = config.get<std::string>("vlm-solver", "type");
  solver.tolerance = config.get<double>("vlm-solver", "tolerance", 1e-15);
  solver.interpolation = config.get<std::string>("vlm-solver", "interpolation_type");
  solver.linearSolver =
      config.get<std::string>("vlm-solver", "linear_solver");
  solver.relaxation = config.get<double>("vlm-solver", "relaxation", 1.0);
  solver.max_iter = config.get<int>("vlm-solver", "max_iter", 100);

  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
  }

}
