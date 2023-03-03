
#include "vlm/input.hpp"
#include "ini/ini.h"
#include <cmath>

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

Vector3d simParam::streamAxis() const {
  return (freeStream().normalized());
}

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

std::tuple<ioParam, simParam, solverParam>
input::importConfigFile(const std::string path) {
  // Looking for config file
  inih::INIReader file = inih::INIReader(path);

  // [SIMULATION]
  simParam sim;
  sim.aoa = file.Get<double>("SIMULATION", "AOA", 5.0);
  sim.sideslip = file.Get<double>("SIMULATION", "SIDESLIP", 0.0);
  sim.vinf = file.Get<double>("SIMULATION", "V_INF", 1.0);
  sim.rho = file.Get<double>("SIMULATION", "DENSITY", 1.0);
  sim.cref = file.Get<double>("SIMULATION", "C_REF", 1.0);
  sim.sref = file.Get<double>("SIMULATION", "S_REF", 1.0);
  auto xref = file.Get<double>("SIMULATION", "X_REF", 0.0);
  auto yref = file.Get<double>("SIMULATION", "Y_REF", 0.0);
  auto zref = file.Get<double>("SIMULATION", "Z_REF", 0.0);
  sim.origin = {xref, yref, zref};
  sim.coreRadius = file.Get<double>("SIMULATION", "LAMB-OSEEN_RADIUS", 0.0);
  sim.databaseFormat = file.Get<std::string>("SIMULATION", "DATABASE_FORMAT",
                                             std::string("NONE"));

  // [IO]
  ioParam io;
  io.baseName = file.Get<std::string>("IO", "BASENAME", std::string("vlm"));
  io.outDir = file.Get<std::string>("IO", "OUTPUT_DIR", std::string("vlm_out"));
  io.meshFile =
      file.Get<std::string>("IO", "MESH_FILE", std::string("mesh.dat"));
  io.databaseFile =
      file.Get<std::string>("IO", "DATABASE_FILE", std::string("database.dat"));

  // [SOLVER]
  solverParam solvP;
  solvP.timeDomain =
      file.Get<std::string>("SOLVER", "TIME_DOMAIN", std::string("STEADY"));
  solvP.type = file.Get<std::string>("SOLVER", "TYPE", std::string("LINEAR"));
  solvP.tolerance = file.Get<double>("SOLVER", "TOLERANCE", 1e-15);
  solvP.interpolation = file.Get<std::string>("SOLVER", "INTERPOLATION_TYPE",
                                              std::string("LAGRANGE"));
  solvP.linearSolver =
      file.Get<std::string>("SOLVER", "LINEAR_SOLVER", std::string("BICGSTAB"));
  solvP.relaxation = file.Get<double>("SOLVER", "RELAXATION", 1.0);
  solvP.max_iter = file.Get<int>("SOLVER", "MAX_ITER", 100);

  return {io, sim, solvP};
}

void input::meshCheck(const meshData &mesh) {
  if (!mesh.nodes.empty() && mesh.nodes.find(0) == mesh.nodes.end()) {
    std::cerr << "\n\033[1;31m ->VLM ERROR: NODE IDs must start at 0 \033[0m"
              << std::endl;
    exit(1);
  }
  if (!mesh.vortexIDs.empty() &&
      mesh.vortexIDs.find(0) == mesh.vortexIDs.end()) {
    std::cerr << "\n\033[1;31m ->VLM ERROR: VORTEX IDs must start at 0 \033[0m"
              << std::endl;
    exit(1);
  }
  if (!mesh.doubletIDs.empty() &&
      mesh.doubletIDs.find(0) == mesh.doubletIDs.end()) {
    std::cerr << "\n\033[1;31m ->VLM ERROR: DOUBLET IDs must start at 0 \033[0m"
              << std::endl;
    exit(1);
  }
  if (!mesh.stationIDs.empty() &&
      mesh.stationIDs.find(0) == mesh.stationIDs.end()) {
    std::cerr
        << "\n\033[1;31m ->VLM ERROR: WINGSTATION IDs must start at 0 \033[0m"
        << std::endl;
    exit(1);
  }
  if (!mesh.wingIDs.empty() && mesh.wingIDs.find(0) == mesh.wingIDs.end()) {
    std::cerr << "\n\033[1;31m ->VLM ERROR: WING IDs must start at 0 \033[0m"
              << std::endl;
    exit(1);
  }
  if (!mesh.patchIDs.empty() && mesh.patchIDs.find(0) == mesh.patchIDs.end()) {
    std::cerr << "\n\033[1;31m ->VLM ERROR: WING IDs must start at 0 \033[0m"
              << std::endl;
    exit(1);
  }
}

meshData input::importMeshFile(const ioParam &names) {
  // Opening file and verifying if already exists
  std::ifstream file(names.meshFile);
  if (!file.is_open()) {
    std::cerr << "\n\033[1;31m ->VLM ERROR: mesh file \"" << names.meshFile
              << "\" not found! \033[0m" << std::endl;
    exit(1);
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
  meshCheck(mesh);
  return (mesh);
}
