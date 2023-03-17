
#include "vlm/input.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <exception>

using namespace vlm;
using namespace input;
using namespace Eigen;

void simParam::set_databaseFormat(const std::string &item) {
  if (!item.compare("NONE")) {
    databaseFormat = 0;
  }
  else if (!item.compare("FILE")) {
    databaseFormat = 1;
  }
}

std::string simParam::get_databaseFormat() const { return databaseFormat_options.at(databaseFormat); }

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

Vector3d simParam::origin() const { return Vector3d(x0,y0,z0); }

// -------------------------------------

void solverParam::set_timeDomain(const std::string &item) {
  if (!item.compare("STEADY")) {
    timeDomain = 0;
  }
  else if (!item.compare("UNSTEADY")) {
    timeDomain = 1;
  }
}

std::string solverParam::get_timeDomain() const {
  return timeDomain_options.at(timeDomain);
}

void solverParam::set_type(const std::string &item) {
  if (!item.compare("LINEAR")) {
    type = 0;
  }
  else if (!item.compare("NONLINEAR")) {
    type = 1;
  }
}

std::string solverParam::get_type() const {
  return type_options.at(type);
}

void solverParam::set_linearSolver(const std::string &item) {
  if (!item.compare("BICGSTAB")) {
    linearSolver = 0;
  }
  else if (!item.compare("DIRECT")) {
    linearSolver = 1;
  }
}

std::string solverParam::get_linearSolver() const {
  return linearSolver_options.at(linearSolver);
}

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

  // [vlm-simulation]
  sim.aoa = config.get<double>("vlm-simulation", "aoa", 5.0);
  sim.sideslip = config.get<double>("vlm-simulation", "sideslip", 0.0);
  sim.vinf = config.get<double>("vlm-simulation", "v_inf", 1.0);
  sim.rho = config.get<double>("vlm-simulation", "density", 1.0);
  sim.cref = config.get<double>("vlm-simulation", "c_ref", 1.0);
  sim.sref = config.get<double>("vlm-simulation", "s_ref", 1.0);
  sim.x0 = config.get<double>("vlm-simulation", "x_ref", 0.0);
  sim.y0 = config.get<double>("vlm-simulation", "y_ref", 0.0);
  sim.z0 = config.get<double>("vlm-simulation", "z_ref", 0.0);
  sim.coreRadius = config.get<double>("vlm-simulation", "lamb-oseen_radius", 0.0);
  sim.set_databaseFormat(config.get<std::string>("vlm-simulation", "database_format"));

  // [vlm-io]
  io.baseName = config.get<std::string>("vlm-io", "basename");
  io.outDir = config.get<std::string>("vlm-io", "output_dir");
  io.meshFile =
      config.get<std::string>("vlm-io", "mesh_file");
  io.databaseFile =
      config.get<std::string>("vlm-io", "database_file");
  io.locationFile =
      config.get<std::string>("vlm-io", "location_file");

  // [vlm-solver]
  solver.set_timeDomain(
      config.get<std::string>("vlm-solver", "time_domain"));
  solver.set_type(config.get<std::string>("vlm-solver", "type"));
  solver.tolerance = config.get<double>("vlm-solver", "tolerance", 1e-15);
  solver.set_linearSolver(
      config.get<std::string>("vlm-solver", "linear_solver"));
  solver.relaxation = config.get<double>("vlm-solver", "relaxation", 1.0);
  solver.max_iter = config.get<int>("vlm-solver", "max_iter", 100);

}

void Settings::export_config_file(tiny::config &config) {

  // Appending sections to config file
  config.sections.push_back("vlm-simulation");
  config.sections.push_back("vlm-io");
  config.sections.push_back("vlm-solver");

  // Adding to config map
  // [vlm-simulation]
  config.config["vlm-simulation"]["aoa"] = std::to_string(sim.aoa);
  config.config["vlm-simulation"]["sideslip"] = std::to_string(sim.sideslip);
  config.config["vlm-simulation"]["v_inf"] = std::to_string(sim.vinf);
  config.config["vlm-simulation"]["density"] = std::to_string(sim.rho);
  config.config["vlm-simulation"]["c_ref"] = std::to_string(sim.cref);
  config.config["vlm-simulation"]["s_ref"] = std::to_string(sim.sref);
  config.config["vlm-simulation"]["x_ref"] = std::to_string(sim.x0);
  config.config["vlm-simulation"]["y_ref"] = std::to_string(sim.y0);
  config.config["vlm-simulation"]["z_ref"] = std::to_string(sim.z0);
  config.config["vlm-simulation"]["lamb-oseen_radius"] = std::to_string(sim.coreRadius);
  config.config["vlm-simulation"]["database_format"] = sim.get_databaseFormat();

  // [vlm-io]
  config.config["vlm-io"]["basename"] = io.baseName;
  config.config["vlm-io"]["output_dir"] = io.outDir;
  config.config["vlm-io"]["database_file"] = io.databaseFile;
  config.config["vlm-io"]["location_file"] = io.locationFile;
  config.config["vlm-io"]["mesh_file"] = io.meshFile;

  // [vlm-solver]
  config.config["vlm-solver"]["time_domain"] = solver.get_timeDomain();
  config.config["vlm-solver"]["type"] = solver.get_type();
  config.config["vlm-solver"]["tolerance"] = std::to_string(solver.tolerance);
  config.config["vlm-solver"]["linear_solver"] = solver.get_linearSolver();
  config.config["vlm-solver"]["relaxation"] = std::to_string(solver.relaxation);
  config.config["vlm-solver"]["max_iter"] = std::to_string(solver.max_iter);

}
