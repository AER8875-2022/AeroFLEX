
#include "vlm/output.hpp"
#include "vtu11/vtu11.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace vlm;

void output::exportForces(const model &object, const int it) {
  auto &io = object.io;

  // Output
  std::string output = "";

  // Converting time iteration to string with padded 0
  std::ostringstream itStream;
  itStream << std::setfill('0') << std::setw(6) << it;

  // File header
  output += "######################################\n";
  output += "# " + io.baseName + "_forces_" + itStream.str() + ".dat" + "\n\n";
  output += "# NOTE: GLOBAL coefficients are based on S_REF and C_REF while \n";
  output += "#       SURFACE coefficients are based on C_REF and their      \n";
  output += "#       surface area.                                          \n";
  output += "######################################\n";
  output += "\n";

  std::ostringstream CL, CD, CY, CMx, CMy, CMz;

  // Global coefficients
  output += "GLOBAL:\n";
  CL << std::setprecision(12) << object.get_cl();
  output += "CL = " + CL.str() + "\n";
  CD << std::setprecision(12) << object.get_cd();
  output += "CD = " + CD.str() + "\n";
  CY << std::setprecision(12) << object.get_cy();
  output += "CY = " + CY.str() + "\n";
  CMx << std::setprecision(12) << object.get_cm()(0);
  output += "CMx = " + CMx.str() + "\n";
  CMy << std::setprecision(12) << object.get_cm()(1);
  output += "CMy = " + CMy.str() + "\n";
  CMz << std::setprecision(12) << object.get_cm()(2);
  output += "CMz = " + CMz.str() + "\n";
  output += "\n";

  // Surfaces coefficients
  int i = 0;
  for (auto &wing : object.wings) {
    std::ostringstream CL, CD, CY, CMx, CMy, CMz;
    output += "SURFACE " + std::to_string(i) + ":\n";
    CL << std::setprecision(12) << wing.get_cl();
    output += "CL = " + CL.str() + "\n";
    CD << std::setprecision(12) << wing.get_cd();
    output += "CD = " + CD.str() + "\n";
    CY << std::setprecision(12) << wing.get_cy();
    output += "CY = " + CY.str() + "\n";
    CMx << std::setprecision(12) << wing.get_cm()(0);
    output += "CMx = " + CMx.str() + "\n";
    CMy << std::setprecision(12) << wing.get_cm()(1);
    output += "CMy = " + CMy.str() + "\n";
    CMz << std::setprecision(12) << wing.get_cm()(2);
    output += "CMz = " + CMz.str() + "\n";
    output += "\n";
    i++;
  }

  // Output to file
  std::filesystem::create_directories(io.outDir);
  std::ofstream file(io.outDir + io.baseName + "_forces_" + itStream.str() +
                     ".dat");
  file << output;
  file.close();
}

void output::exportSurfacesVTU(const model &object, const int it) {
  const auto &io = object.io;
  const auto nPanels = object.vortexRings.size() + object.doubletPanels.size();

  // Converting time iteration to string with padded 0
  std::ostringstream itStream;
  itStream << std::setfill('0') << std::setw(6) << it;

  // mesh vertices
  std::vector<double> points;
  points.reserve(3 * object.nodes.size());
  for (size_t i = 0; i < object.nodes.size(); i++) {
    points.push_back(object.nodes[i](0));
    points.push_back(object.nodes[i](1));
    points.push_back(object.nodes[i](2));
  }
  // Connectivity table
  std::vector<vtu11::VtkIndexType> connectivity;
  connectivity.reserve(4 * nPanels);
  for (auto &vortex : object.vortexRings) {
    for (auto &node : vortex.get_nodeIDs()) {
      connectivity.push_back(node);
    }
  }
  for (auto &doublet : object.doubletPanels) {
    for (auto &node : doublet.get_nodeIDs()) {
      connectivity.push_back(node);
    }
  }
  // Cell types and offsets
  std::vector<vtu11::VtkCellType> types;
  types.reserve(nPanels);
  std::vector<vtu11::VtkIndexType> offsets;
  offsets.reserve(nPanels);
  int count = 0;
  for (size_t i = 0; i != object.vortexRings.size(); i++) {
    types.push_back(9);
    count += 4;
    offsets.push_back(count);
  }
  for (auto &doublet : object.doubletPanels) {
    if (doublet.get_nodeIDs().size() == 3) {
      types.push_back(9);
      count += 4;
      offsets.push_back(count);
    } else if (doublet.get_nodeIDs().size() == 4) {
      types.push_back(5);
      count += 3;
      offsets.push_back(count);
    } else {
      std::cerr << "\033[1;33m==>WARNING: Elements with more than 4 nodes are "
                   "not yet supported\033[0m"
                << std::endl;
      return;
    }
  }
  // Output data
  // STRENGTH
  std::vector<double> strengths;
  strengths.reserve(nPanels);
  for (auto &vortex : object.vortexRings) {
    strengths.push_back(vortex.get_gamma());
  }
  for (auto &doublet : object.doubletPanels) {
    strengths.push_back(doublet.get_sigma());
  }
  // AREA
  std::vector<double> areas;
  areas.reserve(nPanels);
  for (auto &vortex : object.vortexRings) {
    areas.push_back(vortex.get_area());
  }
  for (auto &doublet : object.doubletPanels) {
    areas.push_back(doublet.get_area());
  }
  // CF
  std::vector<double> cf;
  cf.reserve(3 * nPanels);
  for (auto &vortex : object.vortexRings) {
    cf.push_back(vortex.get_cf()(0));
    cf.push_back(vortex.get_cf()(1));
    cf.push_back(vortex.get_cf()(2));
  }
  for (auto &doublet : object.doubletPanels) {
    cf.push_back(0.0);
  }
  // CM
  std::vector<double> cm;
  cm.reserve(3 * nPanels);
  for (auto &vortex : object.vortexRings) {
    cm.push_back(vortex.get_cm()(0));
    cm.push_back(vortex.get_cm()(1));
    cm.push_back(vortex.get_cm()(2));
  }
  for (auto &doublet : object.doubletPanels) {
    cm.push_back(0.0);
    cm.push_back(0.0);
    cm.push_back(0.0);
  }
  // Creating data object
  std::vector<vtu11::DataSetInfo> dataSetInfo{
      {"strength", vtu11::DataSetType::CellData, 1},
      {"area", vtu11::DataSetType::CellData, 1},
      {"cf", vtu11::DataSetType::CellData, 3},
      {"cm", vtu11::DataSetType::CellData, 3}};
  // Creating mesh object
  vtu11::Vtu11UnstructuredMesh mesh{points, connectivity, offsets, types};
  // Writing file
  vtu11::writeVtu(io.outDir + io.baseName + "_surface_" + itStream.str() +
                      ".vtu",
                  mesh, dataSetInfo, {strengths, areas, cf, cm}, "RawBinary");
}

void output::exportWakeVTU(const model &object, const int it) {
  auto &io = object.io;
  // Converting time iteration to string with padded 0
  std::ostringstream itStream;
  itStream << std::setfill('0') << std::setw(6) << it;

  // mesh vertices
  std::vector<double> points;
  points.reserve(3 * object.wakeNodes.size());
  for (size_t i = 0; i < object.wakeNodes.size(); i++) {
    points.push_back(object.wakeNodes[i](0));
    points.push_back(object.wakeNodes[i](1));
    points.push_back(object.wakeNodes[i](2));
  }
  // Connectivity table
  std::vector<vtu11::VtkIndexType> connectivity;
  connectivity.reserve(4 * object.wakePanels.size());
  for (auto &wake : object.wakePanels) {
    for (auto &node : wake.get_nodeIDs()) {
      connectivity.push_back(node);
    }
  }
  // Cell types and offsets
  std::vector<vtu11::VtkCellType> types;
  types.reserve(object.wakePanels.size());
  std::vector<vtu11::VtkIndexType> offsets;
  offsets.reserve(object.wakePanels.size());
  int count = 0;
  for (size_t i = 0; i != object.wakePanels.size(); i++) {
    types.push_back(9);
    count += 4;
    offsets.push_back(count);
  }
  // Output data
  // STRENGTH
  std::vector<double> strengths;
  strengths.reserve(object.wakePanels.size());
  for (auto &wake : object.wakePanels) {
    strengths.push_back(wake.get_gamma());
  }
  // AREA
  std::vector<double> areas;
  areas.reserve(object.wakePanels.size());
  for (auto &wake : object.wakePanels) {
    areas.push_back(wake.get_area());
  }
  // Creating data object
  std::vector<vtu11::DataSetInfo> dataSetInfo{
      {"strength", vtu11::DataSetType::CellData, 1},
      {"area", vtu11::DataSetType::CellData, 1},
  };
  // Creating mesh object
  vtu11::Vtu11UnstructuredMesh mesh{points, connectivity, offsets, types};
  // Writing file
  vtu11::writeVtu(io.outDir + io.baseName + "_wake_" + itStream.str() + ".vtu",
                  mesh, dataSetInfo, {strengths, areas}, "RawBinary");
}
