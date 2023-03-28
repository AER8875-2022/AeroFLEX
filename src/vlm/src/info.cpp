
#include "vlm/info.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace vlm;
using namespace info;

void info::printArtwork() {
  std::cout << "\n";
  std::cout << "██╗   ██╗██╗     ███╗   ███╗"
            << "\n";
  std::cout << "██║   ██║██║     ████╗ ████║"
            << "\n";
  std::cout << "██║   ██║██║     ██╔████╔██║"
            << "\n";
  std::cout << "╚██╗ ██╔╝██║     ██║╚██╔╝██║"
            << "\n";
  std::cout << " ╚████╔╝ ███████╗██║ ╚═╝ ██║"
            << "\n";
  std::cout << "  ╚═══╝  ╚══════╝╚═╝     ╚═╝"
            << "\n";
  std::cout << std::endl;
}

void info::printCaseInfo(const Settings &settings) {
  std::ostringstream aoa, sideslip, vinf, rho, cref, sref, coreRadius;
  aoa << std::setprecision(4) << settings.sim.aoa;
  sideslip << std::setprecision(4) << settings.sim.sideslip;
  vinf << std::setprecision(4) << settings.sim.vinf;
  rho << std::setprecision(4) << settings.sim.rho;
  cref << std::setprecision(4) << settings.sim.cref;
  sref << std::setprecision(4) << settings.sim.sref;
  coreRadius << std::setprecision(4) << settings.sim.coreRadius;

  std::cout << std::endl;
  std::cout << "     \033[1;20mVLM CASE:\033[0m " << settings.io.baseName
            << "\n";
  std::cout << "      \033[1;20mOUT DIR:\033[0m " << settings.io.outDir << "\n";
  std::cout << "    \033[1;20mMESH FILE:\033[0m " << settings.io.meshFile
            << "\n\n";

  std::cout << "  \033[1;20mSOLVER TYPE:\033[0m " << settings.solver.type
            << "\n";
  std::cout << "  \033[1;20mTIME DOMAIN:\033[0m " << settings.solver.timeDomain
            << "\n";
  std::cout << "    \033[1;20mTOLERANCE:\033[0m " << settings.solver.tolerance
            << "\n\n";

  std::cout << "          \033[1;20mAOA:\033[0m " << aoa.str() << "\n";
  std::cout << "     \033[1;20mSIDESLIP:\033[0m " << sideslip.str() << "\n";
  std::cout << "         \033[1;20mVINF:\033[0m " << vinf.str() << "\n";
  std::cout << "          \033[1;20mRHO:\033[0m " << rho.str() << "\n";
  std::cout << "        \033[1;20mC_REF:\033[0m " << cref.str() << "\n";
  std::cout << "        \033[1;20mS_REF:\033[0m " << sref.str() << "\n";
  std::cout << "   \033[1;20mCORERADIUS:\033[0m " << coreRadius.str() << "\n";
  std::cout << std::endl;
}

void info::printMeshInfo(const input::meshData &mesh) {
  std::cout << std::endl;
  std::cout << "             \033[1;20mNUMBER OF NODES:\033[0m "
            << mesh.nodes.size() << "\n";
  std::cout << "      \033[1;20mNUMBER OF VORTEX RINGS:\033[0m "
            << mesh.vortexIDs.size() << "\n";
  std::cout << "            \033[1;20mNUMBER OF PANELS:\033[0m "
            << mesh.doubletIDs.size() << "\n";
  std::cout << "     \033[1;20mNUMBER OF WING STATIONS:\033[0m "
            << mesh.stationIDs.size() << "\n";
  std::cout << "             \033[1;20mNUMBER OF WINGS:\033[0m "
            << mesh.wingIDs.size() << "\n";
  std::cout << "   \033[1;20mNUMBER OF SURFACE PATCHES:\033[0m "
            << mesh.patchIDs.size() << "\n";
  std::cout << std::endl;
}
