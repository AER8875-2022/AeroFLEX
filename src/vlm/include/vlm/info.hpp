
#ifndef VLM_INFO_HPP
#define VLM_INFO_HPP

/** @file info.hpp */

#include "vlm/input.hpp"

namespace vlm {

namespace info {

/** @brief Function that prints the artwork/logo of VLM */
void printArtwork();

/** @brief Function that prints information on the current case to the terminal
 *  @param sim Object holding simulation parameters
 *  @param io Object holding input/output parameters
 *  @param solvP Object holding solver parameters */
void printCaseInfo(const Settings &settings);

/** @brief Function that prints information on the current mesh to the terminal
 *  @param mesh Information on the mesh */
void printMeshInfo(const input::meshData &mesh);

} // namespace info
} // namespace vlm

#endif
