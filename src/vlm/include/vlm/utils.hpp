
#ifndef VLM_UTILS_HPP
#define VLM_UTILS_HPP

/** @file utils.hpp */

#include "vlm/input.hpp"

namespace vlm {

namespace utils {

/** @brief Function that prints the artwork/logo of VLM */
void printArtwork();

/** @brief Function that prints information on the current case to the terminal
 *  @param sim Object holding simulation parameters
 *  @param io Object holding input/output parameters
 *  @param solvP Object holding solver parameters */
void printCaseInfo(const input::simParam &sim, const input::ioParam &io,
                   const input::solverParam &solvP);

/** @brief Function that prints information on the current mesh to the terminal
 *  @param mesh Information on the mesh */
void printMeshInfo(const input::meshData &mesh);

} // namespace utils
} // namespace vlm

#endif
