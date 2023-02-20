
#ifndef VLM_OUTPUT_HPP
#define VLM_OUTPUT_HPP

#include "vlm/input.hpp"
#include "vlm/model.hpp"

/** @file output.hpp */

namespace vlm {

namespace output {

/** @brief Function that exports computed forces to a .dat file
 *  @param object Model holding all the VLM elements
 *  @param it Current iteration number. Useful for unsteady simulations */
void exportForces(const model &object, const int it);

/** @brief Function that exports the surface solution to .vtu format
 *  @param object Model holding all the VLM elements
 *  @param it Current iteration number. Useful for unsteady simulations */
void exportSurfacesVTU(const model &object, const int it);

/** @brief Function that exports wake solutions to .vtu format
 *  @param object Model holding all the VLM elements
 *  @param it Current iteration number. Useful for unsteady simulations */
void exportWakeVTU(const model &object, const int it);

} // namespace output
} // namespace vlm

#endif
