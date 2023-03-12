
#ifndef __VLM__
#define __VLM__

#include "common_aeroflex.hpp"
#include "vlm/input.hpp"
#include "vlm/model.hpp"
#include "vlm/solver.hpp"
#include "vlm/utils.hpp"

namespace vlm {

/** @brief Main struct holding parameters on the current case */
struct Settings {

  /** @brief Object holding case and physics-oriented parameters */
  input::simParam sim;

  /** @brief Object holding input/output information */
  input::ioParam io;

  /** @brief Object holding solver parameters */
  input::solverParam solver;
};

class VLM {

public:
  model object;
  database::table database;
  Settings settings;

  std::vector<double> residuals;
  std::atomic<int> iters = 0;
  GUIHandler &gui;

private:
  bool is_initialized = false;

public:
  VLM(GUIHandler &gui);
  void input();
  void solve();

private:
  void reinitialize();

}; // end class vlm

} // end namespace vlm

#endif
