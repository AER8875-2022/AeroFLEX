
#ifndef __VLM__
#define __VLM__

#include "common_aeroflex.hpp"
#include "vlm/info.hpp"
#include "vlm/input.hpp"
#include "vlm/model.hpp"
#include "vlm/solver.hpp"

namespace vlm {

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
