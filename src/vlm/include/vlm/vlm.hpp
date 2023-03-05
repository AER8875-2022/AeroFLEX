
#ifndef __VLM__
#define __VLM__

#include "common_aeroflex.hpp"
#include "vlm/input.hpp"
#include "vlm/model.hpp"
#include "vlm/solver.hpp"
#include "vlm/utils.hpp"

namespace vlm {

class VLM {

public:
  vlm::model object;
  database::table database;
  input::settings data;

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
