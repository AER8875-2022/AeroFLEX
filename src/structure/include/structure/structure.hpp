
#ifndef __STRUCTURE__
#define __STRUCTURE__

#include "common_aeroflex.hpp"
#include <vector>
#include <atomic>
#include "structure/model.hpp"

namespace structure {

struct Settings {
  // TODO
};

class Structure {

public:
  Settings settings;
  // TODO: STORE SOLUTION

  std::vector<double> residuals;
  std::atomic<int> iters = 0;
  GUIHandler &gui;

public:
  Structure(GUIHandler &gui): gui(gui) {}
  void input() {} // TODO
  void solve() {} // TODO
};

}

#endif
