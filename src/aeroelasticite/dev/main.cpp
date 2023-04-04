//
// Created by bocan on 2023-02-16.
//
#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include "aero/aerotools.h"
#include "vlm/model.hpp"
#include "aero/aeroelasticity.hpp"




int main(int argc, char **argv) {

    GUIHandler gui;
    vlm::VLM vlm(gui);
    structure::Structure structure(gui);
    aero::Aero elast(gui, vlm, structure);

  /* Verifying and handling input arguments*/
     if (argc < 2) {
         std::cerr << "\033[1;31m==>ERROR: expected config file as argument \033[0m"
                  << std::endl;
         std::cerr << "\033[1;36m==>USAGE>> vlm.out path/to/config/file \033[0m"
                   << std::endl;
         return 1;
     } else if (argc > 2) {
         std::cerr << "\033[1;33m==>WARNING: Extra arguments are ignored \033[0m"
              << std::endl;
     std::cerr << "==>\033[1;36mUSAGE>> vlm.out path/to/config/file\033[0m"
              << std::endl;
 }

    // #############################
    // LOADING SIMULATION PARAMETERS
    // #############################
    std::cout << "==>Loading simulation parameters...";

    tiny::config config;

    try {
        config.read(argv[1]);
        elast.vlm.settings.import_config_file(config);
        std::cout << "\033[1;36mDone\033[0m" << std::endl;

        vlm::info::printCaseInfo(elast.vlm.settings);

        elast.structure.settings.import_config_file(config);
    } catch (std::exception &e) {
        std::cerr << "\033[1;31m==>ERROR: " << e.what() << "\033[0m"
                  << std::endl;
        return 1;
    }


    elast.input();
    elast.solve();

    return 0;
}
