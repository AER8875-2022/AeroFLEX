#pragma once

#include <rans/core.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#define __RANS_VERSION__ 0.1
#define __RANS_RELEASE_DATE__ "2022-02-02"


namespace rans {

Settings parse(int argc, char **argv, std::string compiled_file) {

    Settings data;

    // Print headers

    std::cout << "\
    \n\
    \033[0;32mI-RANS Implicit RANS Solver\033[0m\n\
    - Version       : " << __RANS_VERSION__ << "\n\
    - Release Date  : " << __RANS_RELEASE_DATE__ << "\n\
    - Compiled at   : " << __TIMESTAMP__ << "\n\
    - Compiled file : " << compiled_file << "\n";

#ifdef _OPENMP
    std::cout << "\
    - OpenMP enable : true\n\
    - Threads count : " << Eigen::nbThreads() << "\n";
#else
    std::cout << "\
    - OpenMP enable : false\n";
#endif
    std::cout << std::endl;

    // Parse command line inputs

    std::string filename = "";
    for (int i = 0; i < argc-1; ++i) {
        if ((std::string(argv[i]) == "-f")|(std::string(argv[i]) == "--file")) {
            filename = argv[i+1];
        }
    }
    std::string outfilename = "";
    for (int i = 0; i < argc-1; ++i) {
        if ((std::string(argv[i]) == "-o")|(std::string(argv[i]) == "--output")) {
            outfilename = argv[i+1];
        }
    }

    if (filename == "") {
        std::cout << "\033[0;31m";
        std::cout << "Error, no input file name provided. ";
        std::cout << "\033[0mAdd to CLI : -f <input file name>\n";
        data.read_failure = 2;
        return data;
    }
    if (outfilename == "") {
        std::cout << "\033[0;31m";
        std::cout << "Error, no output file name provided. ";
        std::cout << "\033[0mAdd to CLI : -o <output file name>\n";
        data.read_failure = 2;
        return data;
    }

    data.outfilename = outfilename;


    std::cout << "Input file name  : " << filename << "\n";
    std::cout << "Output file name : " << outfilename << "\n";
    std::cout << std::endl;

    // Read input file

    std::ifstream infile;
    infile.open(filename);
    if (infile.fail()) {
        std::cout << "\033[0;31m";
        std::cout << "Error, failed to open input file ";
        std::cout << "\033[0m\n";
        data.read_failure = 2;
        return data;
    }

    std::string tag;
    infile >> tag;
    do {
        if (tag[0] == '[') {
            if (tag == "[gas]") {
                std::string property;
                infile >> property;
                infile >> tag;
                while ((tag != "end") & (!infile.eof())) {
                    double prop_i = std::stod(tag);
                    if (property == "gamma")
                        data.g.gamma = prop_i;
                    else if (property == "R")
                        data.g.R = prop_i;
                    else if (property == "mu_L")
                        data.g.mu_L = prop_i;
                    else if (property == "Pr_L")
                        data.g.Pr_L = prop_i;
                    else if (property == "Pr_T")
                        data.g.Pr_T = prop_i;
                    else if (property == "cp")
                        data.g.cp = prop_i;
                    infile >> property;
                    if (property == "end")
                        tag = "end";
                    else
                        infile >> tag;
                }
            } else if (tag == "[bc]") {
                std::string bc_type;
                std::string bc_name;
                boundary_variables bc_vars;

                std::string property;
                infile >> property;
                infile >> tag;
                while ((tag != "end") & (!infile.eof())) {
                    if (property == "type")
                        bc_type = tag;
                    else if (property == "name")
                        bc_name = tag;
                    else if (property == "T")
                        bc_vars.T = std::stod(tag);
                    else if (property == "mach")
                        bc_vars.mach = std::stod(tag);
                    else if (property == "angle")
                        bc_vars.angle = std::stod(tag);
                    else if (property == "p")
                        bc_vars.p = std::stod(tag);
                    infile >> property;
                    if (property == "end")
                        tag = "end";
                    else
                        infile >> tag;
                }

                boundary_condition bc_i;
                bc_i.bc_type = bc_type;
                bc_i.vars_far = bc_vars;
                data.bcs[bc_name] = bc_i;
            } else if (tag == "[solver]") {
                std::string property;
                infile >> property;
                infile >> tag;
                while ((tag != "end") & (!infile.eof())) {
                    std::string prop_i = tag;
                    if (property == "type")
                        data.set_solver_type(prop_i);
                    else if (property == "second_order") {
                        std::vector<std::string> mylist{"1", "true", "True", "TRUE"};
                        data.second_order = std::find(std::begin(mylist), std::end(mylist), prop_i) != std::end(mylist);
                    } else if (property == "viscosity_model") {
                        data.viscosity_model = prop_i;
                    } else if (property == "relaxation") {
                        data.relaxation = std::stod(prop_i);
                    } else if (property == "tolerance") {
                        data.tolerance = std::stod(prop_i);
                    } else if (property == "start_cfl") {
                        data.start_cfl = std::stod(prop_i);
                    } else if (property == "max_cfl") {
                        data.max_cfl = std::stod(prop_i);
                    } else if (property == "slope_cfl") {
                        data.slope_cfl = std::stod(prop_i);
                    } else if (property == "rhs_iterations") {
                        data.rhs_iterations = (uint) std::stoi(prop_i);
                    } else if (property == "max_iterations") {
                        data.max_iterations = (uint) std::stoi(prop_i);
                    } else if (property == "gradients") {
                        data.gradient_scheme = prop_i;
                    } else if (property == "limiter_k") {
                        data.limiter_k = std::stod(prop_i);
                    } else if (property == "airfoil_name") {
                        data.airfoil_name = prop_i;
                    } else {
                        std::cout << "\033[0;33m";
                        std::cout << "Warning, skipped unknown property " + property;
                        std::cout << "\033[0m\n";
                    }
                    infile >> property;
                    if (property == "end")
                        tag = "end";
                    else
                        infile >> tag;
                }
            } else if (tag == "[mesh]") {
                std::string mesh_name;
                infile >> mesh_name;
                while ((mesh_name != "end") & (!infile.eof())) {
                    data.meshes.push_back(mesh_name);
                    infile >> mesh_name;
                }
            }
        }
        infile >> tag;
    } while (!infile.eof());

    infile.close();

    data.read_failure = 0;

    return data;
}





}


