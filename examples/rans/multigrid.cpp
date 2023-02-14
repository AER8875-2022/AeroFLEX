/*

        ____     ____  ___    _   _______
       /  _/    / __ \/   |  / | / / ___/
       / /_____/ /_/ / /| | /  |/ /\__ \
     _/ /_____/ _, _/ ___ |/ /|  /___/ /
    /___/    /_/ |_/_/  |_/_/ |_//____/

    Implicit RANS Solver
    Description: Example of multigrid solve using three GMSH meshes
    By: Alexis Angers

*/

// Define compiler directives
// #define RANS_DEBUG  // Prints linear solver iterations

// Include rans library
#include <rans.h>



int main() {

    // Load meshes
    rans::mesh m0("naca0012_coarse.msh");  // 4 000 elements
    rans::mesh m1("naca0012_mid.msh");     // 20 000 elements
    rans::mesh m2("naca0012_fine.msh");    // 50 000 elements

    // Define farfield conditions, (T, mach, angle, p)
    rans::farfield_conditions vars_far(1, 0.8, 0, 1);

    // Define gas, (R, mu, Pr_L, Pr_T, gamma)
    rans::gas g(1./1.4, 0.001, 0.72, 0.9, 1.4);

    // Create the multigrid solver
    rans::multigrid<rans::implicitSolver> multi({m0, m1, m2}, g);

    // Solve and return the fine grid converged solution
    rans::solver& s = multi.run(vars_far);

    // Save the solution to VTU file
    save("multigrid.vtu", s);

    return 0;
}

