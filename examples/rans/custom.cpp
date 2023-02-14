/*

        ____     ____  ___    _   _______
       /  _/    / __ \/   |  / | / / ___/
       / /_____/ /_/ / /| | /  |/ /\__ \
     _/ /_____/ _, _/ ___ |/ /|  /___/ /
    /___/    /_/ |_/_/  |_/_/ |_//____/

    Implicit RANS Solver
    Description: Example of a custom solver with custom CFL growth
    By: Alexis Angers

*/

// Use the Sparse LU direct solver
#define RANS_SPARSELU

#include <solver.h>
#include <post.h>


int main() {

    // Load the mesh, create the gas and the solver
    rans::mesh m("naca0012_mid.msh");
    rans::gas g(1./1.4, 0.001, 0.72, 0.9, 1.4);
    rans::explicitSolver s(m, g);
    s.set_second_order(false);

    // Init the solver with farfield conditions
    s.init(rans::farfield_conditions(1, 0.8, 0, 1));

    // Set the boundary conditions in the farfield
    // Important, this must be done after initialization
    s.farfield(rans::farfield_conditions(1, 0.8, 0, 1));

    // Iterations loop
    int i = 0;
    double err = 0;
    double err_last = 0;
    double err_0 = 0;
    do {
        if (i % s.get_print_interval() == 0) std::cout << "Iteration " << i << " " << std::flush;

        // Set CFL
        s.set_cfl(2);

        // Fill system matrices
        s.fill();
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        // Compute system matrices
        int is_error = s.compute();
        if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

        // Check if there's an error with the compute step
        if (is_error) {
            // Set error to -1 and skip solving
            err = -1;
            if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;
        } else {
            // Solve linear system
            err_last = err;

            err = s.solve();

            if (i % s.get_print_interval() == 0) std::cout << "." << std::flush;

            // Evaluate scaled residual
            if (err > 0) {
                if (i == 0) err_0 = err;
                err /= err_0;
            }
        }

        // Print residuals and step up iteration number i
        if (i % s.get_print_interval() == 0) std::cout << " Residual = " << err << std::endl;
        i++;

    } while ((err > 1e-6)&(i < 5000));

    // Save results to vtu file
    save("custom.vtu", s);

    return 0;
}

