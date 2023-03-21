#include <rans/rans.h>
#include <rans/parser.h>

using namespace rans;

template<class T>
void multigrid_run_and_save(std::vector<mesh> &ms, Settings data) {

    multigrid<T> multi(ms, data.g, data.second_order, signal, residuals, iters);

    rans::solver& s = multi.run(data.vars_far, true, data.relaxation);

    save(data.outfilename, s);
    std::cout << "Saved results to file " << data.outfilename << "\n" << std::endl;
}

int main(int argc, char **argv) {

    auto data = parse(argc, argv, __FILE__);

    if (data.read_failure == 1) {
        std::cout << "\033[0;31m";
        std::cout << "Error, unspecified error reading input file.";
        std::cout << "\033[0m\n" << std::endl;
        return 1;
    } else if (data.read_failure == 2) {
        std::cout << std::endl;
        return 1;
    }

    std::vector<mesh> ms;

    for (const auto& mesh_name : data.meshes) {
        ms.push_back(mesh(mesh_name));
    }

    if (data.solver_type() == "implicit") {
        multigrid_run_and_save<implicitSolver>(ms, data);
    } else if (data.solver_type() == "explicit") {
        multigrid_run_and_save<explicitSolver>(ms, data);
    }

    return 0;
}

