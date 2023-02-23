/*

        ____     ____  ___    _   _______
       /  _/    / __ \/   |  / | / / ___/
       / /_____/ /_/ / /| | /  |/ /\__ \
     _/ /_____/ _, _/ ___ |/ /|  /___/ /
    /___/    /_/ |_/_/  |_/_/ |_//____/

    Implicit RANS Solver
    Description: VTU file solution output
    By: Alexis Angers

*/
#pragma once

#include <iomanip>
#include <sstream>
#include <Eigen/Dense>
#include <fstream>

#include <rans/solver.h>



namespace rans {


std::string double2string(const double& x, const int precision) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << x;
    return stream.str();
}



struct wallProfile {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> cp;

    void reserve(const uint& n) {
        x.reserve(n);
        y.reserve(n);
        cp.reserve(n);
    }
};


void save(
    const std::string filename,
    solver& solv
) {
    solution q = solv.get_solution();
    mesh& m = solv.get_mesh();

    // Write solution q to .vtk file filename
    std::string s = "";

    s += "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n";
    s += "  <UnstructuredGrid>\n";
    s += "    <Piece NumberOfPoints=\"" + std::to_string(m.nodesX.size()) + "\" NumberOfCells=\"" + std::to_string(m.nRealCells) + "\">\n";
    s += "      <Points>\n";
    s += "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
    for (uint i=0; i<m.nodesX.size(); ++i) {
        s += "          " + std::to_string(m.nodesX[i]) + " " + std::to_string(m.nodesY[i]) + " 0.0\n";
    }
    s += "        </DataArray>\n";
    s += "      </Points>\n";

    // Figure out the total number of integer tags in the cells
    s += "      <Cells>\n";
    s += "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";

    // Add each cell
    s += "          ";
    for (int fi=0; fi<m.nRealCells; ++fi) {

        std::vector<uint> ns = {
            m.cellsNodes(fi, 0),
            m.cellsNodes(fi, 1),
            m.cellsNodes(fi, 2)
        };
        if (!m.cellsIsTriangle[fi]) {
            ns.push_back(m.cellsNodes(fi, 3));
        }
        for (auto ni : ns) {
            s += std::to_string(ni) + "  ";
        }
        s += "\n          ";
    }
    s += "        </DataArray>\n";

    // Save cell offsets
    s += "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
    s += "          ";
    uint current_offset = 0;
    for (uint i=0; i<m.nRealCells; ++i) {
        current_offset += m.cellsIsTriangle[i] ? 3 : 4;
        s += std::to_string(current_offset) + "  ";
    }
    s += "\n";
    s += "        </DataArray>\n";

    // Write cell types
    s += "        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
    s += "          ";
    for (uint i=0; i<m.nRealCells; ++i) {
        bool isTriangle = m.cellsIsTriangle[i];
        if (isTriangle) {
            // Tri cell
            s += "5  ";
        } else {
            // Quad cell
            s += "9  ";
        }
    }
    s += "\n";
    s += "        </DataArray>\n";
    s += "      </Cells>\n";

    // Save scalar data
    s += "      <CellData Scalars=\"scalars\">\n";

    // Save wall distance
    s += "        <DataArray type=\"Float32\" Name=\"Wall Distance\" Format=\"ascii\">\n";
    for (int j=0; j<m.nRealCells; ++j) {
        s += "          " + double2string(m.wall_dist[j], 16) + "\n";
    }
    s += "        </DataArray>\n";

    // Save variables
    s += "        <DataArray type=\"Float32\" Name=\"Mach\" Format=\"ascii\">\n";
    for (int i=0; i<m.nRealCells; ++i) {
        s += "          " + double2string(q.mach(i), 16) + "\n";
    }
    s += "        </DataArray>\n";
    s += "        <DataArray type=\"Float32\" Name=\"Density\" Format=\"ascii\">\n";
    for (int i=0; i<m.nRealCells; ++i) {
        s += "          " + double2string(q.rho(i), 16) + "\n";
    }
    s += "        </DataArray>\n";
    s += "        <DataArray type=\"Float32\" Name=\"Pressure\" Format=\"ascii\">\n";
    for (int i=0; i<m.nRealCells; ++i) {
        s += "          " + double2string(q.p(i), 16) + "\n";
    }
    s += "        </DataArray>\n";
    s += "        <DataArray type=\"Float32\" Name=\"Temperature\" Format=\"ascii\">\n";
    for (int i=0; i<m.nRealCells; ++i) {
        s += "          " + double2string(q.T(i), 16) + "\n";
    }
    s += "        </DataArray>\n";

    s += "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
    for (int i=0; i<m.nRealCells; ++i) {
        s += "          " + double2string(q.u(i), 16) + " ";
        s += double2string(q.v(i), 16) + " 0.0\n";
    }
    s += "        </DataArray>\n";

    s += "      </CellData>\n";

    s += "    </Piece>\n";
    s += "  </UnstructuredGrid>\n";
    s += "</VTKFile>";

    // Save data to file
    std::ofstream out(filename);
    out << s;
    out.close();

}




wallProfile get_wall_profile(solver& solv, std::string patch_name) {

    wallProfile wall;

    mesh& m = solv.get_mesh();

    auto bcs = solv.get_bcs();

    for (uint b=0; b<m.boundaryEdges.size(); ++b) {
        
    }

    return wall;
}



}


