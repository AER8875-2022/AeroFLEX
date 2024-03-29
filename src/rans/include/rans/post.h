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
#include <limits>
#include <iostream>
#include <atomic>
#include <mutex>

#include <Eigen/Dense>
#include <fstream>
#include "common_aeroflex.hpp"
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

    double cd = 0.;
    double cl = 0.;
    double cm = 0.;

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

class CpProfile {
    public:
    double x_moment;
    double y_moment;
    std::mutex m_mutex;

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> cp;
    std::vector<bool> cp_positive;
    std::vector<double> cp_airfoil_x;
    std::vector<double> cp_airfoil_y;
    double xmax = 0.1;
    double xmin = -0.1;
    double ymax = 0.1;
    double ymin = -0.1;
    bool filled = false;

    void calc_chord_coords(solver &s, std::string& af) {
        std::scoped_lock lock(m_mutex);
        mesh& m = s.get_mesh();
        filled = false;
        xmax = 0.1;
        xmin = -0.1;
        ymax = 0.1;
        ymin = -0.1;
        x.clear();
        y.clear();
        cp.clear();
        cp_positive.clear();
        cp_airfoil_x.clear();
        cp_airfoil_y.clear();
        double xleft = std::numeric_limits<double>::max();
        double xright = std::numeric_limits<double>::min();
        double yleft;
        double yright;

        for (uint b=0; b<m.boundaryEdges.size(); ++b) {
            if (m.boundaryEdgesPhysicals[b] == af) {
                const uint e = m.boundaryEdges[b];
                if (m.edgesCentersX[e] < xleft) {
                    xleft = m.edgesCentersX[e];
                    yleft = m.edgesCentersY[e];
                }
                if (m.edgesCentersX[e] > xright) {
                    xright = m.edgesCentersX[e];
                    yright = m.edgesCentersY[e];
                }
                cp.push_back(0.0);
                cp_airfoil_x.push_back(0.0);
                cp_airfoil_y.push_back(0.0);
                cp_positive.push_back(false);

                x.push_back(m.edgesCentersX[e]);
                y.push_back(m.edgesCentersY[e]);
                xmax = std::max(xmax, m.edgesCentersX[e]);
                ymax = std::max(ymax, m.edgesCentersY[e]);
                xmin = std::min(xmin, m.edgesCentersX[e]);
                ymin = std::min(ymin, m.edgesCentersY[e]);
            }
        }

        y_moment = (yright - yleft)*0.25 + yleft;
        x_moment = (xright - xleft)*0.25 + xleft;
    }

    void calc_cp(solver &s, std::string& af) {
        mesh& m = s.get_mesh();
        boundary_variables vars_far = s.get_boundary_variables();
        const double p_inf = vars_far.p;
        const double mach_inf = vars_far.mach;
        solution sol = s.get_solution();

        int idx = 0;
        int idx_neg = 0;
        int idx_pos = 0;
        std::scoped_lock lock(m_mutex);
        for (uint b=0; b<m.boundaryEdges.size(); ++b) {
            if (m.boundaryEdgesPhysicals[b] == af) {
                const uint e = m.boundaryEdges[b];
                const double x = m.edgesCentersX[e];
                const double y = m.edgesCentersY[e];
                const double nx = m.edgesNormalsX[e];
                const double ny = m.edgesNormalsY[e];
                const uint n0 = m.edgesNodes(e, 0);
                const uint n1 = m.edgesNodes(e, 1);

                const uint cell = m.edgesCells(e, 0);

                // Also calculate cp
                const double rho = sol.rho(cell);
                const double u = sol.u(cell);
                const double v = sol.v(cell);
                const double p = sol.p(cell);
                const double umag = sqrt(u*u + v*v);
                const double cs = sol.c(cell);
                const double mach = umag / cs;
                const double cp_ = 2./(sol.gamma()*mach_inf*mach_inf)*(p/p_inf - 1.);

                cp[idx] = cp_;
                if (cp_ > 0.0) {
                    cp_airfoil_x[idx] = x - nx*cp_*0.1;
                    cp_airfoil_y[idx] = y - ny*cp_*0.1;
                    cp_positive[idx] = true;
                } else {
                    cp_airfoil_x[idx] = x + nx*cp_*0.1;
                    cp_airfoil_y[idx] = y + ny*cp_*0.1;
                }
                xmax = std::max(xmax, cp_airfoil_x[idx]);
                ymax = std::max(ymax, cp_airfoil_y[idx]);
                xmin = std::min(xmin, cp_airfoil_x[idx]);
                ymin = std::min(ymin, cp_airfoil_y[idx]);
                idx++;
            }
        }
        filled = true;
    }
};

wallProfile get_wall_profile(solver& solv, std::string patch_name) {

    // First, get farfield conditions and compute infinity variables
    boundary_variables vars_far = solv.get_boundary_variables();
    const double p_inf = vars_far.p;
    const double mach_inf = vars_far.mach;

    solution sol = solv.get_solution();

    wallProfile wall;

    mesh& m = solv.get_mesh();

    // Find quarter chord position
    double xmin = 0.;
    double xmax = 0.;
    double y_moment = 0.;
    uint n_added = 0;
    bool first_added = true;
    for (uint b=0; b<m.boundaryEdges.size(); ++b) {
        if (m.boundaryEdgesPhysicals[b] == patch_name) {
            const uint e = m.boundaryEdges[b];
            if (first_added) {
                xmin = m.edgesCentersX[e];
                xmax = m.edgesCentersX[e];
                y_moment = m.edgesCentersY[e];
                first_added = false;
                n_added += 1;
            } else {
                xmin = std::min(xmin, m.edgesCentersX[e]);
                xmax = std::max(xmax, m.edgesCentersX[e]);
                y_moment += m.edgesCentersY[e];
                n_added += 1;
            }
        }
    }
    y_moment /= (double) n_added;
    const double x_moment = (xmax - xmin)*0.25 + xmin;
    wall.reserve(n_added);

    // compute cp
    for (uint b=0; b<m.boundaryEdges.size(); ++b) {
        if (m.boundaryEdgesPhysicals[b] == patch_name) {
            const uint e = m.boundaryEdges[b];
            const uint n0 = m.edgesNodes(e, 0);
            const uint n1 = m.edgesNodes(e, 1);

            const uint cell = m.edgesCells(e, 0);

            // Save wall x, y positions
            wall.x.push_back(m.edgesCentersX[e]);
            wall.y.push_back(m.edgesCentersY[e]);

            // Also calculate cp
            const double rho = sol.rho(cell);
            const double u = sol.u(cell);
            const double v = sol.v(cell);
            const double p = sol.p(cell);
            const double umag = sqrt(u*u + v*v);
            const double cs = sol.c(cell);
            const double mach = umag / cs;
            const double cp = 2./(sol.gamma()*mach_inf*mach_inf)*(p/p_inf - 1.);

            wall.cp.push_back(cp);

            // Integrate forces on profile
            const double ecx = m.edgesCentersX[e];
            const double ecy = m.edgesCentersY[e];
            const double fxi = cp * m.edgesNormalsX[e] * m.edgesLengths[e] / (xmax - xmin);
            const double fyi = cp * m.edgesNormalsY[e] * m.edgesLengths[e] / (xmax - xmin);
            const double mi = (ecx - x_moment)/(xmax - xmin)*fyi - (ecy - y_moment)/(xmax - xmin)*fxi;
            wall.cd += fxi;
            wall.cl += fyi;
            wall.cm -= mi;
        }
    }


    const double fx = wall.cd;
    const double fy = wall.cl;
    const double aoa = vars_far.angle;

    wall.cd =  fx*std::cos(aoa) + fy*std::sin(aoa);
    wall.cl = -fx*std::sin(aoa) + fy*std::cos(aoa);

    return wall;
}



void update_wall_profile(wallProfile &wall, solver& solv, std::string patch_name) {

    // First, get farfield conditions and compute infinity variables
    boundary_variables vars_far = solv.get_boundary_variables();
    const double p_inf = vars_far.p;
    const double mach_inf = vars_far.mach;

    solution sol = solv.get_solution();

    mesh& m = solv.get_mesh();

    // Find quarter chord position
    double xmin = 0.;
    double xmax = 0.;
    double y_moment = 0.;
    uint n_added = 0;
    bool first_added = true;
    for (uint b=0; b<m.boundaryEdges.size(); ++b) {
        if (m.boundaryEdgesPhysicals[b] == patch_name) {
            const uint e = m.boundaryEdges[b];
            if (first_added) {
                xmin = m.edgesCentersX[e];
                xmax = m.edgesCentersX[e];
                y_moment = m.edgesCentersY[e];
                first_added = false;
                n_added += 1;
            } else {
                xmin = std::min(xmin, m.edgesCentersX[e]);
                xmax = std::max(xmax, m.edgesCentersX[e]);
                y_moment += m.edgesCentersY[e];
                n_added += 1;
            }
        }
    }
    y_moment /= (double) n_added;
    const double x_moment = (xmax - xmin)*0.25 + xmin;

    // compute cp

    uint k = 0;
    for (uint b=0; b<m.boundaryEdges.size(); ++b) {
        if (m.boundaryEdgesPhysicals[b] == patch_name) {
            const uint e = m.boundaryEdges[b];
            const uint n0 = m.edgesNodes(e, 0);
            const uint n1 = m.edgesNodes(e, 1);

            const uint cell = m.edgesCells(e, 0);

            // Save wall x, y positions
            wall.x[k] = m.edgesCentersX[e];
            wall.y[k] = m.edgesCentersY[e];

            // Also calculate cp
            const double rho = sol.rho(cell);
            const double u = sol.u(cell);
            const double v = sol.v(cell);
            const double p = sol.p(cell);
            const double umag = sqrt(u*u + v*v);
            const double cs = sol.c(cell);
            const double mach = umag / cs;
            const double cp = 2./(sol.gamma()*mach_inf*mach_inf)*(p/p_inf - 1.);

            wall.cp[k] = cp;

            // Integrate forces on profile
            const double ecx = m.edgesCentersX[e];
            const double ecy = m.edgesCentersY[e];
            const double fxi = cp * m.edgesNormalsX[e] * m.edgesLengths[e] / (xmax - xmin);
            const double fyi = cp * m.edgesNormalsY[e] * m.edgesLengths[e] / (xmax - xmin);
            const double mi = (ecx - x_moment)/(xmax - xmin)*fyi - (ecy - y_moment)/(xmax - xmin)*fxi;
            wall.cd += fxi;
            wall.cl += fyi;
            wall.cm -= mi;

            k += 1;
        }
    }


    const double fx = wall.cd;
    const double fy = wall.cl;
    const double aoa = vars_far.angle;

    wall.cd =  fx*std::cos(aoa) + fy*std::sin(aoa);
    wall.cl = -fx*std::sin(aoa) + fy*std::cos(aoa);
}



}


