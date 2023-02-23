/*

        ____     ____  ___    _   _______
       /  _/    / __ \/   |  / | / / ___/
       / /_____/ /_/ / /| | /  |/ /\__ \
     _/ /_____/ _, _/ ___ |/ /|  /___/ /
    /___/    /_/ |_/_/  |_/_/ |_//____/

    Implicit RANS Solver
    Description: Mesh object and gmsh wrappers
    By: Alexis Angers

*/
#pragma once



#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <algorithm>

#include <fstream>
#include <stdexcept>

#include <rans/core.h>

namespace rans {



const uint MESH_EDGE_NULL = 4294967295;


// Cut string by character c
std::vector<std::string> cut_str(const std::string s, const char c) {
	std::vector<std::string> out;
	std::vector<int> I;
	// find zeros
	for (int i=0; i < s.size(); ++i) {
		if (s[i] == c) {
			I.push_back(i);
		}
	}
	// get each integer value
	if (I.size() == 0) {
		out.push_back( s );
	} else {
		out.push_back( s.substr(0, I[0]) );
		for (int j=0; j < (I.size()-1); ++j) {
			out.push_back( s.substr(I[j]+1, I[j+1]-I[j]-1) );
		}
		if (I[I.size()-1] != s.size()-1) {
			out.push_back( s.substr(I[I.size()-1]+1, s.size()-I[I.size()-1]-1) );
		}
	}
	return out;
}

std::vector<uint> str_to_ints(std::string s) {
	std::vector<uint> out;
	std::vector<uint> I;
	// find zeros
	for (uint i=0; i < s.size(); ++i) {
		if (s[i] == ' ') {
			I.push_back(i);
		}
	}
	// get each integer value
	if (I.size() == 0) {
		out.push_back( std::stoi(s) );
	} else {
		out.push_back( std::stoi(s.substr(0, I[0])) );
		for (uint j=0; j < (I.size()-1); ++j) {
			out.push_back( std::stoi(s.substr(I[j]+1, I[j+1]-I[j])) );
		}
		if (I[I.size()-1] != s.size()-1) {
			out.push_back( std::stoi(s.substr(I[I.size()-1]+1, s.size()-I[I.size()-1]-1)) );
		}
	}
	return out;
}

std::vector<double> str_to_floats(const std::string s) {
	std::vector<double> out;
	std::vector<int> I;
	// find zeros
	for (int i=0; i < s.size(); ++i) {
		if (s[i] == ' ') {
			I.push_back(i);
		}
	}
	// get each integer value
	if (I.size() == 0) {
		out.push_back( std::stod(s) );
	} else {
		out.push_back( std::stod(s.substr(0, I[0])) );
		for (int j=0; j < (I.size()-1); ++j) {
			out.push_back( std::stod(s.substr(I[j]+1, I[j+1]-I[j])) );
		}
		if (I[I.size()-1] != s.size()-1) {
			out.push_back( std::stod(s.substr(I[I.size()-1]+1, s.size()-I[I.size()-1]-1)) );
		}
	}
	return out;
}






template<class T>
void move_to_end(std::vector<T>& v, const uint& i) {
    // Move element i in vector to the end
    std::rotate(v.begin() + i,  v.begin() + i + 1, v.end());
}





template<uint N>
class meshArray {
private:
    std::vector<uint> nodes;
    uint n;
public:
    meshArray();
    uint& operator()(const uint& i, const uint& j);
    void push_back(const std::vector<uint> v);
    uint cols() const;
    uint rows() const;
    void dump();
    void swap(const uint& i, const uint& j);
    void move_to_end(const uint& i);
    const std::vector<uint>& get_vector();
};

template<uint N>
meshArray<N>::meshArray() {n=0;}

template<uint N>
uint& meshArray<N>::operator()(const uint& i, const uint& j) {
    return nodes[i*N + j];
}

template<uint N>
void meshArray<N>::push_back(const std::vector<uint> v) {
    if (v.size() != N) {
        throw std::invalid_argument( "size of input array " + std::to_string(v.size()) + " not equal to " + std::to_string(N) );
    }
    for (uint i=0; i<v.size(); ++i) {
        nodes.push_back(v[i]);
    }
    n += 1;
}

template<uint N>
uint meshArray<N>::cols() const {return n;};

template<uint N>
uint meshArray<N>::rows() const {return N;};

template<uint N>
void meshArray<N>::dump() {
    for (uint i=0; i<n; ++i) {
        for (uint j=0; j<N; ++j) {
            std::cout << this->operator()(i, j) << " ";
        }
        std::cout << "\n";
    }
}

template<uint N>
void meshArray<N>::swap(const uint& i, const uint& j) {
    // Swap row i and row j
    for (uint k=0; k<N; ++k) std::iter_swap(nodes.begin()+N*i+k, nodes.begin()+N*j+k);
}

template<uint N>
void meshArray<N>::move_to_end(const uint& i) {
    // Move row i to the end
    std::rotate(nodes.begin() + N*i, nodes.begin() + N*i + N, nodes.end());
}

template<uint N>
const std::vector<uint>& meshArray<N>::get_vector() {
    return nodes;
}



// Class for a partitionned mesh domain
class mesh {
public:

    std::string filename;

    std::map<std::string, std::string> meshFormat;

    std::map<uint, std::string> physicalNames;
    std::map<uint, uint> entityTagToPhysicalTag;
    std::map<std::string, uint> entitiesNumber;

    std::vector<double> nodesX;
    std::vector<double> nodesY;
    std::map<uint, uint> originalNodesRef;

    meshArray<2> edgesNodes;
    meshArray<2> edgesCells;
    std::vector<double> edgesLengths;
    std::vector<double> edgesNormalsX;
    std::vector<double> edgesNormalsY;
    std::vector<double> edgesCentersX;
    std::vector<double> edgesCentersY;

    std::vector<uint> boundaryEdges;
    std::vector<uint> boundaryEdges0;
    std::vector<uint> boundaryEdges1;
    std::vector<uint> boundaryEdgesIntTag;
    std::vector<std::string> boundaryEdgesPhysicals;

    std::map<std::tuple<uint, uint>, uint> edgesRef;

    meshArray<4> cellsNodes;
    meshArray<4> cellsEdges;
    std::vector<bool> cellsIsTriangle;
    std::vector<double> cellsAreas;
    std::vector<double> cellsCentersX;
    std::vector<double> cellsCentersY;
    std::vector<bool> cellsIsGhost;

    std::map<uint, uint> originalToCurrentCells;    // maps original index -> current index
    std::map<uint, uint> currentToOriginalCells;    // maps original index -> current index

    std::vector<uint> ghostCellsOriginalIndices;
    std::vector<uint> ghostCellsCurrentIndices;
    std::vector<uint> ghostCellsOwners;

    std::vector<double> wall_dist;

    uint nRealCells;


    mesh() {};
    mesh(std::string filename) {read_file(filename);}

    void read_entities();
    void read_nodes();
    void read_boundaries();
    void read_elements();
    void read_ghost_elements();

    void add_boundary_cells();

    void read_file(std::string filename);

    uint find_edge_with_nodes(const uint n0, const uint n1);
    uint find_bound_with_nodes(const uint n0, const uint n1);
    bool find_if_edge_in_mesh(const uint n0, const uint n1);

    void convert_node_face_info();
    void compute_mesh();
    void add_cell_edges(uint cell_id);

    void compute_wall_dist();

    void send_mesh_info();
};





}




namespace rans {





uint mesh::find_edge_with_nodes(const uint n0, const uint n1) {
    const uint nmin = std::min(n0, n1);
    const uint nmax = std::max(n0, n1);
    // Also check in bounds
    const auto tp = std::make_tuple(nmin, nmax);
    if (edgesRef.find(tp) != edgesRef.end()) {
        // Found
        return edgesRef.at(tp);
    } else {
        return MESH_EDGE_NULL;
    }
}


bool mesh::find_if_edge_in_mesh(const uint n0, const uint n1) {
    const uint nmin = std::min(n0, n1);
    const uint nmax = std::max(n0, n1);
    const auto tp = std::make_tuple(nmin, nmax);
    if (edgesRef.find(tp) != edgesRef.end()) {
        // Found
        return true;
    }
    return false;
}



void mesh::add_cell_edges(uint cell_id) {
    uint size = cellsIsTriangle[cell_id] ? 3 : 4;

    for (uint i=0; i<size; ++i) {
        uint j =  (i<(size-1)) ? (i+1) : 0;

        if (!find_if_edge_in_mesh(cellsNodes(cell_id, i), cellsNodes(cell_id, j))) {
            edgesNormalsX.push_back(0.);
            edgesNormalsY.push_back(0.);
            edgesCentersX.push_back(0.);
            edgesCentersY.push_back(0.);

            std::vector<uint> efi = {cell_id, 0};
            edgesCells.push_back(efi);

            std::vector<uint> eni = {cellsNodes(cell_id, i), cellsNodes(cell_id, j)};
            edgesNodes.push_back(eni);
            edgesLengths.push_back(0.);

            const uint nmin = std::min(cellsNodes(cell_id, i), cellsNodes(cell_id, j));
            const uint nmax = std::max(cellsNodes(cell_id, i), cellsNodes(cell_id, j));
            edgesRef.insert({ std::make_tuple(nmin, nmax), edgesLengths.size()-1 });
        }
    }
}


void mesh::convert_node_face_info() {
    // Take a mesh input, and convert its info to correct fvm data
    // Currently,
    //      edges contain no face connectivity data
    //      bounds contain no face connectivity data

    // Loop over all faces
    std::vector<uint> cellEdges(4);
    for (int i=0; i<cellsAreas.size(); ++i) {

		// Find edges
        cellEdges[0] = find_edge_with_nodes(cellsNodes(i, 0), cellsNodes(i, 1));
        cellEdges[1] = find_edge_with_nodes(cellsNodes(i, 1), cellsNodes(i, 2));
        if (cellsIsTriangle[i]) {
            cellEdges[2] = find_edge_with_nodes(cellsNodes(i, 2), cellsNodes(i, 0));
            cellEdges[3] = MESH_EDGE_NULL;
        } else {
            cellEdges[2] = find_edge_with_nodes(cellsNodes(i, 2), cellsNodes(i, 3));
            cellEdges[3] = find_edge_with_nodes(cellsNodes(i, 3), cellsNodes(i, 0));
        }
        cellsEdges.push_back(cellEdges);

        // Now update these edges with face connectivity info
        for (int ei : cellEdges) {
            if (ei != MESH_EDGE_NULL) {
                if (i != edgesCells(ei, 0)) {
                    edgesCells(ei, 1) = i;
                }
            }
        }
    }
}


void mesh::compute_mesh() {

	// Compute mesh normals, areas, lengths
    // Loop over each face

    #pragma omp parallel for
	for (int i=0; i<cellsAreas.size(); ++i) {
		// Evaluate face centroid position
        double cellC[2];
        cellC[0] = 0.; cellC[1] = 0.;

        uint n_cells = cellsIsTriangle[i] ? 3 : 4;
        for (uint j=0; j<n_cells; ++j) {
            cellC[0] += nodesX[cellsNodes(i, j)]/((double) n_cells);
            cellC[1] += nodesY[cellsNodes(i, j)]/((double) n_cells);
        }
        cellsCentersX[i] = cellC[0];
        cellsCentersY[i] = cellC[1];
    }

    // Evaluate edge normals and edge centers
	// Loop over each edge

    #pragma omp parallel for
	for (int i=0; i<edgesLengths.size(); ++i) {
        edgesCentersX[i] = (nodesX[edgesNodes(i, 1)] + nodesX[edgesNodes(i, 0)])*0.5;
        edgesCentersY[i] = (nodesY[edgesNodes(i, 1)] + nodesY[edgesNodes(i, 0)])*0.5;

        double dex = nodesX[edgesNodes(i, 1)] - nodesX[edgesNodes(i, 0)];
        double dey = nodesY[edgesNodes(i, 1)] - nodesY[edgesNodes(i, 0)];
        double l = sqrt(dex*dex + dey*dey);
		edgesLengths[i] = l;
		edgesNormalsX[i] = -dey/l;
        edgesNormalsY[i] =  dex/l;

        // Check normal direction
        double dot_f_f =
              edgesNormalsX[i]*(edgesCentersX[i] -  cellsCentersX[edgesCells(i, 0)])
            + edgesNormalsY[i]*(edgesCentersY[i] -  cellsCentersY[edgesCells(i, 0)]);
        if (dot_f_f < 0) {
            edgesNormalsX[i] *= -1.;
            edgesNormalsY[i] *= -1.;
        }
	}

    // Compute face area
    #pragma omp parallel for
    for (int i=0; i<cellsAreas.size(); ++i) {
        double& x1 = nodesX[cellsNodes(i, 0)];
        double& x2 = nodesX[cellsNodes(i, 1)];
        double& x3 = nodesX[cellsNodes(i, 2)];

        double& y1 = nodesY[cellsNodes(i, 0)];
        double& y2 = nodesY[cellsNodes(i, 1)];
        double& y3 = nodesY[cellsNodes(i, 2)];

        if (cellsIsTriangle[i]) {
            cellsAreas[i] = 0.5*abs(
                x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            );
        } else {
            double& x4 = nodesX[cellsNodes(i, 3)];
            double& y4 = nodesY[cellsNodes(i, 3)];

            cellsAreas[i] = 0.5*abs(
                x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            ) + 0.5*abs(
                x1*(y3-y4) + x3*(y4-y1) + x4*(y1-y3)
            );
        }

        if (cellsAreas[i] < 1e-15) {
            std::cout << "Null cell area for node " << i << " ";
        }
    }
}



void mesh::read_entities() {
    std::ifstream infile(filename);
    std::string line;

    uint n = 0;
    uint ns = 0;
    uint nss = 0;

    uint nGhostEntities;

    std::string currentSection = "";

    while (std::getline(infile, line)) {

        if (line[0] == '$') {
            // This line tells us a section number
            currentSection = line.substr(1, line.size());
            // Restart current section line
            ns = 0;
        } else if (currentSection == "MeshFormat") {
            // Read mesh format
            if (ns == 1) {
                meshFormat["format"] = line;
            }
        } else if (currentSection == "PhysicalNames") {
            if (ns > 1) {
                auto lineS = cut_str(line, ' ');
                std::string pnamei = lineS[2];
                pnamei.erase(remove(pnamei.begin(), pnamei.end(), '"'), pnamei.end());
                physicalNames[std::stoi(lineS[1]) + 1000*std::stoi(lineS[0])] = pnamei;
            }
        } else if (currentSection == "Entities") {
            // Entities are nodes on border elements
            if (ns == 1) {
                // First line, contains number of entities
                auto l = str_to_ints(line);
                entitiesNumber["points"] = l[0];
                entitiesNumber["curves"] = l[1] + l[0];
                entitiesNumber["surfaces"] = l[2] + l[1] + l[0];
                entitiesNumber["volumes"] = l[3] + l[2] + l[1] + l[0];

                nss = 0;
            } else if (nss <= entitiesNumber["points"]) {
                // Reading points
            } else if (nss <= entitiesNumber["curves"]) {
                // Reading curves
                auto l = str_to_ints(line);
                entityTagToPhysicalTag[l[0]] = l[8] + 1000*1;
            } else if (nss <= entitiesNumber["surfaces"]) {
                // Reading surfaces
                auto l = str_to_ints(line);
            }
        } else if (currentSection == "Nodes") {
            // Do not read nodes, break
            break;
        }

        n += 1;
        ns += 1;
        nss += 1;
    }

    infile.close();
}



void mesh::read_nodes() {

    std::ifstream infile(filename);
    std::string line;

    uint n = 0;
    uint ns = 0;
    uint nss = 0;

    uint start_offset = 0;
    uint n_in_block = 0;

    std::string currentSection = "";

    std::vector<uint> block_tags;

    while (std::getline(infile, line)) {

        if (line[0] == '$') {
            // This line tells us a section number
            currentSection = line.substr(1, line.size());
            // Restart current section line
            ns = 0;
        } else if (currentSection == "Nodes") {
            // Read nodes
            if (ns == 1) {
                // Read number of node blocks
                start_offset = 1;
                nss = 0;
            } else if (ns > (start_offset + 2*n_in_block)) {
                // Read block
                auto l = str_to_ints(line);
                n_in_block = l[3];
                start_offset = ns;
                block_tags.resize(n_in_block);
                nss = 0;
            } else if (nss <= n_in_block) {
                // Read current node tag
                block_tags[nss-1] = std::stoi(line) - 1;
            } else {
                // Read current node
                uint tag_i = block_tags[nss - 1 - n_in_block];
                auto l = str_to_floats(line);
                originalNodesRef[tag_i] = nodesX.size();
                nodesX.push_back(l[0]);
                nodesY.push_back(l[1]);
            }
        } else if (currentSection == "Elements") {
            break;
        }

        n += 1;
        ns += 1;
        nss += 1;
    }

    infile.close();
}




void mesh::read_boundaries() {

    std::ifstream infile(filename);
    std::string line;

    uint n = 0;
    uint ns = 0;
    uint nss = 0;

    uint start_offset = 0;
    uint n_in_block = 0;

    uint blockDimension = 0;
    uint blockPhysicalTag;

    std::string currentSection = "";

    while (std::getline(infile, line)) {

        if (line[0] == '$') {
            // This line tells us a section number
            currentSection = line.substr(1, line.size());
            // Restart current section line
            ns = 0;
        } else if (currentSection == "Elements") {
            // Read elements
            if (ns == 1) {
                // Read number of node blocks
                start_offset = 1;
                n_in_block = 0;
                nss = 0;
            } else if (ns > (start_offset + n_in_block)) {
                // Read block
                auto l = str_to_ints(line);
                blockDimension = l[0];

                if (blockDimension == 1) {
                    blockPhysicalTag = entityTagToPhysicalTag.at(l[1]);
                }
                n_in_block = l[3];
                start_offset = ns;
                nss = 0;
            } else if (blockDimension == 1) {
                // Read current tag
                auto li = str_to_ints(line);
                uint tag_i = li[0] - 1;
                std::vector<uint> l(li.size()-1);
                for (uint i=1; i<li.size(); ++i) {
                    l[i-1] = li[i] - 1;
                }
                if (l.size() == 2) {
                    // Boundary edge
                    boundaryEdges0.push_back(originalNodesRef.at(l[0]));
                    boundaryEdges1.push_back(originalNodesRef.at(l[1]));
                    boundaryEdgesIntTag.push_back(blockPhysicalTag);
                }
            }
        }

        n += 1;
        ns += 1;
        nss += 1;
    }

    infile.close();
}



void mesh::read_elements() {

    std::ifstream infile(filename);
    std::string line;

    uint n = 0;
    uint ns = 0;
    uint nss = 0;

    uint start_offset = 0;
    uint n_in_block = 0;

    uint blockDimension = 0;
    uint blockPhysicalTag;

    std::string currentSection = "";

    while (std::getline(infile, line)) {

        if (line[0] == '$') {
            // This line tells us a section number
            currentSection = line.substr(1, line.size());
            // Restart current section line
            ns = 0;
        } else if (currentSection == "Elements") {
            // Read elements
            if (ns == 1) {
                // Read number of node blocks
                start_offset = 1;
                n_in_block = 0;
                nss = 0;
            } else if (ns > (start_offset + n_in_block)) {
                // Read block
                auto l = str_to_ints(line);
                blockDimension = l[0];
                if (blockDimension == 1) {
                    blockPhysicalTag = entityTagToPhysicalTag.at(l[1]);
                }
                n_in_block = l[3];
                start_offset = ns;
                nss = 0;
            } else if (blockDimension == 2) {
                // Read current tag
                auto li = str_to_ints(line);
                uint tag_i = li[0] - 1;
                std::vector<uint> l(li.size()-1);
                for (uint i=1; i<li.size(); ++i) {
                    l[i-1] = li[i] - 1;
                }
                if (l.size() != 2) {
                    // Triangle or quad cell
                    currentToOriginalCells[cellsIsTriangle.size()] = tag_i;
                    originalToCurrentCells[tag_i] = cellsIsTriangle.size();

                    uint l_size = l.size();

                    for (uint i=0; i<l.size(); ++i) {
                        l[i] = originalNodesRef.at(l[i]);
                    }

                    if (l_size == 3) {
                        l.push_back(0);
                        cellsIsTriangle.push_back(true);
                    } else {
                        cellsIsTriangle.push_back(false);
                    }

                    cellsNodes.push_back(l);

                    cellsAreas.push_back(0.);
                    cellsCentersX.push_back(0.);
                    cellsCentersY.push_back(0.);
                    cellsIsGhost.push_back(false);
                }
            }
        }

        n += 1;
        ns += 1;
        nss += 1;
    }

    infile.close();
}





void mesh::add_boundary_cells() {
    // Add boundary cells
    for (uint i=0; i<boundaryEdges0.size(); ++i) {

        // Find edge index on boundary
        const uint nmin = std::min(boundaryEdges0[i], boundaryEdges1[i]);
        const uint nmax = std::max(boundaryEdges0[i], boundaryEdges1[i]);
        std::tuple<uint, uint> tp = std::make_tuple(nmin, nmax);
        uint e = 0;
        if (edgesRef.find(tp) != edgesRef.end()) {
            // Found
            e = edgesRef.at(tp);
        } else {
            throw std::invalid_argument( "invalid edge ref (" + std::to_string(nmin) + ", " + std::to_string(nmax) + ")" );
        }

        uint cell = edgesCells(e, 0);

        // Compute virtual cell properties
        double area = cellsAreas[cell];
        double dx = edgesCentersX[e] - cellsCentersX[cell];
        double dy = edgesCentersY[e] - cellsCentersY[cell];
        double dist = sqrt(dx*dx + dy*dy);
        double cx = edgesCentersX[e] + dist*edgesNormalsX[e];
        double cy = edgesCentersY[e] + dist*edgesNormalsY[e];
        std::vector<uint> thisCellNodes = {nmin, nmax, 0, 0};

        // Add virtual cell
        cellsAreas.push_back(area);
        cellsCentersX.push_back(cx);
        cellsCentersY.push_back(cy);
        cellsIsTriangle.push_back(true);
        cellsNodes.push_back(thisCellNodes);

        // Add cell index to edge
        edgesCells(e, 1) = cellsAreas.size() - 1;

        // Add edge index to boundary edges
        boundaryEdges.push_back(e);
        boundaryEdgesPhysicals.push_back(
            physicalNames.at(boundaryEdgesIntTag[i])
        );
    }
}






void mesh::compute_wall_dist() {


    // Compute the distance to the wall
    #pragma omp parallel for
    for (int i=0; i<cellsAreas.size(); ++i) {
        double& cx = cellsCentersX[i];
        double& cy = cellsCentersY[i];

        double mind = 1;
        // Loop over all boundary cells
        for (uint j=0; j<boundaryEdges.size(); ++j) {
            const int& e = boundaryEdges[j];

            double& wx = edgesCentersX[e];
            double& wy = edgesCentersY[e];

            double dx = cx - wx;
            double dy = cy - wy;
            double d = sqrt(dx*dx + dy*dy);

            if (j == 0) {
                mind = d;
            } else {
                mind = std::min(mind, d);
            }
        }
        wall_dist[i] = mind;
    }
}



void mesh::read_file(std::string filename_in) {
    filename = filename_in;

    // Read all contents of the file
    read_entities();
    read_nodes();
    read_boundaries();
    read_elements();

    // Add all cells edges
    for (uint i=0; i<cellsAreas.size(); ++i) {
        add_cell_edges(i);
    }

    // Fix edges part of ghost cells

    #pragma omp parallel for
    for (int i=0; i<cellsAreas.size(); ++i) {
        if (cellsIsGhost[i]) {
            const uint cell_size = cellsIsTriangle[i] ? 3 : 4;
            for (uint j=0; j<cell_size; ++j) {
                uint k = (j==(cell_size-1)) ? 0 : j+1;
                uint n0 = cellsNodes(i, j);
                uint n1 = cellsNodes(i, k);

                const uint nmin = std::min(n0, n1);
                const uint nmax = std::max(n0, n1);
                const auto tp = std::make_tuple(nmin, nmax);
                const uint e = edgesRef.at(tp);

                edgesCells(e, 1) = edgesCells(e, 0);
            }
        }
    }

    nRealCells = cellsAreas.size();

    // Convert nodes and faces info
    convert_node_face_info();

    // Compute mesh metrics
    compute_mesh();

    // Add boundary cells
    add_boundary_cells();

    // Compute wall distance
    wall_dist.resize(cellsAreas.size());
    for (auto& i : wall_dist) i = 0;

}




}





