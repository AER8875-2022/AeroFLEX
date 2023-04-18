#include "geometrie/mesh_vlm.hpp"
#include "geometrie/geometry.hpp"
#include <algorithm>

void remove(std::vector<int> &v) {
    auto end = v.end();
    for (auto it = v.begin(); it != end; ++it) {
        end = std::remove(it + 1, end, *it);
    }
    v.erase(end, v.end());
};

std::vector<double> sum(std::vector<double> a, std::vector<double> b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
};

std::vector<double> substract(std::vector<double> v_A, std::vector<double> v_B) {
    std::vector<double> c_P = {0, 0, 0};
    c_P[0] = v_A[0] - v_B[0];
    c_P[1] = v_A[1] - v_B[1];
    c_P[2] = v_A[2] - v_B[2];
    return c_P;
};

std::vector<std::vector<double>> sum_matrices(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B) {
    std::vector<std::vector<double>> C;
    for (int i=0; i < A.size(); i++) {
        std::vector<double> C_i;
        for (int j=0; j < A[i].size(); j++) {
            C_i.push_back(A[i][j] + B[i][j]);
        };
        C.push_back(C_i);
    };
    return C;
};

std::vector<std::vector<double>> scale_matrices(double s, std::vector<std::vector<double>> A) {
    std::vector<std::vector<double>> C;
    for (int i=0; i < A.size(); i++) {
        std::vector<double> C_i;
        for (int j=0; j < A[i].size(); j++) {
            C_i.push_back(s * A[i][j]);
        };
        C.push_back(C_i);
    };
    return C;
};

std::vector<double> scaleVector(double s, std::vector<double> a) {
    std::vector<double> c;
    for (int i=0; i < a.size(); i++) {c.push_back(s * a[i]);};
    return c;
};

double dotProduct(std::vector<double> v_A, std::vector<double> v_B) {
    double c_P = 0;
    c_P += v_A[0] * v_B[0];
    c_P += v_A[1] * v_B[1];
    c_P += v_A[2] * v_B[2];
    return c_P;
};

std::vector<double> crossProduct(std::vector<double> v_A, std::vector<double> v_B) {
    std::vector<double> c_P = {0, 0, 0};
    c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
    c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
    c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
    return c_P;
};

double sign(double a) {
    if (a > 0) {return 1;}
    else if (a < 0) {return -1;}
    else {return 0;};
};

double signed_tetra_volume(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d) {
    return sign(dotProduct(crossProduct(substract(b, a), substract(c, a)), substract(d, a)) / 6.0);
};

std::vector<double> get_intersection(std::vector<double> a, std::vector<double> b, std::vector<double> x, std::vector<double> y, std::vector<double> z) {
    std::vector<double> n = crossProduct(substract(y, x), substract(z, x));
    double t = dotProduct(substract(x, a), n) / dotProduct(substract(b, a), n);
    return sum(a, scaleVector(t, substract(b, a)));
};

bool check_intersection(std::vector<double> a, std::vector<double> b, std::vector<double> x, std::vector<double> y, std::vector<double> z) {
    if (dotProduct(substract(b, a), crossProduct(substract(y, x), substract(z, x)))!=0) {
        double s3 = signed_tetra_volume(a, b, x, y);
        double s4 = signed_tetra_volume(a, b, y, z);
        double s5 = signed_tetra_volume(a, b, z, x);
        if ((s3 == s4) && (s4 == s5)) {return true;}
        else {return false;};
    } else {return false;};
};

double dis(std::vector<double> m, std::vector<double> n, std::vector<double> o, std::vector<double> p) {
    return (m[0] - n[0])*(o[0] - p[0]) + (m[1] - n[1])*(o[1] - p[1]) + (m[2] - n[2])*(o[2] - p[2]);
};

std::vector<double> lines_intersection(std::vector<double> p1, std::vector<double> p2, std::vector<double> p3, std::vector<double> p4) {
    double tol = 1e-14;
    double tol_2 = tol * tol;
    //mua = ( d1343 d4321 - d1321 d4343 ) / ( d2121 d4343 - d4321 d4321 )
    //mub = ( d1343 + mua d4321 ) / d4343
    double a = dis(p2, p1, p2, p1) * dis(p4, p3, p4, p3) - dis(p4, p3, p2, p1) * dis(p4, p3, p2, p1);
    double b = dis(p4, p3, p4, p3);
    if (a != 0 && b != 0) {
        double mua = (dis(p1, p3, p4, p3) * dis(p4, p3, p2, p1) - dis(p1, p3, p2, p1) * dis(p4, p3, p4, p3)) / a;
        double mub = (dis(p1, p3, p4, p3) + mua * dis(p4, p3, p2, p1)) / b;
        bool are_end_points = (((mua-0)*(mua-0)<=tol_2)||((mua-1)*(mua-1)<=tol_2)) && (((mub-0)*(mub-0)<=tol_2)||((mub-1)*(mub-1)<=tol_2));
        if ((mua >= 0 && mua <= 1 && mub >= 0 && mub <= 1) && !are_end_points) {
            //(mua != 0 && mua != 1) && (mub != 0 && mub != 1);
            std::vector<double> pa = sum(p1, scaleVector(mua, substract(p2, p1)));
            std::vector<double> pb = sum(p3, scaleVector(mub, substract(p4, p3)));
            if (dis(pa, pb, pa, pb) <= tol_2) {
                return pa;
            };
        };
    };
    return {};
};

bool is_point_on_vertex(std::vector<double> a, std::vector<double> p0, std::vector<double> p1, std::vector<double> p2, double tol) {
    if (dotProduct(substract(a, p0), substract(a, p0)) <= tol * tol) {return true;};
    if (dotProduct(substract(a, p1), substract(a, p1)) <= tol * tol) {return true;};
    if (dotProduct(substract(a, p2), substract(a, p2)) <= tol * tol) {return true;};
    return false;
};

Box::Box(std::vector<int> box_id, std::vector<std::vector<double>> list_of_points) {
    this->id = box_id;
    this->points = list_of_points;
    this->min = {1e15, 1e15, 1e15};
    this->max = {-1e15, -1e15, -1e15};
    this->set_boundaries();
};

void Box::set_boundaries() {
    for (int i=0; i < this->points.size(); i++) {
        // Set min
        if (this->points[i][0] < this->min[0]) {this->min[0] = this->points[i][0];};
        if (this->points[i][1] < this->min[1]) {this->min[1] = this->points[i][1];};
        if (this->points[i][2] < this->min[2]) {this->min[2] = this->points[i][2];};
        // Set max
        if (this->points[i][0] > this->max[0]) {this->max[0] = this->points[i][0];};
        if (this->points[i][1] > this->max[1]) {this->max[1] = this->points[i][1];};
        if (this->points[i][2] > this->max[2]) {this->max[2] = this->points[i][2];};
    };
};

std::vector<std::vector<double>> Box::get_points() {return this->points;};
std::vector<int> Box::get_id() {return this->id;};

bool Box::box_touches_box(Box box) {
    bool x_touches = false;
    bool y_touches = false;
    bool z_touches = false;
    if (!(this->min[0] > box.max[0] || this->max[0] < box.min[0])) {x_touches = true;};
    if (!(this->min[1] > box.max[1] || this->max[1] < box.min[1])) {y_touches = true;};
    if (!(this->min[2] > box.max[2] || this->max[2] < box.min[2])) {z_touches = true;};

    if (x_touches && y_touches && z_touches) {
        return true;
    } else {
        return false;
    };
};

bool Box::is_box_inside_box(Box box) {
    bool is_inside = true;
    if (this->min[0] == this->max[0] || this->min[1] == this->max[1] || this->min[2] == this->max[2]) {
        is_inside = false;
    } else {
        for (int i=0; i < box.get_points().size(); i++) {
            double u = (box.get_points()[i][0] - this->min[0]) / (this->max[0] - this->min[0]);
            double v = (box.get_points()[i][1] - this->min[1]) / (this->max[1] - this->min[1]);
            double w = (box.get_points()[i][2] - this->min[2]) / (this->max[2] - this->min[2]);
            if (!((1. >= u >= 0.) && (1. >= v >= 0.) && (1. >= w >= 0.))) {
                is_inside = false;
            };
        };
    };
    return is_inside;
};

/*
    def is_object_inside_box(self, point_object):
        is_inside = True
        if self.min[0] == self.max[0] or self.min[1] == self.max[1] or self.min[2] == self.max[2]:
            is_inside = False
        else:
            for point in point_object.points:
                u = (point[0] - self.min[0]) / (self.max[0] - self.min[0])
                v = (point[1] - self.min[1]) / (self.max[1] - self.min[1])
                w = (point[2] - self.min[2]) / (self.max[2] - self.min[2])
                if (1 >= u >= 0) and (1 >= v >= 0) and (1 >= w >= 0):
                    pass
                else:
                    is_inside = False
                    break
        return is_inside
*/




Vertex::Vertex() {
    this->point = {};
};
Vertex::Vertex(std::vector<double> point) {
    this->point = point;
};
std::vector<double> Vertex::get_point() {return this->point;};


Vortex::Vortex() {
    this->type = "Vortex";
    this->id = -1;
    this->surface = -1;
    this->station = -1;
    this->row = -1;
    // A vortex always has 4 vertices
    this->vertices = {};
};
Vortex::Vortex(std::vector<int> vertices, int surface, int station, int row) {
    this->type = "Vortex";
    // A vortex always has 4 vertices
    this->vertices = vertices;
    this->id = -1;
    this->surface = surface;
    this->station = station;
    this->row = row;
};
std::vector<int> Vortex::get_vertices() {return this->vertices;};
int Vortex::get_id() {return this->id;};
int Vortex::get_surface() {return this->surface;};
int Vortex::get_station() {return this->station;};
int Vortex::get_row() {return this->row;};
void Vortex::set_id(int id) {this->id = id;};
void Vortex::set_vertex(std::vector<int> new_ids) {this->vertices = new_ids;};


Doublet::Doublet() {
    this->type = "Doublet";
    this->vertices = {};
    this->id = -1;
    this->surface = -1;
};
Doublet::Doublet(std::vector<int> vertices, int surface) {
    this->type = "Doublet";
    this->vertices = vertices;
    this->id = -1;
    this->surface = surface;
};
std::vector<int> Doublet::get_vertices() {return this->vertices;};
int Doublet::get_id() {return this->id;};
int Doublet::get_surface() {return this->surface;};
void Doublet::set_id(int id) {this->id = id;};
void Doublet::set_vertex(std::vector<int> new_ids) {this->vertices = new_ids;};


Wingstation::Wingstation() {
    this->elements = {};
    this->id = -1;
    this->surface = -1;
};
Wingstation::Wingstation(int wingstation_id, int surface, std::vector<int> elements) {
    this->elements = elements;
    this->id = wingstation_id;
    this->surface = surface;
};
std::vector<int> Wingstation::get_elements() {return this->elements;};
int Wingstation::get_id() {return this->id;};
int Wingstation::get_surface() {return this->surface;};
void Wingstation::push_vortex(int vortex_id) {this->elements.push_back(vortex_id);}; 


WingSurface::WingSurface() {
    this->id = -1;
    this->type = "Undef";
    this->wingstations = {};
};
WingSurface::WingSurface(int surface_id, std::string surface_type) {
    this->id = surface_id;
    this->type = surface_type;
    this->wingstations = {};
};
int WingSurface::get_id() {return this->id;};
std::string WingSurface::get_type() {return this->type;};
std::vector<int> WingSurface::get_wingstations() {return this->wingstations;};
void WingSurface::push_wingstation(int wingstation_id) {this->wingstations.push_back(wingstation_id);};
void WingSurface::reverse_wingstations() {std::reverse(this->wingstations.begin(), this->wingstations.end());};

Mesh::Mesh() {
    this->vertices = {};
    // this->elements = {};
    this->vortexes = {};
    this->doublets = {};
    this->wingstations = {};
    this->surfaces = {};
    this->vortex_count = 0;
    this->doublet_count = 0;
    this->wingstation_count = 0;
    this->surface_count = 0;  
};

int Mesh::append_vertex(Vertex vertex) {
    if (vertex.get_point().size() == 3) {
        int new_vertex_id = this->vertices.size();
        this->vertices.push_back(vertex);
        return new_vertex_id;
    } else {return -1;};
};

void Mesh::create_vortex(std::vector<int> vertices, int surface, int station, int row) {
    // Assign a unique vortex_id to element, then update mesh vortex_count
    int new_vortex_id = this->vortex_count;
    Vortex new_vortex = Vortex(vertices, surface, station, row);
    new_vortex.set_id(new_vortex_id);
    this->vortex_count += 1;

    // Append Vortex object to mesh elements list
    this->vortexes.push_back(new_vortex);
    //this->vortexes[new_vortex_id] = new_vortex;
    if (this->wingstation_count > station) {
        if (this->wingstations[station].get_surface() == surface) {
                this->wingstations[station].push_vortex(new_vortex_id);
        };
    };
};

void Mesh::create_doublet(std::vector<int> vertices, int surface) {
    // Assign a unique doublet_id to element, then update mesh doublet_count
    int new_doublet_id = this->doublet_count;
    Doublet new_doublet = Doublet(vertices, surface);
    new_doublet.set_id(new_doublet_id);
    this->doublet_count += 1;
    // Append Doublet object to mesh elements list
    this->doublets.push_back(new_doublet);
    //this->doublets[new_doublet_id] = new_doublet;
};

void Mesh::create_new_surface(std::string surface_type) {
    int new_surface_id = this->surface_count;
    WingSurface new_surface = WingSurface(new_surface_id, surface_type);
    this->surfaces.push_back(new_surface);
    this->surface_count += 1;
};

void Mesh::create_new_wingstation(int surface_id, std::vector<int> vortexes) {
    int new_wingstation_id = this->wingstation_count;
    Wingstation new_wingstation = Wingstation(new_wingstation_id, surface_id, vortexes);
    this->wingstations.push_back(new_wingstation);
    this->wingstation_count += 1;
    this->surfaces[surface_id].push_wingstation(new_wingstation_id);
};

void Mesh::push_vortex_to_wingstation(int vortex_id, int wingstation_id) {
    this->wingstations[wingstation_id].push_vortex(vortex_id);
};

std::vector<Vertex> Mesh::get_vertices(std::vector<int> vertex_ids) {
    std::vector<Vertex> vertex_list;
    for (int i=0; i < vertex_ids.size(); i++) {
        vertex_list.push_back(this->vertices[vertex_ids[i]].get_point());
    };
    return vertex_list;
};

std::vector<Vertex> Mesh::get_vertices() {return this->vertices;};
std::vector<WingSurface> Mesh::get_surfaces() {return this->surfaces;};
std::vector<Wingstation> Mesh::get_wingstations() {return this->wingstations;};
//std::unordered_map<int, Vortex> Mesh::get_vortexes() {return this->vortexes;};
//std::unordered_map<int, Doublet> Mesh::get_doublets() {return this->doublets;};
std::vector<Vortex> Mesh::get_vortexes() {return this->vortexes;};
std::vector<Doublet> Mesh::get_doublets() {return this->doublets;};
int Mesh::get_surface_count() {return this->surface_count;};
int Mesh::get_wingstation_count() {return this->wingstation_count;};
void Mesh::reverse_wingstations(int surface_id) {this->surfaces[surface_id].reverse_wingstations();};

void Mesh::change_doublet_vertices(int id, std::vector<int> new_vertices) {this->doublets[id].set_vertex(new_vertices);};
void Mesh::change_vortex_vertices(int id, std::vector<int> new_vertices) {this->vortexes[id].set_vertex(new_vertices);};


std::vector<std::vector<double>> Mesh::get_vortexes_normals() {
    std::vector<std::vector<double>> normals;
    // Get normals of vortexes
    for (int i=0; i < this->vortexes.size(); i++) {
        std::vector<int> indexes = this->vortexes[i].get_vertices();
        std::vector<double> a = this->vertices[indexes[0]].get_point();
        std::vector<double> b = this->vertices[indexes[1]].get_point();
        std::vector<double> c = this->vertices[indexes[2]].get_point();
        std::vector<double> norm = crossProduct(substract(a, b), substract(a, c));
        double norm_length = sqrt(dotProduct(norm, norm));
        if (norm_length > 0) {
            for (int j=0; j < norm.size(); j++) {norm[j] *= 1 / norm_length;};
        } else {norm = {0.0, 0.0, 0.0};};
        normals.push_back(norm);
    };
    return normals;
};

std::vector<std::vector<double>> Mesh::get_doublets_normals() {
    std::vector<std::vector<double>> normals;
    // Get normals of vortexes
    for (int i=0; i < this->doublets.size(); i++) {
        std::vector<int> indexes = this->doublets[i].get_vertices();
        std::vector<double> a = this->vertices[indexes[0]].get_point();
        std::vector<double> b = this->vertices[indexes[1]].get_point();
        std::vector<double> c = this->vertices[indexes[2]].get_point();
        std::vector<double> norm = crossProduct(substract(a, b), substract(a, c));
        double norm_length = sqrt(dotProduct(norm, norm));
        if (norm_length > 0) {
            for (int j=0; j < norm.size(); j++) {norm[j] *= 1 / norm_length;};
        } else {norm = {0.0, 0.0, 0.0};};
        normals.push_back(norm);
    };
    return normals;
};

void Mesh::delete_vortex(int vortex_id) {
    // Remove element
    //this->vortexes.erase(vortex_id);
    this->vortexes.erase(this->vortexes.begin() - 1 + vortex_id);
    // Update vortex count and IDs
    this->vortex_count += -1;
    for (int i=0; i < this->vortexes.size(); i++) {
        int current_vortex_id = this->vortexes[i].get_id();
        if (current_vortex_id >= vortex_id) {this->vortexes[i].set_id(current_vortex_id - 1);};
    };
};

void Mesh::delete_doublet(int doublet_id) {
    // Remove element
    //this->doublets.erase(doublet_id);
    this->doublets.erase(this->doublets.begin() - 1 + doublet_id);
    // Update doublet count and IDs
    this->doublet_count += -1;
    for (int i=0; i < this->doublets.size(); i++) {
        int current_doublet_id = this->doublets[i].get_id();
        if (current_doublet_id >= doublet_id) {this->doublets[i].set_id(current_doublet_id - 1);};
    };
};


std::vector<int> Mesh::get_unused_vertices() {
    // Initialize vector of zeros
    std::vector<int> vertices_usage;
    vertices_usage.resize(this->vertices.size());
    std::fill(vertices_usage.begin(), vertices_usage.end(), 0);
    // Look through vortexes
    for (int i=0; i < this->vortexes.size(); i++) {
        std::vector<int> vortex_vertices = this->vortexes[i].get_vertices();
        for (int j=0; j < vortex_vertices.size(); j++) {vertices_usage[vortex_vertices[j]] += 1;};
    };
    // Look through doublets
    for (int i=0; i < this->doublets.size(); i++) {
        std::vector<int> doublet_vertices = this->doublets[i].get_vertices();
        for (int j=0; j < doublet_vertices.size(); j++) {vertices_usage[doublet_vertices[j]] += 1;};
    };
    // Order 
    std::vector<int> unused_vertices;
    for (int i=0; i < vertices_usage.size(); i++) {if (vertices_usage[i] == 0) {unused_vertices.push_back(i);};};
    return unused_vertices;
};

/*
std::vector<int> Mesh::get_unused_vertices() {
    std::vector<int> vertices_usage;
    for (int i=0; i < this->vertices.size(); i++) {
        vertices_usage.push_back(0);
    };

    for (int i=0; i < this->doublets.size(); i++) {
        for (int j=0; j < this->doublets[i].get_vertices().size(); j++) {
            vertices_usage[this->doublets[i].get_vertices()[j]] += 1;
        };
    };
    std::vector<int> unused_vertices;
    for (int i=0; i < vertices_usage.size(); i++) {
        if (vertices_usage[i] == 0) {
            unused_vertices.push_back(i);
        };
    };
    return unused_vertices;
};
*/

std::vector<int> Mesh::get_transmute_vertices() {
    
    std::vector<int> unused_vertices = this->get_unused_vertices();
    int vertex_number = this->vertices.size();
    std::vector<int> transmute_vertex_to;
    for (int i=0; i < this->vertices.size(); i++) {
        transmute_vertex_to.push_back(i);
    };

    for (int i=0; i < vertex_number; i++) {
        for (int j=0; j < unused_vertices.size(); j++) {
            if (unused_vertices[j] <= i) {
                transmute_vertex_to[i] += -1;
            };
        };
    };
    return transmute_vertex_to;
};

void Mesh::delete_unused_vertices() {
    
    std::vector<int> unused_vertices = this->get_unused_vertices();
    std::vector<int> transmute_vertex_to = this->get_transmute_vertices();
    //int vertex_count = this->vertices.size();

    // re-index vertices to avoid gaps in ordering
    std::vector<int> new_vertices_index;
    for (int i=0; i < this->vertices.size(); i++) {
        new_vertices_index.push_back(i);
    };

    // Remove transmuted elements from vertex_table
    for (int i=0; i < unused_vertices.size(); i++) {
        // this->vertices.pop(unused_vertices[unused_vertices.size() - 1 - i]);
        this->vertices.erase(this->vertices.begin() + unused_vertices[unused_vertices.size() - 1 - i]);
    };

    // Transmute vertices in vortexes
    for (int i=0; i < this->vortexes.size(); i++) {
        std::vector<int> new_vertices = this->vortexes[i].get_vertices();
        for (int j=0; j < this->vortexes[i].get_vertices().size(); j++) {
            new_vertices[j] = transmute_vertex_to[this->vortexes[i].get_vertices()[j]];
        };
        this->vortexes[i].set_vertex(new_vertices);
    };
    // Transmute vertices in doublets
    for (int i=0; i < this->doublets.size(); i++) {
        std::vector<int> new_vertices = this->doublets[i].get_vertices();
        for (int j=0; j < this->doublets[i].get_vertices().size(); j++) {
            new_vertices[j] = transmute_vertex_to[this->doublets[i].get_vertices()[j]];
        };
        this->doublets[i].set_vertex(new_vertices);
    };
};


void Mesh::set_used_vertices_in_containers(int max_subdivisions, double tol) {
    /*
    # Extract list of PointObjects from mesh -> elements
    list_of_boxes = self.convert_used_vertices_to_box_list(tol)
    # Create container structure
    container = generate_containers(None, list_of_boxes, max_subdivisions)
    return container
    */
};


std::vector<Box> Mesh::convert_used_vertices_to_box_list(double tol) {
    std::vector<Box> list_of_boxes;

    // Run through all vortexes
    for (int i=0; i < this->vortexes.size(); i++) {
        for (int j=0; j < vortexes[i].get_vertices().size(); j++) {
            int vertex_id = this->vortexes[i].get_vertices()[j];
            std::vector<double> point = this->vertices[vertex_id].get_point();
            std::vector<std::vector<double>> list_of_points;
            if (tol > 0) {
                std::vector<double> point_min = {point[0] - tol, point[1] - tol, point[2] - tol};
                std::vector<double> point_max = {point[0] + tol, point[1] + tol, point[2] + tol};
                list_of_points = {point_min, point_max};
            } else {
                list_of_points = {point};
            };
            list_of_boxes.push_back(Box({vertex_id, i, j, 1}, list_of_points));
        };
    };
    // Run through all doublets
    for (int i=0; i < this->doublets.size(); i++) {
        for (int j=0; j < doublets[i].get_vertices().size(); j++) {
            int vertex_id = this->doublets[i].get_vertices()[j];
            std::vector<double> point = this->vertices[vertex_id].get_point();
            std::vector<std::vector<double>> list_of_points;
            if (tol > 0) {
                std::vector<double> point_min = {point[0] - tol, point[1] - tol, point[2] - tol};
                std::vector<double> point_max = {point[0] + tol, point[1] + tol, point[2] + tol};
                list_of_points = {point_min, point_max};
            } else {
                list_of_points = {point};
            };
            list_of_boxes.push_back(Box({vertex_id, i, j, 0}, list_of_points));
        };
    };
    return list_of_boxes;
};

/*
void Mesh::stitch_used_vertices(double tol, int nb_subdivisions) {

    // Get container of used vertices
    // Step 1
    // container = this->set_used_vertices_in_containers(nb_subdivisions, 2 * tol);

    // Step 2
    // list_of_collisions = get_container_collisions(container, []);
    std::vector<Box> list_of_boxes = this->convert_used_vertices_to_box_list(tol);
    
    std::vector<std::vector<std::vector<int>>> list_of_collisions;
    for (int i=0; i < list_of_boxes.size(); i++) {
        for (int j=i+1; j < list_of_boxes.size(); j++) {
            if (i != j) {
                list_of_collisions.push_back({list_of_boxes[i].get_id(), list_of_boxes[j].get_id()});
            };
        };
    };
    
    double tol_square = tol * tol;
    for (int i=0; i < list_of_collisions.size(); i++) {

        std::vector<int> ids_1 = list_of_collisions[i][1];
        std::vector<int> ids_2 = list_of_collisions[i][0];

        if (ids_1[0] != ids_2[0]) {
        // if (true) {
            
            std::vector<double> vertex_1 = this->vertices[ids_1[0]].get_point();
            std::vector<double> vertex_2 = this->vertices[ids_2[0]].get_point();

            double distance_square = (vertex_1[0] - vertex_2[0]) * (vertex_1[0] - vertex_2[0]) + (vertex_1[1] - vertex_2[1]) * (vertex_1[1] - vertex_2[1]) + (vertex_1[2] - vertex_2[2]) * (vertex_1[2] - vertex_2[2]);
            // If vertices are the same, set vertex_address_2 equal to vertex_address_1
            if (distance_square <= tol_square) {

                std::cout << ids_1[0] << " " << ids_2[0] << " : " << distance_square << std::endl;

                // ids_i contains [vertex_id, element_id, element_vertex, type_id(0 or 1)]
                if (ids_1[3] == 0) {
                    // Then the element is a doublet

                    //std::cout << doublets[ids_1[1]].get_vertices()[ids_1[2]] << " ";

                    if (ids_2[0] <= ids_1[0]) {
                        std::vector<int> new_vertices = this->doublets[ids_1[1]].get_vertices();
                        new_vertices[ids_1[2]] = ids_2[0];
                        this->doublets[ids_1[1]].set_vertex(new_vertices);
                    } else {
                        std::vector<int> new_vertices = this->doublets[ids_2[1]].get_vertices();
                        new_vertices[ids_2[2]] = ids_1[0];
                        this->doublets[ids_2[1]].set_vertex(new_vertices);
                    };

                    std::cout << doublets[ids_1[1]].get_vertices()[ids_1[2]] << std::endl;

                } else {
                    // Otherwise it is a vortex
                    if (ids_2[0] <= ids_1[0]) {
                        std::vector<int> new_vertices = this->vortexes[ids_1[1]].get_vertices();
                        new_vertices[ids_1[2]] = ids_2[0];
                        this->vortexes[ids_1[1]].set_vertex(new_vertices);
                    } else {
                        std::vector<int> new_vertices = this->vortexes[ids_2[1]].get_vertices();
                        new_vertices[ids_2[2]] = ids_1[0];
                        this->vortexes[ids_2[1]].set_vertex(new_vertices);
                    };
                };
            };
        };
    };
};
*/


void Mesh::stitch_used_vertices(double tol, int nb_subdivisions) {
    // Get list of vertices collisions and transmutation vector associated
    std::vector<std::vector<int>> list_of_collisions = this->get_identical_vertices(tol, nb_subdivisions);
    std::vector<int> transmute_vertices = this->get_identical_vertices_transmute(tol, nb_subdivisions);

    // Deal with vortexes
    for (int i=0; i < this->vortexes.size(); i++) {
        // Initialize updated_vertices
        std::vector<int> updated_vertices = this->vortexes[i].get_vertices();
        for (int j=0; j < this->vortexes[i].get_vertices().size(); j++) {
            updated_vertices[j] = transmute_vertices[updated_vertices[j]];
        };
        // Update vertices in current vortex
        this->vortexes[i].set_vertex(updated_vertices);
    };

    // Deal with doublets
    for (int i=0; i < this->doublets.size(); i++) {
        // Initialize updated_vertices
        std::vector<int> updated_vertices = this->doublets[i].get_vertices();
        for (int j=0; j < this->doublets[i].get_vertices().size(); j++) {
            updated_vertices[j] = transmute_vertices[updated_vertices[j]];
        };
        // Update vertices in current vortex
        this->doublets[i].set_vertex(updated_vertices);
    };
};


std::vector<int> Mesh::get_identical_vertices_transmute(double tol, int nb_subdivisions) {
    std::vector<std::vector<int>> list_of_collisions = this->get_identical_vertices(tol, nb_subdivisions);
    // Get initial transmutation
    std::vector<int> transmute_vertex;
    for (int i=0; i < this->vertices.size(); i++) {transmute_vertex.push_back(i);};
    for (int i=0; i < list_of_collisions.size(); i++) {
        int id_1 = list_of_collisions[i][0];
        int id_2 = list_of_collisions[i][1];
        if (transmute_vertex[id_2] <= transmute_vertex[id_1]) {
            transmute_vertex[id_1] = transmute_vertex[id_2];
        } else {
            transmute_vertex[id_2] = transmute_vertex[id_1];
        };
    };
    return transmute_vertex;
};


std::vector<std::vector<int>> Mesh::get_identical_vertices(double tol, int nb_subdivisions) {
    // Turn all vertices into list of boxes
    std::vector<Box> list_of_boxes;
    for (int i=0; i < this->vertices.size(); i++) {
        std::vector<double> point = this->vertices[i].get_point();
        std::vector<std::vector<double>> list_of_points;
        if (tol > 0) {
            std::vector<double> point_min = {point[0] - 2 * tol, point[1] - 2 * tol, point[2] - 2 * tol};
            std::vector<double> point_max = {point[0] + 2 * tol, point[1] + 2 * tol, point[2] + 2 * tol};
            list_of_points = {point_min, point_max};
        } else {
            list_of_points = {point};
        };
        list_of_boxes.push_back(Box({i}, list_of_points));
    };

    // // Set all vertices in container
    // container = generate_containers(None, list_of_boxes, nb_subdivisions)

    // // Get potential collisions
    // potential_collisions = get_container_collisions(container, [])
    std::vector<std::vector<std::vector<int>>> potential_collisions;
    for (int i=0; i < list_of_boxes.size(); i++) {
        for (int j=i+1; j < list_of_boxes.size(); j++) {
            if (i != j) {
                potential_collisions.push_back({list_of_boxes[i].get_id(), list_of_boxes[j].get_id()});
            };
        };
    };

    double tol_square = tol * tol;
    std::vector<std::vector<int>> list_of_collisions;
    for (int i=0; i < potential_collisions.size(); i++) {
        int id_1 = potential_collisions[i][0][0];
        int id_2 = potential_collisions[i][1][0];
        std::vector<double> vertex_1 = this->vertices[id_1].get_point();
        std::vector<double> vertex_2 = this->vertices[id_2].get_point();
        double distance_square = (vertex_1[0] - vertex_2[0]) * (vertex_1[0] - vertex_2[0]) + (vertex_1[1] - vertex_2[1]) * (vertex_1[1] - vertex_2[1]) + (vertex_1[2] - vertex_2[2]) * (vertex_1[2] - vertex_2[2]);
        // If vertices are the same, set vertex_address_2 equal to vertex_address_1
        if (distance_square <= tol_square) {list_of_collisions.push_back({id_1, id_2});};
    };
    return list_of_collisions;
};

void Mesh::cleanup_doublets() {
    int index;
    for (int i=0; i < this->doublets.size(); i++) {
        index = this->doublets.size() - 1 - i;
        std::vector<int> v = this->doublets[i].get_vertices();
        remove(v);
        if (v.size() > 3) {
            this->doublets[i].set_vertex(v);
        } else if (v.size() == 3) {
            //std::cout << "3 vertices : " << index << std::endl;
            this->doublets[i].set_vertex(v);
        } else {
            //std::cout << "2 or less vertices : " << index << std::endl;
            this->delete_doublet(index);
        };
    };
};

void Mesh::cleanup_vortexes() {
    int index;
    for (int i=0; i < this->vortexes.size(); i++) {
        index = this->vortexes.size() - 1 - i;
        std::vector<int> v = this->vortexes[i].get_vertices();
        remove(v);
        if (v.size() >= 4) {
            this->vortexes[i].set_vertex(v);
        } else {
            this->delete_vortex(index);
        };
    };
};



/*
    # def delete_unused_vertices(self):
    #     vertices_usage = np.zeros([len(self.vertices), 1])
    #     for element in self.elements:
    #         for vertex in element.vertices:
    #             vertices_usage[int(vertex)] += 1

    # def update_surfaces(self):
    #     # Re-initialize element count in mesh surfaces
    #     for surface in self.surfaces:
    #         surface.element_count = 0
    #
    #     for element in self.elements:
    #         is_new_surface = True
    #         for surface in self.surfaces:
    #             if element.surface == surface.id and element.type == surface.type:
    #                 surface.element_count += 1
    #                 # Then surface is already initialized, break out of for loop
    #                 is_new_surface = False
    #                 break
    #
    #         if is_new_surface:
    #             # Create new surface in mesh and initialize element count to 1
    #             self.surfaces.append(Mesh.Surface(element.surface, element.type, 1))
    #             self.surface_count += 1


def convert_doublet_to_mesh_object(mesh, vertex_table, element_table, surface):
    if mesh is None:
        mesh = Mesh()
    for i in range(0, np.shape(vertex_table)[0]):
        # mesh.append_vertex([vertex_table[i, 0], vertex_table[i, 1], vertex_table[i, 2]])
        mesh.append_vertex(vertex_table[i, :])
    for i in range(0, np.shape(element_table)[0]):
        mesh.create_doublet(element_table[i, :], surface)
    return mesh


def convert_vortex_to_mesh_object(mesh, vertex_table, element_table, surface, stations, rows):
    if mesh is None:
        mesh = Mesh()
    for i in range(0, np.shape(vertex_table)[0]):
        # mesh.append_vertex([vertex_table[i, 0], vertex_table[i, 1], vertex_table[i, 2]])
        mesh.append_vertex(vertex_table[i, :])
    for i in range(0, np.shape(element_table)[0]):
        mesh.create_vortex(element_table[i, :], surface, stations[i], rows[i])
    return mesh
*/





double get_triangle_area_square(Mesh &mesh, std::vector<int> vertex_ids) {
    std::vector<double> point_a = mesh.get_vertices(vertex_ids)[0].get_point();
    std::vector<double> point_b = mesh.get_vertices(vertex_ids)[1].get_point();
    std::vector<double> point_c = mesh.get_vertices(vertex_ids)[2].get_point();
    std::vector<double> n = crossProduct(substract(point_b, point_a), substract(point_c, point_a));
    // double area_square = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) / 4;
    return (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) / 4;
};

void mesh_rectangle_to_triangles(Mesh &mesh, std::vector<std::vector<int>> vertices, int surface_id) {
    if ((vertices.size() == 2) && (vertices[0].size() == 2) && (vertices[1].size() == 2)) {
        // Two ways to cut the shape into triangles
        // Arrangement 1
        std::vector<int> triangle_1a = {vertices[0][0], vertices[0][1], vertices[1][1]};
        std::vector<int> triangle_1b = {vertices[0][0], vertices[1][1], vertices[1][0]};
        // Measure the "greatness" of these two triangles
        double area_min_1 = get_triangle_area_square(mesh, triangle_1a);
        area_min_1 = std::min(area_min_1, get_triangle_area_square(mesh, triangle_1b));

        // Arrangement 2
        std::vector<int> triangle_2a = {vertices[1][0], vertices[0][0], vertices[0][1]};
        std::vector<int> triangle_2b = {vertices[1][0], vertices[0][1], vertices[1][1]};
        // Measure the "greatness" of these two triangles
        double area_min_2 = get_triangle_area_square(mesh, triangle_2a);
        area_min_2 = std::min(area_min_2, get_triangle_area_square(mesh, triangle_2b));

        // Choose the best arrangement, currently based on triangles area
        if (area_min_1 < area_min_2) {
            // std::vector<int>::iterator ip;
            // ip = std::unique(triangle_2a.begin(), triangle_2a.end());
            // std::distance(triangle_2a.begin(), ip);


            // Arrangement 2 is superior
            if (std::distance(triangle_2a.begin(), std::unique(triangle_2a.begin(), triangle_2a.end())) == 3) {
                mesh.create_doublet(triangle_2a, surface_id);
            };
            if (std::distance(triangle_2b.begin(), std::unique(triangle_2b.begin(), triangle_2b.end())) == 3) {
                mesh.create_doublet(triangle_2b, surface_id);
            };
        } else {
            // Arrangement 1 is superior
            if (std::distance(triangle_1a.begin(), std::unique(triangle_1a.begin(), triangle_1a.end())) == 3) {
                mesh.create_doublet(triangle_1a, surface_id);
            };
            if (std::distance(triangle_1b.begin(), std::unique(triangle_1b.begin(), triangle_1b.end())) == 3) {
                mesh.create_doublet(triangle_1b, surface_id);
            };
        };
    };
};

void mesh_grid(Mesh &mesh, std::vector<std::vector<int>> vertex_addressing, std::vector<int> grid_offset, int surface_id, bool is_vortex, bool clockwise, bool verify_vertices, std::string mode) {
    int station_count = 0;
    for (int i=0; i < vertex_addressing.size() - 1; i++) {
        int row_count = 0;
        bool wingstation_is_initialized = false;

        for (int j=0; j < vertex_addressing[i].size() - 1; j++) {
            bool do_mesh = true;
            std::vector<std::vector<int>> vertices;
            // Transpose vertices to reverse the order of vertices if not clockwise
            if (!clockwise) {
                vertices = {{vertex_addressing[i][j], vertex_addressing[i][j + 1]}, {vertex_addressing[i + 1][j], vertex_addressing[i + 1][j + 1]}};
            } else {
                vertices = {{vertex_addressing[i][j], vertex_addressing[i + 1][j]}, {vertex_addressing[i][j + 1], vertex_addressing[i + 1][j + 1]}};
            };
            // If specified, make sure the vertices are valid (more than 2 unique vertices)
            if (verify_vertices) {
                std::vector<int> verify_vertices = {vertex_addressing[i][j], vertex_addressing[i][j + 1], vertex_addressing[i + 1][j], vertex_addressing[i + 1][j + 1]};
                if (std::distance(verify_vertices.begin(), std::unique(verify_vertices.begin(), verify_vertices.end())) < 3) {do_mesh = false;};
            };

            if (do_mesh) {
                if (mode == "Structured") {
                    std::vector<int> rectangle = {vertices[0][0], vertices[0][1], vertices[1][1], vertices[1][0]};
                    // Make sure vortexes are not added if they do not have more than 3 distinct nodes

                    if (is_vortex) {
                        if (std::distance(rectangle.begin(), std::unique(rectangle.begin(), rectangle.end())) > 3) {
                            // Add to row count
                            row_count += 1;

                            // On the last pass, initiate new wingstation
                            if (!wingstation_is_initialized) {
                                mesh.create_new_wingstation(surface_id, {});
                                wingstation_is_initialized = true;
                                // Add to station count
                                station_count += 1;
                            };

                            int station = station_count - 1 + grid_offset[0];
                            int row = row_count - 1 + grid_offset[1];
                            mesh.create_vortex(rectangle, surface_id, station, row);

                        };
                    } else {
                        mesh.create_doublet(rectangle, surface_id);
                    };
                } else {
                    mesh_rectangle_to_triangles(mesh, vertices, surface_id);
                };
            };
        };
    };
};

std::vector<std::vector<int>> get_vertices_address(Mesh &mesh, std::vector<std::vector<double>> x, std::vector<std::vector<double>> y, std::vector<std::vector<double>> z) {
    int vertices_nb = mesh.get_vertices().size();
    std::vector<std::vector<int>> vertex_addressing;
    for (int i=0; i < x.size(); i++) {
        std::vector<int> vertex_addressing_i;
        for (int j=0; j < x[i].size(); j++) {
            int vertex_index = i * x[i].size() + j + vertices_nb;
            vertex_addressing_i.push_back(vertex_index);
            Vertex new_vertex = Vertex({x[i][j], y[i][j], z[i][j]});
            mesh.append_vertex(new_vertex);
        };
        vertex_addressing.push_back(vertex_addressing_i);
    };
    return vertex_addressing;
};

std::vector<std::vector<int>> stitch_addressing(Mesh &mesh, std::vector<std::vector<int>> vertex_addressing_1, std::vector<std::vector<int>> vertex_addressing_2, double tol) {
    double tol_square = tol * tol;
    // Run through all addresses in first grid and compare to all elements in second grid
    for (int i1=0; i1 < vertex_addressing_1.size(); i1++) {
        for (int j1=0; j1 < vertex_addressing_1[0].size(); j1++) {
            int address_1 = vertex_addressing_1[i1][j1];
            std::vector<double> vertex_1 = mesh.get_vertices()[address_1].get_point();
            for (int i2=0; i2 < vertex_addressing_2.size(); i2++) {
                for (int j2=0; j2 < vertex_addressing_2[0].size(); j2++) {
                    int address_2 = vertex_addressing_2[i2][j2];
                    std::vector<double> vertex_2 = mesh.get_vertices()[address_2].get_point();
                    // Get distance squared
                    double distance_square = (vertex_1[0] - vertex_2[0]) * (vertex_1[0] - vertex_2[0]) + (vertex_1[1] - vertex_2[1]) * (vertex_1[1] - vertex_2[1]) + (vertex_1[2] - vertex_2[2]) * (vertex_1[2] - vertex_2[2]);
                    // If vertices are the same, set vertex_address_2 equal to vertex_address_1
                    if (distance_square <= tol_square) {
                        vertex_addressing_1[i1][j1] = address_2;
                    };
                };
            };
        };
    };
    return vertex_addressing_1;
};

/*
def convert_addressing_to_box_list(vertex_addressing, addressing_id, tol):
    shape = np.shape(vertex_addressing)
    list_of_boxes = []
    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            vertex = int(vertex_addressing[i, j])
            point = mesh.vertices[vertex]
            if tol > 0:
                point_min = [point[0] - tol, point[1] - tol, point[2] - tol]
                point_max = [point[0] + tol, point[1] + tol, point[2] + tol]
                list_of_points = [point_min, point_max]
            else:
                list_of_points = [point]
            list_of_boxes.append(Box([addressing_id, [i, j]], list_of_points))
    return list_of_boxes

def set_addressing_in_containers(vertex_addressing, addressing_id, max_subdivisions, tol):
    # Extract list of PointObjects from mesh -> elements
    list_of_boxes = convert_addressing_to_box_list(vertex_addressing, addressing_id, tol=tol)
    # Create container structure
    container = generate_containers(None, list_of_boxes, max_subdivisions)
    return container
*/


/*
def stitch_addressing_v1(vertex_addressing_1, vertex_addressing_2, tol, nb_subdivisions):
    # Get container of vertex_addressing_1
    print(" Step 1")
    container_1 = set_addressing_in_containers(vertex_addressing_1, 1, nb_subdivisions, 2 * tol)
    container_2 = set_addressing_in_containers(vertex_addressing_2, 2, nb_subdivisions, 2 * tol)

    print(" Step 2")
    list_of_collisions = get_collisions(container_1, container_2, [])

    tol_square = tol ** 2
    for collision in list_of_collisions:
        ids_1 = collision[1][1]
        ids_2 = collision[0][1]
        address_1 = vertex_addressing_1[ids_1[0], ids_1[1]]
        address_2 = vertex_addressing_2[ids_2[0], ids_2[1]]
        vertex_1 = mesh.vertices[int(address_1)]
        vertex_2 = mesh.vertices[int(address_2)]

        distance_square = (vertex_1[0] - vertex_2[0]) ** 2 + \
                          (vertex_1[1] - vertex_2[1]) ** 2 + (vertex_1[2] - vertex_2[2]) ** 2
        # If vertices are the same, set vertex_address_2 equal to vertex_address_1
        if distance_square <= tol_square:
            vertex_addressing_1[ids_1[0], ids_1[1]] = address_2
    return vertex_addressing_1
*/



void get_surface_mesh_tables(Body body, Mesh &mesh, std::string mode, bool is_vortex) {
    double distance_tol = 1e-14;

    std::string surface_type;
    if (is_vortex) {
        surface_type = "Vortex";
    } else {
        surface_type = "Doublet";
    };
    mesh.create_new_surface(surface_type);
    int surface_id = mesh.get_surfaces()[mesh.get_surface_count() - 1].get_id();
    
    // Get body surfaces
    std::vector<std::vector<std::vector<std::vector<double>>>> surfaces = body.get_body_surfaces();
    std::vector<std::vector<int>> associativity = body.get_associativity();
    std::vector<std::vector<std::vector<int>>> vertex_addressing_list;
    int nb = 20;

    if (is_vortex) {

        for (int i=0; i < associativity.size(); i++) {
            std::vector<int> stitch = associativity[i];
            int id_1 = stitch[0];
            int id_2 = stitch[1];

            std::vector<std::vector<double>> surface_x = scale_matrices(0.5, sum_matrices(surfaces[id_1][0], surfaces[id_2][0]));
            std::vector<std::vector<double>> surface_y = scale_matrices(0.5, sum_matrices(surfaces[id_1][1], surfaces[id_2][1]));
            std::vector<std::vector<double>> surface_z = scale_matrices(0.5, sum_matrices(surfaces[id_1][2], surfaces[id_2][2]));

            std::vector<std::vector<int>> surface_addressing = get_vertices_address(mesh, surface_x, surface_y, surface_z);
            vertex_addressing_list.push_back(surface_addressing);
        };

        // Stitch every surface together to correct for repeated vertices
        for (int ik=0; ik < vertex_addressing_list.size(); ik++) {
            for (int jk=ik+1; jk < vertex_addressing_list.size(); jk++) {
                vertex_addressing_list[ik] = stitch_addressing(mesh, vertex_addressing_list[ik],
                                                               vertex_addressing_list[jk], distance_tol);
                // vertex_addressing_list[ik] = stitch_addressing_v1(vertex_addressing_list[ik],
                //                                                   vertex_addressing_list[jk], distance_tol, nb);
            };
        };

        // for vertex_addressing_array in vertex_addressing_list:
        std::vector<int> current_offset = {mesh.get_wingstation_count(), 0};
        for (int k=0; k < vertex_addressing_list.size(); k++) {
            bool clockwise = true;
            bool verify_vertices = true;
            mesh_grid(mesh, vertex_addressing_list[k], current_offset, surface_id, is_vortex, clockwise, verify_vertices, mode);
            // Update current grid offset. For now, it is as simple as the sum of previous grids shape[1] values in list
            current_offset[0] += vertex_addressing_list[k].size() - 1;
        };

        // Flip wingstations order if body built is marked as flipped
        if (body.get_is_flipped()) {
            mesh.reverse_wingstations(surface_id);
        };

    } else {
        for (int i=0; i < surfaces.size(); i++) {
            // Get basic addressing of surface
            std::vector<std::vector<int>> surface_addressing = get_vertices_address(mesh, surfaces[i][0], surfaces[i][1], surfaces[i][2]);
            vertex_addressing_list.push_back(surface_addressing);
        };

        // Stitch every surface together to correct for repeated vertices
        for (int ik=0; ik < vertex_addressing_list.size(); ik++) {
            for (int jk=ik+1; jk < vertex_addressing_list.size(); jk++) {
                vertex_addressing_list[ik] = stitch_addressing(mesh, vertex_addressing_list[ik],
                                                               vertex_addressing_list[jk], distance_tol);
                // vertex_addressing_list[ik] = stitch_addressing_v1(vertex_addressing_list[ik],
                //                                                   vertex_addressing_list[jk], distance_tol, nb);
            };
        };
        
        /*
        for (int i=0; i < associativity.size(); i++) {
            std::vector<int> stitch = associativity[i];
            std::vector<std::vector<int>> su_addressing = vertex_addressing_list[stitch[0]];
            std::vector<std::vector<int>> sl_addressing = vertex_addressing_list[stitch[1]];

            su_shape = np.shape(su_addressing)
            sl_shape = np.shape(sl_addressing)
            if su_shape[0] == sl_shape[0] and su_shape[1] == sl_shape[1]:
                stitch_vertices = np.vstack([su_addressing[0, :], sl_addressing[0, :]]);
                mesh_grid(stitch_vertices, [], clockwise=True, verify_vertices=True);

                stitch_vertices = np.vstack([su_addressing[su_shape[0] - 1, :], sl_addressing[sl_shape[0] - 1, :]]);
                mesh_grid(stitch_vertices, [], clockwise=True, verify_vertices=True);

                stitch_vertices = np.vstack([su_addressing[:, 0], sl_addressing[:, 0]]);
                mesh_grid(stitch_vertices, [], clockwise=True, verify_vertices=True);

                stitch_vertices = np.vstack([su_addressing[:, su_shape[1] - 1], sl_addressing[:, sl_shape[1] - 1]]);
                mesh_grid(stitch_vertices, [], clockwise=True, verify_vertices=True);
        };
        */

        // for vertex_addressing_array in vertex_addressing_list:
        for (int k=0; k < vertex_addressing_list.size(); k++) {
            bool order_vertices_clockwise = true;
            bool verify_vertices = true;
            if (body.get_surfaces()[k].get_id() == "Lower") {order_vertices_clockwise = false;};
            mesh_grid(mesh, vertex_addressing_list[k], {0, 0}, surface_id, is_vortex, order_vertices_clockwise, verify_vertices, mode);
        };
    };
    
};







std::string write_vertices(std::string file_string, std::vector<Vertex> vertices) {
    int vertex_dim = vertices.size();
    for (int i=0; i < vertex_dim; i++) {
        file_string += "NODE " + std::to_string(i);
        for (int j=0; j < vertices[i].get_point().size(); j++) {

            file_string += " " + std::to_string(vertices[i].get_point()[j]);
        };
        file_string += "\n";
    };
    return file_string;
};

//std::string write_vortexes(std::string file_string, std::unordered_map<int, Vortex> elements) {
std::string write_vortexes(std::string file_string, std::vector<Vortex> elements) {
    // Loop through all elements to print vortexes
    for (int i=0; i < elements.size(); i++) {
        file_string += "VORTEX " + std::to_string(elements[i].get_id());
        for (int j=0; j < elements[i].get_vertices().size(); j++) {
            file_string += " " + std::to_string(elements[i].get_vertices()[j]);
        };
        file_string += "\n";
    };
    return file_string;
};

//std::string write_doublets(std::string file_string, std::unordered_map<int, Doublet> elements) {
std::string write_doublets(std::string file_string, std::vector<Doublet> elements) {
    // Loop through all elements to print doublets
    for (int i=0; i < elements.size(); i++) {
        file_string += "DOUBLET " + std::to_string(elements[i].get_id());
        for (int j=0; j < elements[i].get_vertices().size(); j++) {
            file_string += " " + std::to_string(elements[i].get_vertices()[j]);
        };
        file_string += "\n";
    };
    return file_string;
};

std::string write_wingstations(std::string file_string, std::vector<Wingstation> wingstations) {
    // Loop through all wingstations
    for (int i=0; i < wingstations.size(); i++) {
        file_string += "WINGSTATION " + std::to_string(wingstations[i].get_id());
        // Run through all elements and add if it is part of surface
        for (int j=0; j < wingstations[i].get_elements().size(); j++) {
            file_string += " " + std::to_string(wingstations[i].get_elements()[j]);
        };
        file_string += "\n";
    };
    return file_string;
};

std::string write_wings(std::string file_string, std::vector<WingSurface> surfaces) {
    // Loop through all surfaces and print those of type Vortex
    for (int i=0; i < surfaces.size(); i++) {
        if (surfaces[i].get_type() == "Vortex") {
            file_string += "WING " + std::to_string(surfaces[i].get_id());
            // Run through all elements and add if it is part of surface
            for (int j=0; j < surfaces[i].get_wingstations().size(); j++) {
                file_string += " " + std::to_string(surfaces[i].get_wingstations()[j]);
            };
            file_string += "\n";
        };
    };
    return file_string;
};

//std::string write_patches(std::string file_string, std::vector<WingSurface> surfaces, std::unordered_map<int, Doublet> elements) {
std::string write_patches(std::string file_string, std::vector<WingSurface> surfaces, std::vector<Doublet> elements) {
    // Loop through all surfaces to print surfaces of type Doublet
    for (int i=0; i < surfaces.size(); i++) {
        // Check if surface is of type doublet
        if (surfaces[i].get_type() == "Doublet") {
            file_string += "PATCH " + std::to_string(surfaces[i].get_id());
            // Run through all elements and add if it is part of surface
            for (int j=0; j < elements.size(); j++) {
                if (elements[j].get_surface() == surfaces[i].get_id()) {
                    file_string += " " + std::to_string(elements[j].get_id());
                };
            };
            file_string += "\n";
        };
    };
    return file_string;
};


void export_mesh(std::string filename, Mesh mesh) {

    std::string file_string;
    file_string = write_vertices(file_string, mesh.get_vertices());
    file_string = write_doublets(file_string, mesh.get_doublets());
    file_string = write_vortexes(file_string, mesh.get_vortexes());
    file_string = write_wingstations(file_string, mesh.get_wingstations());
    file_string = write_wings(file_string, mesh.get_surfaces());
    file_string = write_patches(file_string, mesh.get_surfaces(), mesh.get_doublets());
    // Write string to file
    std::ofstream(filename) << file_string;
};


