#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <array>
#include "geometrie/geometry.hpp"


#ifndef mesh_vlm_hpp
#define mesh_vlm_hpp

class Box
{
private:
    std::vector<int> id;
    std::vector<std::vector<double>> points;
    std::vector<double> min;
    std::vector<double> max;
public:
    Box(std::vector<int> box_id, std::vector<std::vector<double>> list_of_points);
    void set_boundaries();
    std::vector<std::vector<double>> get_points();
    bool box_touches_box(Box box);
    bool is_box_inside_box(Box box);
    std::vector<int> get_id();
};


class Vertex
{
private:
    std::vector<double> point;
public:
    Vertex();
    Vertex(std::vector<double> point);
    std::vector<double> get_point();
};

class Vortex
{
private:
    std::string type;
    std::vector<int> vertices;
    int id;
    int surface;
    int station;
    int row;
public:
    Vortex();
    Vortex(std::vector<int> vertices, int surface, int station, int row);
    std::vector<int> get_vertices();
    int get_id();
    int get_surface();
    int get_station();
    int get_row();
    void set_id(int id);
    void set_vertex(std::vector<int> new_ids);
};

class Doublet
{
private:
    std::string type;
    std::vector<int> vertices;
    int id;
    int surface;
public:
    Doublet();
    Doublet(std::vector<int> vertices, int surface);
    std::vector<int> get_vertices();
    int get_id();
    int get_surface();
    void set_id(int id);
    void set_vertex(std::vector<int> new_ids);
};

class Wingstation
{
private:
    std::vector<int> elements;
    int id;
    int surface;
public:
    Wingstation();
    Wingstation(int wingstation_id, int surface, std::vector<int> elements);
    std::vector<int> get_elements();
    int get_id();
    int get_surface();
    void push_vortex(int vortex_id);
};

class WingSurface
{
private:
    int id;
    std::string type;
    std::vector<int> wingstations;
public:
    WingSurface();
    WingSurface(int surface_id, std::string surface_type);
    void push_wingstation(int wingstation_id);
    void reverse_wingstations();
    int get_id();
    std::string get_type();
    std::vector<int> get_wingstations();
};



class Mesh
{
private:
    std::vector<Vertex> vertices;
    std::vector<Vortex> vortexes;
    std::vector<Doublet> doublets;
    std::vector<Wingstation> wingstations;
    std::vector<WingSurface> surfaces;
    int vortex_count;
    int doublet_count;
    int wingstation_count;
    int surface_count;

public:
    Mesh();

    int append_vertex(Vertex vertex);
    void create_vortex(std::vector<int> vertices, int surface, int station, int row);
    void create_doublet(std::vector<int> vertices, int surface);
    void create_new_surface(std::string surface_type);
    void create_new_wingstation(int surface_id, std::vector<int> vortexes);
    void push_vortex_to_wingstation(int vortex_id, int wingstation_id);
    std::vector<Vertex> get_vertices(std::vector<int> vertex_ids);
    std::vector<Vertex> get_vertices();
    std::vector<WingSurface> get_surfaces();
    std::vector<Wingstation> get_wingstations();
    std::vector<Vortex> get_vortexes();
    std::vector<Doublet> get_doublets();
    int get_surface_count();
    int get_wingstation_count();
    void reverse_wingstations(int surface_id);
    void change_doublet_vertices(int id, std::vector<int> new_vertices);
    void change_vortex_vertices(int id, std::vector<int> new_vertices);
    std::vector<std::vector<double>> get_vortexes_normals();
    std::vector<std::vector<double>> get_doublets_normals();
    void delete_vortex(int vortex_id);
    void delete_doublet(int doublet_id);
    std::vector<int> get_unused_vertices();
    std::vector<int> get_transmute_vertices();
    void delete_unused_vertices();
    std::vector<Box> convert_used_vertices_to_box_list(double tol);
    void set_used_vertices_in_containers(int max_subdivisions, double tol);
    void stitch_used_vertices(double tol, int nb_subdivisions);
    std::vector<int> get_identical_vertices_transmute(double tol, int nb_subdivisions);
    std::vector<std::vector<int>> get_identical_vertices(double tol, int nb_subdivisions);
    void cleanup_doublets();
    void cleanup_vortexes();
};

void get_surface_mesh_tables(Body body, Mesh &mesh, std::string mode, bool is_vortex);




#endif /* mesh_vlm_hpp */