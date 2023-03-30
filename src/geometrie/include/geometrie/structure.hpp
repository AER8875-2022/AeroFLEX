#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <array>

#ifndef structure_hpp
#define structure_hpp
#include <stdio.h>


std::vector<std::vector<double>> get_center_point(std::vector<std::vector<double>> X, std::vector<std::vector<double>> Y, std::vector<std::vector<double>> Su, std::vector<std::vector<double>> Sl);
void vecteur_normal(std::array<double, 3> &normal, std::vector<double> p1, std::vector<double> p2);
std::vector<std::vector<std::vector<std::vector<double>>>> get_geometry(class Body);
std::vector<std::tuple<int,std::vector<double>,std::vector<double>,std::vector<double>>> maillage_structure(class Body);
#endif /* structure_hpp */
