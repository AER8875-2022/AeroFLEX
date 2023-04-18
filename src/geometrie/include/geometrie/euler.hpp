#include <cstddef>
#include <fstream>
#include <set>
# include <iostream>
# include <vector>
# include <cmath>
# include <array>
# include <gmsh.h> //ajout SEB

#ifndef euler_hpp
#define euler_hpp
//using namespace std;

#include <stdio.h>


void generer(std::vector<std::vector<double>> xu, std::vector<std::vector<double>> xl, std::vector<std::vector<double>> y, std::vector<std::vector<double>> su, std::vector<std::vector<double>> sl, double disc, std::string file_name, bool RANS);

#endif /* geometrie_hpp */
