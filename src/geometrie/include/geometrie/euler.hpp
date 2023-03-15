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
using namespace std;

#include <stdio.h>


void generer(vector<vector<double>> xu, vector<vector<double>> xl, vector<vector<double>> y, vector<vector<double>> su, vector<vector<double>> sl, double disc, string file_name, bool RANS);

#endif /* geometrie_hpp */
