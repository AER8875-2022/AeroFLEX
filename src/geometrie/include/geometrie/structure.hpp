#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#ifndef structure_hpp
#define structure_hpp
#include <stdio.h>

using namespace std;

vector<vector<double>> get_center_point(vector<vector<double>> X, vector<vector<double>> Y, vector<vector<double>> Su, vector<vector<double>> Sl);
vector<double> vecteur_normal(vector<double> p1, vector<double> p2);
vector<vector<vector<vector<double>>>> get_geometry(class Body);
vector<tuple<int,vector<double>,vector<double>,vector<double>>> maillage_structure(class Body);
#endif /* structure_hpp */
