#pragma once
#ifndef utils_h
#define utils_h
#include <stdio.h>
#include <string>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <sstream>
using namespace std;

int openWithFcntl(const char* filename);
string replace_char(string s, char c1, char c2);
void sort_vect(vector<int>& vect, vector<int>& perm);
void sort_vect(vector<float>& vect, vector<int>& perm);
void sort_vect(vector<vector<int>>& vect, vector<int>& perm);
void sort_vect(vector<vector<float>>& vect, vector<int>& perm);
int get_index(vector<int>& vect, int elem);
int get_index(vector<float>& vect, float elem);
int get_index(vector<string>& vect, string elem);
int vect_max(vector<int>& vect);
int kron_del(int int1, int int2);
int kron_del(float float1, float float2);
int kron_del(int int1, float float2);
int kron_del(float float1, int int2);
float round_to(float val, int digits);
int vect_max(vector<int>& vect);
float vect_max(vector<float>& vect);
float dot_pos(vector<float>& vect1, vector<float>& vect2);
float pos_dist(vector<float> pos1, vector<float> pos2);
vector<int> vect_permut(vector<int>& vect);
vector<int> vect_permut(vector<float>& vect);
vector<int> vect_add(vector<int>& vect1, vector<int>& vect2);
vector<int> scale_vect(vector<int>& vect, int numb);
vector<float> pos_round(vector<float>& pos, int digits);
vector<float> vect_add(vector<float>& vect1, vector<float>& vect2);
vector<float> vect_subtract(vector<float>& vect1, vector<float>& vect2);
vector<float> vect_add(vector<float>& vect1, float vect2[3]);
vector<float> pos_transform(vector<float>& pos, vector<vector<float>>& trans);
vector<float> scale_vect(vector<float>& vect, float numb);
vector<double> vect_add(vector<double>& vect1, vector<float>& vect2);
vector<double> vect_add(vector<double>& vect1, vector<double>& vect2);
vector<string> split(string str, const string delim);
vector<string> split(string str); 
int sgn(float v);
int sign(float v);
bool fcomp(float a, float b, float tol = 0.0001);
bool pos_comp(vector<float>& pos1, vector<float>& pos2, float tol = 0.0001);
vector<float> frac_to_cart(vector<float>& pos, vector<vector<float>>& lat_vec);
vector<float> cart_to_frac(vector<float>& pos, vector<vector<float>>& lat_vec);
#endif
