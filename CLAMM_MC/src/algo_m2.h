#pragma once
#ifndef algo_m2_h
#define algo_m2_h

#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <map>
#include <cmath>
#include <random>
#include <errno.h>

#include "rule.h"
#include "sim_cell.h"
#include "running_stat.h"
#include "session.h"
#include "utils.h"
using namespace std;
#define DONE -1

class Algo_m2 {
public:
	const double Kb = 0.00008617333262; // Boltzmann constant
	const double uB = 0.00005788381806; // Bohr magnaton
	int outfile_count = 0;
	SimCell sim_cell;
	Session session;
	vector<int> chem_list;
	vector<float> spin_list;
	vector<vector<float>> pos_list;
	vector<vector<vector<vector<int>>>> spin_motif_groups;
	vector<vector<vector<vector<int>>>> sro_motif_groups;
	vector<vector<vector<vector<int>>>> chem_motif_groups;
	map <size_t, vector<double>> rule_map_chem;
	map <size_t, vector<double>> rule_map_spin;
	map <size_t, double> rule_map_sro;

	Algo_m2(void);
	Algo_m2(Session& _session, SimCell& _sim_cell);
	void run();
	void fill_CMG(vector<vector<int>>& neigh_ind_list);
	void fill_SMG(vector<vector<int>>& neigh_ind_list);
	void fill_SROMG(vector<vector<int>>& neigh_ind_list);
	void print_state(string contcar_name, int temp);
	bool bc_check(vector<float> check_vect, vector<float>& pos);
	vector<double> eval_lat();
	double eval_lat_chem();
	double eval_lat_spin();
	double eval_site_spin(int site);
	//double eval_site_spin_for_dif(int site);
	vector<double> eval_site_spin_for_dif(int site);
	double eval_site_chem(int site);
	double eval_spin_flip(int site, float old_spin);
	float spec_spin(int site, int spec);
	float spin_dot(int site);
	size_t cust_hash(vector<uint32_t>& vect);
};

#endif
