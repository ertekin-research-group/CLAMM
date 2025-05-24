#pragma once
#ifndef algo1_h
#define algo1_h

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

class Algo1 {
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

	Algo1(void);
	Algo1(Session& _session, SimCell& _sim_cell);
	void run();
	void fill_CMG(vector<vector<int>>& neigh_ind_list);
	void fill_SMG(vector<vector<int>>& neigh_ind_list);
	void fill_SROMG(vector<vector<int>>& neigh_ind_list);
	void print_state(string contcar_name, int temp);
	bool bc_check(vector<float> check_vect, vector<float>& pos);
	vector<double> eval_lat();
	double eval_lat_spin();
	double eval_site_spin(int site);
	double eval_site_chem(int site);
	double eval_spin_flip(int site, float old_spin);
	float spec_spin(int site, int spec);
	float spin_dot(int site);
	size_t cust_hash(vector<uint32_t>& vect);
};

#endif

//#pragma once
//#ifndef algo1_h
//#define algo1_h
//
//#include <chrono>
//#include <vector>
//#include <string>
//#include <fstream>
//#include <iterator>
//#include <iostream>
//#include <sstream>
//#include <algorithm>
//#include <iterator>
//#include <map>
//#include <cmath>
//#include <random>
//#include <fcntl.h>
//#include <errno.h>
//#include <unistd.h>
//#include "rule.h"
//#include "sim_cell.h"
//#include "running_stat.h"
//#include "session.h"
//#include "utils.h"
//using namespace std;
//#define ONE_SPIN 0
//#define NO_SPIN 0
//#define ONE_ATOM 0
//#define METHOD_0 0
//#define METHOD_1 1
//#define DONE -1
//
//class Algo1 {
//public:
//    const double Kb = 0.00008617333262; // Boltzmann constant
//    const double uB = 0.00005788381806; // Bhor magnaton
//    int flip_count = 0;
//    int flip_count2 = 0;
//    int outfile_count = 0;
//    int passes;
//    int numb_vac = 0;
//    float spin_flip = 0.0;
//    double init_spin = 0.0;
//    double e_flip = 0.0;
//    double init_enrg = 0.0;
//    RunningStat rs_E;
//    RunningStat rs_M;
//    SimCell sim_cell;
//    Session session;
//    vector<int> chem_list;
//    vector<float> spin_list;
//    vector<float> site_rule_count_list;
//    vector<double> lat_rule_count_list;
//    vector<double> count_avg;
//    vector<vector<float>> pos_list;
//    vector<vector<vector<vector<int>>>> spin_motif_groups;
//    vector<vector<vector<vector<int>>>> chem_motif_groups;
//    map <size_t, vector<double>> rule_map_chem;
//    map <size_t, double> rule_map_spin;
//    // setup rng for random spin choice and acceptance probability
//    std::mt19937_64 rng;
//    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//    std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
//    std::uniform_real_distribution<double> unif;
//    std::uniform_int_distribution<int> rand_atom;
//    std::uniform_int_distribution<int> rand_method;
//
//    Algo1(void);
//    Algo1(Session& _session, SimCell& _sim_cell);
//    size_t cust_hash(vector<uint32_t>& vect);
//    void run();
//    void fill_CMG(vector<vector<int>>& neigh_ind_list);
//    void fill_SMG(vector<vector<int>>& neigh_ind_list);
//    void print_state(string contcar_name, int temp);
//    void spin_move(int site, int pass, float temp, float new_spin);
//    bool bc_check(vector<float> check_vect, vector<float>& pos);
//    double eval_lat();
//    double eval_lat_spin();
//    double eval_site_spin(int site);
//    double eval_site_chem(int site);
//    double eval_spin_flip(int site, float old_spin);
//};
//
//#endif
