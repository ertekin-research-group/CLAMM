#pragma once
#ifndef session_h
#define session_h
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "utils.h"
#include "rule.h"
using namespace std;

class Session
{
public:
	bool use_poscar = true;
	bool use_states = true;
    bool do_conv_output = false;
	bool write_contcars = true;
	int algo = 0;
	int numb_passes = 1;
	int ta_passes = 1;
	int eq_passes = 1;
	int shape[3] = { 0,0,0 };
	int tot_atoms;
	float start_temp = 1;
	float end_temp = 0;
	float temp_step = 0.5;
	float sro_temp = 0;
	float sro_target;
	float mag_ext = 0;
	float intercept = 0;
	int multiplicity = 1;
	string structure_file; //unused
	string rules_file; //unused
	string sro_file; //unused
	string sim_type; //unused
	string phase_init;
	string spin_init;
	string species_init; //unfinished
	vector<int> atom_numbs;
	vector<int> species_inds;
	vector<int> sro_def;
	vector<float> moments;
	vector<float> unique_dists;
	vector<string> species_str;
	vector<Rule> spin_rule_list;
	vector<Rule> chem_rule_list;
	vector<Rule> sro_rule_list;
	vector<vector<float>> spin_states;
	vector<vector<vector<vector<float>>>> chem_motif_list;
	vector<vector<vector<vector<float>>>> spin_motif_list;
	vector<vector<vector<vector<float>>>> sro_motif_list;


	Session(void);
	Session(Session& _session);
	Session(string input_file);

	void _copy(Session& _session);
	void add_spin_states(string input_file);
	void fill_rule_list(string cluster_file);
	void fill_sro_list(string sro_file);
	void find_unique_dists();
};
#endif
