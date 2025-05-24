#pragma once
#ifndef rule_h
#define rule_h
#include "utils.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;


class Rule {
private:
	float enrg_cont;
	int type; // chem or spin or spect_spin
	int spectator = 2; // spectator species in case of type 2 spin rule (spectating trimer). Default is 2 for no spectator
	int phase;
public:
	vector<vector<float>> motif;
	int clust_ind;
	vector<int> deco;
	vector<float> dists;
	Rule(void);
	Rule(float _enrg_cont, int _type, int _phase, vector<int> _deco, vector<vector<float>> _motif, int clust_ind);
	Rule(float _enrg_cont, vector<int> _type, int _phase, vector<int> _deco, vector<vector<float>> _motif, int clust_ind);

	//Rule(Rule & _rule);
	float GetEnrgCont();
	int GetType();
	int GetSpect();
	int GetPhase();
	int GetLength();
	bool IsRuleSpin(vector<int> test_species, vector<float> test_dists);
	bool IsRuleChem(vector<int> test_species, vector<float> test_dists);
	vector<int> GetDeco();
	vector<float> GetDists();
};

#endif
