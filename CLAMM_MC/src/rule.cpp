#include "rule.h"


// This is the rule class (mc_rules object). It servs as a container for the MC rules and is used to create the rule map
Rule::Rule(void) {
}

Rule::Rule(float _enrg_cont, int _type, int _phase, vector<int> _deco, vector<vector<float>> _motif, int _clust_ind) {
	motif = _motif;
	enrg_cont = _enrg_cont;
	type = _type;
	phase = _phase;
	deco = _deco;
	clust_ind = _clust_ind;
	for (int i = 0; i < motif.size(); i++) {
		for (int j = i + 1; j < motif.size(); j++) {
			float new_dist = 0;
			for (int k = 0; k < 3; k++) { new_dist += pow((motif[i][k] - motif[j][k]), 2.0); }
			new_dist = sqrt(new_dist);
			if (std::find(dists.begin(), dists.end(), new_dist) == dists.end()) {
				dists.push_back(new_dist);
			}
		}
	}
    if ( motif.size() == 1 ) {
        dists.push_back(0);
    }
}

Rule::Rule(float _enrg_cont, vector<int> _type, int _phase, vector<int> _deco, vector<vector<float>> _motif, int _clust_ind) {
	motif = _motif;
	enrg_cont = _enrg_cont;
	type = _type[0];
	spectator = _type[1];
	phase = _phase;
	deco = _deco;
	clust_ind = _clust_ind;
	for (int i = 0; i < motif.size(); i++) {
		for (int j = i + 1; j < motif.size(); j++) {
			float new_dist = 0;
			for (int k = 0; k < 3; k++) { new_dist += pow((motif[i][k] - motif[j][k]), 2.0); }
			new_dist = sqrt(new_dist);
			if (std::find(dists.begin(), dists.end(), new_dist) == dists.end()) {
				dists.push_back(new_dist);
			}
		}
	}
	if (motif.size() == 1) {
		dists.push_back(0);
	}
}

//Rule::Rule(Rule& _rule) {
//	motif = _rule.motif;
//	enrg_cont = _rule.enrg_cont;
//	type = _rule.type;
//	phase = _rule.phase;
//	deco = _rule.deco;
//	clust_ind = _rule.clust_ind;
//	dists = _rule.dists;
//}

int Rule::GetPhase() {
	return phase;
}

float Rule::GetEnrgCont() {
	return enrg_cont;
}

int Rule::GetType() {
	return type;
}

int Rule::GetSpect() {
	return spectator;
}

int Rule::GetLength() {
	return motif.size();
}

vector<int> Rule::GetDeco() {
	return deco;
}

vector<float> Rule::GetDists() {
	return dists;
}

//bool Rule::IsRuleChem(vector<int> test_deco, vector<float> test_dists) {
//	if (type != 0) { return false; }
//	else if (test_dists.size() != dists.size() or test_deco.size() != test_deco.size()) { return false; }
//	else {
//		for (int i = 0; i < test_dists.size(); i++) { test_dists[i] = round(test_dists[i] * 100000) / 100000; }
//		for (int i = 0; i < dists.size(); i++) { dists[i] = round(dists[i] * 100000) / 100000; }
//		if (motif.size() == 1) {
//			if (test_deco[0] != deco[0]) { return false; }
//			else { return true; }
//		}
//		else if (motif.size() == 2) {
//			if (test_dists[0] == dists[0]) {
//				if (test_deco[0] == deco[0] and test_deco[1] == deco[1]) { return true; }
//				else if (test_deco[0] == deco[1] and test_deco[1] == deco[0]) { return true; }
//				else { return false; }
//			}
//			else { return false; }
//		}
//		else if (deco.size() == 3) {
//			if (test_deco[0] == deco[0] and test_deco[1] == deco[1] and test_deco[2] == deco[2] and test_dists[0] == dists[0] and test_dists[1] == dists[1] and test_dists[2] == dists[2]) { return true; }
//			else if (test_deco[0] == deco[0] and test_deco[2] == deco[1] and test_deco[1] == deco[2] and test_dists[1] == dists[0] and test_dists[0] == dists[1] and test_dists[2] == dists[2]) { return true; }
//			else if (test_deco[1] == deco[0] and test_deco[0] == deco[1] and test_deco[2] == deco[2] and test_dists[0] == dists[0] and test_dists[2] == dists[1] and test_dists[1] == dists[2]) { return true; }
//			else if (test_deco[1] == deco[0] and test_deco[2] == deco[1] and test_deco[0] == deco[2] and test_dists[2] == dists[0] and test_dists[0] == dists[1] and test_dists[1] == dists[2]) { return true; }
//			else if (test_deco[2] == deco[0] and test_deco[0] == deco[1] and test_deco[1] == deco[2] and test_dists[1] == dists[0] and test_dists[2] == dists[1] and test_dists[0] == dists[2]) { return true; }
//			else if (test_deco[2] == deco[0] and test_deco[1] == deco[1] and test_deco[0] == deco[2] and test_dists[2] == dists[0] and test_dists[1] == dists[1] and test_dists[0] == dists[2]) { return true; }
//			else { return false; }
//		}
//	}
//}
//
//bool Rule::IsRuleSpin(vector<int> test_deco, vector<float> test_dists) {
//	if (type != 1) { return false; }
//	else if (test_dists.size() != dists.size() or test_deco.size() != test_deco.size()) { return false; }
//	else {
//		if (motif.size() == 1) {
//			if (test_deco[0] != deco[0]) { return false; }
//			else { return true; }
//		}
//		else if (motif.size() == 2) {
//			if (test_dists[0] == dists[0]) {
//				if (test_deco[0] == deco[0] and test_deco[1] == deco[1]) { return true; }
//				else if (test_deco[0] == deco[1] and test_deco[1] == deco[0]) { return true; }
//				else { return false; }
//			}
//			else { return false; }
//		}
//		else if (motif.size() == 3) {
//			if (test_deco[0] == deco[0] and test_deco[1] == deco[1] and test_deco[2] == deco[2] and test_dists[0] == dists[0] and test_dists[1] == dists[1] and test_dists[2] == dists[2]) { return true; }
//			else { return false; }
//		}
//	}
//}
