#include "sim_cell.h"

// This is the sim_cell class. It holds the state of the entire simulation cell. In the MC implementation, it's main components are
//copied to a separate vector to minimize the size or complexity of arguments being passed to certain functions. Still, it does
// serve as a useful tool for initializing the simulation cell and as a reference.
SimCell::SimCell(void) {
	sim_type = "EMPTY";
}

SimCell::SimCell(SimCell& sc_copy) {
	cutoff = sc_copy.cutoff;
	sim_type = sc_copy.sim_type;
	phase_init = sc_copy.phase_init;
	species_types = sc_copy.species_types;
	species_numbs = sc_copy.species_numbs;
	poscar_comp = sc_copy.poscar_comp;
	copy(sc_copy.sup_cell, sc_copy.sup_cell + 3, sup_cell);
	copy(sc_copy.LC, sc_copy.LC + 3, LC);
	numb_atoms = sc_copy.numb_atoms;
	copy(sc_copy.numb_cells, sc_copy.numb_cells + 3, numb_cells);
	copy(sc_copy.unit_LC, sc_copy.unit_LC + 3, unit_LC);
	atom_list = sc_copy.atom_list;
	unit_cell = sc_copy.unit_cell;
	unit_lat_fact = sc_copy.unit_lat_fact;
	unit_lat_vect = sc_copy.unit_lat_vect;
}

void SimCell::_copy(SimCell& sc_copy) {
	cutoff = sc_copy.cutoff;
	sim_type = sc_copy.sim_type;
	phase_init = sc_copy.phase_init;
	species_types = sc_copy.species_types;
	species_numbs = sc_copy.species_numbs;
	poscar_comp = sc_copy.poscar_comp;
	copy(sc_copy.sup_cell, sc_copy.sup_cell + 3, sup_cell);
	copy(sc_copy.LC, sc_copy.LC + 3, LC);
	numb_atoms = sc_copy.numb_atoms;
	copy(sc_copy.numb_cells, sc_copy.numb_cells + 3, numb_cells);
	copy(sc_copy.unit_LC, sc_copy.unit_LC + 3, unit_LC);
	atom_list = sc_copy.atom_list;
	unit_cell = sc_copy.unit_cell;
	unit_lat_fact = sc_copy.unit_lat_fact;
	copy(sc_copy.unit_lat_vect.begin(), sc_copy.unit_lat_vect.end(), back_inserter(unit_lat_vect));
}

// fill the unit cell using the poscar file 
void SimCell::fillUnitCell(string POSCAR_file, Session& sess) {
    int pos_start = 8;
    int comp_start = 6;
    ifstream POS_list;
	string pos_line;
    string type = "frac";
	vector<string> comp_line;
	vector<string> pos_lines;
	vector<string> LCs;
	vector<string> pos_list_s; //single?
	vector<int> pos_species_s;
	vector<vector<float>> pos_list_f; //full?
	vector<vector<int>> pos_species_f;
	vector<int> species_list;
	vector<float> pos;
	// read the poscar file
	POS_list.open(POSCAR_file, ifstream::in);
	if (POS_list.is_open()) {
		while (getline(POS_list, pos_line))
		{
			pos_lines.push_back(pos_line);
		}
		POS_list.close();
		cout << "POSCAR read\n";
	}
	else cout << "Unable to open POSCAR file\n";

	// get the lattice constants vectors and factor
	unit_lat_fact = stof(pos_lines[1]);
	for (int i = 2; i < 5; i++) {
		pos_line = pos_lines[i];
		LCs = split(pos_line);// , " ");
		unit_lat_vect.push_back({ unit_lat_fact * stof(LCs[0]), unit_lat_fact * stof(LCs[1]), unit_lat_fact * stof(LCs[2]) });
		unit_LC[i - 2] = sqrt(pow(unit_lat_vect[i - 2][0], 2) + pow(unit_lat_vect[i - 2][1], 2) + pow(unit_lat_vect[i - 2][2], 2));
	}
    cout << "Lattice constants read\n";

	// get the poscar type
    for (int i = 5; i < 8; i++) {
        if (pos_lines[i].find("artesian") != std::string::npos) {
            pos_start = i+1;
            comp_start = i-1;
            type = "cart";
        }
        else if (pos_lines[i].find("irect") != std::string::npos) {
            pos_start = i+1;
            comp_start = i-1;
            type = "frac";
        }
    }
    
	cout << "Type: " << type << "\n";

    // get the positions of each site and convert to cartesain
	for (int i = pos_start; i < pos_lines.size(); i++) {
		pos_line = pos_lines[i].substr(0, pos_lines[i].find("#"));
		pos_list_s = split(pos_line, " ");
		//for (int j = 0; j < pos_list_s.size(); j++) { cout << pos_list_s[j] << "\n"; }
		for (int j = 0; j < 3; j++) {pos.push_back(stof(pos_list_s[j]));}
		//for (int j = 0; j < pos.size(); j++) { cout << pos[j] << "\n"; }
        if (type == "frac") {
            pos = pos_transform(pos, unit_lat_vect);
        }
		pos_list_f.push_back(pos);
		pos.clear();
		// get the allowed species for each site
		for (int j = 3; j < pos_list_s.size(); j++) {
			string str = pos_list_s[j];
			str.erase(remove(str.begin(), str.end(), '\r'), str.end());
			str.erase(remove(str.begin(), str.end(), ' '), str.end());
			pos_species_s.push_back(get_index(sess.species_str, str));
		}
		if (pos_species_s.size() == 0) {
			pos_species_f.push_back(sess.species_inds);
		}
		else {
			pos_species_f.push_back(pos_species_s);
		}
        pos_species_s.clear();
	}
    
	cout << "Positions read\n";

	// assign the correct atomic species to each atom position
	comp_line = split(pos_lines[comp_start]);
	for (int i = 0; i < comp_line.size(); i++) {
		poscar_comp.push_back(stoi(comp_line[i]));
	}
	int comp_count = 0;
	int comp_running_sum = poscar_comp[comp_count];
	for (int i = 0; i < pos_list_f.size(); i++) {
		if (i < comp_running_sum) {
			species_list.push_back(comp_count);
		}
		else {
			bool inc_cc = true;
			while (inc_cc == true)
			{
				comp_count += 1;
				if (poscar_comp[comp_count] > 0) { inc_cc = false; }
			}
			comp_running_sum += poscar_comp[comp_count];
			species_list.push_back(comp_count);
		}
	}
    
	cout << "Species read\n";

	// create "atom" object for each atom in the unit cell
	for (int i = 0; i < pos_list_f.size(); i++) {
		unit_cell.push_back(Atom(i, species_list[i], 0, 0, pos_list_f[i], pos_species_f[i]));
	}
	cout << "Unit cell read\n";
}

// make a supercell (simulation cell) from the unit cell
void SimCell::make_supercell(Session& sess) {
	int x = sup_cell[0];
	int y = sup_cell[1];
	int z = sup_cell[2];
	vector<float>current_cell{ 0, 0, 0 };
	int index = 0;
	int spin = 0;
	int phase = 0;
//	float spin_rand;
	vector<float> new_atom_pos;
	for (int i = 0; i < 3; i++) { 
		LC[i] = unit_LC[i] * sup_cell[i]; 
		lat_vect.push_back({ unit_lat_vect[i][0] * sup_cell[i], unit_lat_vect[i][1] * sup_cell[i], unit_lat_vect[i][2] * sup_cell[i] });
	}
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	cout << "Looping through supercell deminsions\n";
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			for (int k = 0; k < z; k++) {
				current_cell[0] = i;
				current_cell[1] = j;
				current_cell[2] = k;
				for (int m = 0; m < unit_cell.size(); m++) {
					vector<float> uc_pos; for (int unit_i = 0; unit_i < 3; unit_i++) { uc_pos.push_back(unit_cell[m].pos[unit_i]); }
					vector<float> shifted_pos;
                    shifted_pos = pos_transform(current_cell, unit_lat_vect);
                    new_atom_pos = vect_add(uc_pos, shifted_pos);
					if (sess.spin_init[0] == 'F') { spin = vect_max(sess.spin_states[unit_cell[m].getSpecies()]); }
					else if (sess.spin_init[0] == 'R') {
						if (sess.spin_states.size() == 0) {
							if (unif(rng) <= 0.5) { spin = -1; }
							else { spin = 1; }
						}
						else {
							int rand_int = round(unif(rng) * (sess.spin_states[unit_cell[m].getSpecies()].size() - 1));
							spin = sess.spin_states[unit_cell[m].getSpecies()][rand_int];
						}
					}
					if (sess.phase_init[0] == 'A') { phase = 0; }
					else if (sess.phase_init[0] == 'M') { phase = 1; }
					else if (sess.phase_init[0] == 'R') {
						float phase_rand = unif(rng);
						if (phase_rand <= 1.0 / 3.0) { phase = -1; }
						else if (phase_rand > 1.0 / 3.0 && phase_rand <= 2.0 / 3.0) { phase = 0; }
						else { phase = 1; }
					}
					atom_list.push_back(Atom(index, unit_cell[m].getSpecies(), spin, phase, new_atom_pos, unit_cell[m].allowed_species));
					index += 1;
				}
			}
		}
	}
	numb_atoms = atom_list.size();
	if (sess.use_poscar == true) {
		bool poscar_flag = true;
		if (poscar_comp.size() != species_numbs.size()) { poscar_flag = false; }
		for (int i = 0; i < min(poscar_comp.size(), species_numbs.size()); i++) {
			if (poscar_comp[i] * sup_cell[0] * sup_cell[1] * sup_cell[2] != species_numbs[i]) { poscar_flag = false; }
		}
		if (poscar_flag == false) {
			cout << "The POSCAR and INPUT files have diffrent compositions. Using POSCAR for initialization\n"; //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			species_numbs.clear();
			for (int i = 0; i < poscar_comp.size(); i++) {
				species_numbs.push_back(poscar_comp[i] * sup_cell[0] * sup_cell[0] * sup_cell[2]);
			}
		}
	}
	else {
		vector<bool> reset_spec;
		for (int i = 0; i < numb_atoms; i++) { reset_spec.push_back(true); }
		for (int i = 0; i < species_numbs.size(); i++) {
			int spec_count = 0;
			for (int j = 0; j < numb_atoms; j++) {
				if (spec_count < species_numbs[i]) {
					if (get_index(atom_list[j].allowed_species, i) != -1 && reset_spec[j] == true) {
						atom_list[j].setSpecies(i);
						if (sess.spin_init[0] == 'F') { spin = vect_max(sess.spin_states[i]); }
						else if (sess.spin_init[0] == 'R') {
							if (sess.spin_states.size() == 0) {
								if (unif(rng) <= 0.5) { spin = -1; }
								else { spin = 1; }
							}
							else {
								int rand_int = round(unif(rng) * (sess.spin_states[i].size() - 1));
								spin = sess.spin_states[i][rand_int];
							}
						}
						atom_list[j].setSpin(spin);
						spec_count += 1;
						reset_spec[j] = false;
					}
				}
			}
		}
	}
	setNeighborDists(sess.unique_dists); // add neighbor distance for each atom
}

float SimCell::bc_dist(vector<float>& pos1, vector<float>& pos2) {
	vector<int> dir { -1, 1 };
	vector<float> dists { pos_dist(pos1, pos2) };
	vector<float> shift;
	vector<float> bc_pos { 0, 0, 0 };
	vector<float> lc_shift { 0, 0, 0 };
	
	for (int i = 0; i < 3; i++) { shift.push_back(pos2[i] - pos1[i]); }

	for (int a : dir) {
		for (int i = 0; i < 3; i++) {
			lc_shift = scale_vect(lat_vect[i], a);
			bc_pos = vect_add(pos2, lc_shift);
			dists.push_back(pos_dist(pos1, bc_pos));
		}
		for (int b : dir) {
			lc_shift = scale_vect(lat_vect[0], a);
			bc_pos = vect_add(pos2, lc_shift);
			lc_shift = scale_vect(lat_vect[1], b);
			bc_pos = vect_add(bc_pos, lc_shift);
			dists.push_back(pos_dist(pos1, bc_pos));

			lc_shift = scale_vect(lat_vect[0], a);
			bc_pos = vect_add(pos2, lc_shift);
			lc_shift = scale_vect(lat_vect[2], b);
			bc_pos = vect_add(bc_pos, lc_shift);
			dists.push_back(pos_dist(pos1, bc_pos));

			lc_shift = scale_vect(lat_vect[1], a);
			bc_pos = vect_add(pos2, lc_shift);
			lc_shift = scale_vect(lat_vect[2], b);
			bc_pos = vect_add(bc_pos, lc_shift);
			dists.push_back(pos_dist(pos1, bc_pos));

			for (int c : dir) {
				lc_shift = scale_vect(lat_vect[0], a);
				bc_pos = vect_add(pos2, lc_shift);
				lc_shift = scale_vect(lat_vect[1], b);
				bc_pos = vect_add(bc_pos, lc_shift);
				lc_shift = scale_vect(lat_vect[2], c);
				bc_pos = vect_add(bc_pos, lc_shift);
				dists.push_back(pos_dist(pos1, bc_pos));
			}
		}
	}

	return *min_element(dists.begin(), dists.end());
}

// fill the neighbor matrix using the dist_list (list of relevent neighbors from the mc_rules)
void SimCell::setNeighborDists(vector<float>& dist_list) {
	// declare variables
	vector<float> pos1 { 0, 0, 0 };
	vector<float> pos2 { 0, 0, 0 };
	float dist;
    float max_dist = vect_max(dist_list) + 0.001;
	int added = 0;
	// loop through each atom object
	for (int i = 0; i < atom_list.size(); i++) {
		pos1[0] = atom_list[i].pos[0];
		pos1[1] = atom_list[i].pos[1];
		pos1[2] = atom_list[i].pos[2];
		// loop through all neighbors
		for (int j = 0; j < atom_list.size(); j++) {
			pos2[0] = atom_list[j].pos[0];
			pos2[1] = atom_list[j].pos[1];
			pos2[2] = atom_list[j].pos[2];
			// find straight line distance in BC
			dist = bc_dist(pos1, pos2);
			added = 0;
			// add to neighbor distance list and neighbor index list if the distance is relevent
            if (dist <= max_dist ) {
                atom_list[i].neighbor_dists.push_back(dist);
                atom_list[i].neighbors.push_back(j);
            }
//			for (int k = 0; k < dist_list.size(); k++) {
//				if (fcomp(dist, dist_list[k])) {
//					atom_list[i].neighbor_dists.push_back(dist);
//					atom_list[i].neighbors.push_back(j);
//					added = 1;
//				}
//			}
		}
	}
}

float SimCell::findAtomDists(int atom1, int atom2) {
	auto it = std::find(atom_list[atom1].neighbors.begin(), atom_list[atom1].neighbors.end(), atom2);
	return distance(atom_list[atom1].neighbors.begin(), it);
}

// create the sim cell object and initalize it accordingly
void SimCell::initSimCell(string POSCAR_file, Session& session) {
	int _sup_cell[3] = { session.shape[0], session.shape[1], session.shape[2] };
	vector<int>& _species_numbs = session.atom_numbs;
	float _cutoff = 1; // Double check this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	string _sim_type = session.sim_type;
	string _phase_init = session.phase_init;
	string spin_init = session.spin_init;
	string species_init = session.species_init;
//	bool use_poscar = session.use_poscar;
	sim_type = _sim_type;
	phase_init = _phase_init;
	cutoff = _cutoff;
	for (int i = 0; i < _species_numbs.size(); i++) {
		species_types.push_back(i);
		species_numbs.push_back(_species_numbs[i]);
	}
	cout << "Reading POSCAR file\n";
	fillUnitCell(POSCAR_file, session);
	for (int i = 0; i < 3; i++) {
		sup_cell[i] = _sup_cell[i];
		LC[i] = unit_LC[i] * sup_cell[i];
	}
	vector<vector<float>> _pos_list;
	vector<int> _species_list;
	cout << "Making supercell\n";
	make_supercell(session);
	if (species_init == "RANDOM") {
		cout << "Randomizing species\n";
		randomizeSpecies(session);
	}
	//else if (species_init == "ordered") {
	//	cout << "Ordering species\n";
	//	orderSpecies(session);
	//}
	//else {
	//	cout << "Error: species_init not recognized\n";
	//}
}

// atom object inside SimCell class. used as a container for information about single atoms. Standard constructor and get/set functions are used
SimCell::Atom::Atom(void) {
	species = 0;
	spin = 0;
	phase = 0;
	index = 0;
}

SimCell::Atom::Atom(int _index, int _species, float _spin, int _phase, vector<float> _pos, vector<int> _allowed_species) {
	index = _index;
	species = _species;
	spin = _spin;
	phase = _phase;
	pos[0] = _pos[0];
	pos[1] = _pos[1];
	pos[2] = _pos[2];
	for (int i = 0; i < _allowed_species.size(); i++) {
		allowed_species.push_back(_allowed_species[i]);
	}
}

int SimCell::Atom::getNeighborSpin(int _neighbor, SimCell& sim_cell) {
	int neighbor_index = neighbors[_neighbor];
	int neighbor_spin = sim_cell.atom_list[neighbor_index].getSpin();
	return neighbor_spin;
}

int SimCell::Atom::getNeighborSpecies(int _neighbor, SimCell& sim_cell) {
	int neighbor_index = neighbors[_neighbor];
	int neighbor_species = sim_cell.atom_list[neighbor_index].getSpecies();
	return neighbor_species;
}

int SimCell::Atom::getNeighborPhase(int _neighbor, SimCell& sim_cell) {
	int neighbor_index = neighbors[_neighbor];
	int neighbor_phase = sim_cell.atom_list[neighbor_index].getPhase();
	return neighbor_phase;
}

int SimCell::Atom::getNeighborIndex(int _neighbor, SimCell& sim_cell) {
	return neighbors[_neighbor];
}

float SimCell::Atom::getNeighborDist(int _neighbor, SimCell& sim_cell) {
	return neighbor_dists[_neighbor];
}

int SimCell::Atom::getNumbNeighbors(int _site, SimCell& sim_cell) {
	return sim_cell.atom_list[_site].neighbors.size();
}

float SimCell::Atom::getSpin(void) {
	return spin;
}

int SimCell::Atom::getSpecies(void) {
	return species;
}

int SimCell::Atom::getPhase(void) {
	return phase;
}

void SimCell::Atom::setSpin(float _spin) {
	spin = _spin;
}

void SimCell::Atom::setSpecies(int _species) {
	species = _species;
}

void SimCell::Atom::setPhase(int _phase) {
	phase = _phase;
}

void SimCell::randomizeSpecies(Session& session) {
	for (int i = 0; i < numb_atoms; i++) {
		int rand_index = rand() % numb_atoms;
        rand_index = 0; // to be finished
	}
}
