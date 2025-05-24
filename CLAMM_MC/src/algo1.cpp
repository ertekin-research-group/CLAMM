#include "algo1.h"

Algo1::Algo1(void) {}

Algo1::Algo1(Session& _session, SimCell& _sim_cell) {
    session = _session;
    sim_cell = _sim_cell;
    if (session.numb_passes < 1) {
        cout << "_______________________________________________________________________________" << endl;
        cout << "Error: Algo1 has been given 0 thermal average passes" << endl;
        cout << "This is probably not what you want..." << endl;
        cout << "_______________________________________________________________________________" << endl;
    }
}

size_t Algo1::cust_hash(vector<uint32_t>& vect) {
    std::size_t seed = vect.size();
    for (auto x : vect) {
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = (x >> 16) ^ x;
        seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

double Algo1::eval_site_chem(int site) {
    double enrg = 0;
    map<size_t, vector<double>>::iterator rule_itr;
    vector<uint32_t> rule_info;
    size_t rule_key;
    for (int i = 0; i < chem_motif_groups[site].size(); i++) {
        vector<vector<int>> motif = chem_motif_groups[site][i];
        for (int j = 0; j < motif.size(); j++) {
            rule_info.push_back(0); // chem ind
            rule_info.push_back(i); // clust_ind
            vector<int> group = motif[j];
            for (int k : group) {
                rule_info.push_back(chem_list[k]); // sites ind
            }
            rule_key = cust_hash(rule_info);
            rule_itr = rule_map_chem.find(rule_key);
            if (rule_itr != rule_map_chem.end()) {
                enrg += rule_itr->second[0] / group.size();
            }
            rule_info.clear();
        }
    }
    return enrg;
}

double Algo1::eval_site_spin(int site) {
    double enrg = 0;
    map<size_t, vector<double>>::iterator rule_itr;
    vector<uint32_t> rule_info;
    size_t rule_key;
    for (int i = 0; i < spin_motif_groups[site].size(); i++) {
        vector<vector<int>> motif = spin_motif_groups[site][i];
        for (int j = 0; j < motif.size(); j++) {
            rule_info.push_back(1); // spin ind
            rule_info.push_back(i); // clust_ind
            vector<int> group = motif[j];
            float spin_prod = 1;
            for (int k : group) {
                rule_info.push_back(chem_list[k]); // sites ind
                spin_prod *= spin_list[k];
            }
            rule_key = cust_hash(rule_info);
            rule_itr = rule_map_spin.find(rule_key);
            enrg += (rule_itr != rule_map_spin.end()) ? (rule_itr->second[0] * spin_prod / group.size()) : 0.0;
            rule_info.clear();
        }
    }
    return enrg;
}

double Algo1::eval_spin_flip(int site, float old_spin) {
    double enrg = 0;
    map<size_t, vector<double>>::iterator rule_itr;
    vector<uint32_t> rule_info;
    size_t rule_key;
    for (int i = 0; i < spin_motif_groups[site].size(); i++) {
        vector<vector<int>> motif = spin_motif_groups[site][i];
        for (int j = 0; j < motif.size(); j++) {
            rule_info.push_back(1);
            rule_info.push_back(i);
            vector<int> group = motif[j];
            float spin_prod = 1;
            for (int k : group) {
                rule_info.push_back(chem_list[k]);
                if (k != site) { spin_prod *= spin_list[k]; }
            }
            rule_key = cust_hash(rule_info);
            rule_itr = rule_map_spin.find(rule_key);
            enrg += (rule_itr != rule_map_spin.end()) ? (rule_itr->second[0] * spin_prod) : 0.0;
            rule_info.clear();
        }
    }
    return (enrg * spin_list[site] - enrg * old_spin);
}

vector<double> Algo1::eval_lat() {
    double enrg = 0;
    double chem_enrg = 0;
    double spin_enrg = 0;
    for (int site = 0; site < sim_cell.numb_atoms; site++) {
        chem_enrg += eval_site_chem(site);
        spin_enrg += eval_site_spin(site);
    }
    cout << "Chemical energy is " << chem_enrg / sim_cell.numb_atoms << " per atom\n";
    cout << "Spin energy is " << spin_enrg / sim_cell.numb_atoms << " per atom\n";
    vector<double> enrgs = { chem_enrg,  spin_enrg, session.intercept * sim_cell.numb_atoms };
    return enrgs;
}

double Algo1::eval_lat_spin() {
    double enrg = 0;
    for (int site = 0; site < sim_cell.numb_atoms; site++) {
        enrg += eval_site_spin(site);
    }
    return enrg;
}

bool Algo1::bc_check(vector<float> check_vect, vector<float>& pos) {
    bool bc_test = false;
    vector<int> dir{ -1, 1 };
    vector<float> bc_pos{ 0, 0, 0 };
    vector<float> lc_shift{ 0, 0, 0 };
    vector<vector<float>> bc_shifts{ pos };

    for (int a : dir) {
        for (int i = 0; i < 3; i++) {
            lc_shift = scale_vect(sim_cell.lat_vect[i], a);
            bc_pos = vect_add(pos, lc_shift);
            bc_shifts.push_back(bc_pos);
        }
        for (int b : dir) {
            lc_shift = scale_vect(sim_cell.lat_vect[0], a);
            bc_pos = vect_add(pos, lc_shift);
            lc_shift = scale_vect(sim_cell.lat_vect[1], b);
            bc_pos = vect_add(bc_pos, lc_shift);
            bc_shifts.push_back(bc_pos);

            lc_shift = scale_vect(sim_cell.lat_vect[0], a);
            bc_pos = vect_add(pos, lc_shift);
            lc_shift = scale_vect(sim_cell.lat_vect[2], b);
            bc_pos = vect_add(bc_pos, lc_shift);
            bc_shifts.push_back(bc_pos);

            lc_shift = scale_vect(sim_cell.lat_vect[1], a);
            bc_pos = vect_add(pos, lc_shift);
            lc_shift = scale_vect(sim_cell.lat_vect[2], b);
            bc_pos = vect_add(bc_pos, lc_shift);
            bc_shifts.push_back(bc_pos);

            for (int c : dir) {
                lc_shift = scale_vect(sim_cell.lat_vect[0], a);
                bc_pos = vect_add(pos, lc_shift);
                lc_shift = scale_vect(sim_cell.lat_vect[1], b);
                bc_pos = vect_add(bc_pos, lc_shift);
                lc_shift = scale_vect(sim_cell.lat_vect[2], c);
                bc_pos = vect_add(bc_pos, lc_shift);
                bc_shifts.push_back(bc_pos);
            }
        }
    }
    for (vector<float> check : bc_shifts) {
        if (pos_comp(check, check_vect)) {
            bc_test = true;
            break;
        }
    }
    return bc_test;
}

void Algo1::fill_SMG(vector<vector<int>>& neigh_ind_list) {
    vector<float> new_pos{ 0.0, 0.0, 0.0 };
    vector<float> self_site{ 0.0, 0.0, 0.0 };
    vector<int> sites;
    vector<vector<int>> groups;
    vector<vector<vector<int>>> motifs;
    for (int atom = 0; atom < sim_cell.numb_atoms; atom++) {
        for (vector<vector<vector<float>>> clust : session.spin_motif_list) {
            for (vector<vector<float>> motif : clust) {
                for (vector<float> shift : motif) {
                    if (pos_comp(shift, self_site)) { sites.push_back(atom); }
                    else {
                        new_pos = vect_add(pos_list[atom], shift);
                        for (int neigh : neigh_ind_list[atom]) {
                            if (bc_check(pos_list[neigh], new_pos)) {
                                sites.push_back(neigh);
                            }
                        }
                    }
                }
                groups.push_back(sites);
                sites.clear();
            }
            motifs.push_back(groups);
            groups.clear();
        }
        spin_motif_groups.push_back(motifs);
        motifs.clear();
    }
}

void Algo1::fill_CMG(vector<vector<int>>& neigh_ind_list) {
    vector<float> new_pos{ 0.0, 0.0, 0.0 };
    vector<float> self_site{ 0.0, 0.0, 0.0 };
    vector<int> sites;
    vector<vector<int>> groups;
    vector<vector<vector<int>>> motifs;
    for (int atom = 0; atom < sim_cell.numb_atoms; atom++) {
        for (vector<vector<vector<float>>> clust : session.chem_motif_list) {
            for (vector<vector<float>> motif : clust) {
                for (vector<float> shift : motif) {
                    if (pos_comp(shift, self_site)) { sites.push_back(atom); }
                    else {
                        new_pos = vect_add(pos_list[atom], shift);
                        for (int neigh : neigh_ind_list[atom]) {
                            if (bc_check(pos_list[neigh], new_pos)) {
                                sites.push_back(neigh);
                            }
                        }
                    }
                }
                groups.push_back(sites);
                sites.clear();
            }
            motifs.push_back(groups);
            groups.clear();
        }
        chem_motif_groups.push_back(motifs);
        motifs.clear();
    }
}

void Algo1::fill_SROMG(vector<vector<int>>& neigh_ind_list) {
    vector<float> new_pos{ 0.0, 0.0, 0.0 };
    vector<float> self_site{ 0.0, 0.0, 0.0 };
    vector<int> sites;
    vector<vector<int>> groups;
    vector<vector<vector<int>>> motifs;
    for (int atom = 0; atom < sim_cell.numb_atoms; atom++) {
        for (vector<vector<vector<float>>> clust : session.sro_motif_list) {
            for (vector<vector<float>> motif : clust) {
                for (vector<float> shift : motif) {
                    if (pos_comp(shift, self_site)) { sites.push_back(atom); }
                    else {
                        new_pos = vect_add(pos_list[atom], shift);
                        for (int neigh : neigh_ind_list[atom]) {
                            if (bc_check(pos_list[neigh], new_pos)) {
                                sites.push_back(neigh);
                            }
                        }
                    }
                }
                groups.push_back(sites);
                sites.clear();
            }
            motifs.push_back(groups);
            groups.clear();
        }
        sro_motif_groups.push_back(motifs);
        motifs.clear();
    }
}

void Algo1::print_state(string contcar_name, int temp) {
    vector<int> perm;
    vector<int> temp_spec;
    vector<float> temp_spin;
    vector<vector<float>> temp_pos;
    temp_spec.assign(chem_list.begin(), chem_list.end());
    temp_spin.assign(spin_list.begin(), spin_list.end());
    temp_pos = pos_list;// insert(temp_pos.end(), pos_list.begin(), pos_list.end());
    perm = vect_permut(temp_spin);
    sort_vect(temp_spin, perm);
    sort_vect(temp_pos, perm);
    sort_vect(temp_spec, perm);
    perm = vect_permut(temp_spec);
    sort_vect(temp_spin, perm);
    sort_vect(temp_pos, perm);
    sort_vect(temp_spec, perm);
    ofstream OUT_file;
    string file_name = contcar_name + "_" + to_string(temp);
    if (temp == -1) { file_name = contcar_name; }
    OUT_file.open(file_name);
    if (OUT_file.is_open()) {
        for (string spec : session.species_str) { OUT_file << " " << spec; }
        OUT_file << "\n 1 \n";
        for (int i = 0; i < 3; i++) {
            vector<float>vect = sim_cell.unit_lat_vect[i];
            for (int j = 0; j < 3; j++) { OUT_file << vect[j] * session.shape[i] << " "; }
            OUT_file << "\n";
        }
        for (string spec : session.species_str) { OUT_file << " " << spec; }
        OUT_file << "\n";
        for (int i = 0; i < sim_cell.species_numbs.size(); i++) { OUT_file << sim_cell.species_numbs[i] << " "; }
        OUT_file << "\nCartesian\n";
        for (int i = 0; i < sim_cell.numb_atoms; i++) {
            OUT_file << temp_pos[i][0] << " " << temp_pos[i][1] << " " << temp_pos[i][2] << " ";
            OUT_file << " # " << temp_spec[i] << " " << temp_spin[i] << "\n";
        }
    }
    OUT_file.close();
}

float Algo1::spec_spin(int site, int spec) {
    if (chem_list[site] == spec) {
		return spin_list[site];
	}
    else {
		return 0.0;
	}
}

float Algo1::spin_dot(int site) {
    double dot = 0;
    map<size_t, double>::iterator rule_itr;
    vector<uint32_t> rule_info;
    size_t rule_key;
    for (int i = 0; i < sro_motif_groups[site].size(); i++) {
        vector<vector<int>> motif = sro_motif_groups[site][i];
        for (int j = 0; j < motif.size(); j++) {
            rule_info.push_back(1); // spin ind
            rule_info.push_back(i); // clust_ind
            vector<int> group = motif[j];
            float spin_prod = 1;
            for (int k : group) {
                rule_info.push_back(chem_list[k]); // sites ind
                spin_prod *= spin_list[k];
            }
            rule_key = cust_hash(rule_info);
            rule_itr = rule_map_sro.find(rule_key);
            dot += (rule_itr != rule_map_sro.end()) ? (spin_prod / group.size()) : 0.0;
            rule_info.clear();
        }
    }
    return dot;
}

//void Algo1::run() {
//    // declare variables
//    RunningStat rs_E;
//    RunningStat rs_M;
//    RunningStat rs_M_spec;
//    RunningStat rs_Mdot;
//    double mag_energy = 0.0;
//    double chem_energy = 0.0;
//    double intercept = session.intercept;
//    double init_spin = 0.0;
//    double init_spin_spec = 0.0;
//    double init_spin_dot = 0.0;
//    bool same_spin;
//    int attempts = 0;
//    float rand_spin = 0.0;
//    float sro_target = session.sro_target;
//    float temp1 = session.start_temp;
//    float temp2 = session.end_temp;
//    float temp_step = session.temp_step;
//    float keep_rand;
//    float keep_prob;
//    vector<vector<float>> spin_states = session.spin_states;
//    int passes = session.ta_passes + session.eq_passes;
//    int ta_passes = session.ta_passes;
//    int eq_passes = session.eq_passes;
//    int numb_atoms = sim_cell.numb_atoms;
//    int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
//    vector<vector<int>> neigh_ind_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
//    vector<int> spin_atoms; // atomic species that can have spin
//
//    // setup rng for random spin choice and acceptance probability
//    std::mt19937_64 rng;
//    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//    std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
//    rng.seed(ss);
//    std::uniform_real_distribution<double> unif(0, 1);
//
//    cout << "Making atom list and neighbor index list\n";
//    // Make atom_list more acessable for species and spin and neighbors
//    for (int i = 0; i < sim_cell.numb_atoms; i++) {
//        chem_list.push_back(sim_cell.atom_list[i].getSpecies());
//        spin_list.push_back(sim_cell.atom_list[i].getSpin());
//        pos_list.push_back({ sim_cell.atom_list[i].pos[0], sim_cell.atom_list[i].pos[1], sim_cell.atom_list[i].pos[2] });
//        for (int j = 0; j < numb_neighbors; j++) {
//            neigh_ind_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
//        }
//    }
//
//    cout << "Making rule maps\n";
//    // Make rule_maps for easy lookup
//    float indchem = 0;
//    float indspin = 0;
//    size_t rule_key;
//    vector<uint32_t> rule_info;
//    map<size_t, vector<double>>::iterator rule_itr;
//
//    for (Rule rule : session.chem_rule_list) {
//        rule_info.push_back(rule.GetType());
//        rule_info.push_back(rule.clust_ind);
//        for (int i = 0; i < rule.deco.size(); i++) { rule_info.push_back(rule.deco[i]); }
//        rule_key = cust_hash(rule_info);
//        rule_itr = rule_map_chem.find(rule_key);
//        if (rule_itr == rule_map_chem.end()) {
//            vector<double> enrg_ind_pair = { rule.GetEnrgCont(), indchem };
//            rule_map_chem.insert(pair<size_t, vector<double>>(rule_key, enrg_ind_pair));
//            indchem += 1;
//        }
//        rule_info.clear();
//    }
//    for (Rule rule : session.spin_rule_list) {
//        rule_info.push_back(rule.GetType());
//        rule_info.push_back(rule.clust_ind);
//        for (int i = 0; i < rule.deco.size(); i++) { rule_info.push_back(rule.deco[i]); }
//        rule_key = cust_hash(rule_info);
//        rule_itr = rule_map_spin.find(rule_key);
//        if (rule_itr == rule_map_spin.end()) {
//            vector<double> enrg_ind_pair = { rule.GetEnrgCont(), indspin };
//            rule_map_spin.insert(pair<size_t, vector<double>>(rule_key, enrg_ind_pair));
//            indspin += 1;
//            for (int atom : rule.deco) {
//                if (find(spin_atoms.begin(), spin_atoms.end(), atom) == spin_atoms.end()) { spin_atoms.push_back(atom); } // initialize spin_atoms
//            }
//		}
//        rule_info.clear();
//    }
//    for (Rule rule : session.sro_rule_list) {
//        rule_info.push_back(rule.GetType());
//        rule_info.push_back(rule.clust_ind);
//        for (int i = 0; i < rule.deco.size(); i++) { rule_info.push_back(rule.deco[i]); }
//        rule_key = cust_hash(rule_info);
//        rule_map_sro.insert(pair<size_t, double>(rule_key, rule.GetEnrgCont()));
//        for (int atom : rule.deco) {
//            if (find(spin_atoms.begin(), spin_atoms.end(), atom) == spin_atoms.end()) { spin_atoms.push_back(atom); } // initialize spin_atoms
//        }
//        rule_info.clear();
//    }
//
//    // Fill motif group lists
//    cout << "Making chem motif group lists\n";
//    fill_CMG(neigh_ind_list);
//    cout << "Making spin motif group lists\n";
//    fill_SMG(neigh_ind_list);
//    fill_SROMG(neigh_ind_list);
//    // Create seperate output file to avoid race condition
//    string file_name = "OUTPUT";
//    string contcar_name = "CONTCAR";
//    bool file_exists = true;
//    while (file_exists == true) {
//        const char* c_file = file_name.c_str();
//        int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
//        if (fd < 0) {
//            // file exists or otherwise uncreatable
//            outfile_count += 1;
//            file_name = "OUTPUT" + to_string(outfile_count);
//            contcar_name = "CONTCAR" + to_string(outfile_count);
//        }
//        else {
//            file_exists = false;
//            close(fd);
//        }
//    }
//    const char* c_file = file_name.c_str();
//    ofstream Output;
//    Output.open(c_file);
//    // Output energy and spin for convergence test
//    ofstream Output_converge;
//    if (session.do_conv_output) { Output_converge.open("OUTPUT_CONVERG"); }
//
//    // Begin MC
//    vector<double> init_enrgs = eval_lat();
//    double e0 = init_enrgs[0] + init_enrgs[1] + init_enrgs[2];
//    chem_energy = init_enrgs[0];
//    mag_energy = init_enrgs[1];
//    cout << "Initial chemical energy is " <<  chem_energy / numb_atoms << " per atom\n";
//    cout << "Initial spin energy is " <<  mag_energy / numb_atoms << " per atom\n";
//    cout << "Intercept energy is " <<  init_enrgs[2] / numb_atoms << " per atom\n";
//    cout << "Initial energy energy is " <<  e0/ numb_atoms << " per atom\n";
//    for (int site = 0; site < numb_atoms; site++) {
//        if (find(spin_atoms.begin(), spin_atoms.end(), chem_list[site]) != spin_atoms.end()) {
//            init_spin += spin_list[site];
//            init_spin_spec += spec_spin(site, 1);
//            init_spin_dot += spin_dot(site);
//        }
//    }
//    cout << "Initial spin is " << init_spin / numb_atoms << " per atom\n";
//    Output << "Using Algo1 for spin flip";
//    Output << "\nPhase: " << sim_cell.phase_init;
//    Output << "\nComposition: ";
//    for (int i = 0; i < sim_cell.species_numbs.size(); i++) { Output << sim_cell.species_numbs[i] << " "; }
//    Output << "\nMC passes: " << passes;
//    Output << "\nInitial spin energy per atom: ";
//    Output << init_enrgs[1] / numb_atoms;
//    Output << "\nInitial spin per atom: ";
//    Output << init_spin / numb_atoms;
//    Output << "\nInitial spin spec per atom: ";
//    Output << init_spin_spec / session.atom_numbs[1];
//    Output << "\nInitial spin dot per atom: ";
//    Output << init_spin_dot / numb_atoms;
//    Output << "\ntemp, enrg, mag, mag_spec, mag_dot, var_e, var_spin, Cmag, Xmag, flip_count, flip_count2" << endl;
//
//    //    // initalize system with desired SRO
//    //    Output << "EQ passes: " << session.eq_passes << ", EQ Temp: " << session.sro_temp << "\n";
//    //    Output << "SRO Target: " << session.sro_target << "\n";
//    //    cout << "SRO Target: " << session.sro_target << "\n";
//    //    cout << "Starting Real MC\n";
//
//    // Start MC loop
//    cout << "Entering main loop\n";
//    int temp_count = 0; // index of temperature
//    float inc_dir = 1; // temp increment direction
//    if (signbit(temp2 - temp1) == 1) { inc_dir = -1; }
//    for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_step * inc_dir) {
//        int flip_count = 0;
//        int flip_count2 = 0;
//        for (int pass = 0; pass < passes; pass++) {
//            for (int site = 0; site < numb_atoms; site++) {
//                double e_flip = 0.0;
//                float spin_flip = 0.0;
//                float prod_flip = 0.0;
//                float spec_spin_flip = 0.0;
//                float old_spin = spin_list[site];
//                float new_spin = 0.0;
//                float old_prod = spin_dot(site);
//                float new_prod = 0.0;
//                int state = 1;
//                while (state != DONE) {
//                    switch (state) {
//                    case 0:
//                        state = DONE;
//                        break;
//                    case 1:
//                        if (find(spin_atoms.begin(), spin_atoms.end(), chem_list[site]) == spin_atoms.end()) {
//                            state = 0; //NO_SPIN
//                            break;
//                        }// state set to 0
//                        if (spin_states[chem_list[site]].size() <= 1) {
//                            state = 0; //ONE_SPIN
//                            break;
//                        }
//                        // Flip Spin
//                        attempts = 0;
//                        same_spin = true;
//                        while (same_spin == true and attempts < 20) {
//                            rand_spin = unif(rng);
//                            for (int it_spin_state = 0; it_spin_state < spin_states[chem_list[site]].size(); it_spin_state++) {
//                                if (rand_spin > float(it_spin_state) / float(spin_states[chem_list[site]].size())) {
//                                    new_spin = spin_states[chem_list[site]][it_spin_state];
//                                }
//                            }
//                            if (new_spin != old_spin) { same_spin = false; }
//                            attempts += 1;
//                        }
//                        if (attempts >= 20) {
//                            state = 0;
//                            break;
//                        }
//                        spin_list[site] = new_spin;
//                        e_flip = eval_spin_flip(site, old_spin);
//                        spin_flip = new_spin - old_spin;
//                        new_prod = spin_dot(site);
//                        prod_flip = new_prod - old_prod;
//                        if (e_flip < 0) { flip_count += 1; }
//                        else {
//                            keep_rand = unif(rng);
//                            keep_prob = exp(-1 / (Kb * temp) * (e_flip));
//                            if (keep_rand < keep_prob) { flip_count2 += 1; }
//                            else { spin_list[site] = old_spin; e_flip = 0.0; spin_flip = 0.0; prod_flip = 0.0; }
//                        }
//                        // Record the enrg and spin changes
//                        mag_energy += e_flip;
//                        init_spin += spin_flip;
//                        init_spin_dot += prod_flip;
//                        if (chem_list[site] == 1) { init_spin_spec += spin_flip; }
//                        state = DONE;
//                        break;
//                    }
//                }
//                if (pass >= eq_passes) {
//                    rs_E.Push(mag_energy);
//                    rs_M.Push(init_spin);
//                    rs_Mdot.Push(init_spin_dot);
//                    if (chem_list[site] == 1) { rs_M_spec.Push(init_spin_spec); }
//                }
//            }
//            if (session.do_conv_output) {
//                Output_converge << mag_energy << " " << init_spin << endl;
//            }
//        }
//        double e_avg = (rs_E.Mean() + init_enrgs[0] + init_enrgs[2]) / (numb_atoms);
//        double spin_avg = rs_M.Mean() / numb_atoms;
//        double spin_dot_avg = rs_Mdot.Mean() / numb_atoms;
//        double spin_spec_avg = rs_M_spec.Mean() / session.atom_numbs[1];
//        double var_e = rs_E.Variance();
//        double var_spin = rs_M.Variance();
//        double Cmag = var_e / (Kb * pow(temp, 2));
//        double Xmag = var_spin / (Kb * pow(temp, 2));
//        Output << temp << ", "
//            << e_avg << ", "
//            << spin_avg << ", "
//            << spin_spec_avg << ", "
//            << spin_dot_avg << ", "
//            << var_e << ", "
//            << var_spin << ", "
//            << Cmag << ", "
//            << Xmag << ", "
//            << flip_count << ", "
//            << flip_count2 << endl;
//        rs_E.Clear();
//        rs_M.Clear();
//        rs_Mdot.Clear();
//        rs_M_spec.Clear();
//        if (session.write_contcars == true) {
//            print_state(contcar_name, temp_count);
//        }
//        temp_count += 1;
//    }
//    print_state("CONTCAR_FINAL", -1);
//    Output.close();
//    cout << " MC Finished\n";
//}

void Algo1::run() {
    //declare variables
    RunningStat rs_E;
    RunningStat rs_M;
    RunningStat rs_M_spec;
    RunningStat rs_Mdot;
    double mag_energy = 0.0;
    double chem_energy = 0.0;
    double intercept = session.intercept;
    double init_spin = 0.0;
    double init_spin_spec = 0.0;
    double init_spin_dot = 0.0;
    bool same_spin;
    int attempts = 0;
    float rand_spin = 0.0;
    float sro_target = session.sro_target;
    float temp1 = session.start_temp;
    float temp2 = session.end_temp;
    float temp_step = session.temp_step;
    float keep_rand;
    float keep_prob;
    vector<vector<float>> spin_states = session.spin_states;
    int passes = session.ta_passes + session.eq_passes;
    int ta_passes = session.ta_passes;
    int eq_passes = session.eq_passes;
    int numb_atoms = sim_cell.numb_atoms;
    int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
    vector<vector<int>> neigh_ind_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
    vector<int> spin_atoms; // atomic species that can have spin

    //setup rng for random spin choice and acceptance probability
    std::mt19937_64 rng;
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
    rng.seed(ss);
    std::uniform_real_distribution<double> unif(0, 1);

    cout << "Making atom list and neighbor index list\n";
    //Make atom_list more acessable for species and spin and neighbors
    for (int i = 0; i < sim_cell.numb_atoms; i++) {
        chem_list.push_back(sim_cell.atom_list[i].getSpecies());
        spin_list.push_back(sim_cell.atom_list[i].getSpin());
        pos_list.push_back({ sim_cell.atom_list[i].pos[0], sim_cell.atom_list[i].pos[1], sim_cell.atom_list[i].pos[2] });
        for (int j = 0; j < numb_neighbors; j++) {
            neigh_ind_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
        }
    }

    cout << "Making rule maps\n";
    //Make rule_maps for easy lookup
    float indchem = 0;
    float indspin = 0;
    size_t rule_key;
    vector<uint32_t> rule_info;
    map<size_t, vector<double>>::iterator rule_itr;

    for (Rule rule : session.chem_rule_list) {
        rule_info.push_back(rule.GetType());
        rule_info.push_back(rule.clust_ind);
        for (int i = 0; i < rule.deco.size(); i++) { rule_info.push_back(rule.deco[i]); }
        rule_key = cust_hash(rule_info);
        rule_itr = rule_map_chem.find(rule_key);
        if (rule_itr == rule_map_chem.end()) {
            vector<double> enrg_ind_pair = { rule.GetEnrgCont(), indchem };
            rule_map_chem.insert(pair<size_t, vector<double>>(rule_key, enrg_ind_pair));
            indchem += 1;
        }
        rule_info.clear();
    }
    for (Rule rule : session.spin_rule_list) {
        rule_info.push_back(rule.GetType());
        rule_info.push_back(rule.clust_ind);
        for (int i = 0; i < rule.deco.size(); i++) { rule_info.push_back(rule.deco[i]); }
        rule_key = cust_hash(rule_info);
        rule_itr = rule_map_spin.find(rule_key);
        if (rule_itr == rule_map_spin.end()) {
            vector<double> enrg_ind_pair = { rule.GetEnrgCont(), indspin };
            rule_map_spin.insert(pair<size_t, vector<double>>(rule_key, enrg_ind_pair));
            indspin += 1;
            for (int atom : rule.deco) {
                if (find(spin_atoms.begin(), spin_atoms.end(), atom) == spin_atoms.end()) { spin_atoms.push_back(atom); } // initialize spin_atoms
            }
        }
        rule_info.clear();
    }
    for (Rule rule : session.sro_rule_list) {
        rule_info.push_back(rule.GetType());
        rule_info.push_back(rule.clust_ind);
        for (int i = 0; i < rule.deco.size(); i++) { rule_info.push_back(rule.deco[i]); }
        rule_key = cust_hash(rule_info);
        rule_map_sro.insert(pair<size_t, double>(rule_key, rule.GetEnrgCont()));
        for (int atom : rule.deco) {
            if (find(spin_atoms.begin(), spin_atoms.end(), atom) == spin_atoms.end()) { spin_atoms.push_back(atom); } // initialize spin_atoms
        }
        rule_info.clear();
    }

    //Fill motif group lists
    cout << "Making chem motif group lists\n";
    fill_CMG(neigh_ind_list);
    cout << "Making spin motif group lists\n";
    fill_SMG(neigh_ind_list);
    fill_SROMG(neigh_ind_list);
    //Create seperate output file to avoid race condition
    string file_name = "OUTPUT";
    string contcar_name = "CONTCAR";
    bool file_exists = true;
    while (file_exists == true) {
        const char* c_file = file_name.c_str();
        int fd = openWithFcntl(c_file);
        if (fd < 0) {
            //file exists or otherwise uncreatable
            outfile_count += 1;
            file_name = "OUTPUT" + to_string(outfile_count);
            contcar_name = "CONTCAR" + to_string(outfile_count);
        }
        else {
            file_exists = false;
            //close(fd);
        }
    }
    const char* c_file = file_name.c_str();
    ofstream Output;
    Output.open(c_file);
    //Output energy and spin for convergence test
    ofstream Output_converge;
    if (session.do_conv_output) { Output_converge.open("OUTPUT_CONVERG"); }

    //Begin MC
    vector<double> init_enrgs = eval_lat();
    double e0 = init_enrgs[0] + init_enrgs[1] + init_enrgs[2];
    chem_energy = init_enrgs[0];
    mag_energy = init_enrgs[1];
    cout << "Initial chemical energy is " << chem_energy / numb_atoms << " per atom\n";
    cout << "Initial spin energy is " << mag_energy / numb_atoms << " per atom\n";
    cout << "Intercept energy is " << init_enrgs[2] / numb_atoms << " per atom\n";
    cout << "Initial energy energy is " << e0 / numb_atoms << " per atom\n";
    for (int site = 0; site < numb_atoms; site++) {
        if (find(spin_atoms.begin(), spin_atoms.end(), chem_list[site]) != spin_atoms.end()) {
            init_spin += spin_list[site];
            init_spin_spec += spec_spin(site, 1);
            init_spin_dot += spin_dot(site);
        }
    }
    cout << "Initial spin is " << init_spin / numb_atoms << " per atom\n";
    Output << "Using Algo1 for spin flip";
    Output << "\nPhase: " << sim_cell.phase_init;
    Output << "\nComposition: ";
    for (int i = 0; i < sim_cell.species_numbs.size(); i++) { Output << sim_cell.species_numbs[i] << " "; }
    Output << "\nMC passes: " << passes;
    Output << "\nInitial spin energy per atom: ";
    Output << init_enrgs[1] / numb_atoms;
    Output << "\nInitial spin per atom: ";
    Output << init_spin / numb_atoms;
    Output << "\nInitial spin spec per atom: ";
    Output << init_spin_spec / session.atom_numbs[1];
    Output << "\nInitial spin dot per atom: ";
    Output << init_spin_dot / numb_atoms;
    Output << "\ntemp, enrg, mag, mag_spec, mag_dot, var_e, var_spin, Cmag, Xmag, flip_count, flip_count2" << endl;

        // initalize system with desired SRO
        Output << "EQ passes: " << session.eq_passes << ", EQ Temp: " << session.sro_temp << "\n";
        Output << "SRO Target: " << session.sro_target << "\n";
        cout << "SRO Target: " << session.sro_target << "\n";
        cout << "Starting Real MC\n";

    //Start MC loop
    cout << "Entering main loop\n";
    int temp_count = 0; // index of temperature
    float inc_dir = 1; // temp increment direction
    if (signbit(temp2 - temp1) == 1) { inc_dir = -1; }
    for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_step * inc_dir) {
        int flip_count = 0;
        int flip_count2 = 0;
        for (int pass = 0; pass < passes; pass++) {
            for (int site = 0; site < numb_atoms; site++) {
                double e_flip = 0.0;
                float spin_flip = 0.0;
                float prod_flip = 0.0;
                float spec_spin_flip = 0.0;
                float old_spin = spin_list[site];
                float new_spin = 0.0;
                float old_prod = spin_dot(site);
                float new_prod = 0.0;
                int state = 1;
                while (state != DONE) {
                    switch (state) {
                    case 0:
                        state = DONE;
                        break;
                    case 1:
                        if (find(spin_atoms.begin(), spin_atoms.end(), chem_list[site]) == spin_atoms.end()) {
                            state = 0; //NO_SPIN
                            break;
                        }// state set to 0
                        if (spin_states[chem_list[site]].size() <= 1) {
                            state = 0; //ONE_SPIN
                            break;
                        }
                        double old_enrg = eval_site_spin(site);
                        //Flip Spin
                        attempts = 0;
                        same_spin = true;
                        while (same_spin == true && attempts < 20) {
                            rand_spin = unif(rng);
                            for (int it_spin_state = 0; it_spin_state < spin_states[chem_list[site]].size(); it_spin_state++) {
                                if (rand_spin > float(it_spin_state) / float(spin_states[chem_list[site]].size())) {
                                    new_spin = spin_states[chem_list[site]][it_spin_state];
                                }
                            }
                            if (new_spin != old_spin) { same_spin = false; }
                            attempts += 1;
                        }
                        if (attempts >= 20) {
                            state = 0;
                            break;
                        }
                        spin_list[site] = new_spin;
                        double new_enrg = eval_site_spin(site);
                        e_flip = 2*(new_enrg - old_enrg);
                        e_flip = eval_spin_flip(site, old_spin);
                        spin_flip = new_spin - old_spin;
                        new_prod = spin_dot(site);
                        prod_flip = new_prod - old_prod;
                        if (e_flip < 0) { flip_count += 1; }
                        else {
                            keep_rand = unif(rng);
                            keep_prob = exp(-1 / (Kb * temp) * (e_flip));
                            if (keep_rand < keep_prob) { flip_count2 += 1; }
                            else { spin_list[site] = old_spin; e_flip = 0.0; spin_flip = 0.0; prod_flip = 0.0; }
                        }
                        //Record the enrg and spin changes
                        mag_energy += e_flip;
                        init_spin += spin_flip;
                        init_spin_dot += prod_flip;
                        if (chem_list[site] == 1) { init_spin_spec += spin_flip; }
                        state = DONE;
                        break;
                    }
                }
                if (pass >= eq_passes) {
                    rs_E.Push(mag_energy);
                    rs_M.Push(init_spin);
                    rs_Mdot.Push(init_spin_dot);
                    if (chem_list[site] == 1) { rs_M_spec.Push(init_spin_spec); }
                }
            }
            if (session.do_conv_output) {
                Output_converge << mag_energy << " " << init_spin << endl;
            }
        }
        double e_avg = (rs_E.Mean() + init_enrgs[0] + init_enrgs[2]) / (numb_atoms);
        double spin_avg = rs_M.Mean() / numb_atoms;
        double spin_dot_avg = rs_Mdot.Mean() / numb_atoms;
        double spin_spec_avg = rs_M_spec.Mean() / session.atom_numbs[1];
        double var_e = rs_E.Variance();
        double var_spin = rs_M.Variance();
        double Cmag = var_e / (Kb * pow(temp, 2));
        double Xmag = var_spin / (Kb * pow(temp, 2));
        Output << temp << ", "
            << e_avg << ", "
            << rs_E.Mean() / numb_atoms << ", "
            << spin_avg << ", "
            << spin_spec_avg << ", "
            << spin_dot_avg << ", "
            << var_e << ", "
            << var_spin << ", "
            << Cmag << ", "
            << Xmag << ", "
            << flip_count << ", "
            << flip_count2 << endl;
        rs_E.Clear();
        rs_M.Clear();
        rs_Mdot.Clear();
        rs_M_spec.Clear();
        if (session.write_contcars == true) {
            print_state(contcar_name, temp_count);
        }
        temp_count += 1;
    }
    print_state("CONTCAR_FINAL", -1);
    Output.close();
    cout << " MC Finished\n";
}


//#include "algo1.h"
//
//Algo1::Algo1(void) {}
//
//Algo1::Algo1(Session& _session, SimCell& _sim_cell) {
//    session = _session;
//    sim_cell = _sim_cell;
//    passes = session.ta_passes + session.eq_passes;
//    if (session.ta_passes < 1) {
//        cout << "_______________________________________________________________________________" << endl;
//        cout << "Error: Algo1 has been given 0 thermal average passes" << endl;
//        cout << "This is probably not what you want..." << endl;
//        cout << "_______________________________________________________________________________" << endl;
//    }
//    // setup rng for random spin choice and acceptance probability
//    rng.seed(ss);
//    unif = std::uniform_real_distribution<double>(0.0, 1.0);
//    rand_atom = std::uniform_int_distribution<int>(0, sim_cell.numb_atoms - 1);
//    rand_method = std::uniform_int_distribution<int>(0, passes - 1);
//}
//
//size_t Algo1::cust_hash(vector<uint32_t>& vect) {
//    std::size_t seed = vect.size();
//    for (auto x : vect) {
//        x = ((x >> 16) ^ x) * 0x45d9f3b;
//        x = ((x >> 16) ^ x) * 0x45d9f3b;
//        x = (x >> 16) ^ x;
//        seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//    }
//    return seed;
//}
//
//double Algo1::eval_site_chem(int site) {
//    double enrg = 0;
//    map<size_t, vector<double>>::iterator rule_itr;
//    vector<uint32_t> rule_info;
//    size_t rule_key;
//    if (chem_list[site] == 4) {} //pass eval for vac site
//    else {
//        for (int i = 0; i < chem_motif_groups[site].size(); i++) {
//            vector<vector<int>> motif = chem_motif_groups[site][i];
//            for (int j = 0; j < motif.size(); j++) {
//                int vac_flag = 0;
//                rule_info.push_back(0); // chem ind
//                rule_info.push_back(i); // clust ind
//                vector<int> group = motif[j];
//                for (int k : group) {
//                    if (chem_list[k] == 4) { vac_flag = 1; break; }
//                    rule_info.push_back(chem_list[k]); // site spec ind
//                }
//                if (vac_flag == 1) { enrg += 0; }
//                else if (vac_flag == 0) {
//                    rule_key = cust_hash(rule_info);
//                    rule_itr = rule_map_chem.find(rule_key);
//                    if (rule_itr != rule_map_chem.end()) {
//                        enrg += rule_itr->second[0] / group.size();
//                        lat_rule_count_list[round(rule_itr->second[1])] += 1.0 ;// / group.size();////////////////////////////////////////
//                    }
//                }
//                else { cout << "error in cluster energy evaluation!"; }
//                rule_info.clear();
//            }
//        }
//    }
//    return enrg;
//}
//
//double Algo1::eval_site_spin(int site) {
//    double enrg = 0;
//    map<size_t, double>::iterator rule_itr;
//    vector<uint32_t> rule_info;
//    size_t rule_key;
//    if (chem_list[site] == 4) {} //pass eval for vac site
//    else {
//        for (int i = 0; i < spin_motif_groups[site].size(); i++) {
//            vector<vector<int>> motif = spin_motif_groups[site][i];
//            for (int j = 0; j < motif.size(); j++) {
//                int vac_flag = 0;
//                rule_info.push_back(1); // spin ind
//                rule_info.push_back(i); // clust ind
//                vector<int> group = motif[j];
//                float spin_prod = 1;
//                for (int k : group) {
//                    if (chem_list[k] == 4) { vac_flag = 1; break; }
//                    rule_info.push_back(chem_list[k]); // site spec ind
//                    spin_prod *= spin_list[k];
//                }
//                if (vac_flag == 1) { enrg += 0; }
//                else if (vac_flag == 0) {
//                    rule_key = cust_hash(rule_info);
//                    rule_itr = rule_map_spin.find(rule_key);
//                    enrg += (rule_itr != rule_map_spin.end()) ? (rule_itr->second * spin_prod) : 0.0;// / group.size()) : 0.0;/////////////////////////////////////////
//                }
//                else { cout << "error in cluster energy evaluation!"; }
//                rule_info.clear();
//            }
//        }
//    }
//    return enrg;
//}
//
//double Algo1::eval_spin_flip(int site, float old_spin) {
//    double enrg = 0;
//    map<size_t, double>::iterator rule_itr;
//    vector<uint32_t> rule_info;
//    size_t rule_key;
//    if (chem_list[site] == 4) {} //pass eval for vac site
//    else {
//        for (int i = 0; i < spin_motif_groups[site].size(); i++) {
//            vector<vector<int>> motif = spin_motif_groups[site][i];
//            for (int j = 0; j < motif.size(); j++) {
//                int vac_flag = 0;
//                rule_info.push_back(1);
//                rule_info.push_back(i);
//                vector<int> group = motif[j];
//                float spin_prod = 1;
//                for (int k : group) {
//                    if (chem_list[k] == 4) { vac_flag = 1; break; }
//                    rule_info.push_back(chem_list[k]);
//                    if (k != site) { spin_prod *= spin_list[k]; }
//                }
//                if (vac_flag == 1) { enrg += 0; }
//                else if (vac_flag == 0) {
//                    rule_key = cust_hash(rule_info);
//                    rule_itr = rule_map_spin.find(rule_key);
//                    enrg += (rule_itr != rule_map_spin.end()) ? (rule_itr->second * spin_prod) : 0.0;// / group.size()) : 0.0;///////////////////////////////////////
//                }
//                else { cout << "error in cluster energy evaluation!"; }
//            }
//        }
//    }
//    return (enrg * spin_list[site] - enrg * old_spin);
//}
//
//double Algo1::eval_lat() {
//    double enrg = 0;
//    for (int site = 0; site < sim_cell.numb_atoms; site++) {
//        enrg += eval_site_chem(site);
//        enrg += eval_site_spin(site);
//        if (chem_list[site] == 4) { numb_vac += 1; }
//    }
//    return enrg + session.intercept * (sim_cell.numb_atoms - numb_vac);
//}
//
//double Algo1::eval_lat_spin() {
//    double enrg = 0;
//    for (int site = 0; site < sim_cell.numb_atoms; site++) {
//        enrg += eval_site_spin(site);
//    }
//    return enrg;
//}
//
//bool Algo1::bc_check(vector<float> check_vect, vector<float>& pos) {
//    bool bc_test = false;
//    vector<int> dir{ -1, 1 };
//    vector<float> bc_pos{ 0, 0, 0 };
//    vector<float> lc_shift{ 0, 0, 0 };
//    vector<vector<float>> bc_shifts{ pos };
//
//    for (int a : dir) {
//        for (int i = 0; i < 3; i++) {
//            lc_shift = scale_vect(sim_cell.lat_vect[i], a);
//            bc_pos = vect_add(pos, lc_shift);
//            bc_shifts.push_back(bc_pos);
//        }
//        for (int b : dir) {
//            lc_shift = scale_vect(sim_cell.lat_vect[0], a);
//            bc_pos = vect_add(pos, lc_shift);
//            lc_shift = scale_vect(sim_cell.lat_vect[1], b);
//            bc_pos = vect_add(bc_pos, lc_shift);
//            bc_shifts.push_back(bc_pos);
//
//            lc_shift = scale_vect(sim_cell.lat_vect[0], a);
//            bc_pos = vect_add(pos, lc_shift);
//            lc_shift = scale_vect(sim_cell.lat_vect[2], b);
//            bc_pos = vect_add(bc_pos, lc_shift);
//            bc_shifts.push_back(bc_pos);
//
//            lc_shift = scale_vect(sim_cell.lat_vect[1], a);
//            bc_pos = vect_add(pos, lc_shift);
//            lc_shift = scale_vect(sim_cell.lat_vect[2], b);
//            bc_pos = vect_add(bc_pos, lc_shift);
//            bc_shifts.push_back(bc_pos);
//
//            for (int c : dir) {
//                lc_shift = scale_vect(sim_cell.lat_vect[0], a);
//                bc_pos = vect_add(pos, lc_shift);
//                lc_shift = scale_vect(sim_cell.lat_vect[1], b);
//                bc_pos = vect_add(bc_pos, lc_shift);
//                lc_shift = scale_vect(sim_cell.lat_vect[2], c);
//                bc_pos = vect_add(bc_pos, lc_shift);
//                bc_shifts.push_back(bc_pos);
//            }
//        }
//    }
//    for (vector<float> check : bc_shifts) {
//        if (pos_comp(check, check_vect)) {
//            bc_test = true;
//            break;
//        }
//    }
//    return bc_test;
//}
//
//void Algo1::fill_SMG(vector<vector<int>>& neigh_ind_list) {
//    vector<float> new_pos{ 0.0, 0.0, 0.0 };
//    vector<float> self_site{ 0.0, 0.0, 0.0 };
//    vector<int> sites;
//    vector<vector<int>> groups;
//    vector<vector<vector<int>>> motifs;
//    for (int atom = 0; atom < sim_cell.numb_atoms; atom++) {
//        for (vector<vector<vector<float>>> clust : session.spin_motif_list) {
//            for (vector<vector<float>> motif : clust) {
//                for (vector<float> shift : motif) {
//                    if (pos_comp(shift, self_site)) { sites.push_back(atom); }
//                    else {
//                        new_pos = vect_add(pos_list[atom], shift);
//                        for (int neigh : neigh_ind_list[atom]) {
//                            if (bc_check(pos_list[neigh], new_pos)) {
//                                sites.push_back(neigh);
//                            }
//                        }
//                    }
//                }
//                groups.push_back(sites);
//                sites.clear();
//            }
//            motifs.push_back(groups);
//            groups.clear();
//        }
//        spin_motif_groups.push_back(motifs);
//        motifs.clear();
//    }
//}
//
//void Algo1::fill_CMG(vector<vector<int>>& neigh_ind_list) {
//    vector<float> new_pos{ 0.0, 0.0, 0.0 };
//    vector<float> self_site{ 0.0, 0.0, 0.0 };
//    vector<int> sites;
//    vector<vector<int>> groups;
//    vector<vector<vector<int>>> motifs;
//    for (int atom = 0; atom < sim_cell.numb_atoms; atom++) {
//        for (vector<vector<vector<float>>> clust : session.chem_motif_list) {
//            for (vector<vector<float>> motif : clust) {
//                for (vector<float> shift : motif) {
//                    if (pos_comp(shift, self_site)) { sites.push_back(atom); }
//                    else {
//                        new_pos = vect_add(pos_list[atom], shift);
//                        for (int neigh : neigh_ind_list[atom]) {
//                            if (bc_check(pos_list[neigh], new_pos)) {
//                                sites.push_back(neigh);
//                            }
//                        }
//                    }
//                }
//                groups.push_back(sites);
//                sites.clear();
//            }
//            motifs.push_back(groups);
//            groups.clear();
//        }
//        chem_motif_groups.push_back(motifs);
//        motifs.clear();
//    }
//}
//
//void Algo1::print_state(string contcar_name, int temp) {
//    vector<int> perm;
//    vector<int> temp_spec;
//    vector<float> temp_spin;
//    vector<vector<int>> temp_allowed;
//    vector<vector<float>> temp_pos;
//    temp_spec.assign(chem_list.begin(), chem_list.end());
//    temp_spin.assign(spin_list.begin(), spin_list.end());
//    temp_pos = pos_list;// insert(temp_pos.end(), pos_list.begin(), pos_list.end());
//    perm = vect_permut(temp_spin);
//    for (int i = 0; i < sim_cell.numb_atoms; i++) { temp_allowed.push_back(sim_cell.atom_list[i].allowed_species); }
//    sort_vect(temp_spin, perm);
//    sort_vect(temp_pos, perm);
//    sort_vect(temp_spec, perm);
//    sort_vect(temp_allowed, perm);
//    perm = vect_permut(temp_spec);
//    sort_vect(temp_spin, perm);
//    sort_vect(temp_pos, perm);
//    sort_vect(temp_spec, perm);
//    sort_vect(temp_allowed, perm);
//    ofstream OUT_file;
//    string file_name = contcar_name + "_" + to_string(temp);
//    if (temp == -1) { file_name = contcar_name; }
//    OUT_file.open(file_name);
//    if (OUT_file.is_open()) {
//        for (string spec : session.species_str) { OUT_file << " " << spec; }
//        OUT_file << "\n 1 \n";
//        for (int i = 0; i < 3; i++) {
//            vector<float>vect = sim_cell.unit_lat_vect[i];
//            for (int j = 0; j < 3; j++) { OUT_file << vect[j] * session.shape[i] << " "; }
//            OUT_file << "\n";
//        }
//        for (string spec : session.species_str) { OUT_file << " " << spec; }
//        OUT_file << "\n";
//        for (int i = 0; i < sim_cell.species_numbs.size(); i++) { OUT_file << sim_cell.species_numbs[i] << " "; }
//        OUT_file << "\nCartesian\n";
//        for (int i = 0; i < sim_cell.numb_atoms; i++) {
//            OUT_file << temp_pos[i][0] << " " << temp_pos[i][1] << " " << temp_pos[i][2] << " ";
//            for (int j = 0; j < temp_allowed[i].size(); j++) {
//                OUT_file << session.species_str[temp_allowed[i][j]] << " ";
//            }
//            OUT_file << " # " << temp_spec[i] << " " << temp_spin[i] << "\n";
//        }
//    }
//    OUT_file.close();
//}
//
//void Algo1::spin_move(int site, int pass, float temp, float new_spin) {
//    // Flip Spin
//    float old_spin = spin_list[site];
//    float keep_rand;
//    float keep_prob;
//    spin_list[site] = new_spin;
//    spin_flip += new_spin - old_spin;
//    spin_list[site] = new_spin;
//    e_flip = eval_spin_flip(site, old_spin);
//    spin_flip = new_spin - old_spin;
//    if (e_flip < 0) { flip_count += 1; }
//    else {
//        keep_rand = unif(rng);
//        keep_prob = exp(-1 / (Kb * temp) * (e_flip));
//        if (keep_rand < keep_prob) { flip_count2 += 1; }
//        else { spin_list[site] = old_spin; e_flip = 0.0; spin_flip = 0.0; }
//    }
//    // Record the enrg and spin changes
//    init_enrg += e_flip;
//    init_spin += spin_flip;
//}
//
//void Algo1::run() {
//    // declare variables
//    bool same_spin;
//    bool same_atom;
//    int attempts = 0;
//    int rand_site = 1;
//    float rand_spin = 0.0;
//    float temp1 = session.start_temp;
//    float temp2 = session.end_temp;
//    float temp_step = session.temp_step;
//    vector<vector<float>> spin_states = session.spin_states;
//    int ta_passes = session.ta_passes;
//    int eq_passes = session.eq_passes;
//    int numb_atoms = sim_cell.numb_atoms;
//    int numb_neighbors;
//    vector<vector<int>> neigh_ind_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
//    vector<int> spin_atoms; // atomic species that can have spin
//
//    cout << "Making atom list and neighbor index list\n";
//    // Make atom_list more acessable for species and spin and neighbors
//    for (int i = 0; i < sim_cell.numb_atoms; i++) {
//        chem_list.push_back(sim_cell.atom_list[i].getSpecies());
//        spin_list.push_back(sim_cell.atom_list[i].getSpin());
//        pos_list.push_back({ sim_cell.atom_list[i].pos[0], sim_cell.atom_list[i].pos[1], sim_cell.atom_list[i].pos[2] });
//        numb_neighbors = sim_cell.atom_list[i].getNumbNeighbors(i, sim_cell);
//        for (int j = 0; j < numb_neighbors; j++) {
//            neigh_ind_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
//        }
//    }
//
//    cout << "Making rule maps\n";
//    // Make rule_maps for easy lookup and initialize rule_count_list
//    float ind = 0;
//    size_t rule_key;
//    vector<uint32_t> rule_info;
//    map<size_t, vector<double>>::iterator rule_itr;
//    for (Rule rule : session.chem_rule_list) {
//        rule_info.push_back(rule.GetType());
//        rule_info.push_back(rule.clust_ind);
//        for (int i = 0; i < rule.deco.size(); i++) { rule_info.push_back(rule.deco[i]); }
//        rule_key = cust_hash(rule_info);
//        rule_itr = rule_map_chem.find(rule_key);
//        if (rule_itr == rule_map_chem.end()) {
//            vector<double> enrg_ind_pair = { rule.GetEnrgCont(), ind };
//            rule_map_chem.insert(pair<size_t, vector<double>>(rule_key, enrg_ind_pair));
//            lat_rule_count_list.push_back(0);
//            site_rule_count_list.push_back(0);
//            ind += 1;
//        }
//        rule_info.clear();
//    }
//    for (Rule rule : session.spin_rule_list) {
//        rule_info.push_back(rule.GetType());
//        rule_info.push_back(rule.clust_ind);
//        for (int i = 0; i < rule.deco.size(); i++) { rule_info.push_back(rule.deco[i]); }
//        rule_key = cust_hash(rule_info);
//        rule_map_spin.insert(pair<size_t, double>(rule_key, rule.GetEnrgCont()));
//        for (int atom : rule.deco) {
//            if (find(spin_atoms.begin(), spin_atoms.end(), atom) == spin_atoms.end()) { spin_atoms.push_back(atom); } // initialize spin_atoms
//        }
//        rule_info.clear();
//    }
//
//    // Fill motif group lists
//    cout << "Making chem motif group lists\n";
//    fill_CMG(neigh_ind_list);
//    cout << "Making spin motif group lists\n";
//    fill_SMG(neigh_ind_list);
//
//    // Create seperate output file to avoid race condition
//    string file_name = "OUTPUT";
//    string contcar_name = "CONTCAR";
//    bool file_exists = true;
//    while (file_exists == true) {
//        const char* c_file = file_name.c_str();
//        int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
//        if (fd < 0) {
//            // file exists or otherwise uncreatable
//            outfile_count += 1;
//            file_name = "OUTPUT" + to_string(outfile_count);
//            contcar_name = "CONTCAR" + to_string(outfile_count);
//        }
//        else {
//            file_exists = false;
//            close(fd);
//        }
//    }
//    const char* c_file = file_name.c_str();
//    ofstream Output;
//    Output.open(c_file);
//    // Output energy and spin for convergence test
//    ofstream Output_converge;
//    if (session.do_conv_output) { Output_converge.open("OUTPUT_CONVERG"); }
//
//    // Begin MC
//    init_enrg = eval_lat();
//    cout << "Initial total energy is " << init_enrg / (numb_atoms - numb_vac) << " per atom\n";
//    double init_spin_enrg = eval_lat_spin();
//    cout << "Initial spin energy is " << init_spin_enrg / (numb_atoms - numb_vac) << " per atom\n";
//    for (int site = 0; site < numb_atoms; site++) {
//        if (find(spin_atoms.begin(), spin_atoms.end(), chem_list[site]) != spin_atoms.end()) {
//            init_spin += spin_list[site];
//        }
//    }
//    cout << "Initial spin is " << init_spin / (numb_atoms - numb_vac) << " per atom\n";
//    Output << "Using Algo1 for atom swap and spin flip";
//    Output << "\nPhase: " << sim_cell.phase_init;
//    Output << "\nComposition: ";
//    for (int i = 0; i < sim_cell.species_numbs.size(); i++) { Output << sim_cell.species_numbs[i] << " "; }
//    Output << "\nMC passes: " << passes;
//    Output << "\nInitial total energy per atom: ";
//    Output << init_enrg / (numb_atoms - numb_vac);
//    Output << "\nInitial spin energy per atom: ";
//    Output << init_spin_enrg / (numb_atoms - numb_vac);
//    Output << "\nInitial spin per atom: ";
//    Output << init_spin / (numb_atoms - numb_vac);
//    Output << "\ntemp, enrg, mag, var_e, var_spin, Cp, Xmag, flip_count, flip_count2:" << endl;
//
//    // Start MC loop
//    cout << "Entering main loop\n";
//    int temp_count = 0; // index of temperature
//    float inc_dir = 1; //temp increment direction
//    if (signbit(temp2 - temp1) == 1) { inc_dir = -1; }
//    for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_step * inc_dir) {
//        flip_count = 0;
//        flip_count2 = 0;
//        count_avg = lat_rule_count_list;
//        fill(count_avg.begin(), count_avg.end(), 0.0);
//        for (int pass = 0; pass < passes; pass++) {
//            for (int site = 0; site < numb_atoms; site++) {
//                e_flip = 0.0;
//                spin_flip = 0.0;
//                float new_spin = 0.0;
//                float new_spin1 = 0.0;
//                float new_spin2 = 0.0;
//                int state = METHOD_1;
//                while (state != DONE) {
//                    switch (state) {
//                        //-----------------------------------------------------------
//                    case 0:
//                        state = DONE;
//                        break;
//                        //-----------------------------------------------------------
//                    case METHOD_1:
//                        if (find(spin_atoms.begin(), spin_atoms.end(), chem_list[site]) == spin_atoms.end()) {
//                            state = DONE;
//                            break;
//                        }
//                        if (spin_states[chem_list[site]].size() <= 1) {
//                            state = DONE;
//                            break;
//                        }
//                        same_spin = true;
//                        attempts = 0;
//                        while (same_spin == true and attempts < 20) {
//                            rand_spin = unif(rng);
//                            for (int it_spin_state = 0; it_spin_state < spin_states[chem_list[site]].size(); it_spin_state++) {
//                                if (rand_spin > float(it_spin_state) / float(spin_states[chem_list[site]].size())) {
//                                    new_spin = spin_states[chem_list[site]][it_spin_state];
//                                }
//                            }
//                            if (new_spin != spin_list[site]) { same_spin = false; }
//                            attempts += 1;
//                        }
//                        if (attempts >= 20) {
//                            state = DONE;
//                            break;
//                        }
//                        spin_move(site, pass, temp, new_spin);
//                        state = DONE;
//                        //                        num_method1 += 1;
//                        break;
//                        //-----------------------------------------------------------
//                    }
//                }
//                if (pass >= eq_passes) {
//                    rs_E.Push(init_enrg);
//                    rs_M.Push(init_spin);
//                }
//            }
//            //            if (session.do_conv_output) {
//            //                Output_converge << eval_lat() << "; " << init_enrg << ", " << e_flip << "; " << init_spin << ", " << spin_flip << endl;
//            //            }
//            if (session.do_conv_output) {
//                Output_converge << init_enrg << " " << init_spin << endl;
//            }
//        }
//        //        cout << num_method0 << ", " << num_method1 << ", " << num_method2 << ", " << num_method3 << endl;
//        double scale = 1.0 / (numb_atoms * ta_passes);
//        double e_avg = rs_E.Mean() / (numb_atoms - numb_vac);
//        double spin_avg = rs_M.Mean() / (numb_atoms - numb_vac);
//        for (int i = 0; i < count_avg.size(); i++) { count_avg[i] *= scale; }
//        double var_e = rs_E.Variance();
//        double var_spin = rs_M.Variance();
//        double Cp = var_e / (Kb * pow(temp, 2));
//        double Xmag = var_spin / (Kb * pow(temp, 2));
//        Output << temp << ", "
//            << e_avg << ", "
//            << spin_avg << ", "
//            << var_e << ", "
//            << var_spin << ", "
//            << Cp << ", "
//            << Xmag << ", "
//            << flip_count << ", "
//            << flip_count2 << endl;
//        rs_E.Clear();
//        rs_M.Clear();
//        if (session.write_contcars == true) {
//            print_state(contcar_name, temp_count);
//        }
//        temp_count += 1;
//    }
//    print_state("CONTCAR_FINAL", -1);
//    Output.close();
//    cout << " MC Finished\n";
//}
