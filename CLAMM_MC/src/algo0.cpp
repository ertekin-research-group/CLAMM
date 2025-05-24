#include "algo0.h"

Algo0::Algo0(void) {}

Algo0::Algo0(Session& _session, SimCell& _sim_cell) {
    session = _session;
    sim_cell = _sim_cell;
    passes = session.ta_passes + session.eq_passes;
    for (int i = 0; i < session.sro_def.size(); i++) { numb_sro_atoms += session.atom_numbs[session.sro_def[i]]; }
    if (session.ta_passes < 1) {
        cout << "_______________________________________________________________________________" << endl;
        cout << "Error: Algo0 has been given 0 thermal average passes" << endl;
        cout << "This is probably not what you want..." << endl;
        cout << "_______________________________________________________________________________" << endl;
    }
    // setup rng for random spin choice and acceptance probability
    rng.seed(ss);
    unif = std::uniform_real_distribution<double>(0.0, 1.0);
    rand_atom = std::uniform_int_distribution<int>(0, sim_cell.numb_atoms - 1);
    rand_method = std::uniform_int_distribution<int>(0, passes - 1);
}

size_t Algo0::cust_hash(vector<uint32_t>& vect) {
    std::size_t seed = vect.size();
    for (auto x : vect) {
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = (x >> 16) ^ x;
        seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

double Algo0::eval_site_chem(int site) {
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
                lat_rule_count_list[round(rule_itr->second[1])] += 1.0 / group.size();
            }
            rule_info.clear();
        }
    }
    return enrg;
}

double Algo0::eval_site_spin(int site) {
    double enrg = 0;
    map<size_t, double>::iterator rule_itr;
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
            enrg += (rule_itr != rule_map_spin.end()) ? (rule_itr->second * spin_prod / group.size()) : 0.0;
            rule_info.clear();
        }
    }
    return enrg;
}

double Algo0::eval_spin_flip(int site, float old_spin) {
    double enrg = 0;
    map<size_t, double>::iterator rule_itr;
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
            enrg += (rule_itr != rule_map_spin.end()) ? (rule_itr->second * spin_prod) : 0.0;
            rule_info.clear();
        }
    }
    return (enrg * spin_list[site] - enrg * old_spin);
}

double Algo0::eval_atom_flip(int site) {
    map<size_t, vector<double>>::iterator rule_itr_chem;
    map<size_t, double>::iterator rule_itr_spin;
    vector<uint32_t> rule_info;
    size_t rule_key;
    double enrg = 0.0;
    fill(site_rule_count_list.begin(), site_rule_count_list.end(), 0);
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
            rule_itr_chem = rule_map_chem.find(rule_key);
            if (rule_itr_chem != rule_map_chem.end()) {
                enrg += rule_itr_chem->second[0];
                site_rule_count_list[round(rule_itr_chem->second[1])] += 1.0;
            }
            rule_info.clear();
        }
    }
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
            rule_itr_spin = rule_map_spin.find(rule_key);
            enrg += (rule_itr_spin != rule_map_spin.end()) ? (rule_itr_spin->second * spin_prod) : 0.0;
            rule_info.clear();
        }
    }
    return enrg;
}

double Algo0::eval_lat() {
    double enrg = 0;
    for (int site = 0; site < sim_cell.numb_atoms; site++) {
        enrg += eval_site_chem(site);
        enrg += eval_site_spin(site);
    }
    return enrg + session.intercept * sim_cell.numb_atoms;
}

double Algo0::eval_lat_spin() {
    double enrg = 0;
    for (int site = 0; site < sim_cell.numb_atoms; site++) {
        enrg += eval_site_spin(site);
    }
    return enrg;
}

bool Algo0::bc_check(vector<float> check_vect, vector<float>& pos) {
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

void Algo0::fill_SMG(vector<vector<int>>& neigh_ind_list) {
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

void Algo0::fill_CMG(vector<vector<int>>& neigh_ind_list) {
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

void Algo0::fill_SROMG(vector<vector<int>>& neigh_ind_list) {
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

void Algo0::print_state(string contcar_name, int temp) {
    vector<int> perm;
    vector<int> temp_spec;
    vector<float> temp_spin;
    vector<vector<int>> temp_allowed;
    vector<vector<float>> temp_pos;
    temp_spec.assign(chem_list.begin(), chem_list.end());
    temp_spin.assign(spin_list.begin(), spin_list.end());
    temp_pos = pos_list;// insert(temp_pos.end(), pos_list.begin(), pos_list.end());
    perm = vect_permut(temp_spin);
    for (int i = 0; i < sim_cell.numb_atoms; i++) { temp_allowed.push_back(sim_cell.atom_list[i].allowed_species); }
    sort_vect(temp_spin, perm);
    sort_vect(temp_pos, perm);
    sort_vect(temp_spec, perm);
    sort_vect(temp_allowed, perm);
    perm = vect_permut(temp_spec);
    sort_vect(temp_spin, perm);
    sort_vect(temp_pos, perm);
    sort_vect(temp_spec, perm);
    sort_vect(temp_allowed, perm);
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
            for (int j = 0; j < temp_allowed[i].size(); j++) {
                OUT_file << session.species_str[temp_allowed[i][j]] << " ";
            }
            OUT_file << " # " << temp_spec[i] << "\n";// << " " << temp_spin[i] << "\n";
        }
    }
    OUT_file.close();
}

double Algo0::eval_site_prob(int site) {
    double prob = 0;
    map<size_t, double>::iterator rule_itr;
    vector<uint32_t> rule_info;
    size_t rule_key;
    for (int i = 0; i < sro_motif_groups[site].size(); i++) {
        vector<vector<int>> motif = sro_motif_groups[site][i];
        for (int j = 0; j < motif.size(); j++) {
            rule_info.push_back(0); // chem ind
            rule_info.push_back(i); // clust_ind
            vector<int> group = motif[j];
            for (int k : group) {
                rule_info.push_back(chem_list[k]); // sites ind
            }
            rule_key = cust_hash(rule_info);
            rule_itr = rule_map_sro.find(rule_key);
            prob += (rule_itr != rule_map_sro.end()) ? (rule_itr->second / group.size()) : 0.0;
            rule_info.clear();
        }
    }
    return prob;
}

void Algo0::run() {

    // declare variables
    bool same_spin;
    bool same_atom;
    int attempts = 0;
    int rand_site = 1;
    float rand_spin = 0.0;
    float temp1 = session.start_temp;
    float temp2 = session.end_temp;
    float temp_step = session.temp_step;
    vector<vector<float>> spin_states = session.spin_states;
    int ta_passes = session.ta_passes;
    int eq_passes = session.eq_passes;
    int numb_atoms = sim_cell.numb_atoms;
    int numb_neighbors;
    vector<vector<int>> neigh_ind_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
    vector<int> spin_atoms; // atomic species that can have spin
    cout << "Making atom list and neighbor index list\n";
    // Make atom_list more acessable for species and spin and neighbors
    for (int i = 0; i < sim_cell.numb_atoms; i++) {
        chem_list.push_back(sim_cell.atom_list[i].getSpecies());
        spin_list.push_back(sim_cell.atom_list[i].getSpin());
        pos_list.push_back({ sim_cell.atom_list[i].pos[0], sim_cell.atom_list[i].pos[1], sim_cell.atom_list[i].pos[2] });
        numb_neighbors = sim_cell.atom_list[i].getNumbNeighbors(i, sim_cell);
        for (int j = 0; j < numb_neighbors; j++) {
            neigh_ind_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
        }
    }

    cout << "Making rule maps\n";
    // Make rule_maps for easy lookup and initialize rule_count_list
    float ind = 0;
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
            vector<double> enrg_ind_pair = { rule.GetEnrgCont(), ind };
            rule_map_chem.insert(pair<size_t, vector<double>>(rule_key, enrg_ind_pair));
            lat_rule_count_list.push_back(0);
            site_rule_count_list.push_back(0);
            ind += 1;
        }
        rule_info.clear();
    }
    for (Rule rule : session.spin_rule_list) {
        rule_info.push_back(rule.GetType());
        rule_info.push_back(rule.clust_ind);
        for (int i = 0; i < rule.deco.size(); i++) { rule_info.push_back(rule.deco[i]); }
        rule_key = cust_hash(rule_info);
        rule_map_spin.insert(pair<size_t, double>(rule_key, rule.GetEnrgCont()));
        for (int atom : rule.deco) {
            if (find(spin_atoms.begin(), spin_atoms.end(), atom) == spin_atoms.end()) { spin_atoms.push_back(atom); } // initialize spin_atoms
        }
        rule_info.clear();
    }
    for (Rule rule : session.sro_rule_list) {
        rule_info.push_back(rule.GetType());
        rule_info.push_back(rule.clust_ind);
        for (int i = 0; i < rule.deco.size(); i++) { rule_info.push_back(rule.deco[i]); }
        rule_key = cust_hash(rule_info);
        rule_map_sro.insert(pair<size_t, double>(rule_key, rule.GetEnrgCont()));
        rule_info.clear();
    }


    // Fill motif group lists
    cout << "Making chem motif group lists\n";
    fill_CMG(neigh_ind_list);
    cout << "Making spin motif group lists\n";
    fill_SMG(neigh_ind_list);
    cout << "Making sro motif group lists\n";
    fill_SROMG(neigh_ind_list);

    // Create seperate output file to avoid race condition
    string file_name = "OUTPUT";
    string sro_file_name = "OUTPUT_SRO";
    string contcar_name = "CONTCAR";
    bool file_exists = true;
    while (file_exists == true) {
        const char* c_file = file_name.c_str();
        int fd = openWithFcntl(c_file);
        if (fd < 0) {
            // file exists or otherwise uncreatable
            outfile_count += 1;
            file_name = "OUTPUT" + to_string(outfile_count);
            sro_file_name = "OUTPUT_SRO" + to_string(outfile_count);
            contcar_name = "CONTCAR" + to_string(outfile_count);
        }
        else {
            file_exists = false;
            //close(fd);
        }
    }
    const char* c_file = file_name.c_str();
    const char* sro_file = sro_file_name.c_str();
    ofstream Output;
    ofstream SROout;
    Output.open(c_file);
    SROout.open(sro_file);
    // Output energy and spin for convergence test
    ofstream Output_converge;
    if (session.do_conv_output) { Output_converge.open("OUTPUT_CONVERG"); }

    int old_site_chem;
    int new_site_chem;
    int old_rand_chem;
    int new_rand_chem;
    double prob_site_flip;
    double flip;
    double prob_initial = 0.0;
    double prob = 0.0;
    double sro_initial = 0.0;
    double keep_rand;
    double keep_prob;
    int numb_A = session.atom_numbs[session.sro_def[0]];
    int numb_B = session.atom_numbs[session.sro_def[1]];
    int sro_atom_numbs = 0;
    double prob_target = (1.0 - session.sro_target) * 2.0 * numb_A / numb_atoms * numb_B / numb_atoms;
    cout << "sro_def size: " << session.sro_def.size() << "\n";
    cout <<"sro_motif_groups size: " << sro_motif_groups.size() << "\n";
    cout <<"sro_rule_list size: " << session.sro_rule_list.size() << "\n";
    cout <<"rule_map_sro size: " << rule_map_sro.size() << "\n";
    cout <<"rule_map_chem size: " << rule_map_chem.size() << "\n";
    for (int i = 0; i < session.sro_def.size(); i++) { sro_atom_numbs += session.atom_numbs[session.sro_def[i]]; }
    cout << "sro_atom_numbs: " << sro_atom_numbs << "\n";
    for (int i = 0; i < numb_atoms; i++) {
        prob_initial += eval_site_prob(i);
    }
    prob_initial /= numb_atoms;
    sro_initial = 1.0 - prob_initial / (2.0*numb_A/numb_atoms*numb_B/numb_atoms);
    cout << "Prob Initial: " << prob_initial << " sro_initial: " << sro_initial << "\n";
    Output << "SRO Initial: " << sro_initial << endl;  
    Output << "Temp, SRO, Flips, Flips2\n";    
    // Start MC loop
    cout << "Entering main loop\n";
    int temp_count = 0; // index of temperature
    float inc_dir = 1; //temp increment direction
    double sro_avg = 0.0;
    if (signbit(temp2 - temp1) == 1) { inc_dir = -1; }
    for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_step * inc_dir) {
        flip_count = 0;
        flip_count2 = 0;
        count_avg = lat_rule_count_list;
        fill(count_avg.begin(), count_avg.end(), 0.0);
        for (int pass = 0; pass < passes; pass++) {
            for (int site = 0; site < numb_atoms; site++) {
                e_flip = 0.0;
                double prob_site_old = 0.0;
                double prob_site_new = 0.0;
                double prob_sum_old = 0.0;
                double prob_sum_new = 0.0;
                int state = 0;
                int method_index = rand_method(rng);
                state = METHOD_2;
                while (state != DONE) {
                    switch (state) {
                    case METHOD_0:
                        state = DONE;
                        break;
                    case METHOD_2:
                        if (sim_cell.atom_list[site].allowed_species.size() <= 1) {
                            state = METHOD_0;
                            break;
                        }
                        same_atom = true;
                        attempts = 0;
                        while (same_atom == true && attempts < 100) {
                            rand_site = rand_atom(rng);
                            if (rand_site != site) {
                                if (find(sim_cell.atom_list[rand_site].allowed_species.begin(), sim_cell.atom_list[rand_site].allowed_species.end(), chem_list[site]) != sim_cell.atom_list[rand_site].allowed_species.end() && find(sim_cell.atom_list[site].allowed_species.begin(), sim_cell.atom_list[site].allowed_species.end(), chem_list[rand_site]) != sim_cell.atom_list[site].allowed_species.end()) {
                                    if (chem_list[site] != chem_list[rand_site]) { same_atom = false; }
                                }
                            }
                            attempts += 1;
                        }
                        if (attempts >= 100) {
                            state = METHOD_0;
                            break;
                        }
                        prob_sum_old = eval_site_prob(site) + eval_site_prob(rand_site);
                        prob_site_old = prob_sum_old / 2 * numb_sro_atoms/numb_atoms;
                        old_site_chem = chem_list[site];
                        old_rand_chem = chem_list[rand_site];
                        chem_list[site] = old_rand_chem;
                        chem_list[rand_site] = old_site_chem;
                        prob_sum_new = eval_site_prob(site) + eval_site_prob(rand_site);
                        prob_site_new = prob_sum_new / 2 * numb_sro_atoms/numb_atoms;
                        prob_site_new = (eval_site_prob(site) + eval_site_prob(rand_site))/ 2 * numb_sro_atoms/numb_atoms;
                        flip = (abs(prob_site_new - prob_target) - abs(prob_site_old - prob_target)) * 2 / numb_atoms;
                        if (flip < 0) { prob_initial += (prob_sum_new - prob_sum_old)/numb_atoms; }
                        else {
                            keep_rand = unif(rng);
                            keep_prob = exp(-1 / (Kb * temp) * (flip * 200));
                            cout << "keep_prob: " << keep_prob << " keep_rand: " << keep_rand << "\n";
                            if (keep_rand < keep_prob) { prob_initial += (prob_sum_new - prob_sum_old) / numb_atoms; }
                            else {
                                chem_list[site] = old_site_chem;
                                chem_list[rand_site] = old_rand_chem;
                            }
                        }
                        state = DONE;
                        break;
                    }
                }
				rs_SRO.Push(1 - prob_initial/(2.0 * numb_A/numb_atoms * numb_B/numb_atoms));
                rs_PROB.Push(prob_initial);
            }
            if (session.do_conv_output) {
                Output_converge << init_enrg << " " << init_spin << endl;
            }
        }
        double scale = 1.0 / (numb_atoms * ta_passes);
        double e_avg = rs_SRO.Mean();
        double prob_avg = rs_PROB.Mean();
        for (int i = 0; i < count_avg.size(); i++) { count_avg[i] *= scale; }
        Output << temp << ", "
            << e_avg << ", "
            << prob_avg << ", "
            << flip_count << ", "
            << flip_count2 << endl;
        rs_SRO.Clear();
        rs_PROB.Clear();
        SROout << temp << " ";
        for (float x : count_avg) { SROout << x << " "; }
        SROout << endl;
        if (session.write_contcars == true) {
            print_state(contcar_name, temp_count);
        }
        temp_count += 1;
    }
    
    prob_initial = 0.0;
    for (int i = 0; i < numb_atoms; i++) {
        prob_initial += eval_site_prob(i);
    }
    prob_initial /= numb_atoms;
    sro_initial = 1.0 - prob_initial / (2.0 * numb_A / numb_atoms * numb_B / numb_atoms);
    Output << "SRO Final: " << sro_initial << " Prob Final: " << prob_initial << endl;
    print_state(contcar_name, -1);
    Output.close();
    cout << " MC Finished\n";

}