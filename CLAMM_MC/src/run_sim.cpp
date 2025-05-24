#include "run_sim.h"

Sim::Sim(void) {}

Sim::Sim(Session& _session, SimCell& _sim_cell) {
	sim_cell = _sim_cell;
	session = _session;
}

void Sim::start() {
	cout << "Using algo " << session.algo << "\n";
	if (session.algo == -1) { return; }
	else if (session.algo == -2) { Algo_m2 algo(session, sim_cell); algo.run(); }
	else if (session.algo == 0) { Algo0 algo(session, sim_cell); algo.run(); }
	else if (session.algo == 1) { Algo1 algo(session, sim_cell); algo.run(); }
	else if (session.algo == 2) { Algo2 algo(session, sim_cell); algo.run(); }
    else if (session.algo == 3) { Algo3 algo(session, sim_cell); algo.run(); }
    else if (session.algo == 4) { Algo4 algo(session, sim_cell); algo.run(); }
}
