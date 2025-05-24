#pragma once
#ifndef run_sim_h
#define run_sim_h

#include "session.h"
#include "sim_cell.h"
#include "algo_m2.h"
#include "algo0.h"
#include "algo1.h"
#include "algo2.h"
#include "algo3.h"
#include "algo4.h"

class Sim {
public:
	SimCell sim_cell;
	Session session;

	Sim(void);
	Sim(Session& _session, SimCell& _sim_cell);
	void start();
};

#endif
