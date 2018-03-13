
#pragma once
#include "stdafx.h"

using namespace std;

#ifndef Struct_H
#define Struct_H

struct directedarc {
	int startvertex;
	int endvertex;
	int arcnumber; // A number to uniquely identify each arc. This number should correspond to the position in the initial vector of all arcs.
	float weight;
	float failprob;	
};

struct ndd_donor {
	int bloodtype;
	float failprob;
};

struct patient_donor_pair {
	int donor_bloodtype; // 0 for O, 1 for A, 2 for B, 3 for AB.
	int patient_bloodtype;
	int PRA_level; // 0 for low, 1 for medium, 2 for high.
	float failprob;
};

struct cycle_arcs {
	vector<int> arcs;	// A vector, giving the arc numbers (reference to position in G.arcs)
	vector<int> vertices;	// A vector giving the vertices involved in the cycle in order.
	float weight = 0;		// The weight of the arcs involved in the optimization problem.
};
struct cycle {
	vector<int> vertices;	// A vector giving the vertices involved in the cycle in order.
	float weight = 0;			// The weight of the arcs involved in the optimization problem.
};

struct cycle_arc_reference { // The vector of integers here are references to the individual arcs that make up the cycele
	vector<int> arc_nr;
	float weight;
};

struct chain {
	vector<int> vertices;	// A vector giving the vertices involved in the chain in order, the first is the NDD.
	float weight;			// The weight of the arcs involved in the optimization problem.
};

struct directedgraph { // A directed KEP graph, the first nr_pairs vertices correspond to the donor-patient pairs, the remaining to NDDs.
	int size;
	int nr_pairs;
	int nr_ndd;
	vector<patient_donor_pair> pairs;
	vector<ndd_donor> ndds;
	vector<directedarc> arcs; // This list should be ordered, first based on startvertex, next on endvertex.
	vector<int> first_vertex_arc; // For each vertex, this vector gives the arc number where that vertex first appears as a starvertex.
};

struct matching_result {
	float objective;
	vector<cycle_arcs> cycles;
	vector<chain> chains;
};

struct pre_test_result {
	vector<directedarc> tested_arcs;
	float objective_value;
	float computation_time;
};

struct cycle_variables {
	vector<vector<IloNumVarArray>> Cyclevariable; // First position is graph copy, second is arc position, third refers to arc. 
	vector<vector<vector<int>>> Link_Cyclevar_Arc;
};

struct chain_variables {
	vector<IloNumVarArray> Chainvar;
	vector<vector<int>> Link_Chainvar_Arc;
};

struct arc_variables_uc {
	IloNumVarArray Arcvariable;
	vector<int> Link_Arcvariable;
};

struct configuration {
	int nr_pairs;
	int nr_ndd;
	int nr_scenarios;
	int input_data;
	int cyclelength;
	int chainlength;
	int max_test;
	float failprob;
	float arcprob;
	int subset_size;
	string inputfile;
	string graph_output;
	string testvar_output;
	string solution_input; // The file to read in a previous (pre-test) solution.
	int seed;
	int time_limit;
	int memory_limit;
	int solver;
	int calc_expected_type;
	int failure_type;
	int scen_gen;
	int LW; // Flag for limited world (set number of scenarios).
};



#endif 
