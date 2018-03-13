#pragma once

#include "stdafx.h"
#include "Structures.h"

// Input Function
configuration Readconfig(char configname[], int * error);

// Functions to Generate a new KEP Graph
directedgraph graph_generation(configuration & config);
directedgraph simplegraph(configuration & config);
directedgraph saidman(configuration & config);
directedgraph readgraph(configuration & config, string filename);
directedgraph readgraph2(configuration & config, string filename);

// Solvers for various Models
void Choose_Solver(configuration & config, directedgraph & G);
void pre_test_main(configuration & config, directedgraph G);
void Subset_Recourse(configuration & config, directedgraph G);
pre_test_result Pre_Test(directedgraph G, int chainlength, int cyclelength, int max_tests, int nr_scen, int time_limit, int scen_gen, int failure_type, const configuration & config); // Pre_test outputs a new graph, using only the arcs chosen in the optimization.
matching_result Hybrid_PIEF(directedgraph G, int chainlength, int cyclength, configuration & config);
matching_result PICEF(directedgraph G, configuration & config);


// Pre-Testing Optimization Functions
vector<directedgraph> Generate_Scenarios(const directedgraph &G, int nr_scen);
vector<directedgraph> Generate_Scenarios_Tight(const directedgraph & G, int nr_scen);
vector<directedgraph> Generate_Scenarios_Vertex_Tight(const directedgraph & G, int nr_scen);
cycle_variables Generate_Cycle_Var(IloEnv &env, const directedgraph & G, int cyclelength, int nr_scen);
chain_variables Generate_Chain_Var(IloEnv &env, directedgraph G, int chainlength, int nr_scen);
IloNumVarArray Generate_Testvar(IloEnv &env, directedgraph G);
vector<IloRangeArray> Build_Test_Constraint(IloEnv &env, IloModel model, directedgraph G, const IloNumVarArray & Testvar, const vector<vector<vector<IloNumVarArray>>> & Cyclevar, const vector<vector<vector<vector<int>>>> & cycle_link, const vector<vector<IloNumVarArray>> & Chainvar, const vector<vector<vector<int>>> & Chainvar_arc_link, int nr_scen);
IloRange Build_Max_Test_Constraint(IloEnv & env, IloModel model, const IloNumVarArray & Testvar, int max_test);
void Output_Pre_Test(const pre_test_result & results, const configuration & config);

// HPIEF extra functions.
// Extract Solution
vector<cycle_arcs> Cycle_Solution(const directedgraph &G, const vector<vector<IloNumVarArray>> & cycle_var, const vector<vector<vector<int>>> & cycle_link, IloCplex HPIEF_CPLEX);
vector<chain> Chain_Solution(const directedgraph &G, const vector<IloNumVarArray> & chain_var, const vector<vector<int>> & chain_link, IloCplex HPIEF_CPLEX);
// Variables
cycle_variables Generate_Cycle_Var(IloEnv &env, directedgraph G, int cyclelength);
vector<vector<vector<int>>> cycle_preproces(directedgraph G, const vector<vector<directedarc>> & Acopies, int cyclelength);
vector<vector<int>> distance_calc(const directedgraph & G, const vector<vector<directedarc>> & Acopies,  int cyclelength);
chain_variables Generate_Chain_Var(IloEnv &env, directedgraph G, int chainlength);
vector<vector<int>> chain_preproces(directedgraph G, int chainlength);
// Constraints
IloRangeArray Build_Vertex_Constraint(IloEnv & env, IloModel &model, directedgraph G, vector<vector<IloNumVarArray>> cycle_var, vector<vector<vector<int>>> cycle_link, vector<IloNumVarArray> chain_var, vector<vector<int>> chain_link);
IloRangeArray Build_NDD_Constraint(IloEnv & env, IloModel &model, directedgraph G, vector<IloNumVarArray> chain_var, vector<vector<int>> chain_link);
vector<vector<IloRangeArray>> Build_Vertex_Flow_Constraint(IloEnv & env, IloModel &model, directedgraph G, vector<vector<IloNumVarArray>> cycle_var, vector<vector<vector<int>>> cycle_link, int cyclelength);
vector<IloRangeArray> Build_Vertex_Flow_Chain_Constraint(IloEnv & env, IloModel &model, directedgraph G, vector<IloNumVarArray> chain_var, vector<vector<int>> chain_link, int chainlength);
// Misc
vector<vector<directedarc>> DP_Copy(directedgraph G);
bool arcsort(directedarc a, directedarc b);

// Functions for Subset Recourse
pre_test_result Subset_MIP(const vector<cycle_arcs> & subsets,const directedgraph & G, const configuration & config, double remaining_time, const time_t & start_time);
IloRangeArray Build_Vertex_Constraint_Subset_MIP(IloEnv & env, const directedgraph & G, IloNumVarArray & subsetvar, const vector<cycle_arcs> & subsets);
vector<cycle_arcs> Relevant_Subsets(const directedgraph & G, const configuration & config);
void Relevant_Subsets_Recursion(vector<cycle_arcs> & accepted, cycle_arcs current, queue<cycle_arcs> candidates, int max_subset_size);
void Subset_Arcs(vector<cycle_arcs> & subsets, const directedgraph & G, const configuration & config);
void Subset_Set_Weights_Arc(vector<cycle_arcs> & subsets, const directedgraph & G, const configuration & config);
void Subset_Set_Weights_Vertex(vector<cycle_arcs> & subsets, const directedgraph & G, const configuration & config);
void Subset_Set_Weights_Arc_verbose(vector<cycle_arcs> & subsets, const directedgraph & G, const configuration & config);
void Subset_Set_Weights_Vertex_Verbose(vector<cycle_arcs>& subsets, const directedgraph & G, const configuration & config);
float Subset_Set_Weights_Arc_Root(const directedgraph & G, const configuration & config);
float Subset_Set_Weights_Vertex_Root(const directedgraph & G, const configuration & config);
float Subset_Set_Weights_Arc_Root_verbose(const directedgraph & G, const configuration & config);
float Subset_Set_Weights_Vertex_Root_Verbose(const directedgraph & G, const configuration & config);
float Subset_Set_Weights_Arc_Recursion(IloEnv & env, const directedgraph & G, IloCplex & CPLEX, IloNumVarArray & arcvar, vector<bool> fixed, float Solution, vector<int> arc_solutions, bool succes_fix, int depth);
float Subset_Set_Weights_Vertex_Recursion(IloEnv &env, const directedgraph & G, IloCplex & CPLEX, vector<int> vertex_use, IloRangeArray & vertex_con, vector<bool> fixed, float Solution, bool succes_fix, int depth);
float Subset_Set_Weights_Arc_Recursion_verbose(IloEnv & env, const directedgraph & G, IloCplex & CPLEX, IloNumVarArray & arcvar, vector<bool> fixed, float Solution, vector<int> arc_solutions, bool succes_fix, int depth, deque<bool> & suc_fail);
float Subset_Set_Weights_Vertex_Recursion_Verbose(IloEnv &env, const directedgraph & G, IloCplex & CPLEX, vector<int> vertex_use, IloRangeArray & vertex_con, vector<bool> fixed, float Solution, bool succes_fix, int depth, deque<bool> & suc_fail);
directedgraph Subset_Graph(const cycle_arcs & subsets, const directedgraph & G);
void Vertex_Combine(cycle_arcs & current, const cycle_arcs & candidate);
bool Vertex_Subset_Dominate(const cycle_arcs & parent, const cycle_arcs & candidate);
bool Vertex_Subset_Equal(const cycle_arcs & parent, const cycle_arcs & candidate);
bool Vertex_Disjoint(const cycle_arcs & parent, const cycle_arcs & candidate);
int Vertex_Sum(const cycle_arcs & parent, const cycle_arcs & candidate);
IloRangeArray Build_Vertex_Constraint_SSWR(IloEnv & env, const directedgraph & G, IloNumVarArray & cyclevar, const vector<cycle_arcs> & cycles);
IloRangeArray Build_Cycle_Constraint_SSWR(IloEnv &env, const directedgraph & G, IloNumVarArray & cyclevar, IloNumVarArray & arcvar, const vector<cycle_arcs> & cycles);
void Output_Subset_Recourse(const pre_test_result & results, const configuration & config);

// Functions common to cycle formulations (including robust models)
vector<cycle_arcs> Find_Cycles(directedgraph G, const configuration & config);
vector<int> Find_First_arc(directedgraph & G);
vector<vector<int>> List_Outgoing_Arcs(const directedgraph & G);

// PICEF Functions
// Extract Solution
// Chain_Solution from HPIEF used.
// Variables
// Generate_Chain_Var from HPIEF used.

// Constraints
IloRangeArray Build_Vertex_Constraint_PICEF(IloEnv & env, IloModel &model, directedgraph G, vector<IloNumVarArray> chain_var, vector<vector<int>> chain_link);
// NDD Constraint from HPIEF.
// Chain Flow Constraint from HPIEF.

// Testing Functions
float Calc_Expected_Transplants(configuration & config);
vector<cycle_arcs> Split_SCC(const directedgraph & Tested_G);

// Limited World Comparisons (There is a limited number of scenarios.)
void Limited_World(configuration & config, const directedgraph & G);
void Subset_Set_Weights_LW(vector<cycle_arcs> & subsets, const vector<directedgraph> & Scenarios, const configuration & config);
pre_test_result Subset_MIP_Silent(const vector<cycle_arcs> & subsets, const directedgraph & G, const configuration & config);
directedgraph Subset_Graph_LW(const cycle_arcs & subset, const directedgraph & scenario);
pre_test_result Pre_Test_Scen(directedgraph G, int chainlength, int cyclelength, int max_tests, int nr_scen, int time_limit, int scen_gen, int failure_type, vector<directedgraph> scenarios);
void Output_LW(const pre_test_result & results_SR, const pre_test_result & results_PT, const configuration & config);

// Unlimited Cycle Robust Approach
void UnlimitedCycle_Main(configuration & config, const directedgraph & G);
arc_variables_uc Generate_Arc_Var_UC(IloEnv & env, directedgraph G, int scenario_number);
IloNumVarArray Generate_Test_Var_UC(IloEnv & env, directedgraph G);
IloRangeArray Build_Inflow_Constraint_UC(IloEnv & env, IloModel &model, const directedgraph & G, IloNumVarArray cycle_var, vector<int> cycle_link);
IloRangeArray Build_Flow_Constraint_UC(IloEnv & env, IloModel &model, const directedgraph & G, IloNumVarArray cycle_var, vector<int> cycle_link);
IloRangeArray Build_Test_Constraint_UC(IloEnv & env, IloModel & model, const directedgraph & G, IloNumVarArray cycle_var, vector<int> cycle_link, IloNumVarArray test_var);

// Mathheuristics
void Arc_Heuristic(configuration & config, const directedgraph & G);
pre_test_result Arc_Use_IP(configuration & config, const directedgraph & G, const vector<int> & arc_use, vector<bool> & arc_tested);
IloNumVarArray Generate_Test_Var_Arc_Use_IP(IloEnv &env, directedgraph G);
cycle_variables Generate_Cycle_Var_Arc_Use_IP(IloEnv &env, directedgraph G, int cyclelength, int nr_scen);
chain_variables Generate_Chain_Var_Arc_Use_IP(IloEnv & env, directedgraph G, int chainlength, int nr_scen);
IloRangeArray Arc_Constraints_Arc_Use_IP(IloEnv &env, IloModel &model, const directedgraph G, IloNumVarArray Test_Var, vector<vector<IloNumVarArray>> cycle_variables, vector<vector<vector<int>>> cycle_link, vector<IloNumVarArray> chain_var, vector<vector<int>> chain_link, vector<int> arc_use);
vector<int> Identify_Arcs_To_Delete(const directedgraph & G, vector<int> arc_use, vector<bool> testvar_vector, int nr_scen);
vector<int> Identify_Arcs_To_Delete2(const directedgraph & G, vector<int> arc_use, vector<bool> testvar_vector, int nr_scen);
void delete_from_G(directedgraph & G, const vector<int> & to_delete);
void delete_from_scen(vector<directedgraph> & scen, vector<vector<int>> & arcnumbers, const vector<int> & to_delete);

void Cycle_Heuristic(configuration & config, directedgraph G);
pre_test_result Cycle_Use_IP(configuration & config, const directedgraph & G, const vector<int> & arcu_use, const vector<vector<vector<cycle_arcs>>> cycles_used, vector<bool> & arc_tested);

matching_result Hybrid_PIEF_Heur(directedgraph G, int chainlength, int cyclelength, configuration & config);
vector<directedgraph> Generate_Scenarios_Vertex_Tight_Heur(const directedgraph & G, int nr_scen);