#include "stdafx.h"
#include "Structures.h"
#include "Functions.h"

using namespace std;

void Limited_World(configuration & config, const directedgraph & G)
{
	vector<directedgraph> scenarios;
	if (config.failure_type == 1)
		scenarios =  Generate_Scenarios_Tight(G, config.nr_scenarios);
	else if (config.failure_type == 2)
		scenarios =  Generate_Scenarios_Vertex_Tight(G, config.nr_scenarios);
	else
	{
		cout << "ERROR; Failure type not recognized" << endl;
		cin.get();
	}
	ofstream output;
	output.open("Scenarios.txt");
	for (int i = 0; i < config.nr_scenarios; i++)
	{
		output << "Scenario " << i << endl;
		for (int j = 0; j < scenarios[i].arcs.size(); j++)
		{
			output << "(" << scenarios[i].arcs[j].startvertex << "," << scenarios[i].arcs[j].endvertex << ")" << endl;
		}
		output << endl;
	}
	
	// Build the subsets for Subset Recourse.
	vector<cycle_arcs>subsets = Relevant_Subsets(G, config);

	// Set the weights of the subset based on the scenarios.
	Subset_Set_Weights_LW(subsets, scenarios, config);
	// Solve the Subset Recourse
	time_t now;
	pre_test_result results_SR = Subset_MIP(subsets, G, config, 9999999, now);
	results_SR.objective_value = results_SR.objective_value / config.nr_scenarios;

	// Solve for the Pre-Test
	pre_test_result results_PT = Pre_Test_Scen(G, config.cyclelength, config.chainlength, results_SR.tested_arcs.size(), config.nr_scenarios, config.time_limit, config.scen_gen, config.failure_type, scenarios);
	cout << "No fatal bugs at least" << endl;
	
	Output_LW(results_SR, results_PT, config);
	cin.get();
}

void Subset_Set_Weights_LW(vector<cycle_arcs>& subsets, const vector<directedgraph>& Scenarios, const configuration & config)
{
	for (int i = 0; i < subsets.size(); i++)
	{
		subsets[i].weight = 0;
		//for (int j = 0; j < subsets[i].vertices.size(); j++)
		//	cout << subsets[i].vertices[j] << "\t";
		//cout << endl;
		for (int j = 0; j < Scenarios.size(); j++)
		{
			directedgraph G_sub = Subset_Graph_LW(subsets[i], Scenarios[j]);
			//for (int k = 0; k < G_sub.arcs.size(); k++)
			//	cout << "(" << G_sub.arcs[k].startvertex << "," << G_sub.arcs[k].endvertex << ")" << endl;
			vector<cycle_arcs> cycles =  Find_Cycles(G_sub, config);
			time_t now;
			pre_test_result subset_scen_result = Subset_MIP_Silent(cycles, G_sub, config);
			//cout << "Scen " <<j <<" = " << subset_scen_result.objective_value << endl;
			subsets[i].weight += subset_scen_result.objective_value;
		}
		//cout << "Weight = " << subsets[i].weight << endl;
		//cin.get();
	}
}

pre_test_result Subset_MIP_Silent(const vector<cycle_arcs>& subsets, const directedgraph & G, const configuration & config)
{
	IloEnv env;
	IloModel MIP(env);
	IloNumVarArray subsetvar(env, subsets.size(), subsets.size(), ILOINT);
	for (int i = 0; i < subsets.size(); i++)
		subsetvar[i] = IloNumVar(env, 0, 1, ILOINT);

	IloRangeArray vertex_cons = Build_Vertex_Constraint_Subset_MIP(env, G, subsetvar, subsets);
	MIP.add(vertex_cons);

	IloObjective obj = IloMaximize(env);
	for (int i = 0; i < subsets.size(); i++)
	{
		obj.setLinearCoef(subsetvar[i], subsets[i].weight);
	}
	MIP.add(obj);

	IloCplex Subset_MIP(MIP);
	Subset_MIP.setOut(env.getNullStream());

	Subset_MIP.solve();

	pre_test_result results;
	results.objective_value = Subset_MIP.getObjValue();
	results.computation_time = 0;
	for (int i = 0; i < subsets.size(); i++)
	{
		if (Subset_MIP.getValue(subsetvar[i]) > 0.99)
		{
			for (int j = 0; j < subsets[i].arcs.size(); j++)
				results.tested_arcs.push_back(G.arcs[subsets[i].arcs[j]]);
		}
	}

	env.end();
	return results;
}

directedgraph Subset_Graph_LW(const cycle_arcs & subset, const directedgraph & scenario)
{
	directedgraph G_sub;
	G_sub.size = subset.vertices.size();
	G_sub.nr_pairs = G_sub.size;
	vector<int> vertex_switch(scenario.nr_pairs); // For each vertex in the original graph, we save in which position it occurs, in the list of vertices in the subset.
	for (int i = 0; i < subset.vertices.size(); i++)
	{
		vertex_switch[subset.vertices[i]] = i;
		G_sub.pairs.push_back(scenario.pairs[subset.vertices[i]]);
	}
	for (int i = 0, nr_arc = 0; i < subset.arcs.size(); i++)
	{
		// For each arc in the subset, check whether it also exists in the scenario. If so, add it to the graph.
		for (int j = 0; j < scenario.arcs.size(); j++)
		{
			if (subset.arcs[i] == scenario.arcs[j].arcnumber)
			{
				G_sub.arcs.push_back(scenario.arcs[j]);
				//cout << G_sub.arcs[nr_arc].startvertex << "," << G_sub.arcs[nr_arc].endvertex << endl;
				G_sub.arcs[nr_arc].startvertex = vertex_switch[G_sub.arcs[nr_arc].startvertex];
				G_sub.arcs[nr_arc].endvertex = vertex_switch[G_sub.arcs[nr_arc].endvertex];
				nr_arc++;
			}
		}
	}
	sort(G_sub.arcs.begin(), G_sub.arcs.end(), arcsort);
	Find_First_arc(G_sub);
	return G_sub;
}

pre_test_result Pre_Test_Scen(directedgraph G, int chainlength, int cyclelength, int max_tests, int nr_scen, int time_limit, int scen_gen, int failure_type, vector<directedgraph> Scenarios)
{
	directedgraph Tested_Graph = G;

	time_t start_time;
	time(&start_time);

	IloEnv env;
	IloModel model(env);
	cout << "Generating Variables" << endl;
	vector<vector<vector<IloNumVarArray>>> Cyclevar(nr_scen); // First position is scenario, second index is the Graph Copy, third the position in the Graph, fourth the individual arcs.
	vector<vector<vector<vector<int>>>> Cyclevar_arc_link(nr_scen); // A vector to link the variables to the original arc. Cyclevar_arc_link[i][j][k][l] = m, means that this variable corresponds to the m-th arc in the original arc list.
	for (int i = 0; i < nr_scen; i++)
	{
		cycle_variables cvars = Generate_Cycle_Var(env, Scenarios[i], cyclelength, i);
		Cyclevar[i] = cvars.Cyclevariable;
		Cyclevar_arc_link[i] = cvars.Link_Cyclevar_Arc;
	}
	cout << "Cyclevar Generated" << endl;
	vector<vector<IloNumVarArray>> Chainvar(nr_scen); // First index is the scenario, second index is the position in the graph, third the individual arc.
	vector<vector<vector<int>>> Chainvar_arc_link(nr_scen); // A vector to link the variables to the original arc.
	for (int i = 0; i < nr_scen; i++)
	{
		chain_variables cvars = Generate_Chain_Var(env, Scenarios[i], chainlength, i);
		Chainvar[i] = cvars.Chainvar;
		Chainvar_arc_link[i] = cvars.Link_Chainvar_Arc;
	}
	cout << "Chainvar Generated" << endl;
	cout << "Variables Generated" << endl;
	// Create the Objective Function
	IloObjective obj = IloMaximize(env);
	for (int scen = 0; scen < nr_scen; scen++)
	{
		for (int i = 0; i < G.nr_pairs - 1; i++)
		{
			for (int j = 0; j < cyclelength; j++)
			{
				for (int k = 0; k < Cyclevar[scen][i][j].getSize(); k++)
				{
					obj.setLinearCoef(Cyclevar[scen][i][j][k], G.arcs[Cyclevar_arc_link[scen][i][j][k]].weight);
				}
			}
		}
		for (int i = 0; i < chainlength; i++)
		{
			for (int j = 0; j < Chainvar[scen][i].getSize(); j++)
				obj.setLinearCoef(Chainvar[scen][i][j], G.arcs[Chainvar_arc_link[scen][i][j]].weight);
		}
	}
	model.add(obj);

	// Create Constraints per scenario.
	vector<IloRangeArray> vertex_inflow_cons(nr_scen);
	vector<vector<vector<IloRangeArray>>> vertex_flow_cons(nr_scen);
	vector<vector<IloRangeArray>> vertex_chain_flow_cons(nr_scen);
	vector<IloRangeArray> NDD_Constraint(nr_scen);
	for (int scen = 0; scen < nr_scen; scen++)
	{
		// Max one incoming arc per vertex.
		vertex_inflow_cons[scen] = Build_Vertex_Constraint(env, model, G, Cyclevar[scen], Cyclevar_arc_link[scen], Chainvar[scen], Chainvar_arc_link[scen]);

		// If there is an arc arriving in position i, there should be an outgoing arc in position i+1 in the same copy (excluding the origin vertex in that copy).
		vertex_flow_cons[scen] = Build_Vertex_Flow_Constraint(env, model, G, Cyclevar[scen], Cyclevar_arc_link[scen], cyclelength);

		// If there is an outgoing arc_chain in position i, there should be an incoming arc_chain in position [i-1].
		vertex_chain_flow_cons[scen] = Build_Vertex_Flow_Chain_Constraint(env, model, G, Chainvar[scen], Chainvar_arc_link[scen], chainlength);

		// At most one outgoing arc for each NDD
		NDD_Constraint[scen] = Build_NDD_Constraint(env, model, G, Chainvar[scen], Chainvar_arc_link[scen]);
	}
	// Create Testing Variables and Constraints.
	IloNumVarArray Testvar = Generate_Testvar(env, G);
	vector<IloRangeArray> test_constraint = Build_Test_Constraint(env, model, G, Testvar, Cyclevar, Cyclevar_arc_link, Chainvar, Chainvar_arc_link, nr_scen);
	IloRange Max_Test_Constraint = Build_Max_Test_Constraint(env, model, Testvar, max_tests);

	IloCplex CPLEX(model);
	CPLEX.setParam(IloCplex::TiLim, time_limit);
	CPLEX.solve();

	pre_test_result results;
	results.objective_value = CPLEX.getObjValue() / nr_scen;
	cout << results.objective_value << endl;

	time_t current_time;
	time(&current_time);
	results.computation_time = difftime(current_time, start_time);

	for (int i = 0; i < G.arcs.size(); i++)
	{
		if (CPLEX.getValue(Testvar[i]) > 0.99)
		{
			results.tested_arcs.push_back(G.arcs[i]);
		}
	}


	return results;
}

void Output_LW(const pre_test_result & results_SR, const pre_test_result & results_PT, const configuration & config)
{
	ofstream output;
	output.open(config.testvar_output);
	output << "Nr_Pairs = " << config.nr_pairs << endl;
	output << "Nr_NDD = " << config.nr_ndd << endl << endl;

	output << "Results Pre-Test" << endl;
	output << "Objective Value = " << results_PT.objective_value << endl;
	output << "Computation Time = " << results_PT.computation_time << endl;
	output << "Nr Scenarios = " << config.nr_scenarios << endl;
	output << "Nr Test = " << config.max_test << endl;

	output << "Results Subset Recourse" << endl;
	output << "Objective Value = " << results_SR.objective_value << endl;
	output << "Computation Time = " << results_SR.computation_time << endl;
	output << "Nr Scenarios = " << config.nr_scenarios << endl;
	output << "Nr Test = " << config.max_test << endl;
}
