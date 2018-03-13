#include "stdafx.h"
#include "Structures.h"
#include "Functions.h"

void Arc_Heuristic(configuration & config, const directedgraph & G_original)
{
	directedgraph G = G_original;
	vector<directedgraph> scenarios;
	if (config.failure_type == 2)
	{
		scenarios = Generate_Scenarios_Vertex_Tight(G, config.nr_scenarios);
	}
	else
	{
		scenarios = Generate_Scenarios_Tight(G, config.nr_scenarios);
	}
	
	// Vector which saves how often each arc is used in solutions to the scenarios.
	vector<int> arc_use(G.arcs.size()); 
	// Vector which saves the relation between arcs in the scenarios to the Graph G. Arcnumbers must be changed in the scenarios to avoid things breaking in HPIEF.
	vector<vector<int>> arcnumbers(config.nr_scenarios);

	for (int i = 0; i < config.nr_scenarios; i++)
	{
		// Re-number the arcs, or things break in the HPIEF function. We do save the correct arcnumbers to relate results to the arcs in G.
		arcnumbers[i].resize(scenarios[i].arcs.size());
		for (int j = 0; j < scenarios[i].arcs.size(); j++)
		{
			arcnumbers[i][j] = scenarios[i].arcs[j].arcnumber;
			scenarios[i].arcs[j].arcnumber = j;
		}
		
		matching_result new_result = Hybrid_PIEF(scenarios[i], config.chainlength, config.cyclelength, config);

		// Save which arcs were used in this scenario.
		for (int j = 0; j < new_result.cycles.size(); j++)
		{
			for (int k = 0; k < new_result.cycles[j].arcs.size(); k++)
			{
				arc_use[arcnumbers[i][new_result.cycles[j].arcs[k]]]++;
			}
		}
	}
	for (int i = 0; i < G.arcs.size(); i++)
	{
		cout << "(" << G.arcs[i].startvertex << "," << G.arcs[i].endvertex << ") \t" << arc_use[i] << endl;
	}

	// Set-up the IP to choose arcs based on how often they are used over the different scenarios.
	vector<bool> arc_tested;
	pre_test_result results = Arc_Use_IP(config, G, arc_use, arc_tested);
	
	vector<int> to_delete = Identify_Arcs_To_Delete(G, arc_use, arc_tested, config.nr_scenarios);

	while (to_delete.size() > 0)
	{
		for (int i = 0; i < G.arcs.size(); i++)
		{
			if(arc_tested[i] == 1)
			cout << "(" << G.arcs[i].startvertex << "," << G.arcs[i].endvertex << ") \t" << arc_tested[i] << endl;
		}
		cout << "Size to_delete = " << to_delete.size() << endl;
		delete_from_G(G, to_delete);
		cout << "G deletions finished" << endl;
		delete_from_scen(scenarios, arcnumbers, to_delete);
		cout << "Deletions Finished" << endl;
		arc_use = vector<int> (G.arcs.size(),0);
		for (int i = 0; i < config.nr_scenarios; i++)
		{
			matching_result new_result = Hybrid_PIEF(scenarios[i], config.chainlength, config.cyclelength, config);

			// Save which arcs were used in this scenario.
			for (int j = 0; j < new_result.cycles.size(); j++)
			{
				for (int k = 0; k < new_result.cycles[j].arcs.size(); k++)
				{
					arc_use[arcnumbers[i][new_result.cycles[j].arcs[k]]]++;
				}
			}
		}

		// Set-up the IP to choose arcs based on how often they are used over the different scenarios.
		pre_test_result results = Arc_Use_IP(config, G, arc_use, arc_tested);
		to_delete = Identify_Arcs_To_Delete(G, arc_use, arc_tested, config.nr_scenarios);
	}
	
	Output_Pre_Test(results, config);
	system("PAUSE");
}

void Cycle_Heuristic(configuration & config, directedgraph G)
{
	vector<directedgraph> scenarios;
	if (config.failure_type == 2)
	{
		scenarios = Generate_Scenarios_Vertex_Tight_Heur(G, config.nr_scenarios);
	}
	else
	{
		scenarios = Generate_Scenarios_Tight(G, config.nr_scenarios);
	}
	// Vector which saves the relation between arcs in the scenarios to the Graph G. Arcnumbers must be changed in the scenarios to avoid things breaking in HPIEF.
	vector<vector<int>> arcnumbers(config.nr_scenarios);
	// Re-number the arcs, or things break in the HPIEF function. We do save the correct arcnumbers to relate results to the arcs in G.
	for (int i = 0; i < config.nr_scenarios; i++)
	{
		arcnumbers[i].resize(scenarios[i].arcs.size());
		for (int j = 0; j < scenarios[i].arcs.size(); j++)
		{
			arcnumbers[i][j] = scenarios[i].arcs[j].arcnumber;
			scenarios[i].arcs[j].arcnumber = j;
		}
	}


	// Vector which saves how often each arc is used in solutions to the scenarios.
	vector<int> arc_use(G.arcs.size());
	vector<int> prev_arc_use;
	// Vector which saves which cycles are used in the scenarios. cycles_used[i][j][k] refers to the k-th cycle of length i + 2, that starts in vertex j that is used.
	// cycles_used[i][j][k].weight is how often this cycle is used in the scenarios.
	vector<vector<vector<cycle_arcs>>> cycles_used(config.cyclelength - 1);

	bool run_loop = 1;
	vector<bool> arc_tested;
	vector<bool> prev_arc_test;
	while (run_loop == 1)
	{
		run_loop = 0;
		// Reset the cycles_used vectors.
		arc_use = vector<int>(G.arcs.size(), 0);
		for (int i = 0; i < cycles_used.size(); i++)
		{
			cycles_used[i].clear();
			cycles_used[i].resize(G.nr_pairs);
		}
		for (int i = 0; i < config.nr_scenarios; i++) // Run the scenarios.
		{
			matching_result new_result = Hybrid_PIEF_Heur(scenarios[i], config.chainlength, config.cyclelength, config);
			// Save which cycles and arcs are used in this scenario.
			for (int j = 0; j < new_result.cycles.size(); j++)
			{
				// Save the arc use.
				for (int k = 0; k < new_result.cycles[j].arcs.size(); k++)
				{
					arc_use[arcnumbers[i][new_result.cycles[j].arcs[k]]]++;
				}
				// Check whether the cycle has already been found in an earlier scenario.
				bool exists = 0;
				for (int k = 0; k < cycles_used[new_result.cycles[j].arcs.size() - 2][new_result.cycles[j].vertices[0]].size(); k++)
				{
					if (new_result.cycles[j].vertices == cycles_used[new_result.cycles[j].arcs.size() - 2][new_result.cycles[j].vertices[0]][k].vertices)
					{
						exists = 1;
						cycles_used[new_result.cycles[j].arcs.size() - 2][new_result.cycles[j].vertices[0]][k].weight = cycles_used[new_result.cycles[j].arcs.size() - 2][new_result.cycles[j].vertices[0]][k].weight + new_result.cycles[j].weight;
						break;
					}
				}
				if (exists == 0)
				{
					// Change the arcnumbers to comform to original G.
					for (int k = 0; k < new_result.cycles[j].arcs.size(); k++)
					{
						new_result.cycles[j].arcs[k] = arcnumbers[i][new_result.cycles[j].arcs[k]];
					}
						cycles_used[new_result.cycles[j].arcs.size() - 2][new_result.cycles[j].vertices[0]].push_back(new_result.cycles[j]);
				}
					
				
			}
		}
		for (int i = 0; i < cycles_used.size(); i++) // Illustrative Output
		{
			for (int j = 0; j < cycles_used[i].size(); j++)
			{
				for (int k = 0; k < cycles_used[i][j].size(); k++)
				{
					for (int l = 0; l < cycles_used[i][j][k].vertices.size(); l++)
					{
						cout << cycles_used[i][j][k].vertices[l] << "\t";
					}
					cout << "\t = " << cycles_used[i][j][k].weight << endl;
				}
			}
		}

		pre_test_result result = Cycle_Use_IP(config, G, arc_use, cycles_used, arc_tested);
		Output_Pre_Test(result, config);
		for (int i = 0; i < arc_use.size(); i++)
		{
			if (arc_use[i] > 0)
				cout << "(" << G.arcs[i].startvertex << "," << G.arcs[i].endvertex << ") \t" << arc_use[i] << "\t" << arc_tested[i] << endl;
		}

		vector<int> to_delete = Identify_Arcs_To_Delete2(G, arc_use, arc_tested, config.nr_scenarios);
		if (to_delete.size() > 0)
		{
			run_loop = 1;
			delete_from_G(G, to_delete);
			delete_from_scen(scenarios, arcnumbers, to_delete);
		}
	}
	
}

pre_test_result Cycle_Use_IP(configuration & config, const directedgraph & G, const vector<int>& arc_use, const vector<vector<vector<cycle_arcs>>> cycles_used, vector<bool> & arc_tested)
{
	IloEnv env; 
	IloModel Model(env);

	int total_cycles = 0;
	for (int i = 0; i < cycles_used.size(); i++)
	{
		for (int j = 0; j < cycles_used[i].size(); j++)
		{
			total_cycles += cycles_used[i][j].size(); 
		}
	}

	// Create variables
	IloNumVarArray Test_Var = Generate_Test_Var_Arc_Use_IP(env, G);
	IloNumVarArray Cycle_Var = IloNumVarArray(env, total_cycles, 0, 1, ILOINT);

	IloObjective obj = IloMaximize(env);
	int nr_cycle = 0;
	int number_cycle_arcs = 0; // Summing the number of arcs per cycle over all cycles.
	for (int i = 0; i < cycles_used.size(); i++)
	{
		for (int j = 0; j < cycles_used[i].size(); j++)
		{
			for (int k = 0; k < cycles_used[i][j].size(); k++)
			{
				obj.setLinearCoef(Cycle_Var[nr_cycle], cycles_used[i][j][k].weight);
				nr_cycle++;
				number_cycle_arcs += cycles_used[i][j][k].arcs.size();
			}
		}
	}
	Model.add(obj);
	IloRangeArray constraints(env, number_cycle_arcs);
	int constr_nr = 0;
	int cycle_nr = 0;
	for (int i = 0; i < cycles_used.size(); i++)
	{
		for (int j = 0; j < cycles_used[i].size(); j++)
		{
			for (int k = 0; k < cycles_used[i][j].size(); k++)
			{
				for (int l = 0; l < cycles_used[i][j][k].arcs.size(); l++)
				{
					IloExpr expr(env);
					constraints[constr_nr] = IloRange(expr <= 0);
					constraints[constr_nr].setLinearCoef(Test_Var[cycles_used[i][j][k].arcs[l]], -1);
					constraints[constr_nr].setLinearCoef(Cycle_Var[cycle_nr], 1);
					constr_nr++;
				}
				cycle_nr++;
			}
		}
	}
	Model.add(constraints);

	// Max number of tests
	IloRange Max_Test_Constraint = Build_Max_Test_Constraint(env, Model, Test_Var, config.max_test);

	IloCplex CPLEX(Model);

	CPLEX.solve();
	CPLEX.exportModel("Cycle_Use.lp");
	CPLEX.writeSolution("Cycle_Use.sol");

	pre_test_result results;
	results.objective_value = CPLEX.getObjValue() / config.nr_scenarios;
	cout << results.objective_value << endl;

	arc_tested = vector<bool>(G.arcs.size(), 0);
	for (int i = 0; i < G.arcs.size(); i++)
	{
		if(CPLEX.getValue(Test_Var[i]) > 0.99)
		{
			results.tested_arcs.push_back(G.arcs[i]);
			arc_tested[i] = 1;
		}
		else
		{
			arc_tested[i] = 0;
		}
			
	}

	return results;
}

pre_test_result Arc_Use_IP(configuration & config, const directedgraph & original_G, const vector<int>& arc_use, vector<bool> & arc_tested)
{
	directedgraph G = original_G;
	
	IloEnv env;
	IloModel Arc_Use_IP(env);
	// Create all variables (pre-processing can be implemented here)
	IloNumVarArray Test_Var = Generate_Test_Var_Arc_Use_IP(env, G);

	vector<vector<IloNumVarArray>> Cyclevar; // First index is the Graph Copy, Second the position in the cycle, Third the individual arcs.
	vector<vector<vector<int>>> Cyclevar_arc_link; // A vector to link the variables to the original arc. Cyclevar_arc_link[i][j][k] = l, means that this variable corresponds to the l-th arc in the original arc list. 
	{
		cycle_variables cvars = Generate_Cycle_Var_Arc_Use_IP(env, G, config.cyclelength, config.nr_scenarios);
		Cyclevar = cvars.Cyclevariable;
		Cyclevar_arc_link = cvars.Link_Cyclevar_Arc;
	}
	vector<IloNumVarArray> Chainvar; // First index is the position in the chain, Second the individual arc.
	vector<vector<int>> Chainvar_arc_link; // A vector to link the variables to the original arc.
	{
		chain_variables cvars = Generate_Chain_Var_Arc_Use_IP(env, G, config.chainlength, config.nr_scenarios);
		Chainvar = cvars.Chainvar;
		Chainvar_arc_link = cvars.Link_Chainvar_Arc;
	}

	// Create the Objective Function
	IloObjective obj = IloMaximize(env);
	for (int i = 0; i < G.nr_pairs - 1; i++)
	{
		for (int j = 0; j < config.cyclelength; j++)
		{
			for (int k = 0; k < Cyclevar[i][j].getSize(); k++)
			{
				obj.setLinearCoef(Cyclevar[i][j][k], G.arcs[Cyclevar_arc_link[i][j][k]].weight);
			}
		}
	}
	for (int i = 0; i < config.chainlength; i++)
	{
		for (int j = 0; j < Chainvar[i].getSize(); j++)
			obj.setLinearCoef(Chainvar[i][j], G.arcs[Chainvar_arc_link[i][j]].weight);
	}
	Arc_Use_IP.add(obj);
	
	// Build the model, it is based on the Dickerson et al. (2016) HPIEF, but instead, the upper bound for an arc-variable is how often it is used in the scenarios. 
	IloRangeArray Arc_Constraint = Arc_Constraints_Arc_Use_IP(env, Arc_Use_IP, G, Test_Var, Cyclevar, Cyclevar_arc_link, Chainvar, Chainvar_arc_link, arc_use);

	// Max number of tests
	IloRange Max_Test_Constraint = Build_Max_Test_Constraint(env, Arc_Use_IP, Test_Var, config.max_test);

	// If there is an arc arriving in position i, there should be an outgoing arc in position i+1 in the same copy (excluding the origin vertex in that copy).
	vector<vector<IloRangeArray>> vertex_flow_cons = Build_Vertex_Flow_Constraint(env, Arc_Use_IP, G, Cyclevar, Cyclevar_arc_link, config.cyclelength);

	// If there is an outgoing arc_chain in position i, there should be an incoming arc_chain in position [i-1].
	vector<IloRangeArray> vertex_chain_flow_cons = Build_Vertex_Flow_Chain_Constraint(env, Arc_Use_IP, G, Chainvar, Chainvar_arc_link, config.chainlength);

	IloCplex CPLEX(Arc_Use_IP);
	CPLEX.exportModel("Arc_Use.lp");
	CPLEX.solve();
	CPLEX.writeSolution("Arc_Use.sol");

	pre_test_result results;
	results.objective_value = CPLEX.getObjValue() / config.nr_scenarios;
	cout << results.objective_value << endl;

	arc_tested = vector<bool>(G.arcs.size(), 0);
	for (int i = 0; i < G.arcs.size(); i++)
	{
		if (CPLEX.getValue(Test_Var[i]) > 0.99)
		{
			results.tested_arcs.push_back(G.arcs[i]);
			arc_tested[i] = 1;
		}
		else
			arc_tested[i] = 0;
	}

	return results;
}

IloNumVarArray Generate_Test_Var_Arc_Use_IP(IloEnv &env, directedgraph G)
{
	IloNumVarArray Test_Var(env);
	for (int i = 0; i < G.arcs.size(); i++)
	{
		ostringstream convert;
		convert << "test(" << G.arcs[i].startvertex << "," << G.arcs[i].endvertex << ")";
		string varname = convert.str();
		const char* vname = varname.c_str();
		Test_Var.add(IloNumVar(env, 0, 1, ILOINT, vname));
	}
	return Test_Var;
}

cycle_variables Generate_Cycle_Var_Arc_Use_IP(IloEnv &env, directedgraph G, int cyclelength, int nr_scen)
{
	// Note that we consistently work with nr_pairs - 1. For obvious reasons, the copy corresonding to the last pair can not contain any cycles.
	cycle_variables c;
	c.Cyclevariable.resize(G.nr_pairs - 1);
	c.Link_Cyclevar_Arc.resize(G.nr_pairs - 1);
	// Pre-Processing
	vector<vector<directedarc>> Acopies = DP_Copy(G);
	vector<vector<vector<int>>> copy_pos_arc_possible = cycle_preproces(G, Acopies, cyclelength);

	// Create all variables
	cout << "Variable Creation" << endl;
	for (int i = 0; i < G.nr_pairs - 1; i++)
	{
		c.Link_Cyclevar_Arc[i].resize(cyclelength);
		for (int j = 0; j < cyclelength; j++)
		{
			c.Cyclevariable[i].push_back(IloNumVarArray(env));
			for (int k = 0; k < Acopies[i].size(); k++)
			{
				if (copy_pos_arc_possible[i][j][k] == 1)
				{
					ostringstream convert;
					convert << "x(" << Acopies[i][k].startvertex << "," << Acopies[i][k].endvertex << "," << i << "," << j << ")";
					string varname = convert.str();
					const char* vname = varname.c_str();
					c.Cyclevariable[i][j].add(IloNumVar(env, 0, nr_scen, ILOINT, vname));
					c.Link_Cyclevar_Arc[i][j].push_back(Acopies[i][k].arcnumber);
				}
			}
		}
	}
	return c;
}

chain_variables Generate_Chain_Var_Arc_Use_IP(IloEnv & env, directedgraph G, int chainlength, int nr_scen)
{
	chain_variables c;
	c.Link_Chainvar_Arc.resize(G.nr_pairs - 1);

	// Pre-Processing. For each arc of G, the following function checks whether they can be in a particular position in any feasible solution.
	vector<vector<int>> arc_position_possible = chain_preproces(G, chainlength);
	// The first index is Position, the second is the individual arc.

	// Create all variables
	for (int i = 0; i < chainlength; i++)
	{
		c.Chainvar.push_back(IloNumVarArray(env));
		for (int j = 0; j < G.arcs.size(); j++)
		{
			if (arc_position_possible[i][j] == 1)
			{
				if (i == 0 && G.arcs[j].startvertex >= G.nr_pairs) // If this is true, the starvertex is an NDD and can be placed in the first position of a chain.
				{
					ostringstream convert;
					convert << "y(" << G.arcs[j].startvertex << "," << G.arcs[j].endvertex << "," << i << ")";
					string varname = convert.str();
					const char* vname = varname.c_str();
					c.Chainvar[i].add(IloNumVar(env, 0, nr_scen, ILOINT, vname));
					c.Link_Chainvar_Arc[i].push_back(G.arcs[j].arcnumber);
				}
				else if (i > 0 && G.arcs[j].startvertex < G.nr_pairs) // Any arc between pairs can be used in the following positions
				{
					ostringstream convert;
					convert << "y(" << G.arcs[j].startvertex << "," << G.arcs[j].endvertex << "," << i << ")";
					string varname = convert.str();
					const char* vname = varname.c_str();
					c.Chainvar[i].add(IloNumVar(env, 0, nr_scen, ILOINT, vname));
					c.Link_Chainvar_Arc[i].push_back(G.arcs[j].arcnumber);
				}
			}
		}
	}

	return c;
}

IloRangeArray Arc_Constraints_Arc_Use_IP(IloEnv & env, IloModel &model, const directedgraph G, IloNumVarArray Test_Var, vector<vector<IloNumVarArray>> cycle_variables, vector<vector<vector<int>>> cycle_link, vector<IloNumVarArray> chain_var, vector<vector<int>> chain_link, vector<int> arc_use)
{
	IloRangeArray Arc_Constraint(env, G.arcs.size());
	for (int i = 0; i < G.arcs.size(); i++)
	{
		IloExpr expr(env);
		Arc_Constraint[i] = IloRange(expr <= 0);
		Arc_Constraint[i].setUB(0);
	}

	for (int i = 0; i < G.arcs.size(); i++)
	{
		Arc_Constraint[i].setLinearCoef(Test_Var[i], -arc_use[i]);
	}
	for (int i = 0; i < cycle_variables.size(); i++)
	{
		for (int j = 0; j < cycle_variables[i].size(); j++)
		{
			for (int k = 0; k < cycle_variables[i][j].getSize(); k++)
			{
				Arc_Constraint[cycle_link[i][j][k]].setLinearCoef(cycle_variables[i][j][k], 1);
			}
		}
	}
	for (int i = 0; i < chain_var.size(); i++)
	{
		for (int j = 0; j < chain_var[i].getSize(); j++)
		{
			Arc_Constraint[chain_link[i][j]].setLinearCoef(chain_var[i][j], 1);
		}
	}
	
	model.add(Arc_Constraint);
	return Arc_Constraint;
}

vector<int> Identify_Arcs_To_Delete(const directedgraph & G, vector<int> arc_use, vector<bool> testvar_vector, int nr_scen)
{
	vector<int> to_delete;
	vector<bool> adjacent_arc_deleted(G.size,0); 
	cout << endl;
	int startvert = 0; int arc = 0;
	for (int vertex = 0; vertex < G.size; vertex++)
	{
		int incumbent = -1; int incumbent_use = nr_scen + 1;
		while (G.arcs[arc].startvertex == vertex)
		{
			if (testvar_vector[arc] == 0 && arc_use[arc] < incumbent_use && arc_use[arc] > 0 && adjacent_arc_deleted[G.arcs[arc].startvertex] == 0 && adjacent_arc_deleted[G.arcs[arc].endvertex] == 0)
			{
				incumbent = arc;
				incumbent_use = arc_use[arc];
			}
			arc++;
			if (arc >= G.arcs.size())
				break;
		}
		if (incumbent != -1)
		{
			to_delete.push_back(incumbent);
			adjacent_arc_deleted[G.arcs[incumbent].startvertex] = 1;
			adjacent_arc_deleted[G.arcs[incumbent].endvertex] = 1;
		}
	}
	return to_delete;
}

vector<int> Identify_Arcs_To_Delete2(const directedgraph & G, vector<int> arc_use, vector<bool> testvar_vector, int nr_scen)
{
	vector<int> to_delete;
	vector<bool> adjacent_arc_deleted(G.size, 0);
	vector<int> free_outdegree(G.size, 0);
	int max_outdegree = 0;
	int arc = 0;
	for (int vertex = 0; vertex < G.size; vertex++)
	{
		while (G.arcs[arc].startvertex == vertex && arc < G.arcs.size())
		{
			arc++;
			if(arc_use[arc] > 0)
				free_outdegree[vertex]++;
		}
		if (free_outdegree[vertex] > max_outdegree)
			max_outdegree = free_outdegree[vertex];
	}
	
	int startvert = 0; arc = 0;
	for (int vertex = 0; vertex < G.size; vertex++)
	{
		int incumbent = -1; int incumbent_use = nr_scen + 1;
		while (G.arcs[arc].startvertex == vertex)
		{
			if (testvar_vector[arc] == 0 && arc_use[arc] < incumbent_use && arc_use[arc] > 0 && adjacent_arc_deleted[G.arcs[arc].startvertex] == 0 && adjacent_arc_deleted[G.arcs[arc].endvertex] == 0 && free_outdegree[vertex] == max_outdegree)
			{
				incumbent = arc;
				incumbent_use = arc_use[arc];
			}
			arc++;
			if (arc >= G.arcs.size())
				break;
		}
		if (incumbent != -1)
		{
			to_delete.push_back(incumbent);
			adjacent_arc_deleted[G.arcs[incumbent].startvertex] = 1;
			adjacent_arc_deleted[G.arcs[incumbent].endvertex] = 1;
		}
	}
	cout << "Deleting" << endl;
	for (int i = 0; i < to_delete.size(); i++)
	{
		cout << "(" << G.arcs[to_delete[i]].startvertex << "," << G.arcs[to_delete[i]].endvertex << ") \t" << arc_use[to_delete[i]] << endl;
	}
	return to_delete;
}

void delete_from_G(directedgraph & G, const vector<int>& to_delete)
{
	cout << "Start delete_from_G" << endl;
	if (to_delete.size() == 0)
	{
		cout << "ERROR: Nothing to delete" << endl;
		cin.get();
	}
	cout << to_delete.size() << "\t" << G.arcs.size();
	// In this loop the arcs are actually deleted. Start from the back so we don't need to account for earlier deleted arcs when calculating the index.
	for (int i = to_delete.size() - 1; i >= 0; i--)
	{
		G.arcs.erase(G.arcs.begin() + to_delete[i]);
	}
	for (int i = 0; i < G.arcs.size(); i++)
	{
		G.arcs[i].arcnumber = i;
	}
	cout << "End delete_from_G" << endl;
}

void delete_from_scen(vector<directedgraph>& scenarios, vector<vector<int>> & arcnumbers, const vector<int>& to_delete)
{
	if (to_delete.size() == 0)
	{
		cout << "ERROR: Nothing to delete" << endl;
		cin.get();
	}

	for (int scen = 0; scen < scenarios.size(); scen++)
	{
		//cout << "Scenario " << scen << endl;
		// Delete the unneccesary entries in arcnumbers vector, and delete the arcs from scenarios[scen].
		int j = to_delete.size() -1;
		for (int i = scenarios[scen].arcs.size() -1 ; i >= 0; i--)
		{
			while (arcnumbers[scen][i] < to_delete[j])
				j--;
			if (arcnumbers[scen][i] == to_delete[j])
			{
				arcnumbers[scen].erase(arcnumbers[scen].begin() + i);
				scenarios[scen].arcs.erase(scenarios[scen].arcs.begin() + i);
			}
		}
		// Change the numbers in the arcnumbers vector.
		j = 0;
		for (int i = 0; i < arcnumbers[scen].size(); i++)
		{
			if(j < to_delete.size())
			{
				while (arcnumbers[scen][i] > to_delete[j] && j < to_delete.size())
				{
					j++;
					if (j == to_delete.size())
						break;
				}
			}
			
			arcnumbers[scen][i] = arcnumbers[scen][i] - j;
		}
		// Renumber the arcnumbers in the arcs.
		for (int i = 0; i < scenarios[scen].arcs.size(); i++)
			scenarios[scen].arcs[i].arcnumber = i;

		/*for (int i = 0; i < scenarios[scen].arcs.size(); i++)
		{
			cout << "(" << scenarios[scen].arcs[i].startvertex << "," << scenarios[scen].arcs[i].endvertex << ") \t" << arcnumbers[scen][i] << "\t" << scenarios[scen].arcs[i].arcnumber<< endl;
		}*/
	}
}

matching_result Hybrid_PIEF_Heur(directedgraph G, int chainlength, int cyclelength, configuration & config)
{
	// This function gets the Graph G.
	// First, copies of the graph are made for each patient-donor pair.
	// The copy corresponding to pair [i] (G[i]), only contains vertices [i] and higher.
	// In a particular copy G[i], cycles and chains can only start in vertex [i].
	// "Copy[0]" is the full graph and will be used for the NDDs.
	// We do some pre-processing, figuring out in which position an arc can be in each copy of the graph.
	// Then the HPIEF is constructed.


	// Create all variables (pre-processing can be implemented here)
	IloEnv env;
	IloModel HPIEF(env);
	vector<vector<IloNumVarArray>> Cyclevar; // First index is the Graph Copy, Second the position in the cycle, Third the individual arcs.
	vector<vector<vector<int>>> Cyclevar_arc_link; // A vector to link the variables to the original arc. Cyclevar_arc_link[i][j][k] = l, means that this variable corresponds to the l-th arc in the original arc list. 
	{
		cycle_variables cvars = Generate_Cycle_Var(env, G, cyclelength);
		Cyclevar = cvars.Cyclevariable;
		Cyclevar_arc_link = cvars.Link_Cyclevar_Arc;
	}
	vector<IloNumVarArray> Chainvar; // First index is the position in the chain, Second the individual arc.
	vector<vector<int>> Chainvar_arc_link; // A vector to link the variables to the original arc.
	{
		chain_variables cvars = Generate_Chain_Var(env, G, chainlength);
		Chainvar = cvars.Chainvar;
		Chainvar_arc_link = cvars.Link_Chainvar_Arc;
	}

	// Create the Objective Function
	IloObjective obj = IloMaximize(env);
	for (int i = 0; i < G.nr_pairs - 1; i++)
	{
		for (int j = 0; j < cyclelength; j++)
		{
			for (int k = 0; k < Cyclevar[i][j].getSize(); k++)
			{
				obj.setLinearCoef(Cyclevar[i][j][k], G.arcs[Cyclevar_arc_link[i][j][k]].weight-0.0001*j);
			}
		}
	}
	for (int i = 0; i < chainlength; i++)
	{
		for (int j = 0; j < Chainvar[i].getSize(); j++)
			obj.setLinearCoef(Chainvar[i][j], G.arcs[Chainvar_arc_link[i][j]].weight);
	}
	HPIEF.add(obj);


	// Build the model, references are to the Dickerson et al. 2016 paper.
	// Maximize \sum[arcs, positions, copies] weight[arc] * cycle-x [arc, position, copy] + \sum[arcs, position] weight[arc] * chain-x [arc, position] (7a)
	// Subject to
	// \forall[vertex] \sum[incoming arcs, position, copies] cycle-x[arc,position, copy] + \sum[incoming arcs, position] chain-x[arc, position] <= 1 (7b)
	// \forall[Pair-vertex, position 1..K-1] \sum[incoming arc] cycle-x[arc,position, copy] = \sum[outgoing arc] cycle-x[arc,position +1, copy] (1c)
	// \forall[NDD-Vertex] \sum[outgoing arcs] chain-x[arc, position 1] <= 1 (4c)
	// \forall[Pair-Vertex, position 1..K-1] \sum[incoming arc] chain-x[arc, position] = \sum[outgoing arc] chain-x[arc, position+1] (4d)
	// Integrality Constraints

	// Create Constraints
	// Max one incoming arc per vertex.
	IloRangeArray vertex_inflow_cons = Build_Vertex_Constraint(env, HPIEF, G, Cyclevar, Cyclevar_arc_link, Chainvar, Chainvar_arc_link);

	// If there is an arc arriving in position i, there should be an outgoing arc in position i+1 in the same copy (excluding the origin vertex in that copy).
	vector<vector<IloRangeArray>> vertex_flow_cons = Build_Vertex_Flow_Constraint(env, HPIEF, G, Cyclevar, Cyclevar_arc_link, cyclelength);

	// If there is an outgoing arc_chain in position i, there should be an incoming arc_chain in position [i-1].
	vector<IloRangeArray> vertex_chain_flow_cons = Build_Vertex_Flow_Chain_Constraint(env, HPIEF, G, Chainvar, Chainvar_arc_link, chainlength);

	// At most one outgoing arc for each NDD
	IloRangeArray NDD_Constraint = Build_NDD_Constraint(env, HPIEF, G, Chainvar, Chainvar_arc_link);


	IloCplex HPIEF_CPLEX(HPIEF);
	if (config.solver == 4)
	{
		HPIEF_CPLEX.setOut(env.getNullStream());
	}
	//HPIEF_CPLEX.exportModel("HPIEF.lp");
	HPIEF_CPLEX.solve();
	//HPIEF_CPLEX.writeSolution("HPIEF.sol");

	matching_result results;
	results.objective = HPIEF_CPLEX.getObjValue();

	results.cycles = Cycle_Solution(G, Cyclevar, Cyclevar_arc_link, HPIEF_CPLEX);
	for (int i = 0; i < results.cycles.size(); i++)
	{
		results.cycles[i].weight = results.cycles[i].vertices.size()-1;
	}
	/*cout << "Cycles" << endl;
	for (int i = 0; i < results.cycles.size(); i++)
	{
	for (int j = 0; j < results.cycles[i].vertices.size(); j++)
	{
	cout << results.cycles[i].vertices[j] << "\t";
	}
	cout << endl;
	}
	cout << "Cycles finished" << endl;*/
	results.chains = Chain_Solution(G, Chainvar, Chainvar_arc_link, HPIEF_CPLEX);
	/*cout << "Chains" << endl;
	for (int i = 0; i < results.chains.size(); i++)
	{
	for (int j = 0; j < results.chains[i].vertices.size(); j++)
	{
	cout << results.chains[i].vertices[j] << "\t";
	}
	cout << endl;
	}
	cout << "Chains finished" << endl;*/
	HPIEF.end();
	env.end();
	return results;
}

vector<directedgraph> Generate_Scenarios_Vertex_Tight_Heur(const directedgraph & G, int nr_scen)
{
	ofstream output;
	srand(time(NULL));
	vector<directedgraph> Scenarios(nr_scen);
	// Get the number of successes for each vertex and the number of scenarios.
	vector<int> succes_per_vertex(G.size);
	for (int i = 0; i < G.nr_pairs; i++)
	{
		float expected_succes = nr_scen*(1 - G.pairs[i].failprob);
		float natural;
		float remainder = modf(expected_succes, &natural);
		// Round it to an integer probabilistically.
		if (rand() % 1000 > remainder * 1000)
			remainder = 1;
		else
			remainder = 0;
		succes_per_vertex[i] = natural + remainder;
	}
	for (int i = G.nr_pairs; i < G.size; i++)
	{
		float expected_succes = nr_scen*(1 - G.ndds[i].failprob);
		float natural;
		float remainder = modf(expected_succes, &natural);
		// Round it to an integer probabilistically.
		if (rand() % 1000 > remainder * 1000)
			remainder = 1;
		else
			remainder = 0;
		succes_per_vertex[i] = natural + remainder;
	}
	for (int i = 0; i < nr_scen; i++)
	{
		Scenarios[i] = G;
		Scenarios[i].arcs.resize(0); // Empty out the arcs.
		vector<bool> vertex_succes(G.size);
		for (int j = 0; j < G.size; j++)
		{
			if (rand() % (nr_scen - i) < succes_per_vertex[j])
			{
				succes_per_vertex[j]--;
				vertex_succes[j] = 1;
			}
		}
		int number_arcs = 0;
		for (int j = 0; j < G.arcs.size(); j++)
		{
			if (vertex_succes[G.arcs[j].startvertex] == 1 && vertex_succes[G.arcs[j].endvertex] == 1)
			{
				Scenarios[i].arcs.push_back(G.arcs[j]);
				//Scenarios[i].arcs[number_arcs].weight = 1000 + (rand() % 10);
				number_arcs++;
			}
		}
	}

	return Scenarios;
}