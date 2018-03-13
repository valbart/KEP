#include "stdafx.h"
#include "Structures.h"
#include "Functions.h"

void Pre_Test_Float(directedgraph G, int chainlength, int cyclelength, int max_tests, int nr_scen, int time_limit, int scenario_generator, int failure_type);
cycle_variables Generate_Cycle_Var_Float(IloEnv &env, const directedgraph & G, int cyclelength, int nr_scen);
IloNumVarArray Generate_Testvar_Float(IloEnv & env, directedgraph G);

using namespace std;

void pre_test_main(configuration & config, directedgraph G)
{
	
	pre_test_result results = Pre_Test(G, config.cyclelength, config.chainlength, config.max_test, config.nr_scenarios, config.time_limit, config.scen_gen, config.failure_type, config);
	Output_Pre_Test(results, config);
}

pre_test_result Pre_Test(directedgraph G, int chainlength, int cyclelength, int max_tests, int nr_scen, int time_limit, int scenario_generator, int failure_type, const configuration & config)
{
	directedgraph Tested_Graph = G;

	

	vector<directedgraph> Scenarios;
	if (failure_type == 1)
	{
		cout << "Arcs Fail" << endl;
		if (scenario_generator == 1)
		{
			Scenarios = Generate_Scenarios_Tight(G, nr_scen); cout << "Tight Scen Generator" << endl;
		}
		else
		{
			Scenarios = Generate_Scenarios(G, nr_scen); cout << "Basic Scen Generator" << endl;
		}
			
	}
	else if (failure_type == 2)
	{
		cout << "Vertices Fail" << endl;
		Scenarios = Generate_Scenarios_Vertex_Tight(G, nr_scen);
	}
	cout << "Scenarios Generated" << endl;
	cout << Scenarios.size();

	
	time_t start_time;
	time(&start_time);
	IloEnv env; 
	IloModel model(env);
	cout << "Generating Variables" << endl;
	vector<vector<vector<IloNumVarArray>>> Cyclevar(nr_scen); // First position is scenario, second index is the Graph Copy, third the position in the Graph, fourth the individual arcs.
	vector<vector<vector<vector<int>>>> Cyclevar_arc_link(nr_scen); // A vector to link the variables to the original arc. Cyclevar_arc_link[i][j][k][l] = m, means that this variable corresponds to the m-th arc in the original arc list.
	for(int i = 0; i < nr_scen; i++)
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
	for(int scen = 0; scen < nr_scen; scen++)
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
	CPLEX.setParam(IloCplex::TreLim, config.memory_limit);
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

vector<directedgraph> Generate_Scenarios(const directedgraph & G, int nr_scen)
{
	//ofstream output;
	//output.open("Scenarios.txt");
	
	srand(time(NULL));

	vector<directedgraph> Scenarios;
	for (int i = 0; i < nr_scen; i++)
	{
		//output << "Scenario " << i << endl;
		Scenarios.push_back(G);
		Scenarios[i].arcs.resize(0); // Empty out the arcs.
		for (int j = 0; j < G.arcs.size(); j++)
		{
			if (rand() % 100 >= Scenarios[i].arcs[j].failprob*100)
			{
				Scenarios[i].arcs.push_back(G.arcs[j]);
				//output << "(" << G.arcs[j].startvertex << "," << G.arcs[j].endvertex << ")" << endl;
			}
		}
		//output << endl;
	}
	return Scenarios;
}

vector<directedgraph> Generate_Scenarios_Tight(const directedgraph & G, int nr_scen)
{
	ofstream output;
	srand(time(NULL));
	vector<directedgraph> Scenarios(nr_scen);
	// Get the number of successes for each arc  and the number of scenarios.
	vector<int> succes_per_arc(G.arcs.size());
	for (int i = 0; i < G.arcs.size(); i++)
	{
		float expected_succes = nr_scen*(1 - G.arcs[i].failprob);
		float natural;
		float remainder = modf(expected_succes, & natural);
		// Round it to an integer probabilistically.
		if (rand() % 1000 > remainder * 1000)
			remainder = 1;
		else
			remainder = 0;
		succes_per_arc[i] = natural + remainder;
	}
	for (int i = 0; i < nr_scen; i++)
	{
		Scenarios[i] = G;
		Scenarios[i].arcs.resize(0); // Empty out the arcs.
		for (int j = 0; j < G.arcs.size(); j++)
		{
			if (rand() % (nr_scen - i) < succes_per_arc[j])
			{
				succes_per_arc[j]--;
				Scenarios[i].arcs.push_back(G.arcs[j]);
			}
		}
	}
	return Scenarios;
}

vector<directedgraph> Generate_Scenarios_Vertex_Tight(const directedgraph & G, int nr_scen)
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
		for (int j = 0; j < G.arcs.size(); j++)
		{
			if (vertex_succes[G.arcs[j].startvertex] == 1 && vertex_succes[G.arcs[j].endvertex] == 1)
			{
				Scenarios[i].arcs.push_back(G.arcs[j]);
			}
		}
	}

	return Scenarios;
}

cycle_variables Generate_Cycle_Var(IloEnv &env, const directedgraph & G, int cyclelength, int nr_scen)
{
	// Note that we consistently work with nr_pairs - 1. Since the copy corresponding to the last pair is empty (no arcs), it is useless to include it.
	cycle_variables c;
	c.Cyclevariable.resize(G.nr_pairs - 1);
	c.Link_Cyclevar_Arc.resize(G.nr_pairs - 1);
	// Pre-Processing
	vector<vector<directedarc>> Acopies = DP_Copy(G);
	vector<vector<vector<int>>> copy_pos_arc_possible = cycle_preproces(G, Acopies, cyclelength);

	// Create all variables
	{
		//ostringstream convert;
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
						//convert << "x(" << nr_scen << "," << Acopies[i][k].startvertex << "," << Acopies[i][k].endvertex << "," << i << "," << j << ")";
						//string varname = convert.str();
						//const char* vname = varname.c_str();
						c.Cyclevariable[i][j].add(IloNumVar(env, 0, 1, ILOINT/*, vname*/));
						c.Link_Cyclevar_Arc[i][j].push_back(Acopies[i][k].arcnumber);
						//convert.str("");
						//convert.clear();
					}
				}
			}
		}
	}
	return c;
}

chain_variables Generate_Chain_Var(IloEnv & env, directedgraph G, int chainlength, int nr_scen)
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
					convert << "y(" << nr_scen << "," << G.arcs[j].startvertex << "," << G.arcs[j].endvertex << "," << i << ")";
					string varname = convert.str();
					const char* vname = varname.c_str();
					c.Chainvar[i].add(IloNumVar(env, 0, 1, ILOINT, vname));
					c.Link_Chainvar_Arc[i].push_back(G.arcs[j].arcnumber);
				}
				else if (i > 0 && G.arcs[j].startvertex < G.nr_pairs) // Any arc between pairs can be used in the following positions
				{
					ostringstream convert;
					convert << "y(" << G.arcs[j].startvertex << "," << G.arcs[j].endvertex << "," << i << ")";
					string varname = convert.str();
					const char* vname = varname.c_str();
					c.Chainvar[i].add(IloNumVar(env, 0, 1, ILOINT, vname));
					c.Link_Chainvar_Arc[i].push_back(G.arcs[j].arcnumber);
				}
			}
		}
	}

	return c;
}

IloNumVarArray Generate_Testvar(IloEnv & env, directedgraph G)
{
	IloNumVarArray Testvar(env);
	for (int i = 0; i < G.arcs.size(); i++)
	{
		ostringstream convert;
		convert << "t("  << G.arcs[i].startvertex << "," << G.arcs[i].endvertex << ")";
		string varname = convert.str();
		const char* vname = varname.c_str();
		Testvar.add(IloNumVar(env, 0, 1, ILOINT, vname));
	}
	
	return Testvar;
}

vector<IloRangeArray> Build_Test_Constraint(IloEnv & env, IloModel model, directedgraph G, const IloNumVarArray & Testvar, const vector<vector<vector<IloNumVarArray>>>& Cyclevar, const vector<vector<vector<vector<int>>>>& cycle_link, const vector<vector<IloNumVarArray>>& Chainvar, const vector<vector<vector<int>>>& chain_link, int nr_scen)
{
	vector<IloRangeArray> Test_Constraints(nr_scen);

	for (int scen = 0; scen < nr_scen; scen++)
	{
		Test_Constraints[scen] = IloRangeArray(env, G.arcs.size());
		for (int i = 0; i < G.arcs.size(); i++)
		{
			ostringstream convert;
			convert << "Test(S:" << scen << ",(" << G.arcs[i].startvertex << "," << G.arcs[i].endvertex << "))";
			string varname = convert.str();
			const char* vname = varname.c_str();
			Test_Constraints[scen][i] = IloRange(-Testvar[i] <= 0);
			Test_Constraints[scen][i].setName(vname);
		}
		for (int i = 0; i < Cyclevar[scen].size(); i++) // Scen
		{
			for (int j = 0; j < Cyclevar[scen][i].size(); j++) // Copy
			{
				for (int k = 0; k < Cyclevar[scen][i][j].getSize(); k++) // Position
				{
					Test_Constraints[scen][cycle_link[scen][i][j][k]].setLinearCoef(Cyclevar[scen][i][j][k],1);
				}
			}
		}
		for (int i = 0; i < Chainvar[scen].size(); i++) // Scen
		{
			for (int j = 0; j < Chainvar[scen][i].getSize(); j++) // Position
			{

				Test_Constraints[scen][chain_link[scen][i][j]].setLinearCoef(Chainvar[scen][i][j], 1);
			}
		}
		model.add(Test_Constraints[scen]);
	}
	
	return Test_Constraints;
}

IloRange Build_Max_Test_Constraint(IloEnv & env, IloModel model, const IloNumVarArray & Testvar, int max_test)
{
	IloExpr expr(env);
	IloRange Constraint(expr <= max_test);
	for (int i = 0; i < Testvar.getSize(); i++)
		Constraint.setLinearCoef(Testvar[i], 1);

	model.add(Constraint);
	return Constraint;
}

void Output_Pre_Test(const pre_test_result & results, const configuration & config)
{
	ofstream output;
	output.open(config.testvar_output);
	output << "Nr_Pairs = " << config.nr_pairs << endl;
	output << "Nr_NDD = " << config.nr_ndd << endl;

	for (int i = 0; i < results.tested_arcs.size(); i++)
	{
		output << "(" << results.tested_arcs[i].startvertex << "," << results.tested_arcs[i].endvertex << "), " << results.tested_arcs[i].failprob << ", " << results.tested_arcs[i].weight << endl;
	}

	output << "Objective Value = " << results.objective_value << endl;
	output << "Computation Time = " << results.computation_time << endl;
	output << "Nr Scenarios = " << config.nr_scenarios << endl;
	output << "Nr Test = " << config.max_test << endl;
}

IloNumVarArray Generate_Testvar_Float(IloEnv & env, directedgraph G)
{
	IloNumVarArray Testvar(env);
	for (int i = 0; i < G.arcs.size(); i++)
	{
		ostringstream convert;
		convert << "t(" << G.arcs[i].startvertex << "," << G.arcs[i].endvertex << ")";
		string varname = convert.str();
		const char* vname = varname.c_str();
		Testvar.add(IloNumVar(env, 0, 1, ILOFLOAT, vname));
	}

	return Testvar;
}

cycle_variables Generate_Cycle_Var_Float(IloEnv &env, const directedgraph & G, int cyclelength, int nr_scen)
{
	// Note that we consistently work with nr_pairs - 1. Since the copy corresponding to the last pair is empty (no arcs), it is useless to include it.
	cycle_variables c;
	c.Cyclevariable.resize(G.nr_pairs - 1);
	c.Link_Cyclevar_Arc.resize(G.nr_pairs - 1);
	// Pre-Processing
	vector<vector<directedarc>> Acopies = DP_Copy(G);
	vector<vector<vector<int>>> copy_pos_arc_possible = cycle_preproces(G, Acopies, cyclelength);

	// Create all variables
	{
		//ostringstream convert;
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
						//convert << "x(" << nr_scen << "," << Acopies[i][k].startvertex << "," << Acopies[i][k].endvertex << "," << i << "," << j << ")";
						//string varname = convert.str();
						//const char* vname = varname.c_str();
						c.Cyclevariable[i][j].add(IloNumVar(env, 0, 1, ILOFLOAT/*, vname*/));
						c.Link_Cyclevar_Arc[i][j].push_back(Acopies[i][k].arcnumber);
						//convert.str("");
						//convert.clear();
					}
				}
			}
		}
	}
	return c;
}

void Pre_Test_Float(directedgraph G, int chainlength, int cyclelength, int max_tests, int nr_scen, int time_limit, int scenario_generator, int failure_type)
{
	directedgraph Tested_Graph = G;



	vector<directedgraph> Scenarios;
	if (failure_type == 1)
	{
		cout << "Arcs Fail" << endl;
		if (scenario_generator == 1)
		{
			Scenarios = Generate_Scenarios_Tight(G, nr_scen); cout << "Tight Scen Generator" << endl;
		}
		else
		{
			Scenarios = Generate_Scenarios(G, nr_scen); cout << "Basic Scen Generator" << endl;
		}

	}
	else if (failure_type == 2)
	{
		cout << "Vertices Fail" << endl;
		Scenarios = Generate_Scenarios_Vertex_Tight(G, nr_scen);
	}
	cout << "Scenarios Generated" << endl;
	cout << Scenarios.size();


	time_t start_time;
	time(&start_time);
	IloEnv env;
	IloModel model(env);
	cout << "Generating Variables" << endl;
	vector<vector<vector<IloNumVarArray>>> Cyclevar(nr_scen); // First position is scenario, second index is the Graph Copy, third the position in the Graph, fourth the individual arcs.
	vector<vector<vector<vector<int>>>> Cyclevar_arc_link(nr_scen); // A vector to link the variables to the original arc. Cyclevar_arc_link[i][j][k][l] = m, means that this variable corresponds to the m-th arc in the original arc list.
	for (int i = 0; i < nr_scen; i++)
	{
		cycle_variables cvars = Generate_Cycle_Var_Float(env, Scenarios[i], cyclelength, i);
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
	IloNumVarArray Testvar = Generate_Testvar_Float(env, G);
	vector<IloRangeArray> test_constraint = Build_Test_Constraint(env, model, G, Testvar, Cyclevar, Cyclevar_arc_link, Chainvar, Chainvar_arc_link, nr_scen);
	IloRange Max_Test_Constraint = Build_Max_Test_Constraint(env, model, Testvar, max_tests);

	IloCplex CPLEX(model);
	CPLEX.exportModel("Floatmodel.lp");
	CPLEX.setParam(IloCplex::TiLim, time_limit);
	CPLEX.solve();

	pre_test_result results;
	results.objective_value = CPLEX.getObjValue() / nr_scen;
	cout << results.objective_value << endl;

	time_t current_time;
	time(&current_time);
	results.computation_time = difftime(current_time, start_time);

	ofstream output;
	output.open("Float_LP");
	output << "Nr_Pairs = " << G.nr_pairs << endl;
	output << "Nr_NDD = " << G.nr_ndd << endl;

	for (int i = 0; i < Testvar.getSize(); i++)
	{
		if(CPLEX.getValue(Testvar[i]) > 0)
		output << "(" << G.arcs[i].startvertex << "," << G.arcs[i].endvertex << "), " << CPLEX.getValue(Testvar[i]) << endl;
	}

	output << "Objective Value = " << results.objective_value << endl;
	output << "Computation Time = " << results.computation_time << endl;
	output << "Nr Test = " << max_tests << endl;
	cin.get();
}