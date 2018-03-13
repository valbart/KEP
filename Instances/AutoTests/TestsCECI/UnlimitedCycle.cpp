#include "stdafx.h"
#include "Structures.h"
#include "Functions.h"

using namespace std;

// Unlimited Cycle Formulation
void UnlimitedCycle_Main(configuration & config, const directedgraph & G)
{
	// Build the scenarios
	vector<directedgraph> scenarios;
	if (config.failure_type == 1)
		scenarios = Generate_Scenarios_Tight(G, config.nr_scenarios);
	else if (config.failure_type == 2)
		scenarios = Generate_Scenarios_Vertex_Tight(G, config.nr_scenarios);
	else
	{
		cout << "ERROR; Failure type not recognized" << endl;
		cin.get();
	}

	time_t start_time;
	time(&start_time);

	// Set up the Model
	IloEnv env;
	IloModel model(env); 
	
	// Build Variables per scenario
	vector<IloNumVarArray> Arcvariables(config.nr_scenarios);
	vector<vector<int>> Arcvar_Arc_Link(config.nr_scenarios);
	for (int i = 0; i < config.nr_scenarios; i++)
	{
		arc_variables_uc arcvars = Generate_Arc_Var_UC(env, scenarios[i], i);
		Arcvariables[i] = arcvars.Arcvariable;
		Arcvar_Arc_Link[i] = arcvars.Link_Arcvariable;
	}

	// Build the linking (test) variables.
	IloNumVarArray Testvariables = Generate_Test_Var_UC(env, G);

	// Build Constraints per Scenario.
	vector<IloRangeArray> vertex_inflow_cons(config.nr_scenarios);
	vector<IloRangeArray> vertex_flow_cons(config.nr_scenarios);
	vector<IloRangeArray> arc_test_constraints(config.nr_scenarios);
	for (int scen = 0; scen < config.nr_scenarios; scen++)
	{
		// Max one incoming arc per vertex.
		vertex_inflow_cons[scen] = Build_Inflow_Constraint_UC(env, model, G, Arcvariables[scen], Arcvar_Arc_Link[scen]);

		// If there is an arc arriving in node i, there should be an outgoing arc in node i.
		vertex_flow_cons[scen] = Build_Flow_Constraint_UC(env, model, G, Arcvariables[scen], Arcvar_Arc_Link[scen]);

		// Arcs can only be used if tested.
		arc_test_constraints[scen] = Build_Test_Constraint_UC(env, model, G, Arcvariables[scen], Arcvar_Arc_Link[scen], Testvariables);
	}
	// Constraint for max number of tests.
	IloRange Max_Test_Constraint = Build_Max_Test_Constraint(env, model, Testvariables, config.max_test);

	// Create the Objective Function
	IloObjective obj = IloMaximize(env);
	for (int scen = 0; scen < config.nr_scenarios; scen++)
	{
		for (int j = 0; j < Arcvariables[scen].getSize(); j++)
		{
			obj.setLinearCoef(Arcvariables[scen][j], G.arcs[Arcvar_Arc_Link[scen][j]].weight);
		}
	}
	model.add(obj);

	IloCplex CPLEX(model);
	CPLEX.setParam(IloCplex::TiLim, config.time_limit);
	CPLEX.exportModel("UCmodel.lp");
	CPLEX.solve();

	pre_test_result results;
	time_t current_time;
	time(&current_time);
	results.computation_time = difftime(current_time, start_time);
	results.objective_value = CPLEX.getObjValue() / config.nr_scenarios;
	for (int i = 0; i < G.arcs.size(); i++)
	{
		if (CPLEX.getValue(Testvariables[i]) > 0.99)
		{
			results.tested_arcs.push_back(G.arcs[i]);
		}
	}
	Output_Pre_Test(results, config);
	
	for (int i = 0; i < config.nr_scenarios; i++)
	{
		cout << "Scenario " << i << endl;
		vector<int> outgoing(config.nr_pairs, config.nr_pairs + 1);
		// For each scenario and each vertex, find whether there is and where the ougoing arc goes.
		for (int j = 0; j < scenarios[i].arcs.size(); j++)
		{
			if (CPLEX.getValue(Arcvariables[i][j]) > 0.99)
			{
				outgoing[G.arcs[Arcvar_Arc_Link[i][j]].startvertex] = G.arcs[Arcvar_Arc_Link[i][j]].endvertex;
			}
		}
		for (int j = 0; j < config.nr_pairs; j++)
		{
			while (outgoing[j] < config.nr_pairs + 1)
			{
				cout << j << "\t";
				int next = outgoing[j];
				outgoing[j] = config.nr_pairs + 1;
				j = next;
				if (outgoing[j] == config.nr_pairs + 1)
					cout << endl;
			}
		}
	}
}

arc_variables_uc Generate_Arc_Var_UC(IloEnv & env, directedgraph G, int scenario_number)
{
	arc_variables_uc arcvars;
	arcvars.Arcvariable = IloNumVarArray(env);
	for (int i = 0; i < G.arcs.size(); i++)
	{
		ostringstream convert;
		convert << "x(" << scenario_number << ",(" << G.arcs[i].startvertex << "," << G.arcs[i].endvertex << "))";
		string varname = convert.str();
		const char* vname = varname.c_str();
		arcvars.Arcvariable.add(IloNumVar(env, 0, 1, ILOINT, vname));
		arcvars.Link_Arcvariable.push_back(G.arcs[i].arcnumber);
		convert.str("");
		convert.clear();
	}

	return arcvars;
}

IloNumVarArray Generate_Test_Var_UC(IloEnv & env, directedgraph G)
{
	IloNumVarArray Test_Var(env);
	for (int i = 0; i < G.arcs.size(); i++)
	{
		ostringstream convert;
		convert << "t(" << G.arcs[i].startvertex << "," << G.arcs[i].endvertex << ")";
		string varname = convert.str();
		const char* vname = varname.c_str();
		Test_Var.add(IloNumVar(env, 0, 1, ILOINT, vname));
		convert.str("");
		convert.clear();
	}

	return Test_Var;
}

IloRangeArray Build_Inflow_Constraint_UC(IloEnv &env, IloModel &model, const directedgraph & G, IloNumVarArray cycle_var, vector<int> cycle_link)
{
	IloRangeArray Vertex_Constraint(env, G.nr_pairs);
	for (int i = 0; i < G.nr_pairs; i++)
	{
		IloExpr expr(env);
		Vertex_Constraint[i] = IloRange(expr <= 1);
		Vertex_Constraint[i].setUB(1);
	}

	for (int k = 0; k < cycle_var.getSize(); k++)
	{
		Vertex_Constraint[G.arcs[cycle_link[k]].endvertex].setLinearCoef(cycle_var[k], 1);
	}
	model.add(Vertex_Constraint);
	return Vertex_Constraint;
}

IloRangeArray Build_Flow_Constraint_UC(IloEnv & env, IloModel & model, const directedgraph & G, IloNumVarArray cycle_var, vector<int> cycle_link)
{
	IloRangeArray Vertex_Flow_Constraint(env, G.nr_pairs);
	for (int i = 0; i < G.nr_pairs; i++)
	{
		Vertex_Flow_Constraint[i] = IloRange(env, 0, 0);
	}

	for (int j = 0; j < cycle_var.getSize(); j++)
	{
		Vertex_Flow_Constraint[G.arcs[cycle_link[j]].startvertex].setLinearCoef(cycle_var[j], 1);
		Vertex_Flow_Constraint[G.arcs[cycle_link[j]].endvertex].setLinearCoef(cycle_var[j], -1);
	}
	
	model.add(Vertex_Flow_Constraint);
	return Vertex_Flow_Constraint;
}

IloRangeArray Build_Test_Constraint_UC(IloEnv & env, IloModel & model, const directedgraph & G, IloNumVarArray cycle_var, vector<int> cycle_link, IloNumVarArray test_var)
{
	IloRangeArray Test_Constraint(env, cycle_link.size());
	for (int i = 0; i < cycle_link.size(); i++)
	{
		Test_Constraint[i] = IloRange(env, -1, 0);
		Test_Constraint[i].setLinearCoef(cycle_var[i], 1);
		Test_Constraint[i].setLinearCoef(test_var[cycle_link[i]], -1);
	}

	model.add(Test_Constraint);
	return Test_Constraint;
}
