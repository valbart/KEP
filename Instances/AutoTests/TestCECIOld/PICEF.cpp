#include "stdafx.h"
#include "Structures.h"
#include "Functions.h"

using namespace std;

matching_result PICEF(directedgraph G, configuration & config)
{
	// This model uses position-indexed arc variables for the chains, and variables for each possible cycle.
	IloEnv env;
	IloModel PICEF(env);

	// Create the Variables. (pre-processing can be implemented here)
	vector<IloNumVarArray> Chainvar; // First index is the position in the graph, Second the individual arc.
	vector<vector<int>> Chainvar_arc_link; // A vector to link the variables to the original arc.
	{
		chain_variables cvars = Generate_Chain_Var(env, G, config.chainlength);
		Chainvar = cvars.Chainvar;
		Chainvar_arc_link = cvars.Link_Chainvar_Arc;
	}

	vector<cycle_arcs> cycles = Find_Cycles(G, config);

	// Create the Objective Function
	// Currently an unweigted placeholder.
	IloObjective obj = IloMaximize(env);
	for (int i = 0; i < config.chainlength; i++)
	{
		for (int j = 0; j < Chainvar[i].getSize(); j++)
			obj.setLinearCoef(Chainvar[i][j], 1);
	}
	PICEF.add(obj);


	// Build the model, references are to the Dickerson et al. 2016 paper.
	// Maximize \sum[cycles] cycle_weight[cycle] * cycle-var [cycle] + \sum[arcs, position] weight[arc] * chain-x [arc, position] (4a)
	// Subject to
	// \forall[pair-vertex] \sum[incoming arcs, positions] chain-x[arc, position] + \sum[cycle including vertex] cycle-var. (4b)
	// \forall[NDD-Vertex] \sum[outgoing arcs] chain-x[arc, position 1] <= 1 (4c)
	// \forall[Pair-Vertex, position 1..K-1] \sum[incoming arc] chain-x[arc, position] = \sum[outgoing arc] chain-x[arc, position+1] (4d)
	// Integrality Constraints

	// Create Initial Constraints
	// Max 1 incoming chain var. Later on, add the variables for cycles that use the vertex.
	IloRangeArray vertex_inflow_cons = Build_Vertex_Constraint_PICEF(env, PICEF, G, Chainvar, Chainvar_arc_link);

	// If there is an outgoing arc_chain in position i, there should be an incoming arc_chain in position [i-1].
	vector<IloRangeArray> vertex_chain_flow_cons = Build_Vertex_Flow_Chain_Constraint(env, PICEF, G, Chainvar, Chainvar_arc_link, config.chainlength);

	// Max 1 outgoing chainvar for each NDD.
	IloRangeArray NDD_Constraint = Build_NDD_Constraint(env, PICEF, G, Chainvar, Chainvar_arc_link);


	IloCplex PICEF_CPLEX(PICEF);
	PICEF_CPLEX.exportModel("PICEF.lp");
	PICEF_CPLEX.solve();

	matching_result results;
	results.objective = PICEF_CPLEX.getObjValue();

	results.chains = Chain_Solution(G, Chainvar, Chainvar_arc_link, PICEF_CPLEX);
	cout << "Chains" << endl;
	for (int i = 0; i < results.chains.size(); i++)
	{
		for (int j = 0; j < results.chains[i].vertices.size(); j++)
		{
			cout << results.chains[i].vertices[j] << "\t";
		}
		cout << endl;
	}

	return results;
}


IloRangeArray Build_Vertex_Constraint_PICEF(IloEnv & env, IloModel & model, directedgraph G, vector<IloNumVarArray> chain_var, vector<vector<int>> chain_link)
{
	IloRangeArray Vertex_Constraint(env, G.nr_pairs);
	for (int i = 0; i < G.nr_pairs; i++)
	{
		IloExpr expr(env);
		Vertex_Constraint[i] = IloRange(expr <= 1);
		Vertex_Constraint[i].setUB(1);
	}

	for (int i = 0; i < chain_var.size(); i++)
	{
		for (int j = 0; j < chain_var[i].getSize(); j++)
		{
			Vertex_Constraint[G.arcs[chain_link[i][j]].endvertex].setLinearCoef(chain_var[i][j], 1);
		}
	}

	model.add(Vertex_Constraint);
	return Vertex_Constraint;
	
	return IloRangeArray();
}