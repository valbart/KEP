#include "stdafx.h"
#include "Structures.h"
#include "Functions.h"

using namespace std;

matching_result Hybrid_PIEF(directedgraph G, int chainlength, int cyclelength, configuration & config)
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
	for (int i = 0; i < G.nr_pairs - 1 ; i++)
	{
		for (int j = 0; j < cyclelength; j++)
		{
			for (int k = 0; k < Cyclevar[i][j].getSize(); k++)
			{
				obj.setLinearCoef(Cyclevar[i][j][k], G.arcs[Cyclevar_arc_link[i][j][k]].weight);
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
	if (config.solver == 4 || config.solver == 5)
	{
		HPIEF_CPLEX.setOut(env.getNullStream());
	}
	//HPIEF_CPLEX.exportModel("HPIEF.lp");
	HPIEF_CPLEX.solve();
	//HPIEF_CPLEX.writeSolution("HPIEF.sol");

	matching_result results;
	results.objective = HPIEF_CPLEX.getObjValue();

	results.cycles = Cycle_Solution(G, Cyclevar, Cyclevar_arc_link, HPIEF_CPLEX);
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

vector<cycle_arcs> Cycle_Solution(const directedgraph &G, const vector<vector<IloNumVarArray>> & cycle_var, const vector<vector<vector<int>>> & cycle_link, IloCplex HPIEF_CPLEX)
{
	vector<cycle_arcs> cycles;
	for (int i = 0; i < cycle_var.size(); i++) // Iterate over the graph copies.
	{
		cycle_arcs new_cycle;
		new_cycle.arcs.resize(0);
		// Save all the arcs (in this graph copy) that are used to the 'arcs' vector.
		for (int j = 0; j < cycle_var[i].size(); j++)
		{
			for (int k = 0; k < cycle_var[i][j].getSize(); k++)
			{
				float value = HPIEF_CPLEX.getValue(cycle_var[i][j][k]);
				if (value > 0.9999)
				{
					new_cycle.arcs.push_back(cycle_link[i][j][k]);
				}
					
			}
		}
		if(new_cycle.arcs.size() > 0)
		{
			int j = 0;
			int next_vertex = i; // We want to put the vertices in order. We start with the vertex equal to the graph copy.
			// Add a new cycle
			while(j < new_cycle.arcs.size()+1)
			{
				for (int k = 0; k < new_cycle.arcs.size(); k++)
				{
					if (G.arcs[new_cycle.arcs[k]].startvertex == next_vertex)
					{
						new_cycle.vertices.push_back(next_vertex);
						new_cycle.weight += G.arcs[new_cycle.arcs[k]].weight;
						next_vertex = G.arcs[new_cycle.arcs[k]].endvertex;
						j++;
						break;
					}
				}
			}
			cycles.push_back(new_cycle);
		}
	}
	return cycles;
}

vector<chain> Chain_Solution(const directedgraph & G, const vector<IloNumVarArray>& chain_var, const vector<vector<int>>& chain_link, IloCplex HPIEF_CPLEX)
{
	vector<vector<float>> solution_values(chain_var.size());
	for (int i = 0; i < chain_var.size(); i++) // Extract the solution value to a re-usable vector.
	{
		solution_values[i].resize(chain_var[i].getSize());
		for (int j = 0; j < chain_var[i].getSize(); j++)
		{
			if (HPIEF_CPLEX.getValue(chain_var[i][j]) > 0.99)
			{
				solution_values[i][j] = 1;
			}
		}
	}

	vector<chain> chain_solution;
	for (int i = 0; i < solution_values[0].size(); i++) // We go through all the arcs in position 0, i.e. the possible starts of a chain.
	{
		if (solution_values[0][i] == 1) // If the arc is used, we extract the chain.
		{
			chain new_chain;
			new_chain.vertices.push_back(G.arcs[chain_link[0][i]].startvertex);
			new_chain.vertices.push_back(G.arcs[chain_link[0][i]].endvertex);
			new_chain.weight += G.arcs[chain_link[0][i]].weight;
			for (int j = 1; j < solution_values.size(); j++) // We now go through the following positions.
			{
				for (int k = 0; k < solution_values[j].size(); k++)
				{
					if (solution_values[j][k] == 1 && G.arcs[chain_link[j][k]].startvertex == new_chain.vertices.back()) // If the arc is used and starts in the vertex that is the last in the current chain.
					{
						new_chain.vertices.push_back(G.arcs[chain_link[j][k]].endvertex);
						new_chain.weight += G.arcs[chain_link[j][k]].weight;
					}
				}
			}
			chain_solution.push_back(new_chain);
		}
	}
	
	return chain_solution;
}

cycle_variables Generate_Cycle_Var(IloEnv &env, directedgraph G, int cyclelength)
{
	// Note that we consistently work with nr_pairs - 1. For obvious reasons, the copy corresonding to the last pair can not contain any cycles.
	cycle_variables c;
	c.Cyclevariable.resize(G.nr_pairs - 1);
	c.Link_Cyclevar_Arc.resize(G.nr_pairs-1);
	// Pre-Processing
	vector<vector<directedarc>> Acopies = DP_Copy(G);
	vector<vector<vector<int>>> copy_pos_arc_possible = cycle_preproces(G, Acopies, cyclelength);

	// Create all variables
	for (int i = 0; i < G.nr_pairs-1; i++)
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
					c.Cyclevariable[i][j].add(IloNumVar(env, 0, 1, ILOINT, vname));
					c.Link_Cyclevar_Arc[i][j].push_back(Acopies[i][k].arcnumber);
				}
			}
		}
	}
	return c;
}

vector<vector<vector<int>>> cycle_preproces(directedgraph G, const vector<vector<directedarc>> & Acopies, int cyclelength)
{
	// The pre-proces works as follows. First, for every copy, we calculate the minimum number of steps necessary to get back to the starting
	// vertex of that copy.

	// Next, for each copy and each position, we go through all of the arcs and see whether they can be in that position. We check whether 
	// 1) The arc is valid in that copy (no start or endvertex that falls outside of the copy).
	// 2) In the previous position, a valid arc arrives in the starvertex. (In position 1, only arcs leaving from the "copy-vertex" are valid).
	// 3) Whether there is a short enough path back to the "Copy-vertex" from the arcs endvertex.

	// First index is copy, second is vertex. Distance is G.size+1 if there is no path.
	vector<vector<int>> distance_to_copy_vertex = distance_calc(G, Acopies, cyclelength);
	
	// First index is copy, second is position, third is arc.
	vector<vector<vector<int>>> arc_position_possible(G.nr_pairs - 1);
	for (int i = 0; i < G.nr_pairs - 1; i++)
	{
		if (distance_to_copy_vertex[i][i] > cyclelength) // If there is no possible cycle, set arc_position_possible to 0 for all positions and arcs.
		{
			arc_position_possible[i].resize(cyclelength);
			for (int j = 0; j < cyclelength; j++)
			{
				arc_position_possible[i][j].resize(Acopies[i].size());
			}
		}
		else // Otherwise, we start checking which arcs are possible.
		{
			arc_position_possible[i].resize(cyclelength);
			for(int j = 0; j < cyclelength; j++) // Resize the vectors.
			{
				arc_position_possible[i][j].resize(Acopies[i].size());
			}
			for (int k = 0; k < Acopies[i].size(); k++) // First a loop for the first position, we check whether an arc starts in the copy-vertex and whether a loop is possible through the arc.endvertex.
			{
				if (Acopies[i][k].startvertex == i && distance_to_copy_vertex[i][Acopies[i][k].endvertex] < cyclelength)
					arc_position_possible[i][0][k] = 1;
			}
			for (int j = 1; j < cyclelength; j++) // We start a loop through all further positions.
			{
				for (int k = 0; k < Acopies[i].size(); k++) // We go through all arcs again.
				{
					if (Acopies[i][k].endvertex == i || distance_to_copy_vertex[i][Acopies[i][k].endvertex] + j < cyclelength) // Only look at arcs for which the loop can still be completed (or if arc ends in copy vertex).
					{
						int startvertex = Acopies[i][k].startvertex; // Save the startvertex.
						for (int l = 0; l < Acopies[i].size(); l++) // We go through the arcs again and check whether an arc, valid in the previous positions, ended in the starvertex of the considered arc.
						{
							if (Acopies[i][l].endvertex == Acopies[i][k].startvertex && arc_position_possible[i][j - 1][l] == 1)
							{
								arc_position_possible[i][j][k] = 1; // If this is the case, the arc is possible and added to the model.
								break;
							}
						}
					}
				}
			}
		}
	}

	/*for (int i = 0; i < G.nr_pairs - 1; i++)
	{
		cout << "Copy " << i << endl;
		for (int j = 0; j < cyclelength; j++)
		{
			cout << "Position " << j << endl;
			for (int k = 0; k < Acopies[i].size(); k++)
			{
				cout << "((" << Acopies[i][k].startvertex << "," << Acopies[i][k].endvertex << ")" << arc_position_possible[i][j][k] << ") \t";
			}
			cout << endl;
		}
	}*/
	return arc_position_possible;
}

vector<vector<int>> distance_calc(const directedgraph & G, const vector<vector<directedarc>> & Acopies, int cyclelength)
{
	vector<vector<int>> distance(G.nr_pairs); // Leave out last copy (is empty).

	for (int i = 0; i < G.nr_pairs; i++) // We calculate distances for each graph copy. The distance is the distance back to the vertex = graph copy.
	{
		distance[i].resize(G.nr_pairs, cyclelength+1);
		// First calculate the distance "1" vertices.
		for (int k = 0; k < Acopies[i].size(); k++)
		{
			if (Acopies[i][k].endvertex == i) // If the endingvertex is the "copy-vertex" the distance from the startvertex is 1.
				distance[i][Acopies[i][k].startvertex] = 1;
		}
		for (int j = 2; j <= cyclelength; j++) // We need to consider at most a number of steps equal to the cyclelength. If the path is longer, the exact distance does not matter.
		{
			for (int k = 0; k < Acopies[i].size(); k++)
			{
				if (distance[i][Acopies[i][k].endvertex] == j - 1 && distance[i][Acopies[i][k].startvertex] > j)
					distance[i][Acopies[i][k].startvertex] = j;
			}
		}
		/*for (int j = 0; j < G.nr_pairs; j++)
		{
			cout << distance[i][j] << "\t";
		}
		cout << endl;*/
	}
	
	return distance;
}

chain_variables Generate_Chain_Var(IloEnv & env, directedgraph G, int chainlength)
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
			if(arc_position_possible[i][j] == 1)
			{
			if (i == 0 && G.arcs[j].startvertex >= G.nr_pairs) // If this is true, the starvertex is an NDD and can be placed in the first position of a chain.
			{
				ostringstream convert;
				convert << "y(" << G.arcs[j].startvertex << "," << G.arcs[j].endvertex << "," << i << ")";
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

vector<vector<int>> chain_preproces(directedgraph G, int chainlength)
{
	// First index is position, second is arc. The value is 1 if the arc can be used in this position, 0 otherwise.
	vector<vector<int>> arc_position_possible(chainlength);
	
	for (int j = 0; j < G.arcs.size(); j++)
	{
		if (G.arcs[j].startvertex >= G.nr_pairs) // If the startvertex is an NDD, it is possible to use it in position 0.
			arc_position_possible[0].push_back(1);
		else
			arc_position_possible[0].push_back(0);			
	}
	for (int i = 1; i < chainlength; i++) // For all further positions
	{
		arc_position_possible[i].resize(G.arcs.size());
		for (int j = 0; j < G.arcs.size(); j++) // We go through all arcs and save the startvertex
		{
			int startvertex = G.arcs[j].startvertex;
			for (int k = 0; k < G.arcs.size(); k++) // We check whether an allowed arc in the previous position ends in this startvertex.
													// If so, the current arc is allowed in the current position.
			{
				if (G.arcs[k].endvertex == startvertex && arc_position_possible[i - 1][k] == 1)
				{
					arc_position_possible[i][j] = 1;
					break;
				}
					
			}
		}
	}


	return arc_position_possible;
}

IloRangeArray Build_Vertex_Constraint(IloEnv & env, IloModel &model, directedgraph G, vector<vector<IloNumVarArray>> cycle_var, vector<vector<vector<int>>> cycle_link, vector<IloNumVarArray> chain_var, vector<vector<int>> chain_link)
{
	IloRangeArray Vertex_Constraint(env, G.nr_pairs);
	for (int i = 0; i < G.nr_pairs; i++)
	{
		IloExpr expr(env);
		Vertex_Constraint[i] = IloRange(expr <= 1);
		Vertex_Constraint[i].setUB(1);
	}

	for (int i = 0; i < G.nr_pairs - 1; i++)
	{
		for (int j = 0; j < cycle_var[i].size(); j++)
		{
			for (int k = 0; k < cycle_var[i][j].getSize(); k++)
			{
				Vertex_Constraint[G.arcs[cycle_link[i][j][k]].endvertex].setLinearCoef(cycle_var[i][j][k],1);
			}
		}
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
}

IloRangeArray Build_NDD_Constraint(IloEnv & env, IloModel & model, directedgraph G, vector<IloNumVarArray> chain_var, vector<vector<int>> chain_link)
{
	IloRangeArray NDD_Constraint(env, G.nr_ndd);
	for (int i = 0; i < G.nr_ndd; i++)
	{
		IloExpr expr(env);
		NDD_Constraint[i] = IloRange(expr <= 1);
	}
	for (int i = 0; i < chain_var[0].getSize(); i++)
	{
		NDD_Constraint[G.arcs[chain_link[0][i]].startvertex - G.nr_pairs].setLinearCoef(chain_var[0][i], 1);
	}
	model.add(NDD_Constraint);

	return NDD_Constraint;
}

vector<vector<IloRangeArray>> Build_Vertex_Flow_Constraint(IloEnv & env, IloModel &model, directedgraph G, vector<vector<IloNumVarArray>> cycle_var, vector<vector<vector<int>>> cycle_link, int cyclelength)
{
	// Note, since constraints for a vertex j <= graph copy do not exist, the indexing is as follows, VFC[i][j][k] refers to the "i"th graph copy, "j - i - 1"th vertex and "k"th position. 
	
	vector<vector<IloRangeArray>> Vertex_Flow_Constraint(G.nr_pairs-1); // Last copy is empty.
	// First Graph Copy, then Vertex, then Position.
	for (int i = 0; i < G.nr_pairs - 1; i++)
	{
		Vertex_Flow_Constraint[i].resize(G.nr_pairs - i - 1); // Copy i holds G.nr_pairs - i vertices. No constraints are necessary for the first vertex (-1).
		for (int j = 0; j < G.nr_pairs - i - 1; j++)
		{
			Vertex_Flow_Constraint[i][j] = IloRangeArray(env, cyclelength - 1);
			for (int k = 0; k < cyclelength - 1; k++)
			{
				ostringstream convert;
				convert << "EqFlow(C:" << i << ",V:" << j + i + 1 << ",P:" << k << ")";
				string varname = convert.str();
				const char* vname = varname.c_str();
				Vertex_Flow_Constraint[i][j][k] =  IloRange(env, 0, 0, vname);
			}
		}
	}
	for (int copy = 0; copy < G.nr_pairs - 1; copy++) // Graph Copy, there are no cycle_vars for the final copy, so should be excluded from the loop.
	{
		for (int j = 0; j < cyclelength; j++) // Position in cycle
		{
			for (int k = 0; k < cycle_var[copy][j].getSize(); k++) // We go through all cycle variables at this copy and position.
			{
				if (G.arcs[cycle_link[copy][j][k]].startvertex != copy)
				{
					// This arc leaves the vertex G.arcs[indices].starvertex (=l) in position k. It should be added to VFC[copy][l-1-copy][j-1] with negative coefficient.
					Vertex_Flow_Constraint[copy][G.arcs[cycle_link[copy][j][k]].startvertex - 1 - copy][j - 1].setLinearCoef(cycle_var[copy][j][k], - 1);
				}
				if (G.arcs[cycle_link[copy][j][k]].endvertex != copy)
				{
					// This arc arrives in the vertex G.arcs[indices].endvertex (=l) in position k. It should be added to VFC[copy][l-1-copy][j] with positive coefficnet.
					Vertex_Flow_Constraint[copy][G.arcs[cycle_link[copy][j][k]].endvertex - 1 - copy][j].setLinearCoef(cycle_var[copy][j][k], 1);
				}
			}
		}
	}
	
	for (int i = 0; i < G.nr_pairs - 1; i++)
	{
		for (int j = 0; j < G.nr_pairs - i - 1; j++)
			model.add(Vertex_Flow_Constraint[i][j]);
	}

	return Vertex_Flow_Constraint;
}

vector<IloRangeArray> Build_Vertex_Flow_Chain_Constraint(IloEnv & env, IloModel & model, directedgraph G, vector<IloNumVarArray> chain_var, vector<vector<int>> chain_link, int chainlength)
{
	vector<IloRangeArray> Vertex_Flow_Chain_Constraint(G.nr_pairs); // First Vertex, then Position.
	
	for (int i = 0; i < G.nr_pairs; i++)
	{
		Vertex_Flow_Chain_Constraint[i] = IloRangeArray(env, chainlength-1);
		for (int j = 0; j < chainlength - 1; j++)
		{
			ostringstream convert;
			convert << "EqFlowChain(V:" << i << ",P:" << j << ")";
			string varname = convert.str();
			const char* vname = varname.c_str();
			IloExpr expr(env);
			Vertex_Flow_Chain_Constraint[i][j] = IloRange(expr <= 0);
			Vertex_Flow_Chain_Constraint[i][j].setName(vname);
		}
	}
	for (int i = 0; i < chainlength; i++) // We go through all Chainvar
	{
		for (int j = 0; j < chain_var[i].getSize(); j++)
		{
			if (G.arcs[chain_link[i][j]].startvertex < G.nr_pairs)
			{
				Vertex_Flow_Chain_Constraint[G.arcs[chain_link[i][j]].startvertex][i - 1].setLinearCoef(chain_var[i][j], 1);
			}
			if (G.arcs[chain_link[i][j]].endvertex < G.nr_pairs && i != chainlength-1)
			{
				Vertex_Flow_Chain_Constraint[G.arcs[chain_link[i][j]].endvertex][i].setLinearCoef(chain_var[i][j], -1);
			}
		}
	}

	for (int i = 0; i < G.nr_pairs; i++)
	{
		model.add(Vertex_Flow_Chain_Constraint[i]);
	}

	return Vertex_Flow_Chain_Constraint;
}

vector<vector<directedarc>> DP_Copy(directedgraph G)
{
	// We make copies of the graph (or at least the arc list) for each vertex.
	// The copy associated with a vertex only includes arcs that start and end in vertices equal to or higher than the "copy-number".
	vector<vector<directedarc>> Acopies(G.nr_pairs); // The new copy.
	for (int i = 0; i < G.arcs.size(); i++) // We loop through all of the original arcs.
	{
		int minvertex = min(G.arcs[i].startvertex, G.arcs[i].endvertex); // We look for the lowest vertex in the arc.
		if (G.arcs[i].startvertex < G.nr_pairs) // If the starting vertex is an NDD, the arc does not need to be in the copy (which is exclusively for cycles).
		{
			for (int j = 0; j <= minvertex; j++) // We add the arc to the copy if and only if start and endvertex numbers are higher (or equal) to the graphy copy-number.
			{
				Acopies[j].push_back(G.arcs[i]); // We save the arc to the copy of the graph.
			}
		}
	}
	
	return Acopies;
}