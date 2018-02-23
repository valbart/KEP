#include "stdafx.h"
#include "Structures.h"
#include "Functions.h"
#include <algorithm>
#include <omp.h>

// Recursion to generate relevant subsets.
	// For each cycle in the graph
		// Add another cycle, which overlaps, but is not equal to the current cycle. If this subset is at or below limit, continue adding cycles.

void Subset_Recourse(configuration & config, directedgraph G)
{
	time_t start_time;
	time(&start_time);

	vector<cycle_arcs> subsets = Relevant_Subsets(G, config);

	if (config.failure_type == 1)
		Subset_Set_Weights_Arc(subsets, G, config);
	else if (config.failure_type == 2)
		Subset_Set_Weights_Vertex(subsets, G, config);

	time_t current_time;
	time(&current_time);
	double remaining_time = difftime(current_time, start_time);
	cout << "There are " << subsets.size() << " subsets." << endl;
	pre_test_result results = Subset_MIP(subsets, G, config, remaining_time, start_time);
	Output_Subset_Recourse(results, config);
}

pre_test_result Subset_MIP(const vector<cycle_arcs> & subsets, const directedgraph & G, const configuration & config, double remaining_time, const time_t & start_time)
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

	IloCplex Subset_MIP(MIP) ;
	Subset_MIP.setParam(IloCplex::TiLim, remaining_time);

	Subset_MIP.exportModel("SubsetModel.lp");
	Subset_MIP.solve();

	pre_test_result results;
	results.objective_value = Subset_MIP.getObjValue();
	time_t current_time;
	time(&current_time);
	results.computation_time = difftime(current_time, start_time);
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

void Output_Subset_Recourse(const pre_test_result & results, const configuration & config)
{
	cout << "Writing Output" << endl;
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
	output << "Nr Test = " << results.tested_arcs.size() << endl;
}

IloRangeArray Build_Vertex_Constraint_Subset_MIP(IloEnv & env, const directedgraph & G, IloNumVarArray & subsetvar, const vector<cycle_arcs>& subsets)
{
	IloRangeArray Vertex_Constraint(env, G.nr_pairs);
	for (int i = 0; i < G.nr_pairs; i++)
	{
		IloExpr expr(env);
		Vertex_Constraint[i] = IloRange(expr <= 1);
	}
	for (int i = 0; i < subsets.size(); i++)
	{
		for (int j = 0; j < subsets[i].vertices.size(); j++)
			Vertex_Constraint[subsets[i].vertices[j]].setLinearCoef(subsetvar[i], 1);
	}


	return Vertex_Constraint;
}

vector<cycle_arcs> Relevant_Subsets(const directedgraph & G, const configuration & config)
{
	// Enumeration of relevant subsets, based on Algorithm 2 in Klimentova et al. (2016).
	// Note that K.e.a. does not have any protection from replicating relevant subsets.
	// We add this by ordering the cycles and only adding cycles in ascending order. (Adding 2 to 1 is possible, adding 1 to 2 isn't).
	// One exception is if a lower ordered cycle could not be added previously (through no overlapping vertices). If this happens, we re-order that cycle at the end of the list.

	vector<cycle_arcs> cycles = Find_Cycles(G, config);
	cout << "There are " << cycles.size() << " cycles." << endl;
	vector<cycle_arcs> subsets = cycles;
	for (int i = 0; i < cycles.size(); i++)
	{
		if (cycles[i].vertices.size() < config.cyclelength + config.subset_size) // Only proceed if the max number of vertices not reached.
		{
			queue<cycle_arcs> candidates;
			for (int j = i + 1; j < cycles.size(); j++)
			{
				// If vertices cycle[j] is not a subset of vertices of cycle[i] (and vice versa).
				// Add cycle[j] to candidates
				if (Vertex_Subset_Dominate(cycles[i], cycles[j]) == 0)
				{
					candidates.push(cycles[j]);
				}
			}
			if (candidates.size() > 0)
			{
				Relevant_Subsets_Recursion(subsets, cycles[i], candidates, config.cyclelength + config.subset_size);
			}
			// Go to recursion, return relevant subsets.
			// Input is accepted subsets - current subset - candidate subsets - max subset size.
		}
	}

	// Sort the vertices.
	for (int i = 0; i < subsets.size(); i++)
	{
		sort(subsets[i].vertices.begin(), subsets[i].vertices.end());
	}
	cout << "There are " << subsets.size() << " subsets." << endl;
	// Put the right arcs for each subset.
	Subset_Arcs(subsets, G, config);
	return subsets;

}

void Relevant_Subsets_Recursion(vector<cycle_arcs>& accepted, cycle_arcs current, queue<cycle_arcs> candidates, int max_subset_size)
{
	int nr_candidates = candidates.size();
	for(int i = 0; i < nr_candidates; i++) // At first glance, going until the queue is emptied looks cleaner
											// However, in this way, we can evaluate each candidate once, 
											// while still allowing disjoint candidates to be added back into the queue.
	{
		cycle_arcs parent = current;
		// First check whether combining the two would push the size over the limit.
		// If so, remove the candidate.
		int nr_vertices = Vertex_Sum(parent, candidates.front());
		if (nr_vertices > max_subset_size)
		{
			candidates.pop();
		}
		// Next, check whether parent and candidate are disjoint. If so, the candidate is not removed from candidates, but is moved to the back.
		else if (nr_vertices == parent.vertices.size() + candidates.front().vertices.size())
		{
			candidates.push(candidates.front());
			candidates.pop();
		}
		else // If neither happens, combine the parent and the candidate
			// Add the new subset to the list of accepted subsets (if it is not a duplicate)
			// Go on to find the new set of candidates.
		{
			// If neither of these happen, combine the parent and the candidate.
			Vertex_Combine(parent, candidates.front());
			candidates.pop();
			
			// Add the new current to the set of accepted subsets (if it is not already included.)
			{
				bool included = 0;
				for (int j = 0; j < accepted.size(); j++)
				{

					if (Vertex_Subset_Equal(parent, accepted[j]) == 1)
					{
						included = 1;
						break;
					}
				}
				if (included == 0)
					accepted.push_back(parent);
			}

			// Go through the list of candidates to weed out the dominated ones.
			queue<cycle_arcs> remaining_can = candidates;
			int can_remaining = remaining_can.size();
			for (int j = 0; j < can_remaining; j++)
			{
				if (Vertex_Subset_Dominate(parent, remaining_can.front()) == 0)
				{
					remaining_can.push(remaining_can.front());
					remaining_can.pop();
				}
				else
				{
					remaining_can.pop();
				}
					
			}
			if (remaining_can.size() > 0)
			{
				Relevant_Subsets_Recursion(accepted, parent, remaining_can, max_subset_size);
			}
		}
	}
}

void Subset_Arcs(vector<cycle_arcs>& subsets, const directedgraph & G, const configuration & config)
{
	// This could be improved upon by weeding out the arcs that can not be part of a cycle in the subset.
	
	for (int i = 0; i < subsets.size(); i++)
	{
		subsets[i].arcs.resize(0); // First clean it out.
		for (int j = 0; j < subsets[i].vertices.size(); j++) // For each vertex in the subset
		{
			for (int k = G.first_vertex_arc[subsets[i].vertices[j]]; k < G.first_vertex_arc[subsets[i].vertices[j] + 1]; k++) // Go through the arcs for which it is startvertex
			{
				for (int l = 0; l < subsets[i].vertices.size(); l++) // And for each other vertex in the subset
				{
					if (G.arcs[k].endvertex == subsets[i].vertices[l]) // Check whether it is the endvertex
					{
						subsets[i].arcs.push_back(k); // If it is, add it to the arcs of the subset.
					}
						
				}
			}
		}
	}
}

void Subset_Set_Weights_Arc(vector<cycle_arcs>& subsets, const directedgraph & G, const configuration & config)
{
	int count = 0;
#pragma omp parallel for
	for (int i = 0; i < subsets.size(); i++)
	{
		// Build subset graph
		directedgraph G_sub = Subset_Graph(subsets[i], G);
		// Give Graph to function computing weight for subset
		subsets[i].weight = Subset_Set_Weights_Arc_Root(G_sub, config);
		count++;
		if (count % 100 == 0)
			cout << count << endl;
	}
}

void Subset_Set_Weights_Vertex(vector<cycle_arcs>& subsets, const directedgraph & G, const configuration & config)
{
	int count = 0;
#pragma omp parallel for
	for (int i = 0; i < subsets.size(); i++)
	{
		// Build subset graph
		directedgraph G_sub = Subset_Graph(subsets[i], G);
		// Give Graph to function computing weight for subset
		subsets[i].weight = Subset_Set_Weights_Vertex_Root(G_sub, config);
		count++;
		if (count % 100 == 0)
			cout << count << endl;
	}
}
void Subset_Set_Weights_Vertex_Verbose(vector<cycle_arcs>& subsets, const directedgraph & G, const configuration & config)
{
	int count = 0;
#pragma omp parallel for
	for (int i = 0; i < subsets.size(); i++)
	{
		// Build subset graph
		directedgraph G_sub = Subset_Graph(subsets[i], G);
		// Give Graph to function computing weight for subset
		subsets[i].weight = Subset_Set_Weights_Vertex_Root_Verbose(G_sub, config);
		count++;
		if (count % 100 == 0)
			cout << count << endl;
	}
}

void Subset_Set_Weights_Arc_verbose(vector<cycle_arcs>& subsets, const directedgraph & G, const configuration & config)
{
	int count = 0;
#pragma omp parallel for
	for (int i = 0; i < subsets.size(); i++)
	{
		// Build subset graph
		directedgraph G_sub = Subset_Graph(subsets[i], G);
		// Give Graph to function computing weight for subset
		subsets[i].weight = Subset_Set_Weights_Arc_Root_verbose(G_sub, config);
		count++;
		if (count % 100 == 0)
			cout << count << endl;
	}
}

float Subset_Set_Weights_Arc_Root(const directedgraph & G, const configuration & config)
{
	float weight;
	// We first find all the cycles
	vector<cycle_arcs> cycles = Find_Cycles(G, config);
	// We set up 2 LPs, one for the upper bound (assume all arcs not yet set will succeed) and one for the LB (assume all arcs not yet set will fail).
	IloEnv env;
	IloModel Model(env);
	IloNumVarArray cyclevar(env, cycles.size(), cycles.size(), ILOINT);
	for (int i = 0; i < cycles.size(); i++)
		cyclevar[i] = IloNumVar(env, 0, 1, ILOINT);
	IloNumVarArray arcvar(env, G.arcs.size(), G.arcs.size(), ILOINT);
	for (int i = 0; i < G.arcs.size(); i++)
	{
		arcvar[i] = IloNumVar(env, 0, 1, ILOINT);
	}

	IloRangeArray vertex_cons = Build_Vertex_Constraint_SSWR(env, G, cyclevar, cycles);
	Model.add(vertex_cons);
	// These constraints link whether or not arcs are still under consideration to the cycles.
	IloRangeArray cycle_con = Build_Cycle_Constraint_SSWR(env, G, cyclevar, arcvar, cycles);
	Model.add(cycle_con);

	// Objective function
	IloObjective obj = IloMaximize(env);
	for (int i = 0; i < cycles.size(); i++)
	{
		obj.setLinearCoef(cyclevar[i], cycles[i].weight);
	}
	for (int i = 0; i < G.arcs.size(); i++)
	{
		obj.setLinearCoef(arcvar[i], -0.00001); // Slight negative value so the model only sets needed arcvar to 1.
	}
	Model.add(obj);

	IloCplex CPLEX(Model); CPLEX.setOut(env.getNullStream());
	CPLEX.setWarning(env.getNullStream());
	CPLEX.solve();

	float CPLEX_Solution = ceil(CPLEX.getObjValue()); // Ceiling, to remove the effect of the small negative values of the arcvvars.
	vector<int> arc_sol(G.arcs.size());
	for (int i = 0; i < G.arcs.size(); i++)
	{
		arc_sol[i] = round(CPLEX.getValue(arcvar[i]));
	}

	// We pick one arc used in the current solution and fix it (either fail or success).
	// If success, the number of transplants will not change. We pick a new arc used in the solution and branch again.
	// If failure, we resolve, and pick another arc (in the solution) and branch on it.
	// If there are no arcs used in the solution that are not yet fixed, the branching ends and the current solution is returned.

	vector<bool> fixed(G.arcs.size());
	int arc_fix = 0;
	for (int i = 0; i < G.arcs.size(); i++)
	{
		if (arc_sol[i] == 1 && fixed[i] == 0)
		{
			arc_fix = i;
			fixed[i] = 1;
			break;
		}
	}

	float solution = 0;
	int depth = 0;

	if (CPLEX_Solution > 0)
	{
		// First assume the arc is a succes
		// We go down one level and branch again on another arc used in the solution.
		weight = (1 - G.arcs[arc_fix].failprob)*Subset_Set_Weights_Arc_Recursion(env, G, CPLEX, arcvar, fixed, CPLEX_Solution, arc_sol, 1, depth);
		// Next assume the arc fails.
		arcvar[arc_fix].setUB(0);
		weight = weight + G.arcs[arc_fix].failprob*Subset_Set_Weights_Arc_Recursion(env, G, CPLEX, arcvar, fixed, CPLEX_Solution, arc_sol, 0, depth);
		arcvar[arc_fix].setUB(1);
	}
	env.end();
	return weight;
}
float Subset_Set_Weights_Vertex_Root(const directedgraph & G, const configuration & config)
{
	float weight;
	// We first find all the cycles
	vector<cycle_arcs> cycles = Find_Cycles(G, config);
	// We set up 2 LPs, one for the upper bound (assume all arcs not yet set will succeed) and one for the LB (assume all arcs not yet set will fail).
	IloEnv env;
	IloModel Model(env);
	IloNumVarArray cyclevar(env, cycles.size(), cycles.size(), ILOINT);
	for (int i = 0; i < cycles.size(); i++)
		cyclevar[i] = IloNumVar(env, 0, 1, ILOINT);

	IloRangeArray vertex_cons = Build_Vertex_Constraint_SSWR(env, G, cyclevar, cycles);
	Model.add(vertex_cons);


	// Objective function
	IloObjective obj = IloMaximize(env);
	for (int i = 0; i < cycles.size(); i++)
	{
		obj.setLinearCoef(cyclevar[i], cycles[i].weight);
	}

	Model.add(obj);

	IloCplex CPLEX(Model); CPLEX.setOut(env.getNullStream());
	CPLEX.setWarning(env.getNullStream());
	CPLEX.solve();

	float CPLEX_Solution = ceil(CPLEX.getObjValue());
	vector<int> vertex_use(G.size);
	for (int i = 0; i < G.size; i++)
	{
		vertex_use[i] = 1-round(CPLEX.getSlack(vertex_cons[i]));
	}

	// We pick one vertex used in the current solution and fix it (either fail or success).
	// If success, the number of transplants will not change. We pick a new vertex used in the solution and branch again.
	// If failure, we resolve, and pick another vertex (in the solution) and branch on it.
	// If there are no vertices used in the solution that are not yet fixed, the branching ends and the current solution is returned.

	vector<bool> fixed(G.size);
	int vertex_fix = 0;
	for (int i = 0; i < G.size; i++)
	{
		if (vertex_use[i] == 1 && fixed[i] == 0)
		{
			vertex_fix = i;
			fixed[i] = 1;
			break;
		}
	}

	float solution = 0;
	int depth = 0;
	if (CPLEX_Solution > 0)
	{
		// First assume the vertex is a succes
		// We go down one level and branch again on another arc used in the solution.
		weight = (1 - G.pairs[vertex_fix].failprob)*Subset_Set_Weights_Vertex_Recursion(env, G, CPLEX, vertex_use, vertex_cons, fixed, CPLEX_Solution, 1, depth);
		// Next assume the vertex fails.
		vertex_cons[vertex_fix].setUB(0);
		weight = weight + G.pairs[vertex_fix].failprob*Subset_Set_Weights_Vertex_Recursion(env, G, CPLEX, vertex_use, vertex_cons, fixed, CPLEX_Solution, 0, depth);
		vertex_cons[vertex_fix].setUB(1);
	}
	env.end();
	return weight;
}
float Subset_Set_Weights_Arc_Recursion(IloEnv & env, const directedgraph & G, IloCplex & CPLEX, IloNumVarArray & arcvar, vector<bool> fixed, float CPLEX_Solution, vector<int> arc_sol, bool succes_fix, int depth)
{
	depth++;
	
	float weight = 0;

	if (succes_fix == 0)
	{
		CPLEX.solve();
		CPLEX_Solution = ceil(CPLEX.getObjValue());
		for (int i = 0; i < G.arcs.size(); i++)
		{
			arc_sol[i] = round(CPLEX.getValue(arcvar[i]));
		}
	}
	int arc_fix = 0;
	bool new_fix = 0;
	for (int i = 0; i < G.arcs.size(); i++)
	{
		if (arc_sol[i] == 1 && fixed[i] == 0)
		{
			arc_fix = i;
			fixed[i] = 1;
			new_fix = 1;
			break;
		}
	}
	if (new_fix == 1)
	{
		weight = (1 - G.arcs[arc_fix].failprob)*Subset_Set_Weights_Arc_Recursion(env, G, CPLEX, arcvar, fixed, CPLEX_Solution, arc_sol, 1, depth);
		// Next assume the arc fails.
		arcvar[arc_fix].setUB(0);
		weight = weight + G.arcs[arc_fix].failprob*Subset_Set_Weights_Arc_Recursion(env, G, CPLEX, arcvar, fixed, CPLEX_Solution, arc_sol, 0, depth);
		arcvar[arc_fix].setUB(1);
	}
	else // If there are no more arcs to fix, the weight is equal to the last CPLEX Solution (all arcs used in this solution have been fixed to success).
	{
		weight = CPLEX_Solution;
	}

	return weight;
}
float Subset_Set_Weights_Vertex_Recursion(IloEnv & env, const directedgraph & G, IloCplex & CPLEX, vector<int> vertex_use, IloRangeArray & vertex_cons, vector<bool> fixed, float Solution, bool succes_fix, int depth)
{
	depth++;
	float weight = 0;
	if (succes_fix == 0)
	{
		CPLEX.solve();
		Solution = ceil(CPLEX.getObjValue());
		for (int i = 0; i < G.size; i++)
		{
			vertex_use[i] = 1- round(CPLEX.getSlack(vertex_cons[i]));
		}
	}
	int vertex_fix = 0;
	bool new_fix = 0;
	for (int i = 0; i < G.size; i++)
	{
		if (vertex_use[i] == 1 && fixed[i] == 0)
		{
			vertex_fix = i;
			fixed[i] = 1;
			new_fix = 1;
			break;
		}
	}
	if (new_fix == 1)
	{
		// First assume the vertex is a succes
		// We go down one level and branch again on another arc used in the solution.
		weight = (1 - G.pairs[vertex_fix].failprob)*Subset_Set_Weights_Vertex_Recursion(env, G, CPLEX, vertex_use, vertex_cons, fixed, Solution, 1, depth);
// Next assume the vertex fails.
		vertex_cons[vertex_fix].setUB(0);
		weight = weight + G.pairs[vertex_fix].failprob*Subset_Set_Weights_Vertex_Recursion(env, G, CPLEX, vertex_use, vertex_cons, fixed, Solution, 0, depth);
		vertex_cons[vertex_fix].setUB(1);
	}
	else // If there are no more arcs to fix, the weight is equal to the last CPLEX Solution (all arcs used in this solution have been fixed to success).
	{
		weight = Solution;
	}

	return weight;
}

float Subset_Set_Weights_Arc_Root_verbose(const directedgraph & G, const configuration & config)
{
	float weight;
	// We first find all the cycles
	vector<cycle_arcs> cycles = Find_Cycles(G, config);
	// We set up 2 LPs, one for the upper bound (assume all arcs not yet set will succeed) and one for the LB (assume all arcs not yet set will fail).
	IloEnv env;
	IloModel Model(env);
	IloNumVarArray cyclevar(env, cycles.size(), cycles.size(), ILOINT);
	for (int i = 0; i < cycles.size(); i++)
		cyclevar[i] = IloNumVar(env, 0, 1, ILOINT);
	IloNumVarArray arcvar(env, G.arcs.size(), G.arcs.size(), ILOINT);
	for (int i = 0; i < G.arcs.size(); i++)
	{
		arcvar[i] = IloNumVar(env, 0, 1, ILOINT);
	}

	IloRangeArray vertex_cons = Build_Vertex_Constraint_SSWR(env, G, cyclevar, cycles);
	Model.add(vertex_cons);
	// These constraints link whether or not arcs are still under consideration to the cycles.
	IloRangeArray cycle_con = Build_Cycle_Constraint_SSWR(env, G, cyclevar, arcvar, cycles);
	Model.add(cycle_con);

	// Objective function
	IloObjective obj = IloMaximize(env);
	for (int i = 0; i < cycles.size(); i++)
	{
		obj.setLinearCoef(cyclevar[i], cycles[i].weight);
	}
	for (int i = 0; i < G.arcs.size(); i++)
	{
		obj.setLinearCoef(arcvar[i], -0.00001); // Slight negative value so the model only sets needed arcvar to 1.
	}
	Model.add(obj);

	IloCplex CPLEX(Model); CPLEX.setOut(env.getNullStream());
	CPLEX.setWarning(env.getNullStream());
	CPLEX.solve();

	float CPLEX_Solution = ceil(CPLEX.getObjValue()); // Ceiling, to remove the effect of the small negative values of the arcvvars.
	vector<int> arc_sol(G.arcs.size());
	for (int i = 0; i < G.arcs.size(); i++)
	{
		arc_sol[i] = round(CPLEX.getValue(arcvar[i]));
	}

	// We pick one arc used in the current solution and fix it (either fail or success).
	// If success, the number of transplants will not change. We pick a new arc used in the solution and branch again.
	// If failure, we resolve, and pick another arc (in the solution) and branch on it.
	// If there are no arcs used in the solution that are not yet fixed, the branching ends and the current solution is returned.

	vector<bool> fixed(G.arcs.size());
	deque<bool> suc_fail;
	int arc_fix = 0;
	for (int i = 0; i < G.arcs.size(); i++)
	{
		if (arc_sol[i] == 1 && fixed[i] == 0)
		{
			arc_fix = i;
			fixed[i] = 1;
			break;
		}
	}

	float solution = 0;
	int depth = 0;
	if (CPLEX_Solution > 0)
	{
		// First assume the arc is a succes
		// We go down one level and branch again on another arc used in the solution.
		suc_fail.push_back(1);
		weight = (1 - G.arcs[arc_fix].failprob)*Subset_Set_Weights_Arc_Recursion_verbose(env, G, CPLEX, arcvar, fixed, CPLEX_Solution, arc_sol, 1, depth, suc_fail);
		suc_fail.pop_back();
		suc_fail.push_back(0);
		// Next assume the arc fails.
		arcvar[arc_fix].setUB(0);
		weight = weight + G.arcs[arc_fix].failprob*Subset_Set_Weights_Arc_Recursion_verbose(env, G, CPLEX, arcvar, fixed, CPLEX_Solution, arc_sol, 0, depth, suc_fail);
		suc_fail.pop_back();
		arcvar[arc_fix].setUB(1);
	}
	env.end();
	return weight;
}


bool Is_candidate(int can_id, const vector<int> & vertex_use, const vector<bool> & fixed) {
	return (vertex_use[can_id] == 1 && fixed[can_id] == 0);
}

void Find_next_branch(vector<tuple<int,int>> & vertex_fix, int pos, const vector<int> & vertex_use, vector<bool> & fixed, 
					  list<tuple<int, int>> & branching_order, bool & new_fix) {
	int can_id;

	while (!branching_order.empty()) {
		can_id = get<0>(branching_order.front());

		if (Is_candidate(can_id, vertex_use, fixed)) {
			vertex_fix[pos] = branching_order.front();
			fixed[can_id] = 1;
			new_fix = 1;
			branching_order.pop_front();
			return;
		}

		branching_order.pop_front(); 
	}
	new_fix = 0;
	return;
}


list<tuple<int,int>> recompute_exact_order(const directedgraph & G, const configuration & config, const vector<int> & vertex_use) {
	
	vector<cycle_arcs> cycles = Find_Cycles_Subgraph(G, config, vertex_use);

	vector<int> nr_involved_cycles (G.nr_pairs, 0); // Vertex i is involved in nr_involved_cycles[i] cycles
	for (unsigned int i = 0; i < cycles.size(); ++i) {
		for (unsigned j = 0; j < cycles[i].vertices.size(); ++j) nr_involved_cycles[cycles[i].vertices[j]]++;
	}

	int max_involved_cycles = *(max_element(nr_involved_cycles.begin(), nr_involved_cycles.end())); // Upper bound
	vector<list<tuple<int,int>>> nr_appearance_lists (max_involved_cycles + 1); // nr_appearance_lists[i] is the list of vertex involved in exactly i cycles. 

	for (unsigned int i = 0; i < nr_involved_cycles.size(); ++i) {
		nr_appearance_lists[nr_involved_cycles[i]].push_back(tuple<int,int>(i, nr_involved_cycles[i]));
	}

	list<tuple<int,int>> branching_order;
	for (unsigned int i = 0; i < nr_appearance_lists.size(); ++i) {
		if (!nr_appearance_lists[i].empty()) branching_order.splice(branching_order.begin(), nr_appearance_lists[i]);
	}

	return branching_order;
}

float Subset_Set_Weights_Vertex_Root_Verbose(const directedgraph & G, const configuration & config)
{
	float weight;
	// We first find all the cycles
	vector<cycle_arcs> cycles = Find_Cycles(G, config);

	// Cycles info

	cout << "All cycles:" << endl;
	for (unsigned int i = 0; i < cycles.size(); ++i) {
		for (unsigned int j = 0; j < cycles[i].vertices.size(); ++j) cout << cycles[i].vertices[j] << " ";
		cout << endl;
	}

	// We compute the branching order: first the vertices that are involved in the highest number of cyles

	vector<int> nr_involved_cycles (G.nr_pairs, 0); // Vertex i is involved in nr_involved_cycles[i] cycles
	for (unsigned int i = 0; i < cycles.size(); ++i) {
		for (unsigned j = 0; j < cycles[i].vertices.size(); ++j) nr_involved_cycles[cycles[i].vertices[j]]++;
	}

	int max_involved_cycles = *(max_element(nr_involved_cycles.begin(), nr_involved_cycles.end())); // Upper bound
	vector<list<tuple<int,int>>> nr_appearance_lists (max_involved_cycles + 1); // nr_appearance_lists[i] is the list of vertex involved in exactly i cycles. 

	for (unsigned int i = 0; i < nr_involved_cycles.size(); ++i) {
		nr_appearance_lists[nr_involved_cycles[i]].push_back(tuple<int,int>(i, nr_involved_cycles[i]));
	}

	list<tuple<int,int>> branching_order;
	for (unsigned int i = 0; i < nr_appearance_lists.size(); ++i) {
		if (!nr_appearance_lists[i].empty()) branching_order.splice(branching_order.begin(), nr_appearance_lists[i]);
	}

	// Random order for testing purposes
	/*
	vector<int> random_order (G.size, 0);
	list<tuple<int, int>> branching_order;
	for (int i = 0; i < random_order.size(); ++i) random_order[i] = i;
	srand(time(0));
	random_shuffle(random_order.begin(), random_order.end());
	for (int i = 0; i < random_order.size(); ++i) branching_order.push_back(tuple<int,int>(random_order[i], nr_involved_cycles[random_order[i]]));
	*/
	// We set up 2 LPs, one for the upper bound (assume all arcs not yet set will succeed) and one for the LB (assume all arcs not yet set will fail).
	IloEnv env;
	IloModel Model(env);
	IloNumVarArray cyclevar(env, cycles.size(), cycles.size(), ILOINT);
	for (int i = 0; i < cycles.size(); i++)
		cyclevar[i] = IloNumVar(env, 0, 1, ILOINT);

	IloRangeArray vertex_cons = Build_Vertex_Constraint_SSWR(env, G, cyclevar, cycles);
	Model.add(vertex_cons);


	// Objective function
	IloObjective obj = IloMaximize(env);
	for (int i = 0; i < cycles.size(); i++)
	{
		obj.setLinearCoef(cyclevar[i], cycles[i].weight);
	}

	Model.add(obj);
	cout << "Model Built" << endl;
	IloCplex CPLEX(Model); CPLEX.setOut(env.getNullStream());
	CPLEX.setWarning(env.getNullStream());
	CPLEX.solve();

	float CPLEX_Solution = ceil(CPLEX.getObjValue());
	vector<int> vertex_use(G.size);
	for (int i = 0; i < G.size; i++)
	{
		vertex_use[i] = 1 - round(CPLEX.getSlack(vertex_cons[i]));
	}

	// We pick one vertex used in the current solution and fix it (either fail or success).
	// If success, the number of transplants will not change. We pick a new vertex used in the solution and branch again.
	// If failure, we resolve, and pick another vertex (in the solution) and branch on it.
	// If there are no vertices used in the solution that are not yet fixed, the branching ends and the current solution is returned.

	vector<bool> fixed(G.size, 0);
	vector<tuple<int,int>> vertex_fix (2, tuple<int,int>(-1,-1));


	// Order only with number of outgoing arcs.
	/*
	list<tuple<int,int>> branching_order;
	vector<int> entering_vertex;
	for (unsigned int i = 0; i < vertex_use.size(); ++i) {
		if (vertex_use[i] == 1) entering_vertex.push_back(i);
	}
	update_order(G, branching_order, vertex_use, entering_vertex, fixed);
	branching_order.sort(tuple_compare); 
	*/

	// No branching order, just consider the first unfixed vertex found
	/*
	for (int i = 0; i < G.size; i++)
	{
		if (vertex_use[i] == 1 && fixed[i] == 0)
		{
			vertex_fix = i;
			fixed[i] = 1;
			break;
		}
	}
	*/


	// Test on cut vertices.
	/*if (G.size == 22) {
		cout << "Graph broke" << endl;
		branching_order.push_front(tuple<int,int>(2, nr_involved_cycles[2]));
		branching_order.push_front(tuple<int,int>(0, nr_involved_cycles[0]));
		branching_order.push_front(tuple<int,int>(8, nr_involved_cycles[8]));
		branching_order.push_front(tuple<int,int>(17, nr_involved_cycles[17]));
		branching_order.push_front(tuple<int,int>(7, nr_involved_cycles[7]));
	
		
		branching_order.push_front(tuple<int,int>(15, nr_involved_cycles[15]));
		branching_order.push_front(tuple<int,int>(14, nr_involved_cycles[14]));
		branching_order.push_front(tuple<int,int>(16, nr_involved_cycles[16]));
		branching_order.push_front(tuple<int,int>(21, nr_involved_cycles[21]));
		branching_order.push_front(tuple<int,int>(12, nr_involved_cycles[12]));
	}*/

	cout << "Vertex use info : ";
	for (int i = 0; i < vertex_use.size(); ++i) {
		if (vertex_use[i] == 1) cout << i << " ";
	}
	cout << endl;
	cout << "Branching order info : ";
	for (list<tuple<int,int>>::iterator it = branching_order.begin(); it != branching_order.end(); ++it) cout << get<0>(*it) << ":" << get<1>(*it) << "    ";
	cout << endl;

	bool new_fix = 0;
	bool branch_on_two = 0;

	for(list<tuple<int,int>>::iterator it = branching_order.begin(); it != branching_order.end(); ++it) {
		int id = get<0>(*it);
		if(vertex_use[id] == 1 && fixed[id] == 0) {
			new_fix = 1;
			fixed[id] = 1;
			vertex_fix[0] = *it;
			break;
		}
	}

	// We look for a first candidate
	//Find_next_branch(vertex_fix, 0, vertex_use, fixed, branching_order, new_fix);

	// If success we look for a second one
	//if (new_fix == 1) Find_next_branch(vertex_fix, 1, vertex_use, fixed, branching_order, branch_on_two);

	/*
	for (int i = 0; i < G.size; i++)
	{
		if (vertex_use[i] == 1 && fixed[i] == 0)
		{
			get<0>(vertex_fix[0]) = i;
			fixed[i] = 1;
			break;
		}
	}*/

	cout << "Model Solved, Branching starting" << endl;
	float solution = 0;
	int depth = 0;
	deque<int> suc_fail;
	int nbBranch = 0;

	// TEST: BRANCH ON SEVERAL VERTICES AT ONCE.
	if (CPLEX_Solution > 0 && branch_on_two == 1) {
		depth += 2;
		cout << "Branching on two vertices" << endl;

		// First both success
		suc_fail.push_back(1);
		suc_fail.push_back(1);
		nbBranch++;
		weight = (1 - G.pairs[get<0>(vertex_fix[0])].failprob)*(1 - G.pairs[get<0>(vertex_fix[1])].failprob)*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, CPLEX_Solution, 1, depth, suc_fail, branching_order, nbBranch);
		suc_fail.pop_back();
		suc_fail.pop_back();

		// Next fail success
		vertex_cons[get<0>(vertex_fix[0])].setUB(0);
		nbBranch++;
		suc_fail.push_back(1);
		suc_fail.push_back(0);
		weight = weight + G.pairs[get<0>(vertex_fix[0])].failprob*(1-G.pairs[get<0>(vertex_fix[1])].failprob)*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, CPLEX_Solution, 0, depth, suc_fail, branching_order, nbBranch);
		suc_fail.pop_back();
		suc_fail.pop_back();

		// Next success fail
		vertex_cons[get<0>(vertex_fix[0])].setUB(1);
		vertex_cons[get<0>(vertex_fix[1])].setUB(0);		
		suc_fail.push_back(0);
		suc_fail.push_back(1);
		nbBranch++;
		weight = weight + G.pairs[get<0>(vertex_fix[1])].failprob*(1-G.pairs[get<0>(vertex_fix[0])].failprob)*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, CPLEX_Solution, 0, depth, suc_fail, branching_order, nbBranch);
		suc_fail.pop_back();
		suc_fail.pop_back();

		// Last fail fail
		vertex_cons[get<0>(vertex_fix[0])].setUB(0);
		suc_fail.push_back(0);
		suc_fail.push_back(0);
		nbBranch++;
		weight = weight + G.pairs[get<0>(vertex_fix[0])].failprob*G.pairs[get<0>(vertex_fix[1])].failprob*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, CPLEX_Solution, 0, depth, suc_fail, branching_order, nbBranch);
		suc_fail.pop_back();
		suc_fail.pop_back();
		vertex_cons[get<0>(vertex_fix[0])].setUB(1);
		vertex_cons[get<0>(vertex_fix[1])].setUB(1);
		
	} else if (CPLEX_Solution > 0 && branch_on_two == 0) {
		depth++;
		cout << "Branching on one vertex" << endl;

		// First assume the vertex is a succes
		// We go down one level and branch again on another arc used in the solution.
		suc_fail.push_back(1);
		nbBranch++;
		weight = (1 - G.pairs[get<0>(vertex_fix[0])].failprob)*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, CPLEX_Solution, 1, depth, suc_fail, branching_order, nbBranch);
		suc_fail.pop_back();
		suc_fail.push_back(0);

		// Next assume the vertex fails.
		vertex_cons[get<0>(vertex_fix[0])].setUB(0);
		nbBranch++;
		weight = weight + G.pairs[get<0>(vertex_fix[0])].failprob*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, CPLEX_Solution, 0, depth, suc_fail, branching_order, nbBranch);
		suc_fail.pop_back();
		vertex_cons[get<0>(vertex_fix[0])].setUB(1);
		
	}
	env.end();
	cout << "Number of branching done: " << nbBranch;
	return weight;
}

float Subset_Set_Weights_Arc_Recursion_verbose(IloEnv & env, const directedgraph & G, IloCplex & CPLEX, IloNumVarArray & arcvar, vector<bool> fixed, float CPLEX_Solution, vector<int> arc_sol, bool succes_fix, int depth, deque<bool> & suc_fail)
{
	depth++;

	/*if (depth <15)
	{
		for (int i = 0; i < suc_fail.size(); i++)
			cout << suc_fail[i] << ",";
		cout << endl;
	}*/
	float weight = 0;

	if (succes_fix == 0)
	{
		CPLEX.solve();
		CPLEX_Solution = ceil(CPLEX.getObjValue());
		for (int i = 0; i < G.arcs.size(); i++)
		{
			arc_sol[i] = round(CPLEX.getValue(arcvar[i]));
		}
	}
	int arc_fix = 0;
	bool new_fix = 0;
	for (int i = 0; i < G.arcs.size(); i++)
	{
		if (arc_sol[i] == 1 && fixed[i] == 0)
		{
			arc_fix = i;
			fixed[i] = 1;
			new_fix = 1;
			break;
		}
	}
	if (new_fix == 1)
	{
		// First assume the arc is a succes
		// We go down one level and branch again on another arc used in the solution.
		suc_fail.push_back(1);
		weight = (1 - G.arcs[arc_fix].failprob)*Subset_Set_Weights_Arc_Recursion_verbose(env, G, CPLEX, arcvar, fixed, CPLEX_Solution, arc_sol, 1, depth, suc_fail);
		suc_fail.pop_back();
		suc_fail.push_back(0);
		// Next assume the arc fails.
		arcvar[arc_fix].setUB(0);
		weight = weight + G.arcs[arc_fix].failprob*Subset_Set_Weights_Arc_Recursion_verbose(env, G, CPLEX, arcvar, fixed, CPLEX_Solution, arc_sol, 0, depth, suc_fail);
		suc_fail.pop_back();
		arcvar[arc_fix].setUB(1);
	}
	else // If there are no more arcs to fix, the weight is equal to the last CPLEX Solution (all arcs used in this solution have been fixed to success).
	{
		weight = CPLEX_Solution;
	}

	return weight;
}


bool tuple_compare(tuple<int, int> const & t1, tuple<int,int> const & t2) {
	return get<1>(t1) >= get<1>(t2);
}

list<tuple<int,int>> recompute_order(const directedgraph & G, const vector<int> & vertex_use, const vector<bool> & fixed) {
	vector<int> nr_outgoing_arcs (G.nr_pairs, 0);
	for (unsigned int i = 0; i < G.arcs.size(); ++i) {
		int sv = G.arcs[i].startvertex;
		int ev = G.arcs[i].endvertex;
		if (vertex_use[ev] == 1 && vertex_use[sv] == 1) {
			nr_outgoing_arcs[ev]++;
			nr_outgoing_arcs[sv]++;
		}
	}
	list<tuple<int,int>> order;
	for (unsigned int i = 0; i < G.nr_pairs; ++i) {
		if (vertex_use[i] == 1 && fixed[i] == 0) order.push_back(tuple<int,int>(i, nr_outgoing_arcs[i]));
	}
	order.sort(tuple_compare);
	return order;
}

list<tuple<int,int>> recompute_random_order(const directedgraph & G, const vector<int> & vertex_use, const vector<bool> & fixed) {
	list<tuple<int,int>> order;
	vector<int> random_order(G.size, 0);
	for (int i = 0; i < random_order.size(); ++i) random_order[i] = i;
	random_shuffle(random_order.begin(), random_order.end());
	for (int i = 0; i < random_order.size(); ++i) {
		if (vertex_use[random_order[i]] == 1 && fixed[random_order[i]] == 0) {
			order.push_back(tuple<int,int>(random_order[i], 0));
		}
	}
	return order;
}

void update_order(const directedgraph & G, list<tuple<int,int>> & branching_order, const vector<int> & vertex_use, const vector<int> & entering_vertices, const vector<bool> & fixed) {
	vector<int> nr_outgoing_arcs (G.nr_pairs, 0);
	for (unsigned int i = 0; i < G.arcs.size(); ++i) {
		int sv = G.arcs[i].startvertex;
		int ev = G.arcs[i].endvertex;
		if (vertex_use[ev] == 1 && vertex_use[sv] == 1) {
			nr_outgoing_arcs[ev]++;
			nr_outgoing_arcs[sv]++;
		}
	}
	int current = -1;
	int current_nr_arcs = -1;

	for(list<tuple<int,int>>::iterator it = branching_order.begin(); it != branching_order.end(); ++it) {
		get<1>(*it) = nr_outgoing_arcs[get<0>(*it)];
	}

	for (unsigned int i = 0; i < entering_vertices.size(); ++i) {
		current = entering_vertices[i];
		current_nr_arcs = nr_outgoing_arcs[current];
		tuple<int, int> entering_vertex(current, current_nr_arcs);
		if (fixed[current] != 1) {
			if (branching_order.empty()) branching_order.push_front(entering_vertex);
			else if (current_nr_arcs >= get<1>(branching_order.front())) branching_order.push_front(entering_vertex);
			else branching_order.push_back(entering_vertex); 
		}
	}
	//branching_order.sort(tuple_compare); //slow down the algorithm.
}

float Subset_Set_Weights_Vertex_Recursion_Verbose(const configuration & config, IloEnv & env, const directedgraph & G, IloCplex & CPLEX, vector<int> vertex_use, IloRangeArray & vertex_cons, vector<bool> fixed, float Solution, bool succes_fix, int depth, deque<int> & suc_fail, list<tuple<int, int>> branching_order, int & nbBranch)
{
	// Branching test
	// depth++;
	/*
	if (depth <15)
	{
		for (int i = 0; i < suc_fail.size(); i++)
			cout << suc_fail[i] << ",";
		cout << endl;
	}*/

	float weight = 0;
	if (succes_fix == 0)
	{
		//int old_vertex_use = 0;
		CPLEX.solve();
		Solution = ceil(CPLEX.getObjValue());

		//vector<int> entering_vertices;
		for (int i = 0; i < G.size; i++) 
		{
			//old_vertex_use = vertex_use[i];
			vertex_use[i] = 1 - round(CPLEX.getSlack(vertex_cons[i]));
			//if (old_vertex_use != 1 && vertex_use[i] == 1) {
		//		entering_vertices.push_back(i);
				//branching_order.push_back(tuple<int,int>(i,0));
		//	}
		}

		//if (depth % 1 == 0) {
		//	branching_order = recompute_exact_order(G, config, vertex_use);
		//} else {
		//update_order(G, branching_order, vertex_use, entering_vertices, fixed);
		//}

		//branching_order = recompute_random_order(G, vertex_use, fixed);
		
		/*if (G.size != 22) { 
			update_order(G, branching_order, vertex_use, entering_vertices);
		} else {
			if (depth >= 4) update_order(G, branching_order, vertex_use, entering_vertices);
		}*/
	}


	vector<tuple<int,int>> vertex_fix (2, tuple<int,int>(-1,-1));
	bool new_fix = 0;
	bool branch_on_two = 0;

	// No branching order, just take the first unfixed vertex that you find
	/*
	for (int i = 0; i < G.size; i++)
	{
		if (vertex_use[i] == 1 && fixed[i] == 0)
		{
			get<0>(vertex_fix[0]) = i;
			fixed[i] = 1;
			new_fix = 1;
			break;
		}
	}*/

	for(list<tuple<int,int>>::iterator it = branching_order.begin(); it != branching_order.end(); ++it) {
		int id = get<0>(*it);
		if(vertex_use[id] == 1 && fixed[id] == 0) {
			new_fix = 1;
			fixed[id] = 1;
			vertex_fix[0] = *it;
			break;
		}
	}


	// We look for the first candidate
	//Find_next_branch(vertex_fix, 0, vertex_use, fixed, branching_order, new_fix);
	
	// If success we look for a second one
	// if (new_fix == 1) Find_next_branch(vertex_fix, 1, vertex_use, fixed, branching_order, branch_on_two);

	if (new_fix == 1 && branch_on_two == 0)
	{
		depth++;
		// First assume the vertex is a success
		// We go down one level and branch again on another arc used in the solution.
		suc_fail.push_back(1);
		nbBranch++; 
		weight = (1 - G.pairs[get<0>(vertex_fix[0])].failprob)*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, Solution, 1, depth, suc_fail, branching_order, nbBranch);
		suc_fail.pop_back();


		// Next assume the vertex fails.
		nbBranch++;
		suc_fail.push_back(0);
		vertex_cons[get<0>(vertex_fix[0])].setUB(0);
		weight = weight + G.pairs[get<0>(vertex_fix[0])].failprob*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, Solution, 0, depth, suc_fail, branching_order, nbBranch);
		suc_fail.pop_back();
		vertex_cons[get<0>(vertex_fix[0])].setUB(1);
	}
	else if (new_fix == 1 && branch_on_two == 1) 
	{
		depth += 2;
		// First both success
		suc_fail.push_back(1);
		suc_fail.push_back(1);
		nbBranch++;
		weight = (1 - G.pairs[get<0>(vertex_fix[0])].failprob)*(1 - G.pairs[get<0>(vertex_fix[1])].failprob)*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, Solution, 1, depth, suc_fail, branching_order, nbBranch);
		suc_fail.pop_back();
		suc_fail.pop_back();

		// Next fail success
		vertex_cons[get<0>(vertex_fix[0])].setUB(0);
		nbBranch++;
		suc_fail.push_back(1);
		suc_fail.push_back(0);
		weight = weight + G.pairs[get<0>(vertex_fix[0])].failprob*(1-G.pairs[get<0>(vertex_fix[1])].failprob)*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, Solution, 0, depth, suc_fail, branching_order, nbBranch);
		suc_fail.pop_back();
		suc_fail.pop_back();

		// Next success fail
		vertex_cons[get<0>(vertex_fix[0])].setUB(1);
		vertex_cons[get<0>(vertex_fix[1])].setUB(0);		
		suc_fail.push_back(0);
		suc_fail.push_back(1);
		nbBranch++;
		weight = weight + G.pairs[get<0>(vertex_fix[1])].failprob*(1-G.pairs[get<0>(vertex_fix[0])].failprob)*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, Solution, 0, depth, suc_fail, branching_order, nbBranch);
		suc_fail.pop_back();
		suc_fail.pop_back();

		// Last fail fail
		vertex_cons[get<0>(vertex_fix[0])].setUB(0);
		suc_fail.push_back(0);
		suc_fail.push_back(0);
		nbBranch++;
		weight = weight + G.pairs[get<0>(vertex_fix[0])].failprob*G.pairs[get<0>(vertex_fix[1])].failprob*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, Solution, 0, depth, suc_fail, branching_order, nbBranch);
		suc_fail.pop_back();
		suc_fail.pop_back();
		vertex_cons[get<0>(vertex_fix[0])].setUB(1);
		vertex_cons[get<0>(vertex_fix[1])].setUB(1);
	}
	else // If there are no more arcs to fix, the weight is equal to the last CPLEX Solution (all arcs used in this solution have been fixed to success).
	{
		weight = Solution;
	}

	return weight;
}

IloRangeArray Build_Vertex_Constraint_SSWR(IloEnv & env, const directedgraph & G, IloNumVarArray & cyclevar, const vector<cycle_arcs> & cycles)
{
	IloRangeArray Vertex_Constraint(env, G.nr_pairs);
	for (int i = 0; i < G.nr_pairs; i++)
	{
		IloExpr expr(env);
		Vertex_Constraint[i] = IloRange(expr <= 1);
	}
	for (int i = 0; i < cycles.size(); i++)
	{
		for (int j = 0; j < cycles[i].arcs.size(); j++)
			Vertex_Constraint[G.arcs[cycles[i].arcs[j]].startvertex].setLinearCoef(cyclevar[i], 1);
	}
	
	
	return Vertex_Constraint;
}

IloRangeArray Build_Cycle_Constraint_SSWR(IloEnv & env, const directedgraph & G, IloNumVarArray & cyclevar, IloNumVarArray & arcvar, const vector<cycle_arcs>& cycles)
{
	IloRangeArray Cycle_Constraint(env, cycles.size());
	for (int i = 0; i < cycles.size(); i++)
	{
		IloExpr expr(env);
		int size = cycles[i].arcs.size();
		Cycle_Constraint[i] = IloRange(size*cyclevar[i] <= 0);
		for (int j = 0; j < cycles[i].arcs.size(); j++)
		{
			Cycle_Constraint[i].setLinearCoef(arcvar[cycles[i].arcs[j]], -1);
		}
	}
	
	return Cycle_Constraint;
}

directedgraph Subset_Graph(const cycle_arcs & subset, const directedgraph & G)
{
	directedgraph G_sub;
	G_sub.size = subset.vertices.size();
	G_sub.nr_pairs = G_sub.size;
	vector<int> vertex_switch(G.nr_pairs); // For each vertex in the original graph, we save in which position it occurs, in the list of vertices in the subset.
	for (int i = 0; i < subset.vertices.size(); i++)
	{
		vertex_switch[subset.vertices[i]] = i;
		G_sub.pairs.push_back(G.pairs[subset.vertices[i]]);
	}
	if (subset.vertices.size() == 22) {
		cout << "Vertex 31 switch to " << vertex_switch[31] << endl;
		cout << "Vertex 17 switch to " << vertex_switch[17] << endl;
		cout << "Vertex 45 switch to " << vertex_switch[45] << endl;
		cout << "Vertex 1 switch to " << vertex_switch[1] << endl;
		cout << "Vertex 38 switch to " << vertex_switch[38] << endl;
		cout << "Vertex 29 switch to " << vertex_switch[29] << endl;
		cout << "Vertex 49 switch to " << vertex_switch[49] << endl;
		cout << "Vertex 5 switch to " << vertex_switch[5] << endl;
		cout << "Vertex 19 switch to " << vertex_switch[19] << endl;
		cout << "Vertex 4 switch to " << vertex_switch[4] << endl;	
	}
	for (int i = 0; i < subset.arcs.size(); i++)
	{
		G_sub.arcs.push_back(G.arcs[subset.arcs[i]]);
		G_sub.arcs[i].startvertex = vertex_switch[G_sub.arcs[i].startvertex];
		G_sub.arcs[i].endvertex = vertex_switch[G_sub.arcs[i].endvertex];
	}

	sort(G_sub.arcs.begin(), G_sub.arcs.end(), arcsort);
	Find_First_arc(G_sub);
	return G_sub;
}

void Vertex_Combine(cycle_arcs & current, const cycle_arcs & candidate)
{
	for (int i = 0; i < candidate.vertices.size(); i++)
	{
		bool notinparent = 1;
		for (int j = 0; j < current.vertices.size(); j++)
		{
			if (candidate.vertices[i] == current.vertices[j])
			{
				notinparent = 0;
				break;
			}
		}
		if (notinparent == 1)
			current.vertices.push_back(candidate.vertices[i]);
	}
}

bool Vertex_Subset_Dominate(const cycle_arcs & parent, const cycle_arcs & candidate)
{
	// Functions returns 1 if candidate is a subset of parent or vice versa, 0 otherwise.
	bool CanisSubset = 1;
	for (int i = 0; i < candidate.vertices.size(); i++)
	{
		bool inParent = 0;
		for (int j = 0; j < parent.vertices.size(); j++)
		{
			if (candidate.vertices[i] == parent.vertices[j])
			{
				inParent = 1;
				break;
			}
		}
		if (inParent == 0)
		{
			CanisSubset = 0;
			break;
		}
	}

	bool ParisSubset = 1;
	for (int i = 0; i < parent.vertices.size(); i++)
	{
		bool inCan = 0;
		for (int j = 0; j < candidate.vertices.size(); j++)
		{
			if (parent.vertices[i] == candidate.vertices[j])
			{
				inCan = 1;
				break;
			}
		}
		if (inCan == 0)
		{
			ParisSubset = 0;
			break;
		}
	}

	return CanisSubset + ParisSubset;
}

bool Vertex_Subset_Equal(const cycle_arcs & parent, const cycle_arcs & candidate)
{
	// Functions returns 1 if candidate and parent are equal.
	bool CanisSubset = 1;
	for (int i = 0; i < candidate.vertices.size(); i++)
	{
		bool inParent = 0;
		for (int j = 0; j < parent.vertices.size(); j++)
		{
			if (candidate.vertices[i] == parent.vertices[j])
			{
				inParent = 1;
				break;
			}
		}
		if (inParent == 0)
		{
			CanisSubset = 0;
			break;
		}
	}

	bool ParisSubset = 1;
	if (CanisSubset == 1);
	{
		for (int i = 0; i < parent.vertices.size(); i++)
		{
			bool inCan = 0;
			for (int j = 0; j < candidate.vertices.size(); j++)
			{
				if (parent.vertices[i] == candidate.vertices[j])
				{
					inCan = 1;
					break;
				}
			}
			if (inCan == 0)
			{
				ParisSubset = 0;
				break;
			}
		}
	}
	bool equal = 0;
	if (CanisSubset == 1 && ParisSubset == 1)
		equal = 1;
	return equal;
}

bool Vertex_Disjoint(const cycle_arcs & parent, const cycle_arcs & candidate)
{
	bool disjoint = 1;
	for (int i = 0; i < candidate.vertices.size(); i++)
	{
		for (int j = 0; j < parent.vertices.size(); j++)
		{
			if (candidate.vertices[i] == parent.vertices[j])
			{
				disjoint = 0;
				break;
			}
		}
	}
	
	return disjoint;
}

int Vertex_Sum(const cycle_arcs & parent, const cycle_arcs & candidate)
{
	int nr_vertices = parent.vertices.size(); 
	for (int i = 0; i < candidate.vertices.size(); i++)
	{
		bool notinparent = 1;
		for (int j = 0; j < parent.vertices.size(); j++)
		{
			if (candidate.vertices[i] == parent.vertices[j])
			{
				notinparent = 0;
				break;
			}
		}
		if (notinparent == 1)
			nr_vertices++;
	}
	
	return nr_vertices;
}

bool arcsort(directedarc a, directedarc b)
{
	bool afirst = 0;
	if (a.startvertex < b.startvertex)
		afirst = 1;
	else if (a.startvertex == b.startvertex && a.endvertex < b.endvertex)
		afirst = 1;

	
	return afirst;
}