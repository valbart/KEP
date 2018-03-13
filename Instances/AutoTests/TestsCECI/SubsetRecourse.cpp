#include "stdafx.h"
#include "Structures.h"
#include "Functions.h"
#include <algorithm>
#include <time.h>
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
void Subset_Set_Weights_Vertex_Verbose(vector<cycle_arcs>& subsets, const directedgraph & G, const configuration & config, vector<bool> & fixed)
{
	int count = 0;
//#pragma omp parallel for
	for (int i = 0; i < subsets.size(); i++)
	{
		// Build subset graph
		vector<bool> new_fixed (subsets[i].vertices.size(), 0);
		directedgraph G_sub = Subset_Graph_New(subsets[i], G, new_fixed, fixed);
		// Give Graph to function computing weight for subset
		subsets[i].weight = Subset_Set_Weights_Vertex_Root_Verbose(G_sub, config, new_fixed);
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


// The two following methods are used to branch on the vertices according to an order that is updated
// during the algorithm. Currently they are not used as the order we chose is to simply compute the k.e.a order
// once at the beginning of the algorithm
// Note that the candidates for branching are represented as tuple of two integers: the first one corresponds
// to the candidate ID and the second one can be used to store some data about the current graph (for example
// the current number of cycles this vertex is involved in).

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

// Function to sort the articulations points. They are represented as tuple of 2 int and 1 float. 
// The first integer is the vertex ID, the second integer is the number of SCCs linked to this articulation.
// point. The third number is currently the size of the biggest SCC linked to this point (if removed from graph)
// It has been set as a float so that it can easily be replaced by some other measure that could be more relevant
// to choose the point to branch on. We use currently a lexicographic order that first aim to maximize the number 
// of SCCs and then to minimize the size of the biggest one.
// NB: the order is reversed so that one can simply push these point on the front of the current order
bool Ap_compare(tuple<int, int, float> const & t1, tuple<int, int, float> const & t2) {
	// If same size, we give priority to the smallest one. Not necessary equal as one may 
	// have deleted lone vertices when branching
	if (get<1>(t1) == get<1>(t2)) return get<2>(t1) >= get<2>(t2);
	return get<1>(t1) > get<2>(t1);
}

float Subset_Set_Weights_Vertex_Root_Verbose(const directedgraph & G, const configuration & config, vector<bool> & fixed)
{
	float weight;

	// We first find all the cycles
	vector<cycle_arcs> cycles = Find_Cycles(G, config);
	if (cycles.size() == 0) return 0;
	if (cycles.size() == 1) {
		float tmp_weight = 1;
		for (unsigned int i = 0; i < cycles[0].vertices.size(); ++i) {
			if (fixed[cycles[0].vertices[i]] == 0) tmp_weight *= (1-G.pairs[cycles[0].vertices[i]].failprob);
		}
		tmp_weight *= cycles[0].vertices.size();
		return tmp_weight;
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

	// We concatenate the lists of nr_appearance_lists to obtain the final order
	list<tuple<int,int>> branching_order;
	for (unsigned int i = 0; i < nr_appearance_lists.size(); ++i) {
		if (!nr_appearance_lists[i].empty()) branching_order.splice(branching_order.begin(), nr_appearance_lists[i]);
	}

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

	vector<tuple<int,int>> vertex_fix (2, tuple<int,int>(-1,-1));
	vector<tuple<int, int, float>> a_p = Find_Articulation_Points(G);
	
	// We try to find the best possible articulation point and branch on it first.
	if (a_p.size() > 0) {
		sort(a_p.begin(), a_p.end(), Ap_compare);
		for (int i = a_p.size()-1; i >= 0; --i) {
			if (vertex_use[get<0>(a_p[i])] == 1) {
				branching_order.push_front(tuple<int,int>(get<0>(a_p[i]), 0));
				break;
			}
		}
	}

	bool new_fix = 0;
	// We have to know which vertex have fail in order to be able to compute the potential SCC
	// after each new failure.
	vector<bool> fixed_to_zero(G.nr_pairs, 0);

	for(list<tuple<int,int>>::iterator it = branching_order.begin(); it != branching_order.end(); ++it) {
		int id = get<0>(*it);
		if(vertex_use[id] == 1 && fixed[id] == 0) {
			new_fix = 1;
			fixed[id] = 1;
			vertex_fix[0] = *it;
			break;
		}
	}	
	if (new_fix == 0) return CPLEX_Solution;

	// We look for a first candidate. Use the line below to find the next vertex to branch on
	// with the branching order. We do not have to use it currently as the branching order never changes
	// in the current version.
	// Find_next_branch(vertex_fix, 0, vertex_use, fixed, branching_order, new_fix);

	cout << "Model Solved, Branching starting" << endl;
	float solution = 0;
	int depth = 0;
	deque<int> suc_fail;

	if (CPLEX_Solution > 0) {
		depth++;
		// First assume the vertex is a succes
		// We go down one level and branch again on another arc used in the solution.
		suc_fail.push_back(1);
		weight = (1 - G.pairs[get<0>(vertex_fix[0])].failprob)*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, CPLEX_Solution, 1, depth, suc_fail, branching_order, fixed_to_zero);
		suc_fail.pop_back();
		suc_fail.push_back(0);

		// Next assume the vertex fails.
		vertex_cons[get<0>(vertex_fix[0])].setUB(0);
		fixed_to_zero[get<0>(vertex_fix[0])] = 1;
		weight = weight + G.pairs[get<0>(vertex_fix[0])].failprob*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, CPLEX_Solution, 0, depth, suc_fail, branching_order, fixed_to_zero);
		suc_fail.pop_back();
		vertex_cons[get<0>(vertex_fix[0])].setUB(1);		
	}
	env.end();
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

	if (depth >= 2 && succes_fix == 0)
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

// Function to compare the vertex according to their parameters (nr of cycles or degree, or...)
bool tuple_compare(tuple<int, int> const & t1, tuple<int,int> const & t2) {
	return get<1>(t1) >= get<1>(t2);
}

// Functions to compute a branching order from scratch according to the sum of the in and out degree, 
// sorted in decreasing order.
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

// Function to update the current order. This is the one used for the KEA "heuristic" order. 
// It simply computes the sum of in/out degrees of the entering vertices (the one that are now being
// used after an other vertes has failed). Push it on front if this number is higher than the parameter
// of the vertex currently at the front of the order, else push it back. 
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
}

float Subset_Set_Weights_Vertex_Recursion_Verbose(const configuration & config, IloEnv & env, 
	const directedgraph & G, IloCplex & CPLEX, vector<int> vertex_use, IloRangeArray & vertex_cons, vector<bool> fixed, float Solution, 
	bool succes_fix, int depth, deque<int> & suc_fail, list<tuple<int, int>> branching_order, vector<bool> & fixed_to_zero)
{
	/*
	if (depth <15)
	{
		for (int i = 0; i < suc_fail.size(); i++)
			cout << suc_fail[i] << ",";
		cout << endl;
	}
	*/
	float weight = 0;
	if (succes_fix == 0)
	{
		// First see if graph has been broken into several SCCs.
		// If it is the case, we stop the recursion here and start 
		// again from root with the found SCCs. 
		vector<cycle_arcs> SCC_T = Split_SCC_Tarjan(G, fixed_to_zero);
		if (SCC_T.size() > 1) {
			Subset_Set_Weights_Vertex_Verbose(SCC_T, G, config, fixed);
			for (int i = 0; i < SCC_T.size(); i++) weight = weight + SCC_T[i].weight;
			return weight;
		}

		CPLEX.solve();
		Solution = ceil(CPLEX.getObjValue());

		for (int i = 0; i < G.size; i++) 
		{
			vertex_use[i] = 1 - round(CPLEX.getSlack(vertex_cons[i]));
		}
	}

	vector<tuple<int,int>> vertex_fix (2, tuple<int,int>(-1,-1));
	bool new_fix = 0;

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

	if (new_fix == 1)
	{
		depth++;
		// First assume the vertex is a success
		// We go down one level and branch again on another arc used in the solution.
		suc_fail.push_back(1); 
		weight = (1 - G.pairs[get<0>(vertex_fix[0])].failprob)*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, Solution, 1, depth, suc_fail, branching_order, fixed_to_zero);
		suc_fail.pop_back();

		// Next assume the vertex fails.
		suc_fail.push_back(0);
		vertex_cons[get<0>(vertex_fix[0])].setUB(0);
		weight = weight + G.pairs[get<0>(vertex_fix[0])].failprob*Subset_Set_Weights_Vertex_Recursion_Verbose(config, env, G, CPLEX, vertex_use, vertex_cons, fixed, Solution, 0, depth, suc_fail, branching_order, fixed_to_zero);
		suc_fail.pop_back();
		vertex_cons[get<0>(vertex_fix[0])].setUB(1);
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

// Just a few modification made on the function above. Indeed, when breaking the graph with SCC and starting
// with new tree, the vertices get reindexed. Hence we also have to reindex the vertex that have already been
// fixed so that we do not branch twice on the same vertex.
// Not that we do not have to keep track of the fixed_to_zero vertices has these necessarily disapear 
// when breaking the graph.
directedgraph Subset_Graph_New(const cycle_arcs & subset, const directedgraph & G, vector<bool> & new_fixed, vector<bool> & fixed)
{
	directedgraph G_sub;
	G_sub.size = subset.vertices.size();
	G_sub.nr_pairs = G_sub.size;
	vector<int> vertex_switch(G.nr_pairs); // For each vertex in the original graph, we save in which position it occurs, in the list of vertices in the subset.
	for (int i = 0; i < subset.vertices.size(); i++)
	{
		vertex_switch[subset.vertices[i]] = i;
		new_fixed[i] = fixed[subset.vertices[i]];
		G_sub.pairs.push_back(G.pairs[subset.vertices[i]]);
	}
	if (subset.vertices.size() == 31) {
		cout << "Vertex switch :";
		for (unsigned int i = 0; i < subset.vertices.size(); ++i) cout << subset.vertices[i] << ":" << vertex_switch[subset.vertices[i]] << "  ";
		cout << endl;
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