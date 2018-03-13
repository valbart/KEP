#include "stdafx.h"
#include "Structures.h"
#include "Functions.h"

vector<cycle_arcs> Find_Cycles(directedgraph G, const configuration & config)
{
	// This function identifies all cycles in the graph.
	vector<cycle_arcs> cycles;

	vector<vector<directedarc>> copies = DP_Copy(G);
	vector<vector<int>> distances = distance_calc(G, copies, config.cyclelength);
	vector<int> first_arc = Find_First_arc(G);
	
	for (int i = 0; i < G.nr_pairs; i++) // We loop through all vertices as possible starting points of the cycles.
	{
		queue<vector<int>> paths; // A queue (to allow deleting elements at the front) including all current paths.
		// We intialize the queue by adding all arcs from the initial vertex (if a short enough return path exists).
		for (int j = first_arc[i]; j < first_arc[i+1]; j++)
		{
			if (distances[i][G.arcs[j].endvertex] < config.cyclelength) // Note that due to the distance calc, no arcs to lower numbered vertices will be used in this initialization (since return distance = cyclelength + 1);
			{
				vector<int> init_vector;
				init_vector.push_back(j);
				paths.push(init_vector);
			}
		}
		while (!paths.empty())
		{
			vector<int> path_vector = paths.front(); // Get the first path out and delete it from the queue.
			paths.pop();

			int endvertex = G.arcs[path_vector[path_vector.size() - 1]].endvertex;

			for (int j = first_arc[endvertex]; j < first_arc[endvertex + 1]; j++) // Go through the arcs starting from the current last vertex on the path.
			{
				if(G.arcs[j].endvertex >= i) // We do not include cycles using vertices with number lower than the current starting vertex. These cycles have been added earlier when the starting vertex was lower.
				{
					if (G.arcs[j].endvertex == i) // If the path returns to the intial vertex, save the cycle.
					{
						cycle_arcs temp_cycle;
						temp_cycle.arcs = path_vector;
						temp_cycle.arcs.push_back(j);
						temp_cycle.weight = 0;
						cycles.push_back(temp_cycle);
					}
					else if (distances[i][G.arcs[j].endvertex] + path_vector.size() < config.cyclelength)
					{
						vector<int> new_path = path_vector;
						new_path.push_back(j);
						paths.push(new_path);
					}
				}
			}
		}

	}
	for (int i = 0; i < cycles.size(); i++)
	{
		for (int j = 0; j < cycles[i].arcs.size(); j++)
		{		
			cycles[i].vertices.push_back(G.arcs[cycles[i].arcs[j]].startvertex);
			cycles[i].weight = cycles[i].weight + G.arcs[cycles[i].arcs[j]].weight;
		}
	}

	return cycles;
}

vector<vector<int>> List_Outgoing_Arcs(const directedgraph & G)
{
	// This function lists, for each vertex in the graph, the destinations of the outgoing arcs.
	vector<vector<int>> outgoing_arcs(G.size);

	for (int i = 0; i < G.arcs.size(); i++)
	{
		outgoing_arcs[G.arcs[i].startvertex].push_back(G.arcs[i].endvertex);
	}

	return outgoing_arcs;
}

vector<int> Find_First_arc(directedgraph & G)
{
	// Returns a vector that gives, for each vertex, the first place where it appears as a startvertex in the arc list.
	vector<int> first_arc(G.size+1); // The extra element is to avoid misreferences and is set to the first arc after the pair-pair-arcs are done.

	first_arc[0] = 0;
	int vertex = 1;

	for (int i = 0; i < G.arcs.size(); i++)
	{
		if (G.arcs[i].startvertex >= vertex)
		{
			first_arc[vertex] = i;	// If the if statement is TRUE, this must be the first time G.arcs[i].startvertex >= vertex (since vertex increases
									// in the following line. 
			vertex++;				// We start evaluating the following vertex.
			i--;					// In case one vertex has not outgoing vertices; We lower i by one, so the same arc is evaluated again.
		}
	}

	for (int i = vertex; i <= G.size; i++)
	{
		first_arc[i] = G.arcs.size();
	}
	

	return first_arc;
}
