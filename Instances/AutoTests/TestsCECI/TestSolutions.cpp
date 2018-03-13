#include "stdafx.h"
#include "Structures.h"
#include "Functions.h"

float Exact_Expected_Transplants(configuration & config);
float Scenarios_Expected_Transpants(configuration & config);

// Function to split pre_test_result into strongly connected components.
// Call to Subset_Weights, to find the expected number of transplants in a SCC.

float Exact_Expected_Transplants(configuration & config)
{
	float Expected_Transplants = 0;
	directedgraph G_test;
	if (config.failure_type == 1)
	{
		G_test = readgraph(config, config.solution_input);
	}
	else if (config.failure_type == 2)
	{
		directedgraph G_orig = readgraph2(config, config.inputfile);
		G_test = readgraph(config, config.solution_input);
		G_test.pairs = G_orig.pairs;
		G_test.ndds = G_orig.ndds;
	}
	cout << "Graph Read" << endl;
	vector<bool> fixed_to_zero (G_test.nr_pairs, 0);
	vector<bool> fixed (G_test.nr_pairs, 0);
	vector<cycle_arcs> SCC_T = Split_SCC_Tarjan(G_test, fixed_to_zero);
	if (config.failure_type == 1)
	{
		Subset_Set_Weights_Arc_verbose(SCC_T, G_test, config);
	}
	else if (config.failure_type == 2)
	{
		clock_t start_time = clock();
		Subset_Set_Weights_Vertex_Verbose(SCC_T, G_test, config, fixed);
		clock_t total_time = clock() - start_time;
		cout << endl << "Solution file " << config.solution_input << " PSE time " << ((float)total_time)/CLOCKS_PER_SEC << endl;
	}
	float total_weight = 0;
	for (int i = 0; i < SCC_T.size(); i++)
	{
		cout << "SCC " << i << " weight = " << SCC_T[i].weight << endl;
		total_weight = total_weight + SCC_T[i].weight;
	}
	ofstream output;
	output.open(config.testvar_output, std::ios_base::app);
	output << "Expected Transplants (Weighted) = " << total_weight << endl;
	return Expected_Transplants;
}

float Scenarios_Expected_Transpants(configuration & config)
{
	float Expected_Tranplants = 0;
	directedgraph G_Orig = readgraph2(config, config.inputfile);
	directedgraph G = readgraph(config, config.solution_input);
	G.pairs = G_Orig.pairs;
	G.ndds = G_Orig.ndds;
	vector<directedgraph> Scenarios;

	if (config.failure_type == 1)
	{
		cout << "Arcs Fail" << endl;
		if (config.scen_gen == 1)
		{
			Scenarios = Generate_Scenarios_Tight(G, config.nr_scenarios); cout << "Tight Scen Generator" << endl;
		}
		else
		{
			Scenarios = Generate_Scenarios(G, config.nr_scenarios); cout << "Basic Scen Generator" << endl;
		}

	}
	else if (config.failure_type == 2)
	{
		cout << "Vertices Fail" << endl;
		Scenarios = Generate_Scenarios_Vertex_Tight(G, config.nr_scenarios);
	}

	// Re-number the arcs, or things break in the HPIEF function. We do save the correct arcnumbers to relate results to the arcs in G.
	for (int i = 0; i < Scenarios.size(); i++)
	{
		for (int j = 0; j < Scenarios[i].arcs.size(); j++)
		{
			Scenarios[i].arcs[j].arcnumber = j;
		}
	}

	int total_transplant = 0;
	for(int i = 0; i < config.nr_scenarios; i++)
	{
		matching_result result = Hybrid_PIEF(Scenarios[i], config.chainlength, config.cyclelength, config);
		total_transplant = result.objective + total_transplant;
	}
	cout << total_transplant << endl;
	cout << float (total_transplant) / float (config.nr_scenarios) << endl;
	Expected_Tranplants = float (total_transplant) / float(config.nr_scenarios);
	ofstream output;
	output.open(config.testvar_output, std::ios_base::app);
	output << "Expected Transplants (Weighted, Approximate) = " << Expected_Tranplants << endl;

	return Expected_Tranplants;
}

float Calc_Expected_Transplants(configuration & config)
{
	float Expected_Transplants;
	if (config.calc_expected_type == 1)
		Expected_Transplants = Exact_Expected_Transplants(config);
	else if (config.calc_expected_type == 2)
		Expected_Transplants = Scenarios_Expected_Transpants(config);

	return Expected_Transplants;
}

vector<cycle_arcs> Split_SCC(const directedgraph & tested_G)
{
	// This functions assumes every arc is part of an SCC.
	vector<cycle_arcs> SCC;
	vector<bool> added(tested_G.size);

	for (int i = 0; i < tested_G.size; i++)
	{
		if (added[i] == 0)
		{
			cout << "Starting SCC with " << i << endl;
			added[i] = 1;
			cycle_arcs new_SCC;
			queue<int> to_test;
			to_test.push(i);
			while (to_test.size() > 0)
			{
				for (int j = 0; j < tested_G.arcs.size(); j++)
				{
					if (tested_G.arcs[j].endvertex >= i)
					{
						if (tested_G.arcs[j].startvertex == to_test.front())
						{
							new_SCC.arcs.push_back(j);
						}
						if (tested_G.arcs[j].startvertex == to_test.front() && added[tested_G.arcs[j].endvertex] == 0)
						{
							added[tested_G.arcs[j].endvertex] = 1;
							to_test.push(tested_G.arcs[j].endvertex);
						}
					}
				}
				new_SCC.vertices.push_back(to_test.front());
				to_test.pop();

			}
			if (new_SCC.vertices.size() > 1)
			{
				SCC.push_back(new_SCC);
				for (int j = 0; j < new_SCC.vertices.size(); j++)
				{
					cout << new_SCC.vertices[j] << "\t";
				}
				cout << endl;
				for (int j = 0; j < new_SCC.arcs.size(); j++)
				{
					cout << "(" << tested_G.arcs[new_SCC.arcs[j]].startvertex << "," << tested_G.arcs[new_SCC.arcs[j]].endvertex << ")" << endl;
				}
				cout << endl;
			}

		}
	}

	return SCC;
}

vector<list<int>> G_Adjacency(const directedgraph & tested_G, const vector<bool> & fixed_to_zero) {
	// We keep track of the arcs linked to v and not exactly vertices as we need the
	// index of the arcs to build the SCCs later on. Furthemore the endvertex can be easily
	// gathered through test_G
	vector<list<int>> adj (tested_G.nr_pairs);
	int sv, ev;
	for (unsigned int i = 0; i < tested_G.arcs.size(); ++i) {
		sv = tested_G.arcs[i].startvertex;
		ev = tested_G.arcs[i].endvertex;
		if (fixed_to_zero[sv] == 0 && fixed_to_zero[ev] == 0) adj[sv].push_back(i);
	}
	return adj;
}

vector<cycle_arcs> Split_SCC_Tarjan(const directedgraph & tested_G, const vector<bool> & fixed_to_zero) {
	vector<cycle_arcs> SCC;
	int index = 0;
	vector<int> vertices_index (tested_G.nr_pairs, -1); // -1 means the index of the vertex is undefined.
	vector<int> low_link (tested_G.nr_pairs, -1);
	vector<bool> on_stack(tested_G.nr_pairs, 0);
	vector<list<int>> adj = G_Adjacency(tested_G, fixed_to_zero);
	deque<int> Q;
	vector<vector<int>> SCC_arcs(tested_G.nr_pairs);
	for (unsigned int v = 0; v < tested_G.nr_pairs; ++v) {
		if (vertices_index[v] == -1) Find_SCC_Root(tested_G, adj, v, index, vertices_index, Q, low_link, on_stack, SCC, SCC_arcs);
	}
	return SCC;
}

vector<tuple<int,int,float>> Find_Articulation_Points(const directedgraph & G) {
	vector<tuple<int, int, float>> articulation_points;
	vector<bool> fixed_to_zero;
	vector<cycle_arcs> SCC;
	float score;
	int max_SCC;
	for (unsigned int i = 0; i < G.nr_pairs; ++i) {
		score = 1;
		fixed_to_zero = vector<bool>(G.nr_pairs, 0);
		fixed_to_zero[i] = 1;
		SCC = Split_SCC_Tarjan(G, fixed_to_zero);
		if (SCC.size() > 1) {
			max_SCC = SCC[0].vertices.size();
			for (unsigned int j = 0; j < SCC.size(); j++) {
				if (SCC[j].vertices.size() > max_SCC) max_SCC = SCC[j].vertices.size();
			}
			articulation_points.push_back(tuple<int, int, float>(i, SCC.size(), max_SCC));
		}
	}
	return articulation_points;
}


void Find_SCC_Root(const directedgraph & tested_G, const vector<list<int>> & adj, int v, int & index, vector<int> & vertices_index,
				   deque<int> & Q, vector<int> & low_link, vector<bool> & on_stack, vector<cycle_arcs> & SCC, vector<vector<int>> & SCC_arcs) {
	vertices_index[v] = index;
	low_link[v] = index;
	index++;
	Q.push_front(v);
	on_stack[v] = true;
	for (list<int>::const_iterator it = adj[v].begin(); it != adj[v].end(); ++it) {
		int w = tested_G.arcs[*it].endvertex;
		if (vertices_index[w] == -1) {
			// (*it) has not been visited yet, explore it to find its low link
			Find_SCC_Root(tested_G, adj, w, index, vertices_index, Q, low_link, on_stack, SCC, SCC_arcs);
			low_link[v] = min(low_link[v], low_link[w]);
			if (on_stack[w]) {
				SCC_arcs[v].push_back(*it);
			}
		}
		else if (on_stack[w]) {
			// (*it) has been visited and is still on the stack, hence belongs to the same SCC
			// either rooted in *it (min = index[*it]) or by a vertex above in the DFS tree (min = low_link[v])
			low_link[v] = min(low_link[v], vertices_index[w]);
			SCC_arcs[v].push_back(*it);
		}
	}

	// Then we pop the vertices to build the SCC
	int w = -1;
	if (low_link[v] == vertices_index[v]) {
		cycle_arcs new_SCC;
		do {
			if (Q.empty()) break;
			w = Q.front();
			Q.pop_front();
			on_stack[w] = 0;
			new_SCC.vertices.push_back(w);
		} while(w != v);

		if (new_SCC.vertices.size() > 1) {
			for (int j = 0; j < new_SCC.vertices.size(); j++) {
				// We concatenate the list of arcs that belongs to the current SCC together
				new_SCC.arcs.insert(new_SCC.arcs.end(), SCC_arcs[new_SCC.vertices[j]].begin(), SCC_arcs[new_SCC.vertices[j]].end());
			}
			SCC.push_back(new_SCC);
		}
	}
	return;
}
