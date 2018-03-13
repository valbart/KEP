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
	vector<cycle_arcs> SCC = Split_SCC(G_test);
	if (config.failure_type == 1)
	{
		Subset_Set_Weights_Arc_verbose(SCC, G_test, config);
	}
	else if (config.failure_type == 2)
	{
		clock_t start_time = clock();
		Subset_Set_Weights_Vertex_Verbose(SCC, G_test, config);
		clock_t end_time = clock() - start_time;
		cout << endl << "Solution file : " << config.solution_input << " PSE time OLD : " << ((float)end_time)/CLOCKS_PER_SEC << endl;
	}
	float total_weight = 0;
	for (int i = 0; i < SCC.size(); i++)
	{
		cout << "SSC " << i << " weight = " << SCC[i].weight << endl;
		total_weight = total_weight + SCC[i].weight;
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
			// cout << "Starting SCC with " << i << endl;
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
				/*for (int j = 0; j < new_SCC.vertices.size(); j++)
				{
					cout << new_SCC.vertices[j] << "\t";
				}
				cout << endl;
				for (int j = 0; j < new_SCC.arcs.size(); j++)
				{
					cout << "(" << tested_G.arcs[new_SCC.arcs[j]].startvertex << "," << tested_G.arcs[new_SCC.arcs[j]].endvertex << ")" << endl;
				}
				cout << endl;*/
			}
				
		}
	}

	return SCC;
}
