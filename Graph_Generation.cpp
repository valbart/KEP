
#include "stdafx.h"
#include "Structures.h"
#include "Functions.h"

directedgraph graph_generation(configuration & config)
{
	directedgraph G;
	// Points to the correct function for Graph Generation (or read)
	if (config.input_data == 1)
		G = readgraph(config, config.inputfile);
	else if (config.input_data == 2)
		G = readgraph2(config, config.inputfile);
	else if (config.input_data == 3)
		G = simplegraph(config);
	else if (config.input_data == 4)
		G = saidman(config);
	
	G.first_vertex_arc = Find_First_arc(G);

	return G;
}

directedgraph simplegraph(configuration & config)
{
	cout << "Using simple generator" << endl;
	directedgraph G;
	float arcprob = config.arcprob;
	float failprob = config.failprob;
	G.size = config.nr_ndd + config.nr_pairs; int size = G.size;
	G.nr_ndd = config.nr_ndd;
	G.nr_pairs = config.nr_pairs;
	int count = 0; // The total number of generated arcs.
	ofstream output;
	output.open(config.graph_output);

	output << "Nr_Pairs = " << G.nr_pairs << endl;
	output << "Nr_NDD = " << G.nr_ndd << endl;
	
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < G.nr_pairs; j++)
		{
			float prob = float(rand() % 100 + 1)/100;
			if (prob <= arcprob && i != j)
			{
				directedarc arc;
				arc.startvertex = i; arc.endvertex = j; arc.failprob = failprob; arc.weight = 1; arc.arcnumber = count;
				count++;
				G.arcs.push_back(arc);
				output << "(" << i << "," << j << "), " << arc.failprob << ", " << arc.weight << endl;
			}
		}
	}
	
	return G;
}

directedgraph saidman(configuration & config)
{
	cout << "Using Saidman Generator" << endl;
	directedgraph G;
	G.nr_pairs = config.nr_pairs;
	G.nr_ndd = config.nr_ndd;
	G.size = G.nr_ndd + G.nr_pairs;
	
	
	srand(time(NULL));

	vector<patient_donor_pair> pairs(G.nr_pairs);
	vector<int> NDD(G.nr_ndd);
	// Loop to determine blood types and PRA_level of pairs. 
	for (int i = 0; i < G.nr_pairs; i++)
	{
		bool flag = 1;
		while (flag == 1) // We generate a donor and patient. If they are compatible, the flag stays at 1. If they are compatible, 
			// they would not be added to the KEP pool, so we generate a new patient-donor pair to take their place.
			// If at any point we find they are not compatible, the flag is set to 0.
		{
			int random_nr = rand() % 10000;
			if (random_nr < 4814)
				pairs[i].donor_bloodtype = 0;
			else if (random_nr < 8187)
				pairs[i].donor_bloodtype = 1;
			else if (random_nr < 9615)
				pairs[i].donor_bloodtype = 2;
			else
				pairs[i].donor_bloodtype = 3;

			random_nr = rand() % 10000;
			if (random_nr < 4814)
				pairs[i].patient_bloodtype = 0;
			else if (random_nr < 8187)
				pairs[i].patient_bloodtype = 1;
			else if (random_nr < 9615)
				pairs[i].patient_bloodtype = 2;
			else
				pairs[i].patient_bloodtype = 3;

			random_nr = rand() % 10000;
			if (random_nr < 7019)
				pairs[i].PRA_level = 0; // Low PRA
			else if (random_nr < 9019)
				pairs[i].PRA_level = 1; // Medium PRA
			else
				pairs[i].PRA_level = 2;	// High PRA
			
			// We first check whether the blood types are incompatible. If they are, we set the flag to 0 and break the loop.
			
			if (pairs[i].donor_bloodtype != 0 && pairs[i].patient_bloodtype != 3 && pairs[i].donor_bloodtype != pairs[i].patient_bloodtype)
			{
				flag = 0; 
				break;
				
			}
			{
				// If there is a compatible blood type, there is still the possibility of positive cross-match.
				// Check whether there is a positive cross-match. If there is, break the loop. (low = 95%, medium 45%, high 10%)
				int PRA_odds = rand() % 100;
				if (PRA_odds < 10 || (PRA_odds < 45 && pairs[i].PRA_level == 1) || (PRA_odds < 95 && pairs[i].PRA_level == 2)) // If any of these conditions are met, there is a positive crossmatch.
				{
					flag = 0;
					break;
				}
				// Additional check, because chance of negative cross-match is only 75% of normal given female patient and spouse donor.
				/*if (flag == 1)
				{
					// We first determine the relation of the patient and donor.
					bool female = 0; bool spouse = 0;
					random_nr = rand() % 10000;
					if (random_nr < 4090)
						female = 1;
					if (female == 1)
					{
						random_nr = rand() % 10000;
						if (random_nr < 4897)
							spouse = 1;
					}
					if (female == 1 && spouse == 1)
					{
						PRA_odds = rand() % 100;
						if (PRA_odds < 25)
						{
							flag = 0; break;
						}
					}
				}*/			
			}
		}
		cout << pairs[i].donor_bloodtype << "\t" << pairs[i].patient_bloodtype << "\t" << pairs[i].PRA_level << endl;
	} 
	// Loop to determine blood type of NDD.
	for (int i = 0; i < G.nr_ndd; i++)
	{
		float random_nr = rand() % 10000;
		if (random_nr < 4814)
			NDD[i] = 0;
		else if (random_nr < 8187)
			NDD[i] = 1;
		else if (random_nr < 9615)
			NDD[i] = 2;
		else
			NDD[i] = 3;
	}

	int arcnumber = 0;
	ofstream output;
	output.open(config.graph_output);
	output << "Nr_Pairs = " << G.nr_pairs << endl;
	output << "Nr_NDD = " << G.nr_ndd << endl;

	// Arcs between donor-patient pairs.
	for (int i = 0; i < G.nr_pairs; i++)
	{
		for (int j = 0; j < G.nr_pairs; j++)
		{
			// Determined whether donor of i can donate to patient of j.
			if ((pairs[i].donor_bloodtype == 0 || pairs[j].patient_bloodtype == 3 || pairs[i].donor_bloodtype == pairs[j].patient_bloodtype) && i != j)
			{
				// Check whether PRA allows for it (low = 95%, medium 45%, high 10%)
				int PRA_odds = rand() % 100;
				if (PRA_odds < 10 || (PRA_odds < 45 && pairs[j].PRA_level < 2) || (PRA_odds < 95 && pairs[j].PRA_level < 1))
				{
					directedarc arc;
					arc.startvertex = i; arc.endvertex = j; arc.weight = 1; arc.failprob = config.failprob; arc.arcnumber = arcnumber;
					arcnumber++;
					G.arcs.push_back(arc);
					output << "(" << i << "," << j << "), " << arc.failprob << ", " << arc.weight << endl;
				}
			}
		}
	}
	// Arcs from NDD
	for (int i = 0; i < G.nr_ndd; i++)
	{
		for (int j = 0; j < G.nr_pairs; j++)
		{
			// Determined whether NDD i can donate to patient of j.
			if (NDD[i] == 0 || pairs[j].patient_bloodtype == 3 || NDD[i] == pairs[j].patient_bloodtype)
			{
				// Check whether PRA allows for it (low = 95%, medium 45%, high 10%)
				int PRA_odds = rand() % 100;
				if (PRA_odds < 10 || (PRA_odds < 45 && pairs[j].PRA_level < 2) || (PRA_odds < 95 && pairs[j].PRA_level < 1))
				{
					directedarc arc;
					arc.startvertex = i + G.nr_pairs; arc.endvertex = j; arc.weight = 1; arc.failprob = config.failprob; arc.arcnumber = arcnumber;
					arcnumber++;
					G.arcs.push_back(arc);
					output << "(" << i << "," << j << "), " << arc.failprob << ", " << arc.weight << endl;
				}
			}
		}
	}

	return G;
}

directedgraph readgraph(configuration & config, string filename)
{
	cout << "Reading graph from file" << endl;
	ifstream input(filename);
	if (!input)
	{
		cout << "No input file found." << endl;
	}

	directedgraph G;
	input.ignore(50, '=');
	input >> G.nr_pairs; config.nr_pairs = G.nr_pairs;
	input.ignore(50, '=');
	input >> G.nr_ndd; config.nr_ndd = G.nr_ndd;

	G.size = G.nr_ndd + G.nr_pairs;
	bool end_arcs = 0;
	int nr_arcs = 0;
	string str;
	getline(input, str);
	while (end_arcs == 0)
	{
		string str;
		stringstream stream;
		getline(input, str);
		if (str[0] == '(')
		{
			directedarc arc;
			arc.arcnumber = nr_arcs;
			stream.str(str);
			stream.ignore(1, '(');
			stream >> arc.startvertex;
			stream.ignore(1, ',');
			stream >> arc.endvertex;
			stream.ignore(2, ',');
			stream >> arc.failprob;
			stream.ignore(1, ',');
			stream >> arc.weight;
			nr_arcs++;
			G.arcs.push_back(arc);
		}
		else
			end_arcs = 1;
	}
	return G;
}
directedgraph readgraph2(configuration & config, string filename)
{
	cout << "Reading graph from file" << endl;
	ifstream input(filename);
	if (!input)
	{
		cout << "No input file found." << endl;
	}

	directedgraph G;
	input.ignore(50, '=');
	input >> G.nr_pairs; config.nr_pairs = G.nr_pairs;
	input.ignore(50, '=');
	input >> G.nr_ndd; config.nr_ndd = G.nr_ndd;

	G.size = G.nr_ndd + G.nr_pairs;
	for (int i = 0; i < G.nr_pairs; i++)
	{
		patient_donor_pair new_pair;
		int pair_nr;
		input >> pair_nr;
		input >> new_pair.failprob;
		G.pairs.push_back(new_pair);
	}
	for (int i = 0; i < G.nr_ndd; i++)
	{
		ndd_donor new_donor;
		int ndd_nr;
		input >> ndd_nr;
		input >> new_donor.failprob;
		G.ndds.push_back(new_donor);
	}
	bool end_arcs = 0;
	int nr_arcs = 0;
	string str;
	getline(input, str);
	while (end_arcs == 0)
	{
		string str;
		stringstream stream;
		getline(input, str);
		if (str[0] == '(')
		{
			directedarc arc;
			arc.arcnumber = nr_arcs;
			stream.str(str);
			stream.ignore(1, '(');
			stream >> arc.startvertex;
			stream.ignore(1, ',');
			stream >> arc.endvertex;
			stream.ignore(2, ',');
			stream >> arc.failprob;
			stream.ignore(1, ',');
			stream >> arc.weight;
			nr_arcs++;
			G.arcs.push_back(arc);
		}
		else
			end_arcs = 1;
	}
	return G;
} // Not yet implemented
