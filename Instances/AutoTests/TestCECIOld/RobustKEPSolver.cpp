// RobustKEPSolver.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Structures.h"
#include "Functions.h"

int main(int argc, char  * argv[])
{
	int error = 0;
	
	configuration config = Readconfig(argv[1], &error);
	directedgraph G = graph_generation(config);
	if (config.LW == 0)
		Choose_Solver(config, G);
	else
		Limited_World(config, G);

    return 0;
}

void Choose_Solver(configuration & config, directedgraph & G)
{
	cout << "Choosing Solver: ";
	if (config.solver == 1)
	{
		cout << "Subset Recourse" << endl;
		Subset_Recourse(config, G);
	}
	else if (config.solver == 2)
	{
		cout << "Pre-Test" << endl;
		pre_test_main(config, G);
	}
	else if (config.solver == 3)
	{
		cout << "Unlimited Cycles" << endl;
		UnlimitedCycle_Main(config, G);
	}
	else if (config.solver == 4)
	{
		cout << "Matheuristics" << endl;
		Cycle_Heuristic(config, G);
	}
	else if (config.solver == 5)
	{
		cout << "Expected Transplants" << endl;
		Calc_Expected_Transplants(config);
	}
		
}

configuration Readconfig(char configname[], int * error)

{
	configuration config;
	namespace po = boost::program_options;
	ifstream configfile(configname);
	if (!configfile)
	{
		cout << "No configuration file found." << endl;
		*error = 1;
	}
	po::options_description file_options(configname);
	file_options.add_options()
		("Cyclelength", po::value<int>(&config.cyclelength)->required(), "")
		("Chainlength", po::value<int>(&config.chainlength)->required(), "")
		("Scenarios", po::value<int>(&config.nr_scenarios)->default_value(10), "")
		("Max Test", po::value<int>(&config.max_test)->default_value(10), "")
		("Subset Recourse Size", po::value<int>(&config.subset_size)->default_value(0), "")
		("Input Type", po::value<int>(&config.input_data)->default_value(1), "")
		("Solution Input File", po::value<string>(&config.solution_input)->default_value("Graph.txt"), "")
		("Input File", po::value<string>(&config.inputfile)->default_value("Graph.txt"), "")
		("Graph Output", po::value<string>(&config.graph_output)->default_value("Graph.txt"), "")
		("Testvar Output", po::value<string>(&config.testvar_output)->default_value("Testvar.txt"), "")
		("Time Limit", po::value<int>(&config.time_limit)->default_value(86400), "")
		("Memory Limit", po::value<int>(&config.memory_limit)->default_value(8192),"")
		("Solver", po::value<int>(&config.solver)->required(), "")
		("Expected Type", po::value<int>(&config.calc_expected_type)->default_value(1), "")
		("Failure Type", po::value<int>(&config.failure_type)->default_value(1), "")
		("Scenario Generator", po::value<int>(&config.scen_gen)->default_value(1), "")
		("Limited World", po::value<int>(&config.LW)->default_value(0), "")
		;
	po::variables_map vm;
	ifstream ifs(configname);
	if (!ifs)
	{
		cout << "Failed to open Configuration File" << endl;
		*error = 1;
		return config;
	}
	else
	{
		try {
			store(po::parse_config_file(ifs, file_options), vm);
			notify(vm);
		}
		catch (exception e)
		{
			cout << "Required values not included in config file" << endl; cout << "Press enter to continue" << endl; cin.get(); *error = 1;
		}
	}
	return config;
}

