#ifndef NETWORKPROJECT_RUNNER_MANAGER_H_
#define NETWORKPROJECT_RUNNER_MANAGER_H_

#include <iostream>
#include <string>
#include <map>

class Runner;
typedef Runner* RunnerPtr;

class RunnerManager {
private:
	std::map<std::string, RunnerPtr> _runners;

public:
	// 根据命令返回适用的runner
	RunnerPtr get_runner(const std::string& name) const {
		return _runners.find(name) == _runners.end() ? NULL : _runners.find(name)->second;
	}

	void install(const std::string& name, RunnerPtr p) {
		_runners.insert(std::make_pair(name, p));
	}

	static RunnerManager* instance() {
		static RunnerManager manager;
		return &manager;
	}

	static int help() {
		std::cout << "\nUsage: ./netjob JOBNAME [OPTIONS]\n";
		std::cout << "\tPossible values of JOBNAME are:\n";
		std::cout << "\t* step_source_identification: for source identification. But this version has problems: 1) too slow; 2) different methods use different problem inputs.\n";
		std::cout << "\t* step_source_identification_v2: USE THIS VERSION instead.\n";
		std::cout << "\t* step_get_small_component: [ONLY FOR PEARL] for two-node components and three-node components, find the start part and parts of duration if the part and the duration satisfy some constraints.\n";
		std::cout << "\t* step_calc_static_feature: [ONLY FOR STATIC] calculate network feature\n";
		std::cout << "\t* step_evolve_large_component: [ONLY FOR DYNAMIC] find the best parameters to generate these large components\n";
		std::cout << "\t* step_simulate_evolve: [ONLY FOR DYNAMIC] simulate the evolving process, and calculate features of the networks\n";
		std::cout << "\t* step_reindex_gml: [ONLY FOR STATIC] reindex node id in the gml file so that ids are given continuously from 0\n";
		std::cout << "\t* step_select_nelder_mead_solutions: find the solutions which satisfy constraints.\n";
		std::cout << "\t* step_determine_relationship: determine the relationship between eating graph and friendship.\n";
		std::cout << "For help information of any JOBNAME, use the command:\n";
		std::cout << "\t\t ./netjob JOBNAME -h\n";
		std::cout << "or\n";
		std::cout << "\t\t ./netjob JOBNAME --help\n";
		std::cout << std::endl;
		return 1;
	}

};

#endif // NETWORKPROJECT_RUNNER_MANAGER_H_
