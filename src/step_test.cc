#include "step_test.h"

#include <algorithm>
#include <vector>
#include <string>
#include <memory>

#include "runner.h"
#include "runner_manager.h"
#include "parameter.h"
#include "undirected_graph.h"
#include "network.h"
#include "pearl_network.h"
#include "static_network.h"
#include "dynamic_network.h"
#include "simulator.h"
#include "source_identification.h"
#include "util.h"

StepTest StepTest::_step_test;

StepTest::StepTest(): Runner() {
	_short_options = "";
	_long_options = new struct option[100]{
		{"net_type",	required_argument,	NULL, OPT_NET_TYPE},
		{"net_inroot",	required_argument,	NULL, OPT_NET_INROOT},
		{"net_injson",	required_argument,	NULL, OPT_NET_INJSON},
		{"net_volunteers",required_argument,NULL, OPT_NET_VOLUNTEERS},
		{"disease",		required_argument,	NULL, OPT_DISEASE},
		{"infect_rate",	required_argument,	NULL, OPT_INFECT_RATE},
		{"infect_rate_seconds",	required_argument,	NULL, OPT_INFECT_RATE_SECONDS},
		{"infectious_rate",	required_argument,	NULL, OPT_INFECTIOUS_RATE},
		{"infectious_rate_seconds",	required_argument,	NULL, OPT_INFECTIOUS_RATE_SECONDS},
		{"recover_rate",	required_argument,	NULL, OPT_RECOVER_RATE},
		{"recover_rate_seconds",	required_argument,	NULL, OPT_RECOVER_RATE_SECONDS},
		{"seconds_per_weight",	required_argument, NULL, OPT_SECONDS_PER_WEIGHT},
		{"seconds_per_step",	required_argument, NULL, OPT_SECONDS_PER_STEP},
		{"source_count",	required_argument, NULL, OPT_SOURCE_COUNT},
		{"snapshot_coverage",	required_argument, NULL, OPT_SNAPSHOT_COVERAGE},
		{"max_sim_days",	required_argument, NULL, OPT_MAX_SIM_DAYS},
		{"repeat_times",	required_argument, NULL, OPT_REPEAT_TIMES},
		{NULL,			0,					NULL,  0 } //must end with {0, 0, 0, 0}
	};
	RunnerManager::instance()->install("step_test", this);
}

int StepTest::run(const Parameter& para) {
	std::shared_ptr<Network> net;
	switch (para.get_net_type()) {
		case STATIC: net = std::make_shared<StaticNetwork>(para); break;
		case DYNAMIC: net = std::make_shared<DynamicNetwork>(para); break;
		case PEARL: net = std::make_shared<PearlNetwork>(para); break;
	}
//	UndirectedGraph* g = net->get_merged_graph();
//	g = g->get_sub_graph({"P001", "P002", "P081", "P139"});
//	g->write_graph_gml("test/sub_merge_res.gml");
//	delete g;

//	std::unordered_set<int> initial_infected_indices = Util::gen_rand_indices(net->get_node_size(), 10);
	std::vector<int> all_indices;
	for (int i = 0; i < net->get_node_size(); ++i) {
		all_indices.push_back(i);
	}
	std::random_shuffle(all_indices.begin(), all_indices.end());
	std::vector<int> first_indices;
	for (int i = 0; i < para.get_source_count(); ++i) {
		first_indices.push_back(all_indices[i]);
	}
	std::unordered_set<int> initial_infected_indices = Util::vec2unset(first_indices);
	Simulator simulator(para, net->get_node_size(), initial_infected_indices);
	std::vector<DiseaseStage> res_stages = simulator.get_sim_res(*net, para);
	NodeSet infected_nodes = simulator.get_nodeset_been_infected_from_sim_res(res_stages, para.get_disease(), *net);
	std::cout << "Seeds:";
	std::for_each(initial_infected_indices.begin(), initial_infected_indices.end(), [net] (const int& x) {std::cout << " " << net->get_node_name(x);});
	std::cout << std::endl;
	std::cout << "Infected:";
	std::for_each(infected_nodes.begin(), infected_nodes.end(), [] (const std::string& x) {std::cout << " " << x;});
	std::cout << std::endl;
	std::cout << "Stages: ";
	std::for_each(infected_nodes.begin(), infected_nodes.end(), [net, res_stages] (const std::string& x) { std::cout << " " << x << "--" << res_stages[net->get_node_index(x)]; });
	std::cout << std::endl;
	//std::string source = SourceIdentification::calc_source_RG(*net, infected_nodes, para.get_infect_rate());
	//std::string source = SourceIdentification::calc_source_SSE(*net, infected_nodes);
	//std::string source = SourceIdentification::calc_source_SSEBFS(*net, infected_nodes);
	//std::string source = SourceIdentification::calc_source_SJC(*net, infected_nodes);
	//std::string source = SourceIdentification::calc_source_JCE(*net, infected_nodes);
	//std::string source = SourceIdentification::calc_source_DA(*net, infected_nodes);
	//std::string source = SourceIdentification::calc_source_UB(*net, infected_nodes, para.get_ub_r());
	/*
	int max_sim_days = para.get_max_sim_days();
	if (max_sim_days < 0) max_sim_days = 1000;
	std::string source = SourceIdentification::calc_source_DMP_SIR_unknowntime(*net, res_stages, max_sim_days, para);
	std::string source = SourceIdentification::calc_source_DMP_SIR_knowntime(*net, res_stages, max_sim_days, para);
	*/
	//std::string source = SourceIdentification::calc_source_MCSM_unknowntime(*net, infected_nodes, para.get_max_sim_days(), para);
	//std::string source = SourceIdentification::calc_source_MCSM_knowntime(*net, infected_nodes, para.get_max_sim_days(), para);
	std::string source = SourceIdentification::calc_source_NetSleuth(*net, infected_nodes);
	std::cout << "Estimated source: " << source << std::endl;
	double error_distance = SourceIdentification::calc_error_distance(*net, *infected_nodes.begin(), source);
	std::cout << "Error distance: " << error_distance << std::endl;
	std::cout << "Success!" << std::endl;
	return 0;
}
