#include "step_calc_static_feature.h"

#include <algorithm>
#include <vector>
#include <string>
#include <memory>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "unistd.h"

#include "runner.h"
#include "runner_manager.h"
#include "parameter.h"
#include "undirected_graph.h"
#include "util.h"
#include "util_boost.h"

StepCalcStaticFeature StepCalcStaticFeature::_step_calc_static_feature;

StepCalcStaticFeature::StepCalcStaticFeature(): Runner() {
	_short_options = "h";
	_long_options = new struct option[100]{
		{"help",		no_argument,		NULL, OPT_HELP},
//		{"net_type",	required_argument,	NULL, OPT_NET_TYPE},
		{"net_inroot",	required_argument,	NULL, OPT_NET_INROOT},
//		{"net_injson",	required_argument,	NULL, OPT_NET_INJSON},
//		{"net_volunteers",required_argument,NULL, OPT_NET_VOLUNTEERS},
		{"out_dir",		required_argument,	NULL, OPT_OUT_DIR},
//		{"disease",		required_argument,	NULL, OPT_DISEASE},
//		{"infect_rate",	required_argument,	NULL, OPT_INFECT_RATE},
//		{"infect_rate_seconds",	required_argument,	NULL, OPT_INFECT_RATE_SECONDS},
//		{"infectious_rate",	required_argument,	NULL, OPT_INFECTIOUS_RATE},
//		{"infectious_rate_seconds",	required_argument,	NULL, OPT_INFECTIOUS_RATE_SECONDS},
//		{"recover_rate",	required_argument,	NULL, OPT_RECOVER_RATE},
//		{"recover_rate_seconds",	required_argument,	NULL, OPT_RECOVER_RATE_SECONDS},
//		{"seconds_per_weight",	required_argument, NULL, OPT_SECONDS_PER_WEIGHT},
//		{"seconds_per_step",	required_argument, NULL, OPT_SECONDS_PER_STEP},
//		{"source_count",	required_argument, NULL, OPT_SOURCE_COUNT},
//		{"snapshot_coverage",	required_argument, NULL, OPT_SNAPSHOT_COVERAGE},
//		{"max_sim_days",	required_argument, NULL, OPT_MAX_SIM_DAYS},
//		{"repeat_times",	required_argument, NULL, OPT_REPEAT_TIMES},
//		{"source_identification_method", required_argument, NULL, OPT_SRC_IDN_METHOD},
//		{"source_identification_knowntime", required_argument, NULL, OPT_SRC_IDN_KNOWNTIME},
//		{"start_part",	required_argument,	NULL,	OPT_START_PART},
//		{"end_part",	required_argument,	NULL,	OPT_END_PART},
//		{"last_parts_threshold",	required_argument, NULL,	OPT_LAST_PARTS_THRESHOLD},
		{"calc_edges",		no_argument,		NULL,	OPT_CALC_EDGES},
		{NULL,			0,					NULL,  0 } //must end with {0, 0, 0, 0}
	};
	RunnerManager::instance()->install("step_calc_static_feature", this);
}

void StepCalcStaticFeature::help() {
	std::cout << "\nFunctionality: calc static feature(s) of a graph from gml file.";
	std::cout << "Option list:\n";
	std::cout << "\t* --help (or -h): [ no argument ] print this help information.\n";
	std::cout << "\t* --net_inroot: [ string argument ] gml file of the undirected graph.\n";
	std::cout << "\t* --out_dir: [ string argument ] output file of the results. If this argument is not given, i.e. out_dir.size() == 0, then the results will be shown on screen.\n";
	std::cout << "\t* --calc_edges: [ no argument ] calculate edges or not. In the output, this value follows 'edges' + 1 space.\n";
	std::cout << "\t* To add more features ...\n";
	std::cout << std::endl;

}

int StepCalcStaticFeature::run(const Parameter& para) {
	UGraphPtr g = std::make_shared<UndirectedGraph>(para.get_net_inroot(), NodeSet());
	std::vector<std::string> res = get_calc_res(g, para);
	if (para.get_out_dir().size()) {
		Util::writeLines(res, para.get_out_dir());
	} else {
		Util::writeLines(res, std::cout);
	}
	return 0;
}

std::vector<std::string> StepCalcStaticFeature::get_calc_res(const UGraphPtr& g, const Parameter& para) {
	std::vector<std::string> res;
	if (para.get_calc_edges()) {
		res.push_back("edges " + std::to_string(g->get_edge_size()));
	}
	return res;
}


