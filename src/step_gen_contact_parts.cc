#include "step_gen_contact_parts.h"

#include <algorithm>
#include <vector>
#include <string>
#include <memory>
#include "unistd.h"

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

StepGenContactParts StepGenContactParts::_step_gen_contact_parts;

StepGenContactParts::StepGenContactParts(): Runner() {
	_short_options = "h";
	_long_options = new struct option[100]{
		{"help",		no_argument,		NULL, OPT_HELP},
//		{"net_type",	required_argument,	NULL, OPT_NET_TYPE},
		{"net_inroot",	required_argument,	NULL, OPT_NET_INROOT},
		{"net_injson",	required_argument,	NULL, OPT_NET_INJSON},
		{"net_volunteers",required_argument,NULL, OPT_NET_VOLUNTEERS},
		{"merge_parts",	required_argument,	NULL, OPT_MERGE_PARTS},
		{"out_dir",		required_argument,	NULL, OPT_OUT_DIR},
		{NULL,			0,					NULL,  0 } //must end with {0, 0, 0, 0}
	};
	RunnerManager::instance()->install("step_gen_contact_parts", this);
}

void StepGenContactParts::help() {
	std::cout << "\nFunctionality: generate yes/no of every contact in all parts; net type is restricted as PEARL.\n";
	std::cout << "Option list:\n";
	std::cout << "\t* --help (or -h): [ no argument ] print this help information.\n";
	std::cout << "\t* --net_inroot: [ string argument ] root directory of the networks.\n";
	std::cout << "\t* --net_injson: [ string argument ] json file of the networks.\n";
	std::cout << "\t* --net_volunteers: [ string argument ] volunteer file of the networks.\n";
	std::cout << "\t* --merge_parts: [ string argument ] number of parts to be merged in one element of the resulted matrix.\n";
	std::cout << "\t* --out_dir: [ string argument ] output directory of the results.\n";
	std::cout << std::endl;
}


int StepGenContactParts::run(const Parameter& para) {
	// Read network
	std::shared_ptr<Network> net = std::make_shared<PearlNetwork>(para);
	std::vector<ContactInfo> all_contact_info = net->get_all_contact_info(para.get_merge_parts());
	net->write_all_contact_info(all_contact_info, para.get_out_dir());
	return 0;
}

