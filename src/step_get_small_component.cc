#include "step_get_small_component.h"

#include <algorithm>
#include <vector>
#include <string>
#include <memory>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <set>
#include <tuple>
#include "unistd.h"

#include "runner.h"
#include "runner_manager.h"
#include "parameter.h"
#include "undirected_graph.h"
#include "network.h"
#include "pearl_network.h"
#include "static_network.h"
#include "dynamic_network.h"
//#include "simulator.h"
//#include "source_identification.h"
#include "util.h"
#include "util_boost.h"

StepGetSmallComponent StepGetSmallComponent::_step_get_small_component;

StepGetSmallComponent::StepGetSmallComponent(): Runner() {
	_short_options = "h";
	_long_options = new struct option[100]{
		{"help",		no_argument,		NULL, OPT_HELP},
//		{"net_type",	required_argument,	NULL, OPT_NET_TYPE},
		{"net_inroot",	required_argument,	NULL, OPT_NET_INROOT},
		{"net_injson",	required_argument,	NULL, OPT_NET_INJSON},
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
		{"start_part",	required_argument,	NULL,	OPT_START_PART},
		{"end_part",	required_argument,	NULL,	OPT_END_PART},
		{"last_parts_threshold",	required_argument, NULL,	OPT_LAST_PARTS_THRESHOLD},
		{NULL,			0,					NULL,  0 } //must end with {0, 0, 0, 0}
	};
	RunnerManager::instance()->install("step_get_small_component", this);
}

void StepGetSmallComponent::help() {
	std::cout << "\nFunctionality: find two-node components and three-node components, find the start part and parts of duration if the part and the duration satisfy some constraints.";
	std::cout << "Option list:\n";
	std::cout << "\t* --help (or -h): [ no argument ] print this help information.\n";
	std::cout << "\t* --net_inroot: [ string argument ] root directory of the networks.\n";
	std::cout << "\t* --net_injson: [ string argument ] json file of the networks.\n";
	std::cout << "\t* --out_dir: [ string argument ] output directory of the results.\n";
	std::cout << "\t* --start_part: [ nonnegative integer argument ] start part to consider. 0 denotes the first part.\n";
	std::cout << "\t* --end_part: [ nonnegative integer argument, larger than start_part ] end part to consider. 0 denotes the first part.\n";
	std::cout << "\t* --last_parts_threshold: [ positive integer argument ] part duration threshold, no less than which to be considered.\n";
	std::cout << std::endl;
}

int StepGetSmallComponent::run(const Parameter& para) {
	// Read network and get basic information
	PearlNetwork pn(para);
	int days = pn.get_days();
	int parts = pn.get_parts();
	NodeVec node_vec = pn.get_node_names();

	// deal with node pairs
	typedef std::tuple<int, std::string, std::string> CountConditionDirTuple;
	typedef std::set<CountConditionDirTuple> TupleSet;
	CountConditionDirTuple mytuple;
	TupleSet myset;
//	myset.insert(std::make_tuple(2, "is_line_without_others", "/pairs/"));
	myset.insert(std::make_tuple(3, "is_line_without_others", "/triplets/"));
	myset.insert(std::make_tuple(3, "is_clique_without_others", "/triangles/"));
	for (TupleSet::iterator it = myset.begin(); it != myset.end(); ++it) {
		int count = std::get<0>(*it);
		std::string condition = std::get<1>(*it);
		std::string subdir = std::get<2>(*it);
		// generate all pairs of nodes with pair size(count) == pair_count
		std::vector<NodeVec> node_pairs = get_node_pairs(node_vec, count, condition);
		// generate the directory to save all results
		std::string out_dir = para.get_out_dir() + subdir;
		int _temp = system(std::string("mkdir -p " + out_dir).c_str());
		// For every pair of nodes
		for (int i = 0; i < node_pairs.size(); ++i) {
			// get the parts which satisfy the requirements
			ResVec res_vec = test_component_parts(node_pairs[i], para.get_start_part(), para.get_end_part(), para.get_last_parts_threshold(), pn, condition);
			// write the result into a file only when such parts exist.
			if (res_vec.size() > 0) {
				std::string filename = get_out_filename(node_pairs[i]);
				write_result(out_dir + filename, res_vec);
			}
		}
	}
	return 0;
}

std::vector<Tuple2Str> StepGetSmallComponent::get_node_pairs(const NodeVec& node_vec) {
	std::vector<Tuple2Str> res;
	for (int i = 0; i < node_vec.size(); ++i) {
		for (int j = i + 1; j < node_vec.size(); ++j) {
			res.push_back(std::make_tuple(node_vec[i], node_vec[j]));
		}
	}
	return res;
}

std::vector<Tuple3Str> StepGetSmallComponent::get_node_triplets(const NodeVec& node_vec) {
	std::vector<Tuple3Str> res;
	for (int i = 0; i < node_vec.size(); ++i) {
		for (int j = 0; j < node_vec.size(); ++j) {
			for (int k = j + 1; k < node_vec.size(); ++k) {
				if (j != i && k != i) {
					res.push_back(std::make_tuple(node_vec[j], node_vec[i], node_vec[k]));
				}
			}
		}
	}
	return res;
}


std::vector<NodeVec> StepGetSmallComponent::get_node_pairs(const NodeVec& node_vec, const int& size, const std::string& condition) {
	std::vector<NodeVec> res;
	int n = node_vec.size();
	if (size > n || size == 0) {
		return res;
	} else if (condition.compare("is_line_without_others") == 0) {
		return get_node_lines(node_vec, size);
	} else if (condition.compare("is_clique_without_others") == 0) {
		return get_node_cliques(node_vec, size);
	} else {
		std::cerr << "Error: undefined condition \"" << condition << "\"" << std::endl;
		exit(-1);
	}
}

std::vector<NodeVec> StepGetSmallComponent::get_node_lines(const NodeVec& node_vec, const int& size) {
	std::vector<NodeVec> res;
	int n = node_vec.size();
	if (size == 1) {
		for (int i = 0; i < n; ++i) {
			res.push_back({node_vec[i]});
		}
		return res;
	} else if (size == 2) {
		for (int i = 0; i < n; ++i) {
			for (int j = i + 1; j < n; ++j) {
				res.push_back({node_vec[i], node_vec[j]});
			}
		}
		return res;
	} else if (size == 3) {
		for (int mid = 0; mid < n; ++mid) {
			for (int i = 0; i < n; ++i) {
				for (int j = i+1; j < n; ++j) {
					if (i != mid && j != mid) {
						res.push_back({node_vec[i], node_vec[mid], node_vec[j]});
					}
				}
			}
		}
		return res;
	} else {
		std::cerr << "Haven't implement the cases where pair size > 3." << std::endl;
		exit(-1);
	}
}

std::vector<NodeVec> StepGetSmallComponent::get_node_cliques(const NodeVec& node_vec, const int& size) {
	std::vector<NodeVec> res;
	int n = node_vec.size();
	if (size == 1 || size == 2) {
		return get_node_lines(node_vec, size);
	} else if (size == 3) {
		for (int i = 0; i < n; ++i) {
			for (int j = i + 1; j < n; ++j) {
				for (int k = j + 1; k < n; ++k) {
					res.push_back({node_vec[i], node_vec[j], node_vec[k]});
				}
			}
		}
		return res;
	}
}

StepGetSmallComponent::ResVec StepGetSmallComponent::test_component_parts(const NodeVec& t, const int& start_part, const int& end_part, const int& th, const PearlNetwork& pn, const std::string& condition) {
	ResVec res_vec;
	int days = pn.get_days();
	int parts = pn.get_parts();
	for (int i = 0; i < days; ++i) {
		int start = start_part, last = 0;
		PPVecPtr vec_ptr = std::make_shared<PPVec>(PPVec());
		for (int j = start_part; j <= end_part; ++j) {
			UGraphPtr g = pn.get_undirected_graph(i, j);
//			if (g != nullptr && g->is_line_without_others(t)) {
			if (g != nullptr && check_node_vec(g, t, condition)) {
				if (last == 0) {
					start = j;
					last = 1;
				} else {
					last++;
				}
			} else {
				if (last >= th) {
					vec_ptr->push_back(std::make_pair(start, last));
				}
				last = 0;
			}
		}
		if (last >= th) {
			vec_ptr->push_back(std::make_pair(start, last));
		}
		if (vec_ptr->size() > 0) {
			res_vec.push_back(std::make_pair(i, vec_ptr));
		}
	}
	return res_vec;
}

bool StepGetSmallComponent::check_node_vec(const UGraphPtr& g, const NodeVec& t, const std::string& condition) {
	if (condition.compare("is_line_without_others") == 0) {
		return g->is_line_without_others(t);
	} else if (condition.compare("is_clique_without_others") == 0) {
		return g->is_clique_without_others(t);
	} else {
		return false;
	}
}

std::string StepGetSmallComponent::get_out_filename(const Tuple2Str& t) {
	return std::get<0>(t) + "-" + std::get<1>(t) + ".txt";
}

std::string StepGetSmallComponent::get_out_filename(const Tuple3Str& t) {
	return std::get<0>(t) + "-" + std::get<1>(t) + "-" + std::get<2>(t) + ".txt";
}

std::string StepGetSmallComponent::get_out_filename(const NodeVec& t) {
	std::string res;
	Util::checkTrue(t.size() > 0, "Error: trying to get a filename with no node_name.");
	res = t[0];
	for (int i = 1; i < t.size(); ++i) {
		res.push_back('-');
		res.append(t[i]);
	}
	return res;
}

void StepGetSmallComponent::write_result(const std::string& res_file, const ResVec& res_vec) {
	if (res_vec.size() == 0) {
		return;
	}
	std::ofstream ofs(res_file.c_str());
	Util::checkFalse(ofs.fail(), "Error: write file " + res_file);
	for (int i = 0; i < res_vec.size(); ++i) {
		ofs << res_vec[i].first;
		PPVecPtr vec_ptr = res_vec[i].second;
		for (int j = 0; j < vec_ptr->size(); ++j) {
			ofs << " " << vec_ptr->at(j).first << " " << vec_ptr->at(j).second;
		}
		ofs << std::endl;
	}
	ofs.close();
}


