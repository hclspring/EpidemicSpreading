#include "step_reindex_gml.h"

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
#include "util.h"
#include "util_boost.h"

StepReindexGml StepReindexGml::_step_reindex_gml;

StepReindexGml::StepReindexGml(): Runner() {
	_short_options = "h";
	_long_options = new struct option[100]{
		{"help",		no_argument,		NULL, OPT_HELP},
		{"net_inroot",	required_argument,	NULL, OPT_NET_INROOT},
		{"out_dir",		required_argument,	NULL, OPT_OUT_DIR},
		{NULL,			0,					NULL,  0 } //must end with {0, 0, 0, 0}
	};
	RunnerManager::instance()->install("step_reindex_gml", this);
}

void StepReindexGml::help() {
	std::cout << "\nFunctionality: Reindex the node ids in the gml file such that all nodes are indexed continuously from 0.";
	std::cout << "Option list:\n";
	std::cout << "\t* --help (or -h): [ no argument ] print this help information.\n";
	std::cout << "\t* --net_inroot: [ string argument ] input gml file.\n";
	std::cout << "\t* --out_dir: [ string argument ] output gml file.\n";
	std::cout << std::endl;
}

int StepReindexGml::run(const Parameter& para) {
	std::string infile = para.get_net_inroot();
	std::vector<std::string> lines = Util::readLines(infile);
	std::string outfile = para.get_out_dir();
	reindex_gml(lines);
	Util::writeLines(lines, outfile);
	return 0;
}

void StepReindexGml::reindex_gml(std::vector<std::string>& lines) {
	// start with id --> reindex
	// start with source --> modify id
	// start with target --> modify id
	typedef std::unordered_map<int, int> MapII;
	MapII map;
	int new_id = 0;
	for (int i = 0; i < lines.size(); ++i) {
		int j = get_first_nonblank_pos(lines[i], 0);
		if (starts_with(lines[i], j, std::string("id"))) {
			j = get_first_nonblank_pos(lines[i], j + 2);
			int k = j;
			int old_id = get_first_integer(lines[i], k);
			if (old_id >= 0) {
				lines[i].replace(j, k-j, std::to_string(new_id));
				map.insert(std::make_pair(old_id, new_id));
				++new_id;
			}
		} else if (starts_with(lines[i], j, std::string("source")) ||
				starts_with(lines[i], j, std::string("target"))) {
			j = get_first_nonblank_pos(lines[i], j + 6);
			int k = j;
			int old_id = get_first_integer(lines[i], k);
			Util::checkFalse(map.find(old_id) == map.end(), "Error: node id " + std::to_string(old_id) + " not exist in the map.");
			lines[i].replace(j, k-j, std::to_string(map.find(old_id)->second));
		}
	}
}

int StepReindexGml::get_first_nonblank_pos(const std::string& str, const int& start) {
	int res = start;
	while (res < str.size() && isblank(str[res])) ++res;
	return res;
}

bool StepReindexGml::starts_with(const std::string& str, const int& index, const std::string& str2) {
	int i = 0, j = index;
	if (str2.size() > str.size() - index) {
		return false;
	}
	while (i < str2.size() && j < str.size() && str[j] == str2[i]) {
		++i;
		++j;
	}
	if (i == str2.size()) {
		return true;
	} else {
		return false;
	}
}

int StepReindexGml::get_first_integer(const std::string& str, int& start) {
	int res = 0;
	if (start >= str.size() || !isdigit(str[start])) return -1;
	while (start < str.size() && isdigit(str[start])) {
		res = res * 10 + str[start] - '0';
		++start;
	}
	return res;
}

