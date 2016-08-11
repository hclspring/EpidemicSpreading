#include <iostream>

#include "util.h"
#include "util_boost.h"

PTree UtilBoost::jsonFile2Ptree(const std::string& json_file) {
	PTree pt;
	std::vector<std::string> lines = Util::readLines(json_file);
	std::string content;
	for (int i = 0; i < lines.size(); ++i) {
		content.append(lines[i]);
	}
	std::stringstream sstr(content);
	try {
		BPT::read_json(sstr, pt);
	} catch (BPT::ptree_error& e) {
		std::cerr << "Error parsing the json file: " << json_file << std::endl;
		exit(-1);
	}
	return pt;
}


