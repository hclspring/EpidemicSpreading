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

// This function is used by default.
std::vector<std::string> UtilBoost::parsePtree1layer(const PTree& ptree) {
	std::vector<std::string> res;
	PTree data = ptree.get_child("data");
	for (PTree::iterator it = data.begin(); it != data.end(); ++it) {
		res.push_back(it->second.get<std::string>("file"));
	}
	return res;
}


std::vector<std::string> UtilBoost::parsePtree1layerAI(const PTree& ptree, const int& part_str_length) {
	std::vector<std::string> res;
	std::string prefix = ptree.get<std::string>("prefix");
	std::string range_str = ptree.get<std::string>("range");
	std::vector<int> range = Util::parseIntegers(range_str, ':');
	std::string suffix = ptree.get<std::string>("suffix");
	Util::checkTrue(range.size() == 2, "Error: range must contain two integers splitted by ':'.");
	Util::checkTrue(range[0] < range[1], "Error: range[0] must be less than range[1].");
	for (int i = range[0]; i < range[1]; ++i) {
		res.push_back(prefix + Util::getPartString(i, part_str_length) + suffix);
	}
	return res;
}

std::vector<std::vector<std::string>> UtilBoost::parsePtree2layer(const PTree& ptree) {
	std::vector<std::vector<std::string>> res;
	PTree data = ptree.get_child("data");
	for (PTree::iterator it = data.begin(); it != data.end(); ++it) {
		std::string subdir = it->second.get<std::string>("dir");
		PTree subdata = it->second.get_child("data");
		std::vector<std::string> temp_res;
		for (PTree::iterator jt = subdata.begin(); jt != subdata.end(); ++jt) {
			std::string str = jt->second.get<std::string>("file");
			temp_res.push_back(str);
		}
		res.push_back(Util::addSamePrefix(subdir + "/", temp_res));
	}
	return res;
}

// This function is used for default.
std::vector<std::vector<std::string>> UtilBoost::parsePtree2layerAI(const PTree& ptree, const int& part_str_length) {
	std::vector<std::vector<std::string>> res;
	PTree data = ptree.get_child("data");
	for (PTree::iterator it = data.begin(); it != data.end(); ++it) {
		std::string subdir = it->second.get<std::string>("dir");
		std::string prefix = it->second.get<std::string>("prefix");
		std::string range_str = it->second.get<std::string>("range");
		std::vector<int> range = Util::parseIntegers(range_str, ':');
		std::string suffix = it->second.get<std::string>("suffix");
		Util::checkTrue(range.size() == 2, "Error: range must contain two integers splitted by ':'.");
		Util::checkTrue(range[0] < range[1], "Error: range[0] must be less than range[1].");
		std::vector<std::string> temp_res;
		for (int i = range[0]; i < range[1]; ++i) {
			temp_res.push_back(subdir + "/" + prefix + Util::getPartString(i, part_str_length) + suffix);
		}
		res.push_back(temp_res);
	}
	return res;
}


