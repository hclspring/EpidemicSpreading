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
//		std::vector<int> range = Util::parseIntegers(range_str, ':');
		std::string suffix = it->second.get<std::string>("suffix");
//		Util::checkTrue(range.size() == 2, "Error: range must contain two integers splitted by ':'.");
//		Util::checkTrue(range[0] < range[1], "Error: range[0] must be less than range[1].");
		std::vector<int> range_ints = Util::expandInterval(std::vector<std::string>{range_str}, ':');
		std::vector<std::string> temp_res;
		for (int i = 0; i < range_ints.size(); ++i) {
			temp_res.push_back(subdir + "/" + prefix + Util::getPartString(range_ints[i], part_str_length) + suffix);
		}
		res.push_back(temp_res);
	}
	return res;
}

std::vector<std::vector<std::string>> UtilBoost::parsePtree2layerAI2(const PTree& ptree, const int& part_str_length) {
	std::vector<std::vector<std::string>> res;
	PTree ptree_dates = ptree.get_child("dates");
	PTree ptree_parts = ptree.get_child("parts");
//	PTree ptree_pattern = ptree.get_child("pattern");
	std::vector<std::string> date_strs, raw_parts;
	std::vector<int> part_ints;
	for (PTree::iterator it = ptree_dates.begin(); it != ptree_dates.end(); ++it) {
		date_strs.push_back(it->second.get<std::string>(""));
	}
	for (PTree::iterator it = ptree_parts.begin(); it != ptree_parts.end(); ++it) {
		raw_parts.push_back(it->second.get<std::string>(""));
	}
	part_ints = Util::expandInterval(raw_parts, ':');
	std::string pattern_str = ptree.get<std::string>("pattern");
//	std::cout << pattern_str << std::endl;
	res.resize(date_strs.size(), std::vector<std::string>(part_ints.size()));
	for (int i = 0; i < res.size(); ++i) {
		for (int j = 0; j < res[i].size(); ++j) {
			res[i][j] = pattern_str;
			Util::replace_all(res[i][j], "${date}", date_strs[i]);
			Util::replace_all(res[i][j], "${part}", Util::getPartString(part_ints[i], part_str_length));
		}
	}
	return res;
}


