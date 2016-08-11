#ifndef NETWORKPROJECT_UTILBOOST_H_
#define NETWORKPROJECT_UTILBOOST_H_

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace BPT = boost::property_tree;
typedef BPT::ptree PTree;

class UtilBoost {
public:
	static PTree jsonFile2Ptree(const std::string&);
};

#endif // NETWORKPROJECT_UTILBOOST_H_
