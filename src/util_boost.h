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

//#define BOOST_NO_CXX11_SCOPED_ENUMS
//#include <boost/filesystem.hpp>
//#undef BOOST_NO_CXX11_SCOPED_ENUMS

namespace BPT = boost::property_tree;
typedef BPT::ptree PTree;

class UtilBoost {
public:
	// About json
	static PTree jsonFile2Ptree(const std::string&);
	static std::vector<std::string> parsePtree1layer(const PTree& ptree); // used for default
	static std::vector<std::string> parsePtree1layerAI(const PTree& ptree, const int& part_str_length);
	static std::vector<std::vector<std::string>> parsePtree2layer(const PTree& ptree);
	static std::vector<std::vector<std::string>> parsePtree2layerAI(const PTree& ptree, const int& part_str_length); // used for default
	static std::vector<std::vector<std::string>> parsePtree2layerAI2(const PTree& ptree, const int& part_str_length);

};

#endif // NETWORKPROJECT_UTILBOOST_H_
