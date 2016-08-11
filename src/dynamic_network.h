#ifndef NETWORKPROJECT_DYNAMICNETWORK_H_
#define NETWORKPROJECT_DYNAMICNETWORK_H_

#include <vector>
#include <string>
#include <memory>
#include <igraph.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "network.h"

class Parameter;
class UndirectedGraph;
class NeighborInfo;

typedef std::list<NeighborInfo> NeighborList;
typedef boost::property_tree::ptree PTree;

class DynamicNetwork : public Network {
private:
	std::vector<std::shared_ptr<UndirectedGraph>> _ugraphs;

public:
	DynamicNetwork(const Parameter& para);
	~DynamicNetwork();

	virtual NeighborList get_neighbor_list(const std::string& nodename, const int& day_index, const int& part_index) const override;
	virtual std::shared_ptr<UndirectedGraph> get_merged_graph() override;
};


#endif // NETWORKPROJECT_DYNAMICNETWORK_H_

