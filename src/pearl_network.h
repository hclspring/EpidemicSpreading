#ifndef NETWORKPROJECT_PEARLNETWORK_H_
#define NETWORKPROJECT_PEARLNETWORK_H_

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

class PearlNetwork : public Network {
private:
	std::vector<std::vector<std::shared_ptr<UndirectedGraph>>> _ugraphss;

public:
	PearlNetwork(const Parameter& para);
	~PearlNetwork();

	virtual NeighborList get_neighbor_list(const std::string& nodename, const int& day_index, const int& part_index) const override;
	virtual std::shared_ptr<UndirectedGraph> get_merged_graph() override;
	virtual std::shared_ptr<UndirectedGraph> get_undirected_graph(const int& day_index, const int& part_index) const override;
};


#endif // NETWORKPROJECT_PearlNETWORK_H_

