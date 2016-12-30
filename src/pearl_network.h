#ifndef NETWORKPROJECT_PEARLNETWORK_H_
#define NETWORKPROJECT_PEARLNETWORK_H_

#include <vector>
#include <string>
#include <memory>
#include <tuple>
#include <igraph.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "network.h"

class Parameter;
class UndirectedGraph;
class NeighborInfo;

typedef std::list<NeighborInfo> NeighborList;
typedef boost::property_tree::ptree PTree;
typedef std::vector<std::vector<bool>> BoolVecVecType;
typedef std::tuple<std::string, std::string, BoolVecVecType> ContactInfo;

class PearlNetwork : public Network {
private:
	std::vector<std::vector<std::shared_ptr<UndirectedGraph>>> _ugraphss;

public:
	PearlNetwork(const Parameter& para);
	~PearlNetwork();

	virtual NeighborList get_neighbor_list(const std::string& nodename, const int& day_index, const int& part_index) const override;
	virtual std::shared_ptr<UndirectedGraph> get_merged_graph() override;
	virtual std::shared_ptr<UndirectedGraph> get_undirected_graph(const int& day_index, const int& part_index) const override;
	ContactInfo get_contact_info(const std::string& node1, const std::string& node2) const;
	std::vector<ContactInfo> get_contact_info(const std::string& nodename);
	std::vector<ContactInfo> get_all_contact_info();
};


#endif // NETWORKPROJECT_PearlNETWORK_H_

