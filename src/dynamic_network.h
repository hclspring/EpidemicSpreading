#ifndef NETWORKPROJECT_DYNAMICNETWORK_H_
#define NETWORKPROJECT_DYNAMICNETWORK_H_

#include <vector>
#include <string>
#include <memory>
#include <cstddef>
#include <igraph.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "network.h"

class Parameter;
class UndirectedGraph;
class NeighborInfo;

typedef std::list<NeighborInfo> NeighborList;
typedef boost::property_tree::ptree PTree;
typedef std::unordered_set<std::string> NodeSet;
typedef std::shared_ptr<UndirectedGraph> UGraphPtr;
typedef std::vector<UGraphPtr> UGPVec;

class DynamicNetwork : public Network {
//	typedef double (*CalFd) (const int&, const std::vector<double>&);
	typedef double (*CalProb) (const std::string& type, const UGraphPtr&, const std::string&, const std::string&, const int&, const std::vector<double>&);
private:
	std::vector<std::shared_ptr<UndirectedGraph>> _ugraphs;

public:
	DynamicNetwork(const Parameter& para);
	DynamicNetwork(const std::vector<std::shared_ptr<UndirectedGraph>>&);
	~DynamicNetwork();

	virtual NeighborList get_neighbor_list(const std::string& nodename, const int& day_index, const int& part_index) const override;
	virtual std::shared_ptr<UndirectedGraph> get_merged_graph() override;
	virtual std::shared_ptr<UndirectedGraph> get_undirected_graph(const int& day_index, const int& part_index) const override;

	size_t get_size() const;
	std::shared_ptr<UndirectedGraph>& operator[](const size_t& index);
	std::shared_ptr<UndirectedGraph>& at(const size_t& index);
	std::list<NodeSet> get_diff_nodes() const;
	UGPVec get_accumulating_graphs() const;

	// About network evolvement
	std::shared_ptr<DynamicNetwork> get_evolved_networks_end_with_largest_component();
	std::shared_ptr<DynamicNetwork> get_simulated_evolving_networks(const std::string& type, const UGraphPtr& friendship_graph, CalProb cal_prob, const std::vector<double>& x);
	double get_avg_frobenius_norm_of_diff(DynamicNetwork& dn);

	void write_gml(const std::string& output_dir);
};


#endif // NETWORKPROJECT_DYNAMICNETWORK_H_

