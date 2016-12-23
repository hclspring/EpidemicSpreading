#ifndef NETWORKPROJECT_STATICNETWORK_H_
#define NETWORKPROJECT_STATICNETWORK_H_

#include <list>
#include <string>
#include <memory>
#include <igraph.h>

#include "network.h"

class Parameter;
class UndirectedGraph;
class NeighborInfo;

typedef std::list<NeighborInfo> NeighborList;

class StaticNetwork : public Network {
private:
	std::shared_ptr<UndirectedGraph> _ugraph;

public:
	StaticNetwork(const Parameter& para);
	~StaticNetwork();

	virtual NeighborList get_neighbor_list(const std::string& nodename, const int& day_index, const int& part_index) const override;
	virtual std::shared_ptr<UndirectedGraph> get_merged_graph() override;
	virtual std::shared_ptr<UndirectedGraph> get_undirected_graph(const int& day_index, const int& part_index) const override;
};


#endif // NETWORKPROJECT_STATICNETWORK_H_

