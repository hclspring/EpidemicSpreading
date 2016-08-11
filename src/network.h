#ifndef NETWORKPROJECT_NETWORK_H_
#define NETWORKPROJECT_NETWORK_H_

/*
 * Class Network is the basis class of StaticNetwork, DynamicNetwork, PearlNetwork
 */

#include <vector>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <list>
#include <memory>

class UndirectedGraph;
class Parameter;
class NeighborInfo;

typedef std::unordered_set<std::string> NodeSet;
typedef std::unordered_set<int> IndexSet;
typedef std::vector<std::string> NodeVec;
typedef std::unordered_map<std::string, int> NodeNameIndexMap;
typedef std::list<NeighborInfo> NeighborList;

class Network {
protected:
	NodeVec _nodes;
	NodeNameIndexMap _node_name2index;
	int _days;
	int _parts;
	std::shared_ptr<UndirectedGraph> _merged_graph = NULL;

public:
	Network(const Parameter& para);
	virtual ~Network();

	// get basic information of the graph
	int get_days() const;
	int get_parts() const;
	int get_node_size() const;
	NodeVec get_node_names() const;
	NodeSet get_node_names(const IndexSet& index_set) const;
	std::string get_node_name(const int& i) const;
	bool contains_node(const std::string& nodename) const;
	int get_node_index(const std::string& nodename) const; // NOTE: -1 for non-exist
	bool check_day_index(const int& day_index) const;
	bool check_part_index(const int& day_index) const;

	// get information of the structure
	virtual NeighborList get_neighbor_list(const std::string& nodename, const int& day_index, const int& part_index) const = 0;
	virtual std::shared_ptr<UndirectedGraph> get_merged_graph() = 0;
//	virtual std::shared_ptr<UndirectedGraph> get_undirected_graph(const int& day_index, const int& part_index) = 0;
	
protected:
	void update_node_map();
private:
	void read_volunteers(const std::string& vol_file);
};


#endif // NETWORKPROJECT_NETWORK_H_
