#ifndef NETWORKPROJECT_UNDIRECTEDGRAPH_H_
#define NETWORKPROJECT_UNDIRECTEDGRAPH_H_

#include <tuple>
#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>
#include <memory>

#include <igraph.h>

#include "neighbor_info.h"

class NeighborInfo;

typedef std::unordered_map<std::string, int> NodeNameIndexMap;
typedef std::unordered_set<std::string> NodeSet;
typedef std::vector<std::string> NodeVec;
typedef std::unordered_set<int> IndexSet;
typedef std::list<NeighborInfo> NeighborList;
typedef std::tuple<std::string, std::string, double> EdgeInfo;
typedef std::list<EdgeInfo> EdgeList;
typedef std::vector<int> IndexVec;
typedef std::list<std::string> ChildrenList;
typedef std::vector<ChildrenList> AllChildrenVec;
typedef std::vector<double> AdjacencyRow;
typedef std::vector<AdjacencyRow> AdjacencyMatrix;
typedef std::vector<double> LaplacianRow;
typedef std::vector<LaplacianRow> LaplacianMatrix;

/* 
 * IMPORTANT: 
 * You MUST write the following sentence in the very beginning of the main function.
 * igraph_i_set_attribute_table(&igraph_cattribute_table);
 */

class UndirectedGraph {
	typedef std::tuple<bool, NeighborList> NeighborsTuple;

protected:
	igraph_t _graph;
	
	bool _has_weight = false;
	NodeNameIndexMap _node_name2index;
	std::vector<NeighborsTuple> _neighbors_vec;
	std::string _node_name_attr = "label";
	std::string _edge_weight_attr = "weight";

public:
	UndirectedGraph();
	UndirectedGraph(const UndirectedGraph& original);
	UndirectedGraph(const std::string& filename_gml, const NodeSet& volunteers);
	~UndirectedGraph();


	// get basic information of the graph
	int get_node_size() const;
	int get_edge_size() const;
	NodeSet get_node_names() const;
	NodeVec get_node_names_vec() const;
	NodeSet get_node_names(const IndexSet& index_set) const;
	std::string get_node_name(const int& index) const;
	bool contains_node(const std::string& nodename) const;
	int get_node_index(const std::string& nodename) const; // NOTE: -1 for non-exist

	// get information about the structure
	NeighborList get_neighbor_list(const std::string& nodename); // NOTE: nodename MUST exist.
	int get_degree(const std::string& nodename); // NOTE: nodename MUST exist;
	bool is_connected(const std::string& nodename1, const std::string& nodename2) const;
	bool is_connected() const;
	double get_edge_weight(const std::string& nodename1, const std::string& nodename2) const; // NOTE: nodename1 and nodename2 MUST be connected.
	NodeSet get_same_component(const std::string& nodename) const;
	NodeVec get_bfs_nodes(const std::string& start_node);
	EdgeList get_bfs_edges(const std::string& start_node, NodeVec& fathers, AllChildrenVec& children); // Father of start_node is itself.
	EdgeList get_bfs_edges_with_necessary_extra_nodes(const std::string& start_node, const NodeSet& necessary_nodes, NodeNameIndexMap& bfs_node2index, NodeVec& fathers, AllChildrenVec& children); //Father of start_node is itself.
	bool is_tree() const;
	AdjacencyMatrix get_adjacency_matrix();
	LaplacianMatrix get_laplacian_matrix();

	// get structure properties of the graph
	int get_shortest_distance(const std::string& nodename1, const std::string& nodename2);
	std::vector<std::vector<int>> get_shortest_distance(const NodeVec& nodes1, const NodeVec& nodes2);
	double get_betweenness(const std::string& nodename);
	std::vector<double> get_betweenness_all();
	std::vector<double> get_betweenness(const NodeVec& nodenames);

public:
	// about generating graph
	static std::shared_ptr<UndirectedGraph> merge(const std::vector<std::shared_ptr<UndirectedGraph>>& ugraphs);
	std::shared_ptr<UndirectedGraph> get_sub_graph(const NodeSet& nodenames);
	void write_graph_gml(const std::string& filename_gml) const;

protected:
	// about reading file
	void read_graph_gml(const std::string& filename_gml);

	// about update graph_nodes
	void delete_nodes(const NodeSet& nodenames); // NOTE: nodenames MUST exist.
public:
	void add_nodes(const NodeSet& nodenames); // NOTE: nodenames MUST NOT exist.
protected:
	void add_edge(const std::string& nodename1, const std::string& nodename2, const double& weight); // NOTE: nodenames MUST exist and the edge MUST NOT exist, and the weight MUST be positive.
	void add_edges(const EdgeList& edges); // NOTE: same requirements as the function add_edge(*, *, *).

	// about update on class members after updating graph
	void update_members_with_nodes_changed();
	void update_node_name2index();
	void init_neighbors_vec();

	void update_has_weight();

	bool check_to_add_edge(const std::string& nodename1, const std::string& nodename2, const double& weight);
//	igraph_vs_t get_vertex_selectors(const NodeSet& nodenames);
};


#endif // NETWORKPROJECT_UNDIRECTEDGRAPH_H_
