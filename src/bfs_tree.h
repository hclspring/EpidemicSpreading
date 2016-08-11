#ifndef NETWORKPROJECT_BFSTREE_H_
#define NETWORKPROJECT_BFSTREE_H_

#include <vector>
#include <string>
#include <list>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <tuple>
#include <igraph.h>

#include "undirected_graph.h"

typedef std::vector<int> IndexVec;
typedef std::list<std::string> ChildrenList;
typedef std::vector<std::string> NodeVec;
typedef std::vector<ChildrenList> AllChildrenVec;
typedef std::unordered_set<std::string> NodeSet;
typedef std::tuple<std::string, std::string, double> EdgeInfo;
typedef std::list<EdgeInfo> EdgeList;

/*
 * BFSTree is the BFS tree of an undirected graph.
 * If the graph is not connected, the tree only contains the connected component with the root node.
 */

class BFSTree : public UndirectedGraph {
private:
	std::string _root_name;
	NodeVec _parents; // Parent of root is itself.
	AllChildrenVec _children;
public:
	BFSTree(UndirectedGraph& graph, const std::string& root_name);
	BFSTree(UndirectedGraph& graph, const std::string& root_name, const NodeSet& necessary_nodes); // Generate a BFS tree from root_name with all "necessary_nodes" and only those necessary extra nodes.
	~BFSTree();

	NodeSet get_children(const std::string& nodename) const;
	std::string get_parent(const std::string& nodename) const; // parent of root is itself.
	std::string get_root() const;
	bool is_root(const std::string& nodename) const;
	bool is_leaf(const std::string& nodename) const;
	NodeSet get_all_descendants(const std::string& nodename) const;
	NodeSet get_subtree_nodes(const std::string& nodename) const;

	// calculate properties of the structure
//	std::vector<double> calc_rumor_centralities() const;
	std::vector<int> calc_depths() const; // Depth is the distance to root.
	std::vector<int> calc_heights() const; // Height is the distance to the bottom.
	void calc_height(const std::string& node, std::vector<int>& heights) const;

	// change the structure
	void rewire(const std::string& rewired_node, const std::string& old_neighbor, const std::string& new_neighbor);


};

#endif // NETWORKPROJECT_BFSTREE_H_
