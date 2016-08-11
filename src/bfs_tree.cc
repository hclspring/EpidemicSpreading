#include "bfs_tree.h"

#include "util.h"

BFSTree::BFSTree(UndirectedGraph& graph, const std::string& root_name) : UndirectedGraph() {
	NodeSet nodes_in_same_components = graph.get_same_component(root_name);
	add_nodes(nodes_in_same_components);
	// add the edges
	NodeVec parents;
	AllChildrenVec children;
	EdgeList edge_list = graph.get_bfs_edges(root_name, parents, children);
	add_edges(edge_list);
	// set _root_name, _parents, _children;
	_root_name = root_name;
	int n = get_node_size();
	_parents.resize(n);
	_children.resize(n);
	for (int i = 0; i < n; ++i) {
		std::string node_name = get_node_name(i);
		int index_in_graph = graph.get_node_index(node_name);
		_parents[i] = parents[index_in_graph];
		_children[i] = children[index_in_graph];
	}
}

BFSTree::BFSTree(UndirectedGraph& graph, const std::string& root_name, const NodeSet& necessary_nodes) : UndirectedGraph() {
	NodeNameIndexMap bfs_name2index;
	NodeVec fathers;
	AllChildrenVec children;
	EdgeList edge_list = graph.get_bfs_edges_with_necessary_extra_nodes(root_name, necessary_nodes, bfs_name2index, fathers, children);
	// add the nodes
	NodeSet nodes_to_add;
	for (NodeNameIndexMap::iterator it = bfs_name2index.begin(); it != bfs_name2index.end(); ++it) {
		nodes_to_add.insert(it->first);
	}
	add_nodes(nodes_to_add);
	// add the edges
	add_edges(edge_list);
	// set _root_name, _parents, _children
	_root_name = root_name;
	int n = get_node_size();
	for (int i = 0; i < n; ++i) {
		std::string node_name = get_node_name(i);
		int index_in_map = bfs_name2index.find(node_name)->second;
		_parents.push_back(fathers[index_in_map]);
		_children.push_back(children[index_in_map]);
	}
}

BFSTree::~BFSTree() {
	igraph_destroy(&_graph);
}

NodeSet BFSTree::get_children(const std::string& nodename) const {
	NodeSet res;
	int index = get_node_index(nodename);
	if (index < 0) {
		std::cerr << "Error for node \"" << nodename << "\" does not exist." << std::endl;
		exit(-1);
	}
	ChildrenList children = _children[get_node_index(nodename)];
	for_each (children.begin(), children.end(), 
			[&res, this] (const std::string& x) { res.insert(x); });
	return res;
}

std::string BFSTree::get_parent(const std::string& nodename) const {
	int index = get_node_index(nodename);
	if (index < 0) {
		std::cerr << "Error for node \"" << nodename << "\" does not exist." << std::endl;
		exit(-1);
	}
	return _parents[index];
}

std::string BFSTree::get_root() const {
	return _root_name;
}

bool BFSTree::is_root(const std::string& nodename) const {
	return nodename.compare(_root_name) == 0;
}

bool BFSTree::is_leaf(const std::string& nodename) const {
	return _children[get_node_index(nodename)].size() == 0;
}

NodeSet BFSTree::get_all_descendants(const std::string& nodename) const {
	NodeSet res;
	std::queue<std::string> q;
	q.push(nodename);
//	res.insert(nodename);
	while (!q.empty()) {
		std::string name = q.front();
		q.pop();
		NodeSet children = get_children(name);
		for (NodeSet::iterator it = children.begin(); it != children.end(); ++it) {
			q.push(*it);
			res.insert(*it);
		}
	}
	return res;
}

NodeSet BFSTree::get_subtree_nodes(const std::string& nodename) const {
	NodeSet res = get_all_descendants(nodename);
	res.insert(nodename);
	return res;
}

std::vector<int> BFSTree::calc_depths() const {
	int n = get_node_size();
	std::vector<int> depths(n, -1);
	std::queue<int> q;
	int root_index = get_node_index(_root_name);
	q.push(root_index);
	depths[root_index] = 0;
	while (!q.empty()) {
		int index = q.front();
		int depth = depths[index];
		q.pop();
		ChildrenList children = _children[index];
		for (ChildrenList::iterator it = children.begin(); it != children.end(); ++it) {
			std::string child = *it;
			int index = get_node_index(child);
			q.push(index);
			depths[index] = depth + 1;
		}
	}
	return depths;
}

std::vector<int> BFSTree::calc_heights() const {
	int n = get_node_size();
	std::vector<int> heights(n, -1);
	calc_height(get_root(), heights);
	return heights;
}

void BFSTree::calc_height(const std::string& node, std::vector<int>& heights) const {
	int index = get_node_index(node);
	if (index < 0) {
		std::cerr << "Error for no node named \"" << node << "\"." << std::endl;
		exit(-1);
	}
	if (is_leaf(node)) {
		heights[index] = 0;
		return;
	}
	ChildrenList children = _children[index];
	int max_height = 0;
	for (ChildrenList::iterator it = children.begin(); it != children.end(); ++it) {
		calc_height(*it, heights);
		int height = heights[get_node_index(*it)];
		if (height > max_height) {
			max_height = height;
		}
	}
	heights[index] = max_height + 1;
	return;
}

void BFSTree::rewire(const std::string& rewired_node, const std::string& old_neighbor, const std::string& new_neighbor) {
	// check whether rewired_node is connected to old_neighbor currently.
	if (is_connected(rewired_node, old_neighbor) == false) {
		std::cerr << "Error for node \"" << rewired_node << "\" and \"" << old_neighbor << "\" not connected." << std::endl;
		exit(-1);
	}
	// check whether old_neighbor and new_neighbor are same.
	if (old_neighbor.compare(new_neighbor) == 0) {
		return;
	}
	// check whether rewired_node is not connected to new_neighbor currently.
	if (is_connected(rewired_node, new_neighbor)) {
		std::cerr << "Error for node \"" << rewired_node << "\" and \"" << new_neighbor << "\" connected currently." << std::endl;
		exit(-1);
	}
	// check whether old_neighbor is the parent of rewired_node (otherwise after rewiring, the graph is not a tree).
	if (old_neighbor.compare(get_parent(rewired_node)) != 0) {
		std::cerr << "Error for node \"" << old_neighbor << "\" not being the parent of node \"" << rewired_node << "\"." << std::endl;
		exit(-1);
	}
	int rewired_index = get_node_index(rewired_node);
	int old_index = get_node_index(old_neighbor);
	int new_index = get_node_index(new_neighbor);
	igraph_add_edge(&_graph, rewired_index, new_index);
	igraph_es_t edges_to_delete;
	igraph_es_pairs_small(&edges_to_delete, false, rewired_index, old_index, -1);
	igraph_delete_edges(&_graph, edges_to_delete);
	igraph_es_destroy(&edges_to_delete);
	init_neighbors_vec();
	_parents[rewired_index] = new_neighbor;
	_children[new_index].push_back(rewired_node);
	_children[old_index].remove(rewired_node);
	return;
}



