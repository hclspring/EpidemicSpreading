#include "undirected_graph.h"

#include <igraph.h>
#include "util.h"
#include "util_igraph.h"
#include "neighbor_info.h"

UndirectedGraph::UndirectedGraph() {
	igraph_empty(&_graph, 0, false);
}

UndirectedGraph::UndirectedGraph(const UndirectedGraph& original) {
	igraph_copy(&_graph, &(original._graph));
	_has_weight = original._has_weight;
	_node_name2index = original._node_name2index;
	_neighbors_vec = original._neighbors_vec;
	_node_name_attr = original._node_name_attr;
	_edge_weight_attr = original._edge_weight_attr;
}

//UndirectedGraph::UndirectedGraph(const std::string& filename_gml);

UndirectedGraph::UndirectedGraph(
		const std::string& filename_gml, 
		const NodeSet& volunteers
		) {
	read_graph_gml(filename_gml);
	if (volunteers.size() > 0) {
		NodeSet current_nodenames = get_node_names();
		delete_nodes(Util::getDiff(current_nodenames, volunteers));
		add_nodes(Util::getDiff(volunteers, current_nodenames));
	} else {
		update_node_name2index();
	}
	update_has_weight();
}

UndirectedGraph::~UndirectedGraph() {
	igraph_destroy(&_graph);
}

void UndirectedGraph::read_graph_gml(const std::string& filename_gml) {
	FILE* ifile = fopen(filename_gml.c_str(), "r");
	if (!ifile) {
		std::cerr << "Error reading file: " << filename_gml << std::endl;
		exit(-1);
	}
	igraph_read_graph_gml(&_graph, ifile);
	fclose(ifile);
	update_has_weight();
	update_members_with_nodes_changed();
}

int UndirectedGraph::get_node_size() const {
	return igraph_vcount(&_graph);
}

int UndirectedGraph::get_edge_size() const {
	return igraph_ecount(&_graph);
}

NodeSet UndirectedGraph::get_node_names() const {
	int n = get_node_size();
	NodeSet res;
	for (int i = 0; i < n; ++i) {
		res.insert(get_node_name(i));
	}
	return res;
}

NodeVec UndirectedGraph::get_node_names_vec() const {
	NodeVec res;
	for (int i = 0; i < get_node_size(); ++i) {
		res.push_back(get_node_name(i));
	}
	return res;
}

NodeSet UndirectedGraph::get_node_names(const IndexSet& index_set) const {
	NodeSet res;
	for (IndexSet::const_iterator it = index_set.begin(); it != index_set.end(); ++it) {
		res.insert(get_node_name(*it));
	}
	return res;
}

std::string UndirectedGraph::get_node_name(const int& index) const {
	if (igraph_cattribute_has_attr(&_graph, IGRAPH_ATTRIBUTE_VERTEX, _node_name_attr.c_str())) {
		std::string res = VAS(&_graph, _node_name_attr.c_str(), index);
		return res;
	} else {
		std::cerr << "Error: no vertex attribute called \"" << _node_name_attr << "\"." << std::endl;
	}
}

bool UndirectedGraph::contains_node(const std::string& nodename) const {
	NodeNameIndexMap::const_iterator it = _node_name2index.find(nodename);
	return (it != _node_name2index.end());
}

bool UndirectedGraph::contains_all_nodes(const NodeVec& node_vec) const {
	bool res = true;
	for (int i = 0; i < node_vec.size(); ++i) {
		res = res && contains_node(node_vec[i]);
	}
	return res;
}

int UndirectedGraph::get_node_index(const std::string& nodename) const {
	NodeNameIndexMap::const_iterator it = _node_name2index.find(nodename);
	if (it == _node_name2index.end()) {
		return -1;
	} else {
		return it->second;
	}
}

NeighborList UndirectedGraph::get_neighbor_list(const std::string& nodename) {
	int index = get_node_index(nodename);
	if (index < 0) { // Cannot find the node
		std::cerr << "Error: no vertex labeled \"" << nodename << "\"" << std::endl;
		exit(-1);
	} else if (std::get<0>(_neighbors_vec[index])) { // Information has been recorded.
		return std::get<1>(_neighbors_vec[index]);
	} else {
		igraph_vector_t neighbors;
		igraph_vector_init(&neighbors, 0);
		igraph_neighbors(&_graph, &neighbors, index, IGRAPH_ALL);
		int neighbors_size = igraph_vector_size(&neighbors);
		NeighborList nl;
		for (int i = 0; i < neighbors_size; ++i) {
			int neighbor_index = VECTOR(neighbors)[i];
			std::string neighbor_name = get_node_name(neighbor_index);
			double neighbor_weight = get_edge_weight(nodename, neighbor_name);
			NeighborInfo ni(neighbor_name, neighbor_weight);
			nl.push_back(ni);
		}
		std::get<0>(_neighbors_vec[index]) = true;
		std::get<1>(_neighbors_vec[index]) = nl;
		igraph_vector_destroy(&neighbors);
		return nl;
	}
}

int UndirectedGraph::get_degree(const std::string& nodename) {
	return get_neighbor_list(nodename).size();
}

std::unordered_map<std::string, int> UndirectedGraph::get_degree_map() {
	std::unordered_map<std::string, int> res;
	for (int i = 0; i < get_node_size(); ++i) {
		std::string nodename = get_node_name(i);
		res.insert(std::make_pair(nodename, get_degree(nodename)));
	}
	return res;
}

std::vector<int> UndirectedGraph::get_degree_vector() {
	std::vector<int> res(get_node_size());
	for (int i = 0; i < get_node_size(); ++i) {
		res[i] = get_degree(get_node_name(i));
	}
	return res;
}

bool UndirectedGraph::is_connected(const std::string& nodename1, const std::string& nodename2) const {
	// Time complexity of this implementation is O(log(d)).
	int nodeindex1 = get_node_index(nodename1);
	int nodeindex2 = get_node_index(nodename2);
	int eid;
	igraph_get_eid(&_graph, &eid, nodeindex1, nodeindex2, false, false);
	if (eid >= 0) {
		return true;
	} else {
		return false;
	}
}

bool UndirectedGraph::is_connected() const {
	igraph_bool_t res;
	igraph_is_connected(&_graph, &res, IGRAPH_STRONG);
	return res;
}

double UndirectedGraph::get_edge_weight(const std::string& nodename1, const std::string& nodename2) const {
	// Time complexity of this implementation is O(log(d)).
	int nodeindex1 = get_node_index(nodename1);
	int nodeindex2 = get_node_index(nodename2);
	int eid;
	igraph_get_eid(&_graph, &eid, nodeindex1, nodeindex2, false, false);
	if (eid < 0) {
		std::cerr << "Error: vertices \"" << nodename1 << "\" and \"" << nodename2 << "\" are not connected." << std::endl;
		exit(-1);
	} else if (_has_weight) {
		return static_cast<double>(EAN(&_graph, _edge_weight_attr.c_str(), eid));
	} else {
		return 1;
	}
}

NodeSet UndirectedGraph::get_same_component(const std::string& nodename) const {
	igraph_vector_t vt;
	if (!contains_node(nodename)) {
		std::cerr << "Error finding node \"" << nodename << "\" in the graph." << std::endl;
		exit(-1);
	}
	igraph_vector_init(&vt, 0);
	igraph_subcomponent(&_graph, &vt, get_node_index(nodename), IGRAPH_ALL);
	int n = igraph_vector_size(&vt);
	NodeSet res;
	for (int i = 0; i < n; ++i) {
		res.insert(get_node_name(VECTOR(vt)[i]));
	}
	igraph_vector_destroy(&vt);
	return res;
}

NodeVec UndirectedGraph::get_bfs_nodes(const std::string& start_node) {
	int n = get_node_size();
	NodeVec fathers(n, "");
	NodeVec res;
	std::queue<std::string> visit_q;
	int root_index = get_node_index(start_node);
	fathers[root_index] = get_node_name(root_index); // father of root is itself
	visit_q.push(start_node);
	res.push_back(start_node);
	while (!visit_q.empty()) {
		std::string cur_name = visit_q.front();
		visit_q.pop();
		NeighborList nl = get_neighbor_list(cur_name);
		for (NeighborList::iterator it = nl.begin(); it != nl.end(); ++it) {
			std::string neighbor_name = it->get_name();
			double weight = it->get_weight();
			int neighbor_index = get_node_index(neighbor_name);
			if (fathers[neighbor_index].size() <= 0) { // no parent, means haven't been visited.
				fathers[neighbor_index] = cur_name;
				visit_q.push(neighbor_name);
				res.push_back(neighbor_name);
			}
		}
	}
	return res;
}

EdgeList UndirectedGraph::get_bfs_edges(const std::string& start_node, NodeVec& fathers, AllChildrenVec& children) {
	fathers.clear();
	children.clear();
	int n = get_node_size();
	fathers.resize(n, "");
	children.resize(n, ChildrenList());
	EdgeList edge_list;
	std::queue<std::string> visit_q;
	int root_index = get_node_index(start_node);
	fathers[root_index] = get_node_name(root_index); // father of root is itself
	visit_q.push(start_node);
	while (!visit_q.empty()) {
		std::string cur_name = visit_q.front();
		int cur_index = get_node_index(cur_name);
		visit_q.pop();
		NeighborList nl = get_neighbor_list(cur_name);
		for (NeighborList::iterator it = nl.begin(); it != nl.end(); ++it) {
			std::string neighbor_name = it->get_name();
			double weight = it->get_weight();
			int neighbor_index = get_node_index(neighbor_name);
			if (fathers[neighbor_index].size() <= 0) { // no parent, means haven't been visited.
				fathers[neighbor_index] = cur_name;
				children[cur_index].push_back(neighbor_name);
				visit_q.push(neighbor_name);
				edge_list.push_back(std::make_tuple(cur_name, neighbor_name, weight));
			}
		}
	}
	return edge_list;
}

EdgeList UndirectedGraph::get_bfs_edges_with_necessary_extra_nodes(const std::string& start_node, const NodeSet& necessary_nodes, NodeNameIndexMap& bfs_node2index, NodeVec& fathers, AllChildrenVec& children) {
	NodeSet same_component_nodes = get_same_component(start_node);
	NodeSet diff_nodes = Util::getDiff(necessary_nodes, same_component_nodes);
	if (diff_nodes.size() > 0) {
		std::cerr << "Error for finding these necessary nodes not in the same component with the root: " << std::endl;
		std::for_each(diff_nodes.begin(), diff_nodes.end(),
			[] (const std::string& node_name) { std::cerr << node_name << std::endl; });
		exit(-1);
	}

	bfs_node2index.clear();
	fathers.clear();
	children.clear();

	NodeVec bfs_total_nodes, bfs_total_fathers;
	AllChildrenVec bfs_total_children;
	get_bfs_edges(start_node, bfs_total_fathers, bfs_total_children);

	EdgeList edge_list;
	NodeVec bfs_nodes;
	for (NodeSet::const_iterator it = necessary_nodes.begin(); it != necessary_nodes.end(); ++it) {
		std::string cur_node = *it;
		// Visit every nodes' father: if its father hasn't been visited, then insert the node and the edge.
		while (bfs_node2index.find(cur_node) == bfs_node2index.end()) {
			bfs_node2index.insert(std::make_pair(cur_node, bfs_nodes.size()));
			bfs_nodes.push_back(cur_node);
			// If cur_node is not the root, then add the edge between this node and its father.
			if (cur_node.compare(start_node) != 0) {
				std::string father = bfs_total_fathers[get_node_index(cur_node)];
				double weight = get_edge_weight(cur_node, father);
				edge_list.push_back(std::make_tuple(cur_node, father, weight));
				cur_node = father;
			}
		}
	}
	fathers.resize(bfs_nodes.size());
	children.resize(bfs_nodes.size());
	for (int i = 0; i < bfs_nodes.size(); ++i) {
		std::string name = bfs_nodes[i];
		int index_in_graph = get_node_index(name);
		fathers[i] = bfs_total_fathers[index_in_graph];
		ChildrenList total_children_list = bfs_total_children[index_in_graph];
		for (ChildrenList::iterator it = total_children_list.begin(); it != total_children_list.end(); ++it) {
			if (bfs_node2index.find(*it) != bfs_node2index.end()) {
				children[i].push_back(*it);
			}
		}
	}
	return edge_list;
}

bool UndirectedGraph::is_tree() const {
	int root_index = Util::gen_rand_int(get_node_size());
	if (is_connected() && 
			igraph_ecount(&_graph) == igraph_vcount(&_graph) - 1) {
		return true;
	} else {
		return false;
	}
}

AdjacencyMatrix UndirectedGraph::get_adjacency_matrix(const bool& weighted) {
	int n = get_node_size();
	AdjacencyMatrix res(n, AdjacencyRow(n, 0));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == j) {
				res[i][j] = 0;
			} else {
				std::string nodei = get_node_name(i);
				std::string nodej = get_node_name(j);
				if (is_connected(nodei, nodej)) {
					if (weighted) {
						res[i][j] = get_edge_weight(nodei, nodej);
					} else {
						res[i][j] = 1;
					}
				} else {
					res[i][j] = 0;
				}
			}
		}
	}
	return res;
}

AdjacencyMatrix UndirectedGraph::get_adjacency_matrix(const bool& weighted, const NodeVec& node_vec) {
	int n = node_vec.size();
	AdjacencyMatrix res(n, AdjacencyRow(n, 0));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == j) {
				res[i][j] = 0;
			} else if (!contains_node(node_vec[i]) || !contains_node(node_vec[j])) {
				res[i][j] = 0;
			} else {
				std::string nodei = node_vec[i];
				std::string nodej = node_vec[j];
				if (nodei.compare(nodej) != 0 && is_connected(nodei, nodej)) {
					if (weighted) {
						res[i][j] = get_edge_weight(nodei, nodej);
					} else {
						res[i][j] = 1;
					}
				} else {
					res[i][j] = 0;
				}
			}
		}
	}
	return res;
}

AdjacencyMatrix UndirectedGraph::get_unweighted_adjacency_matrix() {
	return get_adjacency_matrix(false);
}

AdjacencyMatrix UndirectedGraph::get_unweighted_adjacency_matrix(const NodeVec& node_vec) {
	return get_adjacency_matrix(false, node_vec);
}


AdjacencyMatrix UndirectedGraph::get_weighted_adjacency_matrix() {
	return get_adjacency_matrix(true);
}

AdjacencyMatrix UndirectedGraph::get_weighted_adjacency_matrix(const NodeVec& node_vec) {
	return get_adjacency_matrix(true, node_vec);
}

LaplacianMatrix UndirectedGraph::get_laplacian_matrix() {
	LaplacianMatrix res = get_adjacency_matrix(true);
	for (int i = 0; i < res.size(); ++i) {
		for (int j = 0; j < res[i].size(); ++j) {
			if (i == j) {
				res[i][j] = get_degree(get_node_name(i)) - res[i][j];
			} else {
				res[i][j] = -res[i][j];
			}
		}
	}
	return res;
}

int UndirectedGraph::calc_components(IndexVec& components) const {
	components.resize(get_node_size());
	igraph_vector_t membership;
	igraph_vector_init(&membership, get_node_size());
	int number_of_components;
	int igraph_res = igraph_clusters(&_graph, &membership, NULL, 
			&number_of_components, IGRAPH_STRONG);
	Util::checkFalse(igraph_res == IGRAPH_EINVAL, "IGRAPH: Invalid mode argument.");
	for (int i = 0; i < components.size(); ++i) {
		components[i] = VECTOR(membership)[i];
	}
	igraph_vector_destroy(&membership);
	return number_of_components;
}

std::vector<UGraphPtr> UndirectedGraph::get_components() {
	IndexVec component_ids;
	int number_of_components = calc_components(component_ids);
	typedef std::unordered_map<int, std::shared_ptr<IndexVec>> IndexVecPtrMap;
	IndexVecPtrMap index_ids_map = Util::get_vecptrmap(component_ids);
	IndexVecPtrMap::iterator it;
	std::vector<UGraphPtr> res;
	for (it = index_ids_map.begin(); it != index_ids_map.end(); ++it) {
		std::shared_ptr<IndexVec> indexVecPtr = it->second;
		NodeSet nodenames;
		for (size_t i = 0; i < indexVecPtr->size(); ++i) {
			nodenames.insert(get_node_name(indexVecPtr->at(i)));
		}
		res.push_back(get_sub_graph(nodenames));
	}
	return res;
}

UGraphPtr UndirectedGraph::get_largest_component() {
	IndexVec component_ids;
	int number_of_components = calc_components(component_ids);
	typedef std::unordered_map<int, std::shared_ptr<IndexVec>> IndexVecPtrMap;
	IndexVecPtrMap index_ids_map = Util::get_vecptrmap(component_ids);
	IndexVecPtrMap::iterator it, max_it;
	max_it = index_ids_map.begin();
	it = max_it;
	int max_v = max_it->second->size();
	for (++it; it != index_ids_map.end(); ++it) {
		if (it->second->size() > max_v) {
			max_v = it->second->size();
			max_it = it;
		}
	}
	std::shared_ptr<IndexVec> indexVecPtr = max_it->second;
	NodeSet nodenames;
	for (size_t i = 0; i < indexVecPtr->size(); ++i) {
		nodenames.insert(get_node_name(indexVecPtr->at(i)));
	}
	return get_sub_graph(nodenames);
}

bool UndirectedGraph::is_pair_component(const std::string& p1, const std::string& p2) {
	return (contains_node(p1) && contains_node(p2) &&
			is_connected(p1, p2) && 
			get_degree(p1) == 1 && get_degree(p2) == 1);
}

bool UndirectedGraph::is_triplet_component(const std::string& p1, const std::string& pmid, const std::string& p2) {
	return (contains_node(p1) && contains_node(p2) && contains_node(pmid) &&
			is_connected(p1, pmid) && is_connected(p2, pmid) &&
			get_degree(p1) == 1 && get_degree(p2) == 1 && get_degree(pmid) == 2);
}

bool UndirectedGraph::is_line_without_others(const NodeVec& node_vec) {
	if (!contains_all_nodes(node_vec)) {
		return false;
	}
	if (node_vec.size() == 0) {
		return false;
	} else if (node_vec.size() == 1) { 
		return get_degree(node_vec[0]) == 0;
	} else {
		int n = node_vec.size();
		bool res = (get_degree(node_vec[0]) == 1) && (get_degree(node_vec[n-1]) == 1);
		for (int i = 1; i < n - 1; ++i) {
			res = res && (get_degree(node_vec[i]) == 2) && is_connected(node_vec[i-1], node_vec[i]);
		}
		res = res && is_connected(node_vec[n-2], node_vec[n-1]);
		return res;
	}
}

bool UndirectedGraph::is_clique(const NodeVec& node_vec) {
	if (!contains_all_nodes(node_vec)) {
		return false;
	} else if (node_vec.size() == 0) {
		return false;
	} else {
		for (int i = 0; i < node_vec.size(); ++i) {
			for (int j = i + 1; j < node_vec.size(); ++j) {
				if (!is_connected(node_vec[i], node_vec[j])) {
					return false;
				}
			}
		}
		return true;
	}
}

bool UndirectedGraph::is_clique_without_others(const NodeVec& node_vec) {
	return is_clique(node_vec) &&
		get_same_component(node_vec[0]).size() == node_vec.size();
}

int UndirectedGraph::get_shortest_distance(const std::string& nodename1, const std::string& nodename2) {
	igraph_vs_t vertex_set1, vertex_set2;
	igraph_vs_1(&vertex_set1, get_node_index(nodename1));
	igraph_vs_1(&vertex_set2, get_node_index(nodename2));
	igraph_matrix_t imt;
	igraph_matrix_init(&imt, 1, 1);
	igraph_shortest_paths(&_graph, &imt, vertex_set1, vertex_set2, IGRAPH_ALL);
	int result = MATRIX(imt, 0, 0);
	igraph_matrix_destroy(&imt);
	igraph_vs_destroy(&vertex_set1);
	igraph_vs_destroy(&vertex_set2);
	return result;
}

std::vector<std::vector<int>> UndirectedGraph::get_shortest_distance(const NodeVec& nodes1, const NodeVec& nodes2) {
	int n1 = nodes1.size(), n2 = nodes2.size();
	std::vector<std::vector<int>> res(n1, std::vector<int>(n2, 0));
	igraph_vector_t vec1, vec2;
	igraph_vector_init(&vec1, n1);
	igraph_vector_init(&vec2, n2);
	for (int i = 0; i < n1; ++i) {
		VECTOR(vec1)[i] = get_node_index(nodes1[i]);
	}
	for (int i = 0; i < n2; ++i) {
		VECTOR(vec2)[i] = get_node_index(nodes2[i]);
	}
	igraph_vs_t vs1, vs2;
	igraph_vs_vector(&vs1, &vec1);
	igraph_vs_vector(&vs2, &vec2);
	igraph_matrix_t im;
	igraph_matrix_init(&im, n1, n2);
	igraph_shortest_paths(&_graph, &im, vs1, vs2, IGRAPH_ALL);
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			res[i][j] = MATRIX(im, i, j);
		}
	}
	igraph_vs_destroy(&vs1);
	igraph_vs_destroy(&vs2);
	igraph_vector_destroy(&vec1);
	igraph_vector_destroy(&vec2);
	igraph_matrix_destroy(&im);
	return res;
}

double UndirectedGraph::get_betweenness(const std::string& nodename) {
	int index = get_node_index(nodename);
	if (index < 0) {
		std::cerr << "Error for node \"" << nodename << "\" not existing." << std::endl;
		exit(-1);
	}
	igraph_vs_t vs;
	igraph_vs_1(&vs, index);
	igraph_vector_t bet;
	igraph_vector_init(&bet, 0);
	igraph_betweenness(&_graph, &bet, vs, false, NULL, true);
	double res = VECTOR(bet)[0];
	igraph_vs_destroy(&vs);
	igraph_vector_destroy(&bet);
	return res;
}

double UndirectedGraph::get_average_distance() {
	igraph_real_t res;
	igraph_average_path_length(&_graph, &res, IGRAPH_UNDIRECTED, true);
	return res;
}

double UndirectedGraph::get_clustering_coefficient() {
	igraph_real_t res;
	igraph_transitivity_undirected(&_graph, &res, IGRAPH_TRANSITIVITY_ZERO);
	return res;
}

double UndirectedGraph::get_edge_density() {
	int n = get_node_size();
	return get_edge_size() * 2.0 / (n * (n - 1)) ;
}

std::vector<int> UndirectedGraph::get_all_degree_int() {
	std::vector<int> res;
	for (int i = 0; i < get_node_size(); ++i) {
		res.push_back(get_degree(get_node_name(i)));
	}
	return res;
}

std::vector<double> UndirectedGraph::get_all_degree_double() {
	std::vector<double> res;
	for (int i = 0; i < get_node_size(); ++i) {
		res.push_back(get_degree(get_node_name(i)));
	}
	return res;
}

std::vector<double> UndirectedGraph::get_all_closeness() {
	int n = get_node_size();
	igraph_vector_t closenesses, indices;
	igraph_vs_t vids;
	igraph_vector_init(&closenesses, n);
	igraph_vector_init(&indices, n);
	for (int i = 0; i < n; ++i) {
		VECTOR(indices)[i] = i;
	}
	igraph_vs_vector(&vids, &indices);
	igraph_closeness(&_graph, &closenesses, vids, IGRAPH_ALL, NULL, false);
	std::vector<double> res(n, 0);
	for (int i = 0; i < n; ++i) {
		res[i] = VECTOR(closenesses)[i];
	}
	igraph_vs_destroy(&vids);
	igraph_vector_destroy(&closenesses);
	igraph_vector_destroy(&indices);
	return res;
}

std::vector<double> UndirectedGraph::get_all_clustering_coefficient() {
	int n = get_node_size();
	igraph_vector_t ccs, indices;
	igraph_vs_t vids;
	igraph_vector_init(&ccs, n);
	igraph_vector_init(&indices, n);
	for (int i = 0; i < n; ++i) {
		VECTOR(indices)[i] = i;
	}
	igraph_vs_vector(&vids, &indices);
	igraph_transitivity_local_undirected(&_graph, &ccs, vids, IGRAPH_TRANSITIVITY_ZERO);
	std::vector<double> res(n, 0);
	for (int i = 0; i < n; ++i) {
		res[i] = VECTOR(ccs)[i];
	}
	igraph_vs_destroy(&vids);
	igraph_vector_destroy(&ccs);
	igraph_vector_destroy(&indices);
	return res;
}


std::vector<double> UndirectedGraph::get_betweenness_all() {
	std::vector<double> res;
	igraph_vector_t bet;
	igraph_vector_init(&bet, 0);
	igraph_betweenness(&_graph, &bet, igraph_vss_all(), false, NULL, true);
	int n = igraph_vector_size(&bet);
	res.resize(n);
	for (int i = 0; i < n; ++i) {
		res[i] = VECTOR(bet)[i];
	}
	igraph_vector_destroy(&bet);
	return res;
}

std::vector<double> UndirectedGraph::get_betweenness(const NodeVec& nodenames) {
	igraph_vs_t vs;
	igraph_vector_t indices, bets;
	igraph_vector_init(&indices, nodenames.size());
	igraph_vector_init(&bets, nodenames.size());
	for (int i = 0; i < nodenames.size(); ++i) {
		int index = get_node_index(nodenames[i]);
		if (index < 0) {
			std::cerr << "Error for node \"" << nodenames[i] << "\" not existing." << std::endl;
			exit(-1);
		}
		VECTOR(indices)[i] = index;
	}
	igraph_vs_vector(&vs, &indices);
	igraph_betweenness(&_graph, &bets, vs, false, NULL, true);
	std::vector<double> res(nodenames.size(), 0);
	for (int i = 0; i < nodenames.size(); ++i) {
		res[i] = VECTOR(bets)[i];
	}
	igraph_vs_destroy(&vs);
	igraph_vector_destroy(&indices);
	igraph_vector_destroy(&bets);
	return res;
}

std::shared_ptr<UndirectedGraph> UndirectedGraph::merge(const std::vector<std::shared_ptr<UndirectedGraph>>& ugraphs, const bool& with_weight) {
	std::shared_ptr<UndirectedGraph> res = std::make_shared<UndirectedGraph>();
	NodeSet node_set;
	for (int i = 0; i < ugraphs.size(); ++i) {
		node_set = Util::getUnion(node_set, ugraphs[i]->get_node_names());
	}
	res->add_nodes(node_set);

//	res->update_members_with_nodes_changed();
	int n = res->get_node_size();
	std::vector<NeighborList> neighbors(n, NeighborList());
	for (int t = 0; t < ugraphs.size(); ++t) {
		// In graph t
		for (int i = 0; i < n; ++i) {
			// For node i in result graph
			std::string nodename = res->get_node_name(i);
			if (ugraphs[t]->contains_node(nodename)) {
				NeighborList nl = ugraphs[t]->get_neighbor_list(nodename); // get NeighborList of nodename in graph t
				for (NeighborList::iterator jt = nl.begin(); jt != nl.end(); ++jt) {
					std::string neighbor_name = jt->get_name();
					double weight = (with_weight ? jt->get_weight() : 1.0);
					NeighborList::iterator it;
					for (it = neighbors[i].begin(); it != neighbors[i].end(); ++it) {
						if (it->get_name().compare(neighbor_name) == 0) {
							if (with_weight) {
								it->set_weight(it->get_weight() + weight);
							}
							break;
						}
					}
					if (it == neighbors[i].end()) {
						neighbors[i].push_back(NeighborInfo(neighbor_name, weight));
					}
				}
			}
		}
	}
	std::list<EdgeInfo> edges;
	for (int i = 0; i < n; ++i) {
		std::string nodename1 = res->get_node_name(i);
		NeighborList nl = neighbors[i];
		for (NeighborList::iterator it = nl.begin(); it != nl.end(); ++it) {
			std::string nodename2 = it->get_name();
			double weight = it->get_weight();
			int index2 = res->get_node_index(nodename2);
			if (index2 > i) {
				edges.push_back(std::make_tuple(nodename1, nodename2, weight));
			}
		}
	}
	res->update_has_weight();
	res->add_edges(edges);
	return res;
}

std::shared_ptr<UndirectedGraph> UndirectedGraph::merge(const std::shared_ptr<UndirectedGraph>& ugraph1, const std::shared_ptr<UndirectedGraph>& ugraph2, const bool& with_weight) {
	std::vector<UGraphPtr> x {ugraph1, ugraph2};
	return merge(x, with_weight);
}

std::shared_ptr<UndirectedGraph> UndirectedGraph::get_sub_graph(const NodeSet& nodenames) {
	if (Util::getDiff(nodenames, get_node_names()).size() > 0) {
		std::cerr << "Error for extra nodes in the subgraph but not in the original graph." << std::endl;
		exit(-1);
	}
	std::shared_ptr<UndirectedGraph> res = std::make_shared<UndirectedGraph>();
	res->add_nodes(nodenames);
	int n = nodenames.size();
	std::list<EdgeInfo> edges;
	for (int i = 0; i < n; ++i) {
		std::string nodename1 = res->get_node_name(i);
		NeighborList nl = get_neighbor_list(nodename1);
		for (NeighborList::iterator it = nl.begin(); it != nl.end(); ++it) {
			std::string nodename2 = it->get_name();
			double weight = it->get_weight();
			int index2 = res->get_node_index(nodename2);
			if (index2 > i) {
				edges.push_back(std::make_tuple(nodename1, nodename2, weight));
			}
		}
	}
	res->update_has_weight();
	res->add_edges(edges);
	return res;
}

void UndirectedGraph::write_graph_gml(const std::string& filename_gml) const {
	FILE* ofile = fopen(filename_gml.c_str(), "w");
	if (!ofile) {
		std::cerr << "Error writting file: " << filename_gml << std::endl;
		exit(-1);
	}
	igraph_write_graph_gml(&_graph, ofile, NULL, NULL);
	fclose(ofile);
}

void UndirectedGraph::delete_nodes(const NodeSet& nodenames) {
	igraph_vector_t iv;
	int n = nodenames.size();
	igraph_vector_init(&iv, n);
	int i = 0;
	for (NodeSet::const_iterator it = nodenames.begin(); it != nodenames.end(); ++it) {
		int index = get_node_index(*it);
		if (index < 0) {
			std::cerr << "Error trying to delete a node named \"" << *it << "\" which does not exist." << std::endl;
			exit(-1);
		}
		VECTOR(iv)[i] = get_node_index(*it);
		++i;
	}
	igraph_vs_t selectors = igraph_vss_vector(&iv);

	igraph_delete_vertices(&_graph, selectors);
	igraph_vector_destroy(&iv);
	igraph_vs_destroy(&selectors);

	update_members_with_nodes_changed();
}

void UndirectedGraph::add_nodes(const NodeSet& nodenames) {
	int n = get_node_size();
	int delta = nodenames.size();
	igraph_add_vertices(&_graph, delta, NULL);
	NodeSet::const_iterator it;
	int i;
	for (i = 0, it = nodenames.begin(); it != nodenames.end(); ++i, ++it) {
		if (contains_node(*it)) {
			std::cerr << "Error trying add a node named \"" << *it << "\" which already exits." << std::endl;
			exit(-1);
		}
		SETVAS(&_graph, _node_name_attr.c_str(), n + i, (*it).c_str());
	}
	update_members_with_nodes_changed();
}

void UndirectedGraph::add_edge(const std::string& nodename1, const std::string& nodename2, const double& weight) {
	if (check_to_add_edge(nodename1, nodename2, weight)) {
		int nodeindex1 = get_node_index(nodename1);
		int nodeindex2 = get_node_index(nodename2);
		igraph_add_edge(&_graph, nodeindex1, nodeindex2);
		int eid;
		igraph_get_eid(&_graph, &eid, nodeindex1, nodeindex2, false, false);
		if (_has_weight) {
			SETEAN(&_graph, _edge_weight_attr.c_str(), eid, weight);
		}
		init_neighbors_vec();
	}
}

void UndirectedGraph::add_edges(const EdgeList& edges) {
	igraph_vector_t ie;
	igraph_vector_init(&ie, 2 * edges.size());
	typedef std::tuple<int, int, double> TempTuple;
	std::vector<TempTuple> temp_tuple;
	int i = 0;
	for (EdgeList::const_iterator it = edges.begin(); it != edges.end(); ++it, ++i) {
		std::string node_from = std::get<0>(*it);
		std::string node_to = std::get<1>(*it);
		double weight = std::get<2>(*it);
		if (check_to_add_edge(node_from, node_to, weight)) {
			int index_from = get_node_index(node_from);
			int index_to = get_node_index(node_to);
			VECTOR(ie)[2 * i] = index_from;
			VECTOR(ie)[2 * i + 1] = index_to;
			temp_tuple.push_back(std::make_tuple(index_from, index_to, weight));
		}
	}
	igraph_add_edges(&_graph, &ie, NULL);
	igraph_vector_destroy(&ie);
	if (_has_weight) {
		int eid;
		for (i = 0; i < temp_tuple.size(); ++i) {
			igraph_get_eid(&_graph, &eid, std::get<0>(temp_tuple[i]), std::get<1>(temp_tuple[i]), false, false);
			SETEAN(&_graph, _edge_weight_attr.c_str(), eid, std::get<2>(temp_tuple[i]));
		}
	}
	//init_neighbors_vec();
}

void UndirectedGraph::update_members_with_nodes_changed() {
	update_node_name2index();
	init_neighbors_vec();
}

void UndirectedGraph::update_node_name2index() {
	_node_name2index.clear();
	int n = get_node_size();
	for (int i = 0; i < n; ++i) {
		_node_name2index.insert(std::make_pair(get_node_name(i), i));
	}
}

void UndirectedGraph::init_neighbors_vec() {
	_neighbors_vec.clear();
	_neighbors_vec.resize(get_node_size(), std::make_tuple(false, NeighborList()));
}

void UndirectedGraph::update_has_weight() {
	_has_weight = igraph_cattribute_has_attr(&_graph, IGRAPH_ATTRIBUTE_EDGE, _edge_weight_attr.c_str());
}

bool UndirectedGraph::check_to_add_edge(const std::string& nodename1, const std::string& nodename2, const double& weight) {
	if (!contains_node(nodename1)) {
		std::cerr << "Error trying to add an edge from a node named \"" << nodename1 << "\" which does not exist." << std::endl;
		exit(-1);
	} else if (!contains_node(nodename2)) {
		std::cerr << "Error trying to add an edge to a node named \"" << nodename2 << "\" which does not exist." << std::endl;
		exit(-1);
	} else if (is_connected(nodename1, nodename2)) {
		std::cerr << "Error trying to add an edge between two nodes named \"" << nodename1 << "\" and \"" << nodename2 << "\" which are already connected." << std::endl;
		exit(-1);
	} else if (weight <= 0) {
		std::cerr << "Error trying to add an edge with weight = " << weight << "." << std::endl;
		exit(-1);
	} else {
		return true;
	}
}

