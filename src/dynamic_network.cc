#include "dynamic_network.h"

#include "network.h"
#include "parameter.h"
#include "undirected_graph.h"
#include "neighbor_info.h"
#include "util.h"
#include "util_boost.h"

DynamicNetwork::DynamicNetwork(const Parameter& para) : Network(para) {
	std::vector<std::string> filenames = UtilBoost::parsePtree1layer(
			UtilBoost::jsonFile2Ptree(para.get_net_injson()));
	for (int i = 0; i < filenames.size(); ++i) {
		_ugraphs.push_back(std::make_shared<UndirectedGraph>(
					para.get_net_inroot() + "/" + filenames[i], Util::vec2unset(_nodes)));
	}
	_days = _ugraphs.size();
	_parts = 1;
	if (_nodes.size() == 0) {
		NodeSet all_nodes;
		for (int i = 0; i < _ugraphs.size(); ++i) {
			all_nodes = Util::getUnion(all_nodes, _ugraphs[i]->get_node_names());
		}
		_nodes = Util::unset2vec(all_nodes);
		update_node_map();
		for (int i = 0; i < _ugraphs.size(); ++i) {
			_ugraphs[i]->add_nodes(Util::getDiff(all_nodes, _ugraphs[i]->get_node_names()));
		}
	}
}

DynamicNetwork::DynamicNetwork(const std::vector<std::shared_ptr<UndirectedGraph>>& x) : Network() {
	_ugraphs = x;
	get_merged_graph();
	_nodes = _merged_graph->get_node_names_vec();
	update_node_map();
	_days = _ugraphs.size();
	_parts = 1;
}

DynamicNetwork::~DynamicNetwork() {
}

NeighborList DynamicNetwork::get_neighbor_list(const std::string& nodename, const int& day_index, const int& part_index) const {
	if (check_day_index(day_index) && check_part_index(part_index)) {
		return _ugraphs[day_index]->get_neighbor_list(nodename);
	} else {
		std::cerr << "Error: day_index = " << day_index << " and part_index = " << part_index << "." << std::endl;
		exit(-1);
	}
}

std::shared_ptr<UndirectedGraph> DynamicNetwork::get_merged_graph() {
	if (!_merged_graph) {
		_merged_graph = UndirectedGraph::merge(_ugraphs, true); // with_weight = true
	}
	return _merged_graph;
}

std::shared_ptr<UndirectedGraph> DynamicNetwork::get_undirected_graph(const int& day_index, const int& part_index) const {
	if (check_day_index(day_index) && check_part_index(part_index)) {
		return _ugraphs[day_index];
	} else {
		return nullptr;
	}
}

size_t DynamicNetwork::get_size() const {
	return _ugraphs.size();
}

std::shared_ptr<UndirectedGraph>& DynamicNetwork::operator[](const size_t& index) {
	return this->at(index);
}

std::shared_ptr<UndirectedGraph>& DynamicNetwork::at(const size_t& index) {
	Util::checkTrue(index < _ugraphs.size(), "Error for accessing a graph with illegal index " + std::to_string(index) + ".");
	return _ugraphs[index];
}
	
std::list<NodeSet> DynamicNetwork::get_diff_nodes() const {
	Util::checkTrue(get_size() > 1, "Error: the dynamic network contains " + std::to_string(get_size()) + " graphs.");
	std::list<NodeSet> res;
	NodeSet current_set = _ugraphs[0]->get_node_names();
	for (int i = 1; i < get_size(); ++i) {
		NodeSet delta = Util::getDiff(_ugraphs[i]->get_node_names(), current_set);
		res.push_back(delta);
		for (NodeSet::iterator it = delta.begin(); it != delta.end(); ++it) {
			current_set.insert(*it);
		}
	}
	return res;
}

UGPVec DynamicNetwork::get_accumulating_graphs() const {
	UGPVec res;
	switch (get_size()) {
		case 0: return res;
		default: {
			res.push_back(std::make_shared<UndirectedGraph>(*(_ugraphs[0])));
			for (int i = 1; i < get_size(); ++i) {
				res.push_back(UndirectedGraph::merge(res[i-1], _ugraphs[i], true)); //with_weight = true
			}
			return res;
		}
	}
}

std::shared_ptr<DynamicNetwork> DynamicNetwork::get_evolved_networks_end_with_largest_component() {
	int n = get_size();
	Util::checkTrue(n > 0, "Error: DynamicNetwork contains no graph.");
	std::vector<UGraphPtr> ugp_vec(n);
	bool is_connected = at(n-1)->is_connected();
	if (is_connected == false) {
		ugp_vec[n-1] = at(n-1)->get_largest_component();
	} else {
		ugp_vec[n-1] = at(n-1);
	}
	NodeSet nodenames = ugp_vec[n-1]->get_node_names();
	ugp_vec[0] = at(0)->get_sub_graph(nodenames)->get_largest_component();
	Util::checkTrue(ugp_vec[0]->get_node_size() > 0, "Error for finding first graph no edges.");
	for (int j = 1; j < n; ++j) {
		NodeSet last_ugp_nodeset = ugp_vec[j - 1]->get_node_names();
		std::vector<UGraphPtr> cps = at(j)->get_sub_graph(nodenames)->get_components();
		std::vector<int> joint_size(cps.size());
		for (int k = 0; k < cps.size(); ++k) {
			joint_size[k] = Util::getIntersection(cps[k]->get_node_names(), last_ugp_nodeset).size();
		}
		int max_index = Util::getMaxIndex(joint_size);
		ugp_vec[j] = cps[max_index];
	}
	return std::make_shared<DynamicNetwork>(ugp_vec);
}


std::shared_ptr<DynamicNetwork> DynamicNetwork::get_simulated_evolving_networks(const std::string& type, const UGraphPtr& friendship_graph, CalProb cal_prob, const std::vector<double>& x) {
	int n = get_size();
	Util::checkTrue(n > 1, "Error: the dynamic network contains less than 2 graphs.");
	Util::checkTrue(at(n - 1)->is_connected(), "Error: the dynamic networks end with an unconnected graph.");
	typedef std::list<NodeSet> NodeSetList;
	NodeSetList diff_nodes = get_diff_nodes();
	UGraphPtr real_merged_graph = get_merged_graph();
	UGraphPtr my_merged_graph = std::make_shared<UndirectedGraph>(*(at(0)));

	std::vector<int> degree_vector = at(0)->get_degree_vector();
	NodeVec node_vector = at(0)->get_node_names_vec();

	std::vector<UGraphPtr> res (1, my_merged_graph);

	for (NodeSetList::iterator it = diff_nodes.begin(); 
			it != diff_nodes.end(); 
			++it) {
		if (it->size() > 0) {
			// record previous graph size
			int pre_size = degree_vector.size();
			my_merged_graph->add_nodes(*it);
			// determine the "actual" sequence of nodes to add
			NodeVec nodes_to_add = Util::unset2vec(*it);
			//std::cout << nodes_to_add.size() << std::endl;
			std::random_shuffle(nodes_to_add.begin(), nodes_to_add.end());
			// add nodes
			for (int j = 0; j < nodes_to_add.size(); ++j) {
				node_vector.push_back(nodes_to_add[j]);
			}
			// resize the degree vector
			degree_vector.resize(node_vector.size(), 0);

			// initialize edges to add
			EdgeList edge_list;

			// for every new node, determine the edge to add
			for (int j = 0; j < nodes_to_add.size(); ++j) {
				std::string new_node = nodes_to_add[j];
				// node j should connect to previous (pre_size + j) nodes
				for (int k = 0; k < pre_size + j; ++k) {
					std::string old_node = node_vector[k];
					double prob = cal_prob(type, friendship_graph, new_node, old_node, degree_vector[k], x);
					if (Util::gen_rand_double() < prob) {
						edge_list.push_back(std::make_tuple(old_node, new_node, 1));
						degree_vector[k]++;
						degree_vector[pre_size + j]++;
					}
				}
			}
			// add the edges
			my_merged_graph->add_edges(edge_list);
			//my_merged_graph->init_neighbors_vec();
		}
		res.push_back(my_merged_graph);
	}
	return std::make_shared<DynamicNetwork>(res);
}

double DynamicNetwork::get_avg_frobenius_norm_of_diff(DynamicNetwork& dn) {
	int n = get_size();
	Util::checkFalse(n == 0, "Error: this dynamic network contains no graph.");
	Util::checkTrue(dn.get_size() == n, "Error: dn and this contain different number of graphs.");
//	int m = at(n-1)->get_node_size();
//	Util::checkTrue(dn.at(n-1)->get_node_size() == m, "Error: dn and this contain different number of nodes in the last graph.");
	NodeSet this_node_set = get_merged_graph()->get_node_names();
	NodeSet dn_node_set = dn.get_merged_graph()->get_node_names();
//	NodeSet intersection_set = Util::getIntersection(this_node_set, dn_node_set);
	NodeSet union_set = Util::getUnion(this_node_set, dn_node_set);
	NodeVec node_vec = Util::unset2vec(union_set);
	double res = 0.0;
	for (int i = 0; i < n; ++i) {
		res += Util::getFrobeniusNorm(Util::getDiffMatrix(
					at(i)->get_unweighted_adjacency_matrix(node_vec),
					dn[i]->get_unweighted_adjacency_matrix(node_vec)));
	}
	res /= (n * node_vec.size());
	return res;
}

void DynamicNetwork::write_gml(const std::string& output_dir) {
	int a = system(std::string("mkdir -p " + output_dir).c_str());
	for (int i = 0; i < get_size(); ++i) {
		at(i)->write_graph_gml(output_dir + "/index_" + std::to_string(i) + ".gml");
	}
}

