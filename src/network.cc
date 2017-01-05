#include "network.h"

#include "util.h"
#include "parameter.h"
#include "neighbor_info.h"
#include "undirected_graph.h"

Network::Network() {
}

Network::Network(const Parameter& para) {
	if (para.get_net_volunteers().size() > 0) {
		read_volunteers(para.get_net_volunteers());
	}
}

Network::~Network() {
	_nodes.clear();
	_node_name2index.clear();
	_days = 0;
	_parts = 0;
	_merged_graph = nullptr;
}

int Network::get_days() const {
	return _days;
}

int Network::get_parts() const {
	return _parts;
}

int Network::get_node_size() const {
	return _nodes.size();
}

NodeVec Network::get_node_names() const {
	return _nodes;
}

NodeSet Network::get_node_names(const IndexSet& index_set) const {
	NodeSet res;
	for (IndexSet::const_iterator it = index_set.begin(); it != index_set.end(); ++it) {
		res.insert(get_node_name(*it));
	}
	return res;
}

std::string Network::get_node_name(const int& i) const {
	return _nodes[i];
}

bool Network::contains_node(const std::string& nodename) const {
	NodeNameIndexMap::const_iterator it = _node_name2index.find(nodename);
	return (it != _node_name2index.end());
}

int Network::get_node_index(const std::string& nodename) const {
	NodeNameIndexMap::const_iterator it = _node_name2index.find(nodename);
	if (it == _node_name2index.end()) {
		return -1;
	} else {
		return it->second;
	}
}

bool Network::check_day_index(const int& day_index) const {
	return day_index >= 0 && day_index < _days;
}

bool Network::check_part_index(const int& part_index) const {
	return part_index >= 0 && part_index < _parts;
}

ContactInfo Network::get_contact_info(const std::string& node1, const std::string& node2) {
	ContactInfo res;
	Util::checkTrue(get_merged_graph()->contains_node(node1), "Error: node " + node1 + " does not exist.");
	Util::checkTrue(get_merged_graph()->contains_node(node2), "Error: node " + node2 + " does not exist.");
	std::get<0>(res) = node1;
	std::get<1>(res) = node2;
	DoubleVecVecType vec_vec(get_days(), std::vector<double>(get_parts()));
	for (int i = 0; i < vec_vec.size(); ++i) {
		for (int j = 0; j < vec_vec[i].size(); ++j) {
			vec_vec[i][j] = get_undirected_graph(i, j)->is_connected(node1, node2) ? get_undirected_graph(i, j)->get_edge_weight(node1, node2) : 0;
		}
	}
	std::get<2>(res) = std::move(vec_vec);
	return res;
}

std::vector<ContactInfo> Network::get_contact_info(const std::string& nodename) {
	Util::checkTrue(get_merged_graph()->contains_node(nodename), "Error: node " + nodename + " does not exist.");
	NeighborList neighbor_list = get_merged_graph()->get_neighbor_list(nodename);
	std::vector<ContactInfo> res;
	for (NeighborList::iterator it = neighbor_list.begin(); it != neighbor_list.end(); ++it) {
		res.push_back(get_contact_info(nodename, it->get_name()));
	}
	return res;
}

std::vector<ContactInfo> Network::get_all_contact_info() {
	NodeVec all_nodes_vec = get_merged_graph()->get_node_names_vec();
	std::sort(all_nodes_vec.begin(), all_nodes_vec.end());
	std::vector<ContactInfo> res;
	for (int i = 0; i < all_nodes_vec.size(); ++i) {
		NeighborList neighbor_list = get_merged_graph()->get_neighbor_list(all_nodes_vec[i]);
		for (NeighborList::iterator it = neighbor_list.begin(); it != neighbor_list.end(); ++it) {
			if (all_nodes_vec[i].compare(it->get_name()) > 0) {
				res.push_back(get_contact_info(all_nodes_vec[i], it->get_name()));
			}
		}
	}
	return res;
}

void Network::write_all_contact_info(const std::vector<ContactInfo>& all_contact_info, const std::string& output) const {
	std::string mkdir_command = "mkdir -p " + output;
	int temp = system(mkdir_command.c_str());
	for (int i = 0; i < all_contact_info.size(); ++i) {
		std::string filename = std::get<0>(all_contact_info[i]) + "-" + std::get<1>(all_contact_info[i]);
		Util::writeLines(Util::getLines(std::get<2>(all_contact_info[i]), " "), output + "/" + filename);
	}
}


void Network::update_node_map() {
	_node_name2index.clear();
	for (int i = 0; i < _nodes.size(); ++i) {
		_node_name2index.insert(std::make_pair(_nodes[i], i));
	}
}

void Network::read_volunteers(const std::string& vol_file) {
	_nodes = Util::readLines(vol_file);
	update_node_map();
}

