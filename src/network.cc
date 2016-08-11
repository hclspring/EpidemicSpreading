#include "network.h"

#include "util.h"
#include "parameter.h"
#include "neighbor_info.h"
#include "undirected_graph.h"

Network::Network(const Parameter& para) {
	if (para.get_net_volunteers().size() > 0) {
		read_volunteers(para.get_net_volunteers());
	}
}

Network::~Network() {
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

