#include "static_network.h"

#include <cstdio>
#include "network.h"
#include "parameter.h"
#include "undirected_graph.h"
#include "neighbor_info.h"

#include "util.h"

StaticNetwork::StaticNetwork(const Parameter& para) : Network(para) {
	_ugraph = std::make_shared<UndirectedGraph>(para.get_net_inroot(), Util::vec2unset(_nodes));
	_days = 1;
	_parts = 1;
	if (_nodes.size() == 0) {
		_nodes = Util::unset2vec(_ugraph->get_node_names());
		update_node_map();
	}
}

StaticNetwork::~StaticNetwork() {
}

NeighborList StaticNetwork::get_neighbor_list(const std::string& nodename, const int& day_index, const int& part_index) const {
	if (check_day_index(day_index) && check_part_index(part_index)) {
		return _ugraph->get_neighbor_list(nodename);
	} else {
		std::cerr << "Error: day_index = " << day_index << " and part_index = " << part_index << "." << std::endl;
		exit(-1);
	}
}

std::shared_ptr<UndirectedGraph> StaticNetwork::get_merged_graph() {
	if (!_merged_graph) {
		_merged_graph = std::make_shared<UndirectedGraph>(*_ugraph);
	}
	return _merged_graph;
}
