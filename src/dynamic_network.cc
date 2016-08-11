#include "dynamic_network.h"

#include "network.h"
#include "parameter.h"
#include "undirected_graph.h"
#include "neighbor_info.h"
#include "util.h"
#include "util_boost.h"

DynamicNetwork::DynamicNetwork(const Parameter& para) : Network(para) {
	PTree pt = UtilBoost::jsonFile2Ptree(para.get_net_injson());
	PTree data = pt.get_child("data");
	for (PTree::iterator it = data.begin(); it != data.end(); ++it) {
		std::string filename = para.get_net_inroot();
		filename.append("/").append(it->second.get<std::string>("file"));
		_ugraphs.push_back(std::make_shared<UndirectedGraph>(filename, Util::vec2unset(_nodes)));
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
		_merged_graph = UndirectedGraph::merge(_ugraphs);
	}
	return _merged_graph;
}
