#include "pearl_network.h"

#include "network.h"
#include "parameter.h"
#include "undirected_graph.h"
#include "neighbor_info.h"
#include "util.h"
#include "util_boost.h"

PearlNetwork::PearlNetwork(const Parameter& para) : Network(para) {
	PTree pt = UtilBoost::jsonFile2Ptree(para.get_net_injson());
	PTree data = pt.get_child("data");
	for (PTree::iterator it = data.begin(); it != data.end(); ++it) {
		std::string subdir = it->second.get<std::string>("dir");
//		std::cout << "subdir = " << subdir << std::endl;
		PTree subdata = it->second.get_child("data");
		std::vector<std::shared_ptr<UndirectedGraph>> gs;
		for (PTree::iterator jt = subdata.begin(); jt != subdata.end(); ++jt) {
			std::string filename = para.get_net_inroot();
			filename.append("/").append(subdir).append("/").append(jt->second.get<std::string>("file"));
			gs.push_back(std::make_shared<UndirectedGraph>(filename, Util::vec2unset(_nodes)));
		}
		_ugraphss.push_back(gs);
	}
	_days = _ugraphss.size();
	if (_days <= 0) {
		std::cerr << "Error with _days = " << _days << "." << std::endl;
		exit(-1);
	} else {
		_parts = _ugraphss[0].size();
		for (int i = 1; i < _ugraphss.size(); ++i) {
			if (_ugraphss[i].size() != _parts) {
				std::cerr << "Error with _parts = " << _parts << " but _ugraphss[" << i << "].size() = " << _ugraphss[i].size() << "." << std::endl;
				exit(-1);
			}
		}
	}
	if (_nodes.size() == 0) {
		NodeSet all_nodes;
		for (int i = 0; i < _ugraphss.size(); ++i) {
			for (int j = 0; j < _ugraphss[i].size(); ++j) {
				all_nodes = Util::getUnion(all_nodes, _ugraphss[i][j]->get_node_names());
			}
		}
		_nodes = Util::unset2vec(all_nodes);
		update_node_map();
		for (int i = 0; i < _ugraphss.size(); ++i) {
			for (int j = 0; j < _ugraphss[i].size(); ++j) {
				_ugraphss[i][j]->add_nodes(Util::getDiff(all_nodes, _ugraphss[i][j]->get_node_names()));
			}
		}
	}
}

PearlNetwork::~PearlNetwork() {
}

NeighborList PearlNetwork::get_neighbor_list(const std::string& nodename, const int& day_index, const int& part_index) const {
	if (check_day_index(day_index) && check_part_index(part_index)) {
		return _ugraphss[day_index][part_index]->get_neighbor_list(nodename);
	} else {
		std::cerr << "Error: day_index = " << day_index << " and part_index = " << part_index << "." << std::endl;
		exit(-1);
	}
}

std::shared_ptr<UndirectedGraph> PearlNetwork::get_merged_graph() {
	if (!_merged_graph) {
		std::vector<std::shared_ptr<UndirectedGraph>> temp(_days * _parts, NULL);
		for (int i = 0; i < _ugraphss.size(); ++i) {
			for (int j = 0; j < _ugraphss[i].size(); ++j) {
				temp[i * _parts + j] = _ugraphss[i][j];
			}
		}
		_merged_graph = UndirectedGraph::merge(temp);
	}
	return _merged_graph;
}
