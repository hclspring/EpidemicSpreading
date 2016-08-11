#include "simulator.h"

#include <unordered_set>
#include <climits>

#include "network.h"
#include "static_network.h"
#include "dynamic_network.h"
#include "pearl_network.h"
#include "parameter.h"
#include "neighbor_info.h"
#include "util_constant.h"
#include "util.h"
#include "disease_dynamics.h"

Simulator::Simulator(const Parameter& para, const int& population, const IndexSet& source_set) {
	initialize_stages(para.get_disease(), population, source_set);
}

Simulator::~Simulator() {}

std::vector<DiseaseStage> Simulator::get_sim_res(const Network& contact_network, const Parameter& para) {
	_network_days = contact_network.get_days();
	_network_parts = contact_network.get_parts();
	_stages = _initial_stages;
	int max_sim_days = para.get_max_sim_days();
	if (max_sim_days < 0) {
		max_sim_days = INT_MAX;
	}
	for (int day = 0; day < max_sim_days; ++day) {
		for (int part = 0; part < _network_parts; ++part) {
			if (is_unstable(para.get_disease())) {
				_stages = spread(_stages, contact_network, para, day % _network_days, part % _network_parts);
			}
		}
	}
	return _stages;
}

NodeSet Simulator::get_nodeset_been_infected_from_sim_res(const std::vector<DiseaseStage>& sim_res, const DiseaseModel& disease, const Network& network) {
	NodeSet res;
	for (int i = 0; i < sim_res.size(); ++i) {
		if (UtilConstant::hasBeenInfected(disease, sim_res[i])) {
			res.insert(network.get_node_name(i));
		}
	}
	return res;
}

void Simulator::initialize_stages(const DiseaseModel& disease, const int& population, const IndexSet& source_set) {
	_initial_stages.clear();
	_initial_stages.resize(population, SUSCEPTIBLE);
	for (IndexSet::const_iterator it = source_set.begin(); it != source_set.end(); ++it) {
		_initial_stages[*it] = UtilConstant::getInitialInfectedStage(disease);
	}
	_stages = _initial_stages;
	_unstable_stages = UtilConstant::getUnstableStages(disease);
}

bool Simulator::is_unstable(const DiseaseModel& disease) const {
	return !UtilConstant::isStable(disease, _stages, _unstable_stages);
}

std::vector<DiseaseStage> Simulator::spread(const std::vector<DiseaseStage>& previous_stages, const Network& contact_network, const Parameter& para, const int& day, const int& part) {
	std::vector<DiseaseStage> next_stages(previous_stages.size());
	for (int i = 0; i < previous_stages.size(); ++i) {
		std::string nodename = contact_network.get_node_name(i);
		if (UtilConstant::canBeInfected(previous_stages[i])) {
			std::list<double> infect_seconds = get_infect_seconds(i, contact_network, para, day, part);
			next_stages[i] = DiseaseDynamics::get_next_status_bycontact(para, previous_stages[i], infect_seconds);
		} else {
			next_stages[i] = DiseaseDynamics::get_next_status_byself(para, previous_stages[i], para.get_seconds_per_step());
		}
	}
	return next_stages;
}

std::list<double> Simulator::get_infect_seconds(const int& node_index, const Network& contact_network, const Parameter& para, const int& day, const int& part) const {
	std::list<double> res;
	std::string nodename = contact_network.get_node_name(node_index);
	NeighborList nl = contact_network.get_neighbor_list(nodename, day, part);
	for (NeighborList::const_iterator it = nl.begin(); it != nl.end(); ++it) {
		std::string neighbor_name = it->get_name();
		double weight = it->get_weight();
		if (UtilConstant::isInfectious(para.get_disease(), _stages[contact_network.get_node_index(neighbor_name)])) {
			res.push_back(para.get_seconds_per_weight() * weight);
		}
	}
	return res;
}
