#ifndef NETWORKPROJECT_SIMULATOR_H_
#define NETWORKPROJECT_SIMULATOR_H_

#include <vector>
#include <list>
#include <unordered_set>
#include <algorithm>

#include "constant.h"
#include "util_constant.h"

class Network;
class Parameter;
class NeighborInfo;

typedef std::unordered_set<std::string> NodeSet;
typedef std::unordered_set<int> IndexSet;
typedef std::list<NeighborInfo> NeighborList;

class Simulator {
private:
	std::vector<DiseaseStage> _initial_stages;
	std::vector<DiseaseStage> _stages;
	std::vector<DiseaseStage> _unstable_stages;

	int _network_days;
	int _network_parts;
public:
	Simulator(const Parameter& para, const int& population, const IndexSet& source_set);
	~Simulator();

	std::vector<DiseaseStage> get_sim_res(const Network& contact_network, const Parameter& para);
	static NodeSet get_nodeset_been_infected_from_sim_res(const std::vector<DiseaseStage>& sim_res, const DiseaseModel& disease, const Network& network);
private:
//	void initialize_stages(const DiseaseModel& disease, const int& population, const int& source_count);
	void initialize_stages(const DiseaseModel& disease, const int& population, const IndexSet& source_set);
	bool is_unstable(const DiseaseModel& disease) const;

	std::vector<DiseaseStage> spread(const std::vector<DiseaseStage>& previous_stages, const Network& contact_network, const Parameter& para, const int& day, const int& part);
	std::list<double> get_infect_seconds(const int& node_index, const Network& contact_network, const Parameter& para, const int& day, const int& part) const;
};


#endif // NETWORKPROJECT_SIMULATOR_H_

