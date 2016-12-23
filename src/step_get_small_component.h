#ifndef NETWORKPROJECT_STEP_GET_SMALL_COMPONENT_H_
#define NETWORKPROJECT_STEP_GET_SMALL_COMPONENT_H_

#include <unistd.h>
#include <getopt.h>
#include <map>
#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <tuple>

#include "runner.h"

class RunnerManager;
class Parameter;
class PearlNetwork;
class UndirectedGraph;

/*
struct SourceIdentificationRes {
	std::vector<double> running_times;
	std::vector<double> error_distances;
	double running_time_mean;
	double running_time_sigma;
	double error_distance_mean;
	double error_distance_sigma;
	double detection_rate;
	std::string memory_peak_usage;
};
*/

typedef std::tuple<std::string, std::string> Tuple2Str;
typedef std::tuple<std::string, std::string, std::string> Tuple3Str;

class StepGetSmallComponent: public Runner {
	typedef std::vector<std::string> NodeVec;
	typedef std::pair<int, int> StartLastPartPair;
	typedef std::vector<StartLastPartPair> PPVec;
	typedef std::shared_ptr<PPVec> PPVecPtr;
	typedef std::pair<int, PPVecPtr> PartsRes;
	typedef std::vector<PartsRes> ResVec;
	typedef std::shared_ptr<UndirectedGraph> UGraphPtr;

private:
	static StepGetSmallComponent _step_get_small_component;

	StepGetSmallComponent();

public:
	virtual void help();
	virtual int run(const Parameter& para);

private:
	std::vector<Tuple2Str> get_node_pairs(const NodeVec& node_vec);
	std::vector<Tuple3Str> get_node_triplets(const NodeVec& node_vec);

	std::vector<NodeVec> get_node_pairs(const NodeVec& node_vec, const int& size, const std::string& condition);
	std::vector<NodeVec> get_node_lines(const NodeVec& node_vec, const int& size);
	std::vector<NodeVec> get_node_cliques(const NodeVec& node_vec, const int& size);
	// legal conditions are: 
	//	* is_line_without_others
	//	* is_clique_without_others
	ResVec test_component_parts(const NodeVec& t, const int& start_part, const int& end_part, const int& th, const PearlNetwork& pn, const std::string& condition); 
	bool check_node_vec(const UGraphPtr& g, const NodeVec& t, const std::string& condition);

	std::string get_out_filename(const Tuple2Str& t);
	std::string get_out_filename(const Tuple3Str& t);
	std::string get_out_filename(const NodeVec& t);

	void write_result(const std::string& res_file, const ResVec& res_vec);
};


#endif // NETWORKPROJECT_STEP_GET_SMALL_COMPONENT_H_
