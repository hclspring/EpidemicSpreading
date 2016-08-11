#include "step_source_identification.h"

#include <algorithm>
#include <vector>
#include <string>
#include <memory>
#include "unistd.h"

#include "runner.h"
#include "runner_manager.h"
#include "parameter.h"
#include "undirected_graph.h"
#include "network.h"
#include "pearl_network.h"
#include "static_network.h"
#include "dynamic_network.h"
#include "simulator.h"
#include "source_identification.h"
#include "util.h"

StepSourceIdentification StepSourceIdentification::_step_source_identification;

StepSourceIdentification::StepSourceIdentification(): Runner() {
	_short_options = "";
	_long_options = new struct option[100]{
//		{"net_type",	required_argument,	NULL, OPT_NET_TYPE},
		{"net_inroot",	required_argument,	NULL, OPT_NET_INROOT},
//		{"net_injson",	required_argument,	NULL, OPT_NET_INJSON},
//		{"net_volunteers",required_argument,NULL, OPT_NET_VOLUNTEERS},
		{"out_dir",		required_argument,	NULL, OPT_OUT_DIR},
		{"disease",		required_argument,	NULL, OPT_DISEASE},
		{"infect_rate",	required_argument,	NULL, OPT_INFECT_RATE},
		{"infect_rate_seconds",	required_argument,	NULL, OPT_INFECT_RATE_SECONDS},
//		{"infectious_rate",	required_argument,	NULL, OPT_INFECTIOUS_RATE},
//		{"infectious_rate_seconds",	required_argument,	NULL, OPT_INFECTIOUS_RATE_SECONDS},
		{"recover_rate",	required_argument,	NULL, OPT_RECOVER_RATE},
		{"recover_rate_seconds",	required_argument,	NULL, OPT_RECOVER_RATE_SECONDS},
		{"seconds_per_weight",	required_argument, NULL, OPT_SECONDS_PER_WEIGHT},
		{"seconds_per_step",	required_argument, NULL, OPT_SECONDS_PER_STEP},
//		{"source_count",	required_argument, NULL, OPT_SOURCE_COUNT},
//		{"snapshot_coverage",	required_argument, NULL, OPT_SNAPSHOT_COVERAGE},
		{"max_sim_days",	required_argument, NULL, OPT_MAX_SIM_DAYS},
		{"repeat_times",	required_argument, NULL, OPT_REPEAT_TIMES},
		{"source_identification_method", required_argument, NULL, OPT_SRC_IDN_METHOD},
		{"source_identification_knowntime", required_argument, NULL, OPT_SRC_IDN_KNOWNTIME},
		{NULL,			0,					NULL,  0 } //must end with {0, 0, 0, 0}
	};
	RunnerManager::instance()->install("step_source_identification", this);
}

int StepSourceIdentification::run(const Parameter& para) {
	// Read network
	std::shared_ptr<Network> net = std::make_shared<StaticNetwork>(para);
	int n= net->get_node_size();
	// Define clocks
	clock_t start_clock, end_clock;
	// Repeat simulation and source identification
	std::vector<double> detection_rates;
	SourceIdentificationRes res_struct;
	for (int i = 0; i < para.get_repeat_times(); ++i) {
		// Init seed
		int true_seed_index = Util::gen_rand_int(n);
		std::string true_seed_node = net->get_node_name(true_seed_index);
		IndexSet seed_set{true_seed_index};
		// Init simulator
		Simulator simulator(para, n, seed_set);
		// Simulate the disease spreading
		std::vector<DiseaseStage> sim_res = simulator.get_sim_res(*net, para);
		// Infer seed and record time usage
		start_clock = clock();
		std::string inferred_seed_node = SourceIdentification::calc_source(*net, sim_res, para, para.get_source_identification_method(), para.get_source_identification_knowntime());
		end_clock = clock();
		// Calculate error distance, detection rate, and running time
		double error_distance = SourceIdentification::calc_error_distance(*net, true_seed_node, inferred_seed_node);
		double detection_rate = SourceIdentification::calc_detection_rate(true_seed_node, inferred_seed_node);
		double running_time = 1000.0 * (end_clock - start_clock) / CLOCKS_PER_SEC;
		res_struct.error_distances.push_back(error_distance);
		detection_rates.push_back(detection_rate);
		res_struct.running_times.push_back(running_time);
	}
	// Calculate error distance, detection rate, running time, and memory peak usage.
	res_struct.error_distance_mean = Util::getMean(res_struct.error_distances);
	res_struct.error_distance_sigma = Util::getDeviation(res_struct.error_distances, res_struct.error_distance_mean);
	res_struct.detection_rate = Util::getMean(detection_rates);
	res_struct.running_time_mean = Util::getMean(res_struct.running_times);
	res_struct.running_time_sigma = Util::getDeviation(res_struct.running_times, res_struct.running_time_mean);
	res_struct.memory_peak_usage = Util::getMemoryPeakUsage();
	// Write results to files.
	write_result(para.get_out_dir(), res_struct);
	return 0;
}

void StepSourceIdentification::write_result(const std::string& res_dir, const SourceIdentificationRes& res_struct) {
	std::string details = res_dir + "/details.txt", summary = res_dir + "/summary.txt";
	// write details
	std::ofstream ofs(details.c_str());
	if (ofs.fail()) {
		std::cerr << "Error writing file: " << details << std::endl;
		exit(-1);
	}
	ofs << "#index error_distance running_time(ms)" << std::endl;
	for (int i = 0; i < res_struct.error_distances.size(); ++i) {
		ofs << i + 1 << " " << (res_struct.error_distances)[i] << " " << (res_struct.running_times)[i] << std::endl;
	}
	ofs.close();
	ofs.open(summary.c_str());
	if (ofs.fail()) {
		std::cerr << "Error writing file: " << summary << std::endl;
		exit(-1);
	}
	ofs << "error_distance_mean " << res_struct.error_distance_mean << std::endl;
	ofs << "error_distance_sigma " << res_struct.error_distance_sigma << std::endl;
	ofs << "detection_rate " << res_struct.detection_rate << std::endl;
	ofs << "running_time_mean " << res_struct.running_time_mean << std::endl;
	ofs << "running_time_sigma " << res_struct.running_time_sigma << std::endl;
	ofs << "memory_peak_usage " << res_struct.memory_peak_usage << std::endl;
	ofs.close();
}
