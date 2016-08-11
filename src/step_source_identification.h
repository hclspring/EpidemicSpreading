#ifndef NETWORKPROJECT_STEP_SOURCE_IDENTIFICATION_H_
#define NETWORKPROJECT_STEP_SOURCE_IDENTIFICATION_H_

#include <unistd.h>
#include <getopt.h>
#include <map>
#include <iostream>
#include <string>

#include "runner.h"

class RunnerManager;
class Parameter;

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

class StepSourceIdentification: public Runner {
private:
	static StepSourceIdentification _step_source_identification;

	StepSourceIdentification();

public:
	virtual int run(const Parameter& para);

private:
	void write_result(const std::string& res_dir, const SourceIdentificationRes& res_struct);
};


#endif // NETWORKPROJECT_STEP_SOURCE_IDENTIFICATION_H_
