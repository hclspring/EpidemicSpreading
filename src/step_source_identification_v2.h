#ifndef NETWORKPROJECT_STEP_SOURCE_IDENTIFICATION_V2_H_
#define NETWORKPROJECT_STEP_SOURCE_IDENTIFICATION_V2_H_

#include <unistd.h>
#include <getopt.h>
#include <map>
#include <iostream>
#include <string>

#include "runner.h"

class RunnerManager;
class Parameter;

struct SourceIdentificationResV2 {
	std::string method_name;
	std::vector<double> running_times;
	std::vector<double> detection_rates;
	std::vector<double> error_distances;
	double running_time_mean;
	double running_time_sigma;
	double error_distance_mean;
	double error_distance_sigma;
	double detection_rate;
};

class StepSourceIdentificationV2: public Runner {
private:
	static StepSourceIdentificationV2 _step_source_identification_v2;

	StepSourceIdentificationV2();

public:
	virtual void help();
	virtual int run(const Parameter& para);

private:
	void write_result(const std::string& res_dir, const std::vector<SourceIdentificationResV2>& res_structs);
};


#endif // NETWORKPROJECT_STEP_SOURCE_IDENTIFICATION_V2_H_
