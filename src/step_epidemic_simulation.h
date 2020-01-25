#ifndef NETWORKPROJECT_STEP_EPIDEMIC_SIMULATION_H_
#define NETWORKPROJECT_STEP_EPIDEMIC_SIMULATION_H_

#include <unistd.h>
#include <getopt.h>
#include <map>
#include <iostream>
#include <string>

#include "runner.h"

class RunnerManager;
class Parameter;

struct EpidemicSimulationRes { //TODO: Design a reasonable simulation result structure
	std::vector<double> running_times; // how much time used in each simulation
	double running_time_mean;
	double running_time_sigma;
};

class StepEpidemicSimulation: public Runner {
private:
	static StepEpidemicSimulation _step_epidemic_simulation;

	StepEpidemicSimulation();

public:
	virtual void help();
	virtual int run(const Parameter& para);

private:
	void write_result(const std::string& res_dir, const EpidemicSimulationRes& res_struct);
};


#endif // NETWORKPROJECT_STEP_EPIDEMIC_SIMULATION_H_
