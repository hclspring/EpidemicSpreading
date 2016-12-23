#ifndef NETWORKPROJECT_STEP_CALC_STATIC_FEATURE_H_
#define NETWORKPROJECT_STEP_CALC_STATIC_FEATURE_H_

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
class UndirectedGraph;

typedef std::shared_ptr<UndirectedGraph> UGraphPtr;

class StepCalcStaticFeature: public Runner {
private:
	static StepCalcStaticFeature _step_calc_static_feature;

	StepCalcStaticFeature();

public:
	virtual void help();
	virtual int run(const Parameter& para);

private:
	std::vector<std::string> get_calc_res(const UGraphPtr& g, const Parameter& para);
};


#endif // NETWORKPROJECT_STEP_CALC_STATIC_FEATURE_H_
