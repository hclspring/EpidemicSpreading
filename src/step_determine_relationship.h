#ifndef NETWORKPROJECT_STEP_DETERMINE_RELATIONSHIP_H_
#define NETWORKPROJECT_STEP_DETERMINE_RELATIONSHIP_H_

#include <unistd.h>
#include <getopt.h>
#include <map>
#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <unordered_set>
#include <utility>
#include <tuple>

#include "runner.h"

class RunnerManager;
class Parameter;
class UndirectedGraph;

class StepDetermineRelationship: public Runner {
	typedef std::vector<std::string> NodeVec;
	typedef std::unordered_set<std::string> NodeSet;
	typedef std::shared_ptr<UndirectedGraph> UGraphPtr;

private:
	static StepDetermineRelationship _step_determine_relationship;

	StepDetermineRelationship();

public:
	virtual void help();
	virtual int run(const Parameter& para);

private:
	int get_bound(const int& N, const int& K, const int& n, const double& alpha);
	int get_bound(const std::vector<double>& x, const double& alpha);
	std::vector<double> get_hypergeometric_prob(const int& N, const int& K, const int& n);
	long long choose(const int& N, const int& n);
};


#endif // NETWORKPROJECT_STEP_DETERMINE_RELATIONSHIP_H_
