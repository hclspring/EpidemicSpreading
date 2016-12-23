#ifndef NETWORKPROJECT_STEP_REINDEX_GML_H_
#define NETWORKPROJECT_STEP_REINDEX_GML_H_

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

class StepReindexGml: public Runner {

private:
	static StepReindexGml _step_reindex_gml;

	StepReindexGml();

public:
	virtual void help();
	virtual int run(const Parameter& para);

private:
	void reindex_gml(std::vector<std::string>& lines);

	int get_first_nonblank_pos(const std::string& str, const int& start);
	bool starts_with(const std::string& str, const int& index, const std::string& str2);
	int get_first_integer(const std::string& str, int& start);
};


#endif // NETWORKPROJECT_STEP_REINDEX_GML_H_
