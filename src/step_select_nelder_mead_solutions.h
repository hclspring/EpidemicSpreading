#ifndef NETWORKPROJECT_STEP_SELECT_NELDER_MEAD_SOLUTIONS_H_
#define NETWORKPROJECT_STEP_SELECT_NELDER_MEAD_SOLUTIONS_H_

#include <unistd.h>
#include <getopt.h>
#include <map>
#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <tuple>

#include "util_gsl.h"
#include "runner.h"

class RunnerManager;
class Parameter;

class StepSelectNelderMeadSolutions: public Runner {
private:
	static StepSelectNelderMeadSolutions _step_select_nelder_mead_solutions;

	StepSelectNelderMeadSolutions();

public:
	virtual void help();
	virtual int run(const Parameter& para);

private:
	static bool is_satisfied_solution(const std::string& line, const std::string& func, double& value);
	static double get_x(std::string& str);
	static double satisfy_constraint(const std::vector<double>& vecx, const std::string& func);

};


#endif // NETWORKPROJECT_STEP_SELECT_NELDER_MEAD_SOLUTIONS_H_

