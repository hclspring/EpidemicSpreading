#include "step_select_nelder_mead_solutions.h"

#include <cstdlib>
#include <limits>
#include <algorithm>
#include <vector>
#include <string>
#include <memory>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <set>
#include <tuple>
#include <functional>
#include <unistd.h>

#include "runner.h"
#include "runner_manager.h"
#include "parameter.h"
#include "undirected_graph.h"
#include "network.h"
#include "pearl_network.h"
#include "static_network.h"
#include "dynamic_network.h"
//#include "simulator.h"
//#include "source_identification.h"
#include "util.h"
#include "util_boost.h"
#include "util_gsl.h"

StepSelectNelderMeadSolutions StepSelectNelderMeadSolutions::_step_select_nelder_mead_solutions;

StepSelectNelderMeadSolutions::StepSelectNelderMeadSolutions(): Runner() {
	_short_options = "h";
	_long_options = new struct option[100]{
		{"help",		no_argument,		NULL, OPT_HELP},
		{"net_inroot",	required_argument,	NULL, OPT_NET_INROOT},
//		{"net_injson",	required_argument,	NULL, OPT_NET_INJSON},
		{"out_dir",		required_argument,	NULL, OPT_OUT_DIR},
		{"fd_func",		required_argument,	NULL,	OPT_FD_FUNC},
		{NULL,			0,					NULL,  0 } //must end with {0, 0, 0, 0}
	};
	RunnerManager::instance()->install("step_select_nelder_mead_solutions", this);
}

void StepSelectNelderMeadSolutions::help() {
	std::cout << "\nFunctionality: find the solutions which satisfy constraints.\n";
	std::cout << "Option list:\n";
	std::cout << "\t* --help (or -h): [ no argument ] print this help information.\n";
	std::cout << "\t* --net_inroot: [ string argument ] input file of the previous solutions.\n";
//	std::cout << "\t* --net_injson: [ string argument ] json file of the networks.\n";
	std::cout << "\t* --fd_func: [ string argument ] name of the function to calculate f(d) in the equation p = a * s + f(d).\n";
	std::cout << "\t\t\t Possible values of fd_func: order1, order2, fraction, fraction2, fraction3, exponential1 exponential2.\n";
	std::cout << "\t* --out_dir: [ string argument ] output file of the screened results.\n";
	std::cout << std::endl;
}

int StepSelectNelderMeadSolutions::run(const Parameter& para) {
	std::vector<std::string> lines = Util::readLines(para.get_net_inroot());
	std::vector<std::string> out_lines;
	double best_value = (std::numeric_limits<double>::max)();
	double value;
	int best_index = -1;
	for (int i = 0; i < lines.size(); ++i) {
		if (is_satisfied_solution(lines[i], para.get_fd_func(), value)) {
			if (value < best_value) {
				best_value = value;
				best_index = out_lines.size();
			}
			out_lines.push_back(lines[i]);
		}
	}
	if (best_index >= 0) {
		out_lines.push_back("Best solution:");
		out_lines.push_back(out_lines[best_index]);
	}
	Util::writeLines(out_lines, para.get_out_dir());
	return 0;
}


bool StepSelectNelderMeadSolutions::is_satisfied_solution(const std::string& line, const std::string& func, double& value) {
	std::stringstream sstr(line);
	std::string str;
	double x;
	std::vector<double> vecx;
	double y;
	sstr >> str;
	if (str.compare("Trial") == 0) {
		sstr >> str >> str >> str >> str;
		vecx.push_back(get_x(str));
		sstr >> str;
		vecx.push_back(get_x(str));
		if (func.find_first_of("fraction") != 0) {
			sstr >> str;
			vecx.push_back(get_x(str));
		}
		sstr >> str >> str >> y;
		value = y;
		return satisfy_constraint(vecx, func);
	} else {
		return false;
	}
}

double StepSelectNelderMeadSolutions::get_x(std::string& str) {
	str[str.length() - 1] = '\0';
	return std::stod(str);
}

double StepSelectNelderMeadSolutions::satisfy_constraint(const std::vector<double>& vecx, const std::string& func) {
	if (vecx[0] < 0 || vecx[0] > 1) {
		return false;
	} else if (func.compare("order2") == 0) {
		if (vecx[1] < 0 || vecx[1] > 1 - vecx[0]) {
			return false;
		} else if (vecx[2] < 0 || vecx[2] > vecx[0] / (27 * 27)) {
			return false;
		} else {
			return true;
		}
	} else if (func.compare("order1") == 0) {
		if (vecx[1] < 0 || vecx[1] > 1 - vecx[0]) {
			return false;
		} else if (vecx[2] < 0 || vecx[2] > vecx[0] / 27) {
			return false;
		} else {
			return true;
		}
	} else if (func.find_first_of("fraction") == 0) {
		if (vecx[1] < 0 || vecx[1] > 1 - vecx[0]) {
			return false;
		} else {
			return true;
		}
	} else if (func.find_first_of("exponential") == 0) {
		if ((func.compare("exponential2") == 0 && (vecx[1] < 0 || vecx[1] > 10))  ||
				(func.compare("exponential1") == 0 && (vecx[1] < 1 || vecx[1] > 2))) {
			return false;
		} else if (vecx[2] < 0 || vecx[2] > 1 - vecx[0]) {
			return false;
		} else {
			return true;
		}
	} else {
		return true;
	}
}


