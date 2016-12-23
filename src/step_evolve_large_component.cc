#include "step_evolve_large_component.h"

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

//#include <boost/thread/thread.hpp>
//#include <boost/thread/locks.hpp>
//#include <boost/bind.hpp>

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

StepEvolveLargeComponent StepEvolveLargeComponent::_step_evolve_large_component;
//boost::mutex StepEvolveLargeComponent::_mtx;

StepEvolveLargeComponent::StepEvolveLargeComponent(): Runner() {
	_short_options = "h";
	_long_options = new struct option[100]{
		{"help",		no_argument,		NULL, OPT_HELP},
		{"net_inroot",	required_argument,	NULL, OPT_NET_INROOT},
		{"net_injson",	required_argument,	NULL, OPT_NET_INJSON},
//		{"out_dir",		required_argument,	NULL, OPT_OUT_DIR},
		{"net_friendship",	required_argument,	NULL,	OPT_NET_FRIENDSHIP},
		{"fd_func",		required_argument,	NULL,	OPT_FD_FUNC},
		{NULL,			0,					NULL,  0 } //must end with {0, 0, 0, 0}
	};
	RunnerManager::instance()->install("step_evolve_large_component", this);
}

void StepEvolveLargeComponent::help() {
	std::cout << "\nFunctionality: find two-node components and three-node components, find the start part and parts of duration if the part and the duration satisfy some constraints.\n";
	std::cout << "Option list:\n";
	std::cout << "\t* --help (or -h): [ no argument ] print this help information.\n";
	std::cout << "\t* --net_inroot: [ string argument ] root directory of the networks.\n";
	std::cout << "\t* --net_injson: [ string argument ] json file of the networks.\n";
	std::cout << "\t* --net_friendship: [ string argument ] gml file of the friendship network.\n";
	std::cout << "\t* --fd_func: [ string argument ] name of the function to calculate f(d) in the equation p = a * s + f(d).\n";
	std::cout << "\t\t\t Possible values of fd_func: order1, order2, fraction, fraction2, fraction3, exponential1, exponential2.\n";
//	std::cout << "\t* --out_dir: [ string argument ] output directory of the results.\n";
	std::cout << std::endl;
}

int StepEvolveLargeComponent::run(const Parameter& para) {
	std::string real_inroot = para.get_net_inroot();
	std::vector<std::string> paths = Util::getAllDirPaths(real_inroot, 2);
	DyNetPtrVec dynamic_networks = get_dynamic_networks(paths, para.get_net_injson());
	dynamic_networks = get_evolved(dynamic_networks);
	UGraphPtr friendship_graph = std::make_shared<UndirectedGraph>(para.get_net_friendship(), NodeSet());
//	CalFd fd_func = get_fd_func(para.get_fd_func());
	void** gsl_nelder_mead_para = get_gsl_nelder_mead_para(
			&target_function_for_gsl, 
			get_target_function_para(dynamic_networks, friendship_graph, para.get_fd_func()));

	int iter_max = 100;
	double simplex_size_max = 0.0001;
	bool detailed_iteration_wanted = false;
	std::vector<double> temp_res, best_res;
	double min_cost = (std::numeric_limits<double>::max)();
	for (int i = 0; i < 200; ++i) {
		std::vector<double> x;
		std::vector<double> step;
		x.push_back(Util::gen_rand_double());
		if (para.get_fd_func().compare("order2") == 0) {
			// For order 2
			x.push_back(Util::gen_rand_double() * (1 - x[0]));
			x.push_back(Util::gen_rand_double() * x[0] / (27 * 27));
			step = Util::gen_rand_vector_double(3, 0, 0.2);
		} else if (para.get_fd_func().compare("order1") == 0) {
			x.push_back(Util::gen_rand_double() * (1 - x[0]));
			x.push_back(Util::gen_rand_double() * x[0] / 27);
			step = Util::gen_rand_vector_double(3, 0, 0.2);
		} else if (para.get_fd_func().compare("fraction") == 0 ||
				para.get_fd_func().compare("fraction2") == 0 ||
				para.get_fd_func().compare("fraction3") == 0) {
			x.push_back(Util::gen_rand_double() * (1 - x[0]));
			step = Util::gen_rand_vector_double(2, 0, 0.2);
		} else if (para.get_fd_func().compare("exponential1") == 0) {
			x.push_back(Util::gen_rand_double() + 1);
			x.push_back(Util::gen_rand_double() * (1 - x[0]));
			step = Util::gen_rand_vector_double(3, 0, 0.2);
		} else if (para.get_fd_func().compare("exponential2") == 0) {
			x.push_back(Util::gen_rand_double() * 10);
			x.push_back(Util::gen_rand_double() * (1 - x[0]));
			step = Util::gen_rand_vector_double(3, 0, 0.2);
		} else {
			Util::checkTrue(true, "Error: unknown function f(d) named " + para.get_fd_func());
			//x = Util::gen_rand_vector_double(3, 0, 0.5);
			//step = Util::gen_rand_vector_double(2, 0, 0.2);
		}
		double temp_cost = UtilGsl::getMinWithNelderMead(x, step, iter_max, simplex_size_max, gsl_nelder_mead_para, detailed_iteration_wanted, temp_res);
		print_solution(temp_res, temp_cost, i);
		if (temp_cost < min_cost) {
			min_cost = temp_cost;
			best_res = temp_res;
			std::cout << "Best solution is updated." << std::endl;
		}
	}
	return 0;
}

double StepEvolveLargeComponent::test(const std::vector<double>& x) {
	return x[0] * x[0] + 0.1 * x[1] * x[1];
}


DyNetPtrVec StepEvolveLargeComponent::get_dynamic_networks(const std::vector<std::string>& paths, const std::string& json_file) {
	Parameter para;
	para.set_net_injson(json_file);
	DyNetPtrVec res(paths.size());
	for (int i = 0; i < paths.size(); ++i) {
		para.set_net_inroot(paths[i]);
		res[i] = std::make_shared<DynamicNetwork>(para);
//		std::cout << res[i]->get_node_size() << std::endl;
	}
	return res;
}

DyNetPtrVec StepEvolveLargeComponent::get_evolved(const DyNetPtrVec& dynamic_networks) {
	DyNetPtrVec res(dynamic_networks.size());
	for (int i = 0; i < res.size(); ++i) {
		res[i] = dynamic_networks[i]->get_evolved_networks_end_with_largest_component();
	}
	return res;
}

CalFd StepEvolveLargeComponent::get_fd_func(const std::string& str) {
	if (str.compare("order1") == 0) {
		return &calc_fd_order1;
	} else if (str.compare("order2") == 0) {
		return &calc_fd_order2;
	} else if (str.compare("fraction") == 0) {
		return &calc_fd_fraction;
	} else if (str.compare("fraction2") == 0) {
		return &calc_fd_fraction2;
	} else if (str.compare("fraction3") == 0) {
		return &calc_fd_fraction3;
	} else if (str.compare("exponential1") == 0) {
		return &calc_fd_exponential1;
	} else if (str.compare("exponential2") == 0) {
		return &calc_fd_exponential2;
	} else {
		Util::checkTrue(true, "Error: unknown function f(d) with name " + str + ".");
		return &calc_fd_default;
	}
}

void** StepEvolveLargeComponent::get_gsl_nelder_mead_para(UtilGsl::FuncBeforePtr func_ptr, void* target_function_para) {
//	void** res = new (void*)[2];
	void** res = (void**) malloc (sizeof(void*) * 2);
	res[0] = (void*) func_ptr;
	res[1] = target_function_para;
	return res;
}


void* StepEvolveLargeComponent::get_target_function_para(const DyNetPtrVec& dynamic_networks, const UGraphPtr& friendship_graph, const std::string& type) {
	TargetFunctionPara* res = new TargetFunctionPara(dynamic_networks, friendship_graph, type);
	return (void*) res;
}

void StepEvolveLargeComponent::print_solution(const std::vector<double>& x, const double& y, const int& trial) {
	std::cout << "Trial\t" << trial << ": x = ";
	for (int i = 0; i < x.size(); ++i) {
		std::cout << x[i] << ", ";
	}
	std::cout << "y = " << y << std::endl;
}

double StepEvolveLargeComponent::target_function_for_gsl(const std::vector<double>& x, void* para) {
	TargetFunctionPara* tf_para = (TargetFunctionPara*) para;
	return target_function_for_all_nets(x, tf_para->_dn, tf_para->_fg, tf_para->_type);
}

double StepEvolveLargeComponent::target_function_for_all_nets(const std::vector<double>& x, const DyNetPtrVec& dynamic_networks, const UGraphPtr& friendship_graph, const std::string& type) {
	double sum = 0, temp;
	for (int j = 0; j < 100; ++j) {
		for (int i = 0; i < dynamic_networks.size(); ++i) {
			sum += target_function_for_one_nets(x, *(dynamic_networks[i]), friendship_graph, type, temp);
		}
	}
	return sum / (100 * dynamic_networks.size());
}

double StepEvolveLargeComponent::target_function_for_one_nets(const std::vector<double>& x, DynamicNetwork& dn, const UGraphPtr& friendship_graph, const std::string& type, double &re) {
	std::shared_ptr<DynamicNetwork> simulated_evolving_networks = dn.get_simulated_evolving_networks(type, friendship_graph, calc_conn_prob, x);
	return re = dn.get_avg_frobenius_norm_of_diff(*simulated_evolving_networks);
}


std::vector<double> StepEvolveLargeComponent::get_fd_x(const std::vector<double>& x) {
	std::vector<double> res(x.size() - 1);
	for (int i = 1; i < x.size(); ++i) {
		res[i - 1] = x[i];
	}
	return res;
}

double StepEvolveLargeComponent::calc_conn_prob(const std::string& type, const UGraphPtr& fg, const std::string& new_node, const std::string& old_node, const int& old_node_degree, const std::vector<double>& x) {
	Util::checkTrue(x.size() > 0, "Error: x has only " + std::to_string(x.size()) + " variables.");
	CalFd calc_fd = get_fd_func(type);
	double alpha = x[0];
	std::vector<double> fd_x = get_fd_x(x);
	bool isf = is_friend(fg, new_node, old_node);
	return alpha * isf + calc_fd(old_node_degree, fd_x);
}

/*
double StepEvolveLargeComponent::calc_conn_prob(const bool& is_friend, const int& degree, const double& alpha, CalFd calc_fd_func, const std::vector<double>& fd_x) {
	return alpha * is_friend + calc_fd_func(degree, fd_x);
}
*/

bool StepEvolveLargeComponent::is_friend(const UGraphPtr& fg, const std::string& node1, const std::string& node2) {
	return fg->contains_node(node1) && fg->contains_node(node2) && fg->is_connected(node1, node2);
}

bool StepEvolveLargeComponent::do_with_prob(const double& prob) {
	return Util::gen_rand_double() <= prob;
}

double StepEvolveLargeComponent::calc_fd_order1(const int& degree, const std::vector<double>& x) {
	Util::checkTrue(x.size() >= 2, "Error for not enough input x (calc_fd_order1).");
	return x[0] - x[1] * degree;
}

double StepEvolveLargeComponent::calc_fd_order2(const int& degree, const std::vector<double>& x) {
	Util::checkTrue(x.size() >= 2, "Error for not enough input x (calc_fd_order1).");
	return x[0] - x[1] * degree * degree;
}

double StepEvolveLargeComponent::calc_fd_fraction(const int& degree, const std::vector<double>& x) {
	Util::checkTrue(x.size() >= 1, "Error for not enough input x (calc_fd_order1).");
	if (degree <= 0) {
		return 1e10;
	} else {
		return x[0] / (1.0 * degree);
	}
}

double StepEvolveLargeComponent::calc_fd_fraction2(const int& degree, const std::vector<double>& x) {
	Util::checkTrue(x.size() >= 1, "Error for not enough input x (calc_fd_order1).");
	if (degree <= 0) {
		return 1e10;
	} else {
		return x[0] / (1.0 * degree * degree);
	}
}

double StepEvolveLargeComponent::calc_fd_fraction3(const int& degree, const std::vector<double>& x) {
	Util::checkTrue(x.size() >= 1, "Error for not enough input x (calc_fd_order1).");
	if (degree <= 0) {
		return 1e10;
	} else {
		return x[0] / (1.0 * degree * degree * degree);
	}
}

double StepEvolveLargeComponent::calc_fd_exponential1(const int& degree, const std::vector<double>& x) {
	Util::checkTrue(x.size() >= 2, "Error for not enough input x (calc_fd_order1).");
	return x[1] * pow(x[0], -degree);
}

double StepEvolveLargeComponent::calc_fd_exponential2(const int& degree, const std::vector<double>& x) {
	Util::checkTrue(x.size() >= 2, "Error for not enough input x (calc_fd_order1).");
	return x[1] * pow(degree, -x[0]);
}


double StepEvolveLargeComponent::calc_fd_default(const int& degree, const std::vector<double>& x) {
	return 0;
}
