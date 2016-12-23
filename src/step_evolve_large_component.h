#ifndef NETWORKPROJECT_STEP_EVOLVE_LARGE_COMPONENT_H_
#define NETWORKPROJECT_STEP_EVOLVE_LARGE_COMPONENT_H_

#include <unistd.h>
#include <getopt.h>
#include <map>
#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <tuple>

//#include <boost/thread/mutex.hpp>

#include "util_gsl.h"
#include "runner.h"

class RunnerManager;
class Parameter;
class PearlNetwork;
class UndirectedGraph;
class DynamicNetwork;
//class StepSimulateEvolve;

typedef std::vector<std::string> NodeVec;
typedef std::shared_ptr<UndirectedGraph> UGraphPtr;
typedef std::shared_ptr<DynamicNetwork> DyNetPtr;
typedef std::vector<DyNetPtr> DyNetPtrVec;
typedef std::vector<DyNetPtrVec> DyNetPtrVec2;

typedef double (*CalFd)(const int& degree, const std::vector<double>& para);
//typedef double (*CalProb) (CalFd, const UGraphPtr&, const std::string&, const std::string&, const int&, const std::vector<double>&);
typedef double (*CalProb) (const std::string&, const UGraphPtr&, const std::string&, const std::string&, const int&, const std::vector<double>&);

class StepEvolveLargeComponent: public Runner {
private:
	static StepEvolveLargeComponent _step_evolve_large_component;

	struct TargetFunctionPara {
		DyNetPtrVec _dn;
		UGraphPtr _fg;
//		CalFd fd_func;
		std::string _type;
		
		TargetFunctionPara
			(const DyNetPtrVec& o_dn, const UGraphPtr& o_fg, const std::string& o_type) :
			_dn(o_dn), _fg(o_fg), _type(o_type) {}
	};
//	TargetFunctionPara* target_function_para;

	StepEvolveLargeComponent();
//	static boost::mutex _mtx;

public:
	virtual void help();
	virtual int run(const Parameter& para);

	static double test(const std::vector<double>& x);
//	static double function_for_nelder_mead(const std::vector<double>& x, void* para);

//	friend StepSimulateEvolve::run(const Parameter&);

public:
	static DyNetPtrVec get_dynamic_networks(const std::vector<std::string>& paths, const std::string& json_file);
	static DyNetPtrVec get_evolved(const DyNetPtrVec& dynamic_networks);
	static double calc_conn_prob(const std::string& type, const UGraphPtr& fg, const std::string& new_node, const std::string& old_node, const int& old_node_degree, const std::vector<double>& x);

private:
	static CalFd get_fd_func(const std::string&);
	static void** get_gsl_nelder_mead_para(UtilGsl::FuncBeforePtr func_ptr, void* target_function_para);
	static void* get_target_function_para(const DyNetPtrVec& dynamic_networks, const UGraphPtr& friendship_graph, const std::string& type);
	static void print_solution(const std::vector<double>& x, const double& y, const int& trial);

	static double target_function_for_gsl(const std::vector<double>& x, void* para);
	static double target_function_for_all_nets(const std::vector<double>& x, const DyNetPtrVec& dynamic_networks, const UGraphPtr& friendship_graph, const std::string& type);
	static double target_function_for_one_nets(const std::vector<double>& x, DynamicNetwork& dn, const UGraphPtr& friendship_graph, const std::string& type, double &re);

	static std::vector<double> get_fd_x(const std::vector<double>& x);
//	static double calc_conn_prob(const bool& is_friend, const int& degree, const double& alpha, CalFd calc_fd_func, const std::vector<double>& fd_x);
	static bool is_friend(const UGraphPtr& fg, const std::string& node1, const std::string& node2);
	static bool do_with_prob(const double& prob);

	static double calc_fd_order1(const int& degree, const std::vector<double>& x);
	static double calc_fd_order2(const int& degree, const std::vector<double>& x);
	static double calc_fd_fraction(const int& degree, const std::vector<double>& x);
	static double calc_fd_fraction2(const int& degree, const std::vector<double>& x);
	static double calc_fd_fraction3(const int& degree, const std::vector<double>& x);
	static double calc_fd_exponential1(const int& degree, const std::vector<double>& x);
	static double calc_fd_exponential2(const int& degree, const std::vector<double>& x);
	static double calc_fd_default(const int& degree, const std::vector<double>& x);

};


#endif // NETWORKPROJECT_STEP_EVOLVE_LARGE_COMPONENT_H_

