#ifndef NETWORKPROJECT_STEP_SIMULATE_EVOLVE_H_
#define NETWORKPROJECT_STEP_SIMULATE_EVOLVE_H_

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

typedef std::vector<std::string> NodeVec;
typedef std::shared_ptr<UndirectedGraph> UGraphPtr;
typedef std::shared_ptr<DynamicNetwork> DyNetPtr;
typedef std::vector<DyNetPtr> DyNetPtrVec;
typedef std::vector<DyNetPtrVec> DyNetPtrVec2;

typedef double (*CalFd)(const int& degree, const std::vector<double>& para);
//typedef double (*CalProb) (CalFd, const UGraphPtr&, const std::string&, const std::string&, const int&, const std::vector<double>&);
typedef double (*CalProb) (const std::string&, const UGraphPtr&, const std::string&, const std::string&, const int&, const std::vector<double>&);

class StepSimulateEvolve: public Runner {
private:
	static StepSimulateEvolve _step_simulate_evolve;
/*
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
*/
	StepSimulateEvolve();
//	static boost::mutex _mtx;

public:
	virtual void help();
	virtual int run(const Parameter& para);

//	static double test(const std::vector<double>& x);
//	static double function_for_nelder_mead(const std::vector<double>& x, void* para);
private:
	static std::vector<double> get_sim_para(const Parameter& para);
	static DyNetPtrVec get_simulate_res(const DyNetPtrVec& dn_vec, const std::string& type, const UGraphPtr& fg, CalProb calc_prob, const std::vector<double>& x);
	static void calc_diff_features(const DyNetPtrVec& ori_dns, const DyNetPtrVec& sim_dns, std::vector<double>& diff_ccs, std::vector<double>& diff_avg_dists, std::vector<double>& diff_densities, std::vector<double>& diff_fnorms);
	static void calc_self_features(const DyNetPtrVec& dns, std::vector<std::vector<double>>& all_ccs, std::vector<std::vector<double>>& all_closenesses, std::vector<std::vector<double>>& all_degrees);
	static void write_dns(const std::string& out_dir, const DyNetPtrVec& dns, const std::vector<std::string>& subpaths);
	static void write_features(const std::string& out_dir, const std::vector<double>& diff_ccs, const std::vector<double>& diff_avg_dists, const std::vector<double>& diff_densities, const std::vector<double>& diff_fnorms);
	static void write_self_features(const std::string& out_dir, const std::vector<std::vector<double>>& all_ccs, const std::vector<std::vector<double>>& all_closenesses, const std::vector<std::vector<double>>& all_degrees, const std::vector<std::string>& subpaths);
};


#endif // NETWORKPROJECT_STEP_SIMULATE_EVOLVE_H_

