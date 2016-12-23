#include "step_simulate_evolve.h"

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

#include "step_evolve_large_component.h"


StepSimulateEvolve StepSimulateEvolve::_step_simulate_evolve;
//boost::mutex StepSimulateEvolve::_mtx;

StepSimulateEvolve::StepSimulateEvolve(): Runner() {
	_short_options = "h";
	_long_options = new struct option[100]{
		{"help",		no_argument,		NULL, OPT_HELP},
		{"net_inroot",	required_argument,	NULL, OPT_NET_INROOT},
		{"net_injson",	required_argument,	NULL, OPT_NET_INJSON},
		{"out_dir",		required_argument,	NULL, OPT_OUT_DIR},
		{"net_friendship",	required_argument,	NULL,	OPT_NET_FRIENDSHIP},
		{"fd_func",		required_argument,	NULL,	OPT_FD_FUNC},
		{"evolve_para_alpha",	required_argument,	NULL,	OPT_EVOLVE_PARA_ALPHA},
		{"evolve_para_a",	required_argument,	NULL,	OPT_EVOLVE_PARA_A},
		{"evolve_para_b",	required_argument,	NULL,	OPT_EVOLVE_PARA_B},
		{NULL,			0,					NULL,  0 } //must end with {0, 0, 0, 0}
	};
	RunnerManager::instance()->install("step_simulate_evolve", this);
}

void StepSimulateEvolve::help() {
	std::cout << "\nFunctionality: simulate the evolving process, and calculate features of the networks.\n";
	std::cout << "Option list:\n";
	std::cout << "\t* --help (or -h): [ no argument ] print this help information.\n";
	std::cout << "\t* --net_inroot: [ string argument ] root directory of the networks.\n";
	std::cout << "\t* --net_injson: [ string argument ] json file of the networks.\n";
	std::cout << "\t* --out_dir: [ string argument ] root directory of the output.\n";
	std::cout << "\t* --net_friendship: [ string argument ] gml file of the friendship network.\n";
	std::cout << "\t* --fd_func: [ string argument ] name of the function to calculate f(d) in the equation p = a * s + f(d).\n";
	std::cout << "\t* --evolve_para_alpha: [ double argument ] parameter alpha of the probability.\n";
	std::cout << "\t* --evolve_para_a: [ double argument ] parameter a of the probability.\n";
	std::cout << "\t* --evolve_para_b: [ double argument ] parameter b of the probability.\n";
	std::cout << "\t\t\t Possible values of fd_func: order1, order2, fraction, fraction2, fraction3, exponential1, exponential2.\n";
//	std::cout << "\t* --out_dir: [ string argument ] output directory of the results.\n";
	std::cout << std::endl;
}

int StepSimulateEvolve::run(const Parameter& para) {
	std::string real_inroot = para.get_net_inroot();
	std::vector<std::string> paths = Util::getAllDirPaths(real_inroot, 2);
	std::vector<std::string> subpaths = Util::getAllSubDirPaths(real_inroot, 2);

	DyNetPtrVec ori_dns = StepEvolveLargeComponent::get_dynamic_networks(paths, para.get_net_injson());
	/*
	for (int i = 0; i < ori_dns.size(); ++i) {
		ori_dns[i]->write_gml(para.get_out_dir() + "/ori_dns_old/" + std::to_string(i) + "/");
	}
	*/
	ori_dns = StepEvolveLargeComponent::get_evolved(ori_dns);
	std::vector<std::vector<double>> ori_ccs, ori_closenesses, ori_degrees;
	calc_self_features(ori_dns, ori_ccs, ori_closenesses, ori_degrees);
	std::string output_dir = para.get_out_dir() + "/original/";
	int a = std::system(std::string("mkdir -p " + output_dir).c_str());
	write_self_features(output_dir, ori_ccs, ori_closenesses, ori_degrees, subpaths);
	/*
	for (int i = 0; i < ori_dns.size(); ++i) {
		ori_dns[i]->write_gml(para.get_out_dir() + "/ori_dns_new/" + std::to_string(i) + "/");
	}
	*/
	UGraphPtr fg = std::make_shared<UndirectedGraph>(para.get_net_friendship(), NodeSet());
	std::string type = para.get_fd_func();
	std::vector<double> sim_para = get_sim_para(para);
	for (int i = 0; i < 100; ++i) {
		DyNetPtrVec sim_dns = get_simulate_res(ori_dns, type, fg, StepEvolveLargeComponent::calc_conn_prob, sim_para);
		std::vector<double> diff_ccs, diff_avg_dists, diff_densities, diff_fnorms;
		calc_diff_features(ori_dns, sim_dns, diff_ccs, diff_avg_dists, diff_densities, diff_fnorms);
		std::vector<std::vector<double>> all_ccs, all_closenesses, all_degrees;
		calc_self_features(sim_dns, all_ccs, all_closenesses, all_degrees);
		std::string output_dir = para.get_out_dir() + "/trial_" + std::to_string(i) + "/";
		int a = std::system(std::string("mkdir -p " + output_dir).c_str());
		write_dns(output_dir, sim_dns, subpaths);
		write_features(output_dir, diff_ccs, diff_avg_dists, diff_densities, diff_fnorms);
		write_self_features(output_dir, all_ccs, all_closenesses, all_degrees, subpaths);
	}
	return 0;
}

std::vector<double> StepSimulateEvolve::get_sim_para(const Parameter& para) {
	std::vector<double> res {para.get_evolve_para_alpha(), para.get_evolve_para_a()};
	if (para.get_fd_func().find_first_of("fraction") != 0) {
		res.push_back(para.get_evolve_para_b());
	}
	return res;
}


DyNetPtrVec StepSimulateEvolve::get_simulate_res(const DyNetPtrVec& dn_vec, const std::string& type, const UGraphPtr& fg, CalProb calc_prob, const std::vector<double>& x) {
	int n = dn_vec.size();
	DyNetPtrVec res(n);
	for (int i = 0; i < n; ++i) {
		DyNetPtr ori_dnp = dn_vec[i];
		res[i] = ori_dnp->get_simulated_evolving_networks(type, fg, calc_prob, x);
	}
	return res;
}

void StepSimulateEvolve::calc_diff_features(const DyNetPtrVec& ori_dns, const DyNetPtrVec& sim_dns, std::vector<double>& diff_ccs, std::vector<double>& diff_avg_dists, std::vector<double>& diff_densities, std::vector<double>& diff_fnorms) {
	int m = ori_dns.size();
	diff_ccs.resize(m);
	diff_avg_dists.resize(m);
	diff_densities.resize(m);
	diff_fnorms.resize(m);

	for (int i = 0; i < m; ++i) {
		DyNetPtr ori_dnp = ori_dns[i];
		DyNetPtr sim_dnp = sim_dns[i];
		int n = ori_dnp->get_size();
		UGraphPtr ori_g = ori_dnp->at(n - 1);
		UGraphPtr sim_g = sim_dnp->at(n - 1);

		double ori_cc = ori_g->get_clustering_coefficient();
		double ori_avg_dist = ori_g->get_average_distance();
		double ori_density = ori_g->get_edge_density();

		double sim_cc = sim_g->get_clustering_coefficient();
		double sim_avg_dist = sim_g->get_average_distance();
		double sim_density = sim_g->get_edge_density();

		diff_ccs[i] = fabs(ori_cc - sim_cc);
		diff_avg_dists[i] = fabs(ori_avg_dist - sim_avg_dist);
		diff_densities[i] = fabs(ori_density - sim_density);

		double diff_fnorm = ori_dnp->get_avg_frobenius_norm_of_diff(*sim_dnp);
		/*
		double diff_fnorm = Util::getFrobeniusNorm(Util::getDiffMatrix(
			ori_g->get_unweighted_adjacency_matrix(), 
			sim_g->get_unweighted_adjacency_matrix()))
		/ ori_g->get_node_size();
		*/
		diff_fnorms[i] = diff_fnorm;
	}
	return;
}

void StepSimulateEvolve::calc_self_features(const DyNetPtrVec& dns, std::vector<std::vector<double>>& all_ccs, std::vector<std::vector<double>>& all_closenesses, std::vector<std::vector<double>>& all_degrees) {
	int m = dns.size();
	all_ccs.resize(m);
	all_closenesses.resize(m);
	all_degrees.resize(m);

	for (int i = 0; i < m; ++i) {
		DyNetPtr dnp = dns[i];
		int n = dnp->get_size();
		UGraphPtr g = dnp->at(n - 1);

		all_ccs[i] = g->get_all_clustering_coefficient();
		all_closenesses[i] = g->get_all_closeness();
		all_degrees[i] = g->get_all_degree_double();
	}
}

void StepSimulateEvolve::write_dns(const std::string& out_root_dir, const DyNetPtrVec& dns, const std::vector<std::string>& subpaths) {
	Util::checkTrue(dns.size() == subpaths.size(), "Error: Number of DyNets and paths are different.");
	int n = dns.size();
	for (int i = 0; i < n; ++i) {
		std::string output_dir = out_root_dir + "/" + subpaths[i];
		int _temp = system(std::string("mkdir -p " + output_dir).c_str());
		int m = dns[i]->get_size();
		for (int j = 0; j < m; ++j) {
			std::string filename = output_dir + "/index_" + std::to_string(j) + ".gml";
			dns[i]->at(j)->write_graph_gml(filename);
		}
	}
}

void StepSimulateEvolve::write_features(const std::string& out_dir, const std::vector<double>& diff_ccs, const std::vector<double>& diff_avg_dists, const std::vector<double>& diff_densities, const std::vector<double>& diff_fnorms) {
	std::string filename = out_dir + "/diff_features.txt";
	int n = diff_ccs.size();
	Util::checkTrue(diff_avg_dists.size() == n, "Error: sizes of diff_ccs and diff_avg_dists are different.");
	Util::checkTrue(diff_densities.size() == n, "Error: sizes of diff_ccs and diff_densities are different.");
	Util::checkTrue(diff_fnorms.size() == n, "Error: sizes of diff_ccs and diff_fnorms are different.");
	std::vector<std::string> lines (n + 1);
	lines[0] = "# CC avg_dist density fnorm";
	for (int i = 0; i < n; ++i) {
		lines[i + 1] = std::to_string(diff_ccs[i]) + " " + std::to_string(diff_avg_dists[i]) + " "
			+ std::to_string(diff_densities[i]) + " " + std::to_string(diff_fnorms[i]);
	}
	Util::writeLines(lines, filename);
}

void StepSimulateEvolve::write_self_features(const std::string& out_root_dir, const std::vector<std::vector<double>>& all_ccs, const std::vector<std::vector<double>>& all_closenesses, const std::vector<std::vector<double>>& all_degrees, const std::vector<std::string>& subpaths) {
	int m = subpaths.size();
	Util::checkTrue(all_ccs.size() == m, "Error: number of cc files is not equal to the number of subpaths.");
	Util::checkTrue(all_closenesses.size() == m, "Error: number of closeness files is not equal to the number of subpaths.");
	Util::checkTrue(all_degrees.size() == m, "Error: number of degree files is not equal to the number of subpaths.");
	for (int i = 0; i < m; ++i) {
		std::string out_dir = out_root_dir + "/" + subpaths[i] + "/";
		int a = system(std::string("mkdir -p " + out_dir).c_str());
		Util::writeVector(all_ccs[i], out_dir + "cc.txt");
		Util::writeVector(all_closenesses[i], out_dir + "closeness.txt");
		Util::writeVector(all_degrees[i], out_dir + "degree.txt");
	}
}


