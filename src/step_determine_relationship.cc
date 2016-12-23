#include "step_determine_relationship.h"

#include <algorithm>
#include <vector>
#include <string>
#include <memory>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <set>
#include <tuple>
#include "unistd.h"

#include "runner.h"
#include "runner_manager.h"
#include "parameter.h"
#include "undirected_graph.h"
#include "util.h"
#include "util_boost.h"

StepDetermineRelationship StepDetermineRelationship::_step_determine_relationship;

StepDetermineRelationship::StepDetermineRelationship(): Runner() {
	_short_options = "h";
	_long_options = new struct option[100]{
		{"help",		no_argument,		NULL, OPT_HELP},
		{"if",			required_argument,	NULL, OPT_IF},
		{"ie",			required_argument,	NULL, OPT_IE},
		{NULL,			0,					NULL,  0 } //must end with {0, 0, 0, 0}
	};
	RunnerManager::instance()->install("step_determine_relationship", this);
}

void StepDetermineRelationship::help() {
	std::cout << "\nFunctionality: Based on contact network on eating time and the friendship network, determine the relationship between eating graph and friendship.";
	std::cout << "Option list:\n";
	std::cout << "\t* --help (or -h): [ no argument ] print this help information.\n";
	std::cout << "\t* --if: [ no argument ] input friendship gml file.\n";
	std::cout << "\t* --ie: [ no argument ] input eating time gml file.\n";
	std::cout << std::endl;
}

int StepDetermineRelationship::run(const Parameter& para) {
	UndirectedGraph gf(para.get_if());
	UndirectedGraph ge(para.get_ie());
	NodeSet nsf = gf.get_node_names();
	NodeSet nse = ge.get_node_names();
	NodeVec nv = Util::unset2vec(Util::getIntersection(nsf, nse));
	int n = nv.size();
	AdjacencyMatrix matf = gf.get_unweighted_adjacency_matrix(nv);
	AdjacencyMatrix mate = ge.get_unweighted_adjacency_matrix(nv);
	AdjacencyMatrix mat_min = Util::getMin(matf, mate);
	int ef = Util::count(matf, 1.0, 0.000001) / 2;
	int ee = Util::count(mate, 1.0, 0.000001) / 2;
	int em = Util::count(mat_min, 1.0, 0.000001) / 2;
	int max = n * (n - 1) / 2;
	std::cout << "ef = " << ef << std::endl;
	std::cout << "ee = " << ee << std::endl;
	std::cout << "em = " << em << std::endl;
	std::cout << "max = " << max << std::endl;
	double r = 1.0 * em / ee;
	double r_lb = 1.0 * get_bound(max, ee, ef, 0.05) / ee;
//	double r_lb = 1.0 * get_bound(get_hypergeometric_prob(max, ee, ef), 0.05) / ee;
	std::cout << "r = " << r << std::endl;
	std::cout << "lower bound = " << r_lb << std::endl;
	return 0;
}

int StepDetermineRelationship::get_bound(const int& N, const int& K, const int& n, const double& alpha) {
	double p = 1.0 * K / N;
	double t = sqrt(log(alpha) / (-2 * n));
	return (p + t) * n;
}

int StepDetermineRelationship::get_bound(const std::vector<double>& x, const double& alpha) {
	Util::checkTrue(x.size() > 0, "Error: distribution has no elements.");
	Util::checkTrue(alpha > 0 && alpha <= 1, "Error: illegal alpha.");
	std::cout << Util::getSum(x) << std::endl;
	Util::checkTrue(fabs(Util::getSum(x) - 1) < 0.0000001, "Error: probabilities do not sum to 1.");
	Util::checkTrue(std::all_of(x.begin(), x.end(), [] (const double& y) { return y >= 0 && y <= 1; }), "Error: contains illegal probability.");
	double sum = 0;
	for (int i = x.size(); i > 0; --i) {
		sum += x[i - 1];
		if (sum > alpha) {
			return i;
		}
	}
}

std::vector<double> StepDetermineRelationship::get_hypergeometric_prob(const int& N, const int& K, const int& n) {
	std::cout << "N = " << N << ", K = " << K << ", n = " << n << std::endl;
	Util::checkTrue(K > 0 && N >= K && n > 0 && N >= n, "Error: illegal hypergeometric parameters.");
	std::vector<double> res(n + 1, 0);
	res[0] = 1.0 * Util::choose(N - K, n) / Util::choose(N, n);
	std::cout << Util::choose(N - K, n) << " " << Util::choose(N, n) << " " << res[0] << std::endl;
	for (int k = 1; k < std::min(n + 1, K + 1); ++k) {
		double f1 = 1.0 * (K - k + 1) / k;
		double f2 = 1.0 * (n - k + 1) / (N + k - n - K);
		res[k] = res[k-1] * f1 * f2;
		std::cout << res[k] << std::endl;
	}
	for (int k = std::min(n + 1, K + 1); k < n + 1; ++k) {
		res[k] = 0;
	}
	return res;
}


