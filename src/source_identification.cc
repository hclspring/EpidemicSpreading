#include "source_identification.h"

#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cfloat>
#include <vector>

#include "undirected_graph.h"
#include "network.h"
#include "static_network.h"
#include "dynamic_network.h"
#include "pearl_network.h"
#include "disease_dynamics.h"
#include "bfs_tree.h"
#include "util.h"
#include "util_gsl.h"
#include "parameter.h"
#include "simulator.h"

/*
SourceIdentification::SourceIdentification(const Network& contact_network, const NodeSet& infected_nodes) {
	_infection_graph = contact_network.get_merged_graph() -> get_sub_graph(infected_nodes);
}

SourceIdentification::~SourceIdentification() {
}
*/

double SourceIdentification::calc_error_distance(Network& network, const std::string& real_source, const std::string& inferred_source) {
	std::shared_ptr<UndirectedGraph> merged_graph = network.get_merged_graph();
	double res = merged_graph->get_shortest_distance(real_source, inferred_source);
	return res;
}

double SourceIdentification::calc_error_distance(Network& network, const NodeSet& real_source, const NodeSet& inferred_source, const double& diff_count_penalty) {
	std::shared_ptr<UndirectedGraph> merged_graph = network.get_merged_graph();
	std::vector<std::vector<int>> distancess = merged_graph->get_shortest_distance(Util::getVector(real_source), Util::getVector(inferred_source));
	int n_real = real_source.size();
	int n_inferred = inferred_source.size();
	double res = 0.0;
	res += diff_count_penalty * abs(n_real - n_inferred);
	for (int i = 0; i < std::min(n_real, n_inferred); ++i) {
		res += Util::getMin(distancess[i]);
	}
	res /= n_real;
	return res;
}

double SourceIdentification::calc_detection_rate(const std::string& real_source, const std::string& inferred_source) {
	if (real_source.compare(inferred_source) == 0) {
		return 1;
	} else {
		return 0;
	}
}

/*
double SourceIdentification::calc_detection_rate(const NodeSet& real_source, const NodeSet& inferred_source) {
	//TODO
	return 0;
}
*/

std::string SourceIdentification::calc_source(Network& contact_network, const std::vector<DiseaseStage>& sim_res, const Parameter& para, const SrcIdnMethod& method, const bool& known_sim_days) {
	NodeSet nodes_infected = Simulator::get_nodeset_been_infected_from_sim_res(sim_res, para.get_disease(), contact_network);
	const int sim_days = para.get_max_sim_days();
	switch (method) {
		case SSE:
		case SSEBFS:
		case SJC:
		case JCE:
		case RG:
		case DA:
		case UB:
		case SLEUTH:
			return calc_source_with_infection_graph(contact_network, nodes_infected, para, method);
		case DMP:
			if (known_sim_days) {
				return calc_source_DMP_SIR_knowntime(contact_network, sim_res, sim_days, para);
			} else {
				return calc_source_DMP_SIR_unknowntime(contact_network, sim_res, 10 * sim_days, para);
			}
		case MCSM:
			if (known_sim_days) {
				return calc_source_MCSM_knowntime(contact_network, nodes_infected, sim_days, para);
			} else {
				return calc_source_MCSM_unknowntime(contact_network, nodes_infected, 10 * sim_days, para);
			}
		default:
			std::cerr << "Error for illegal method." << std::endl;
			exit(-1);
	}
}

std::string SourceIdentification::calc_source_with_infection_graph(Network& contact_network, const NodeSet& nodes_infected, const Parameter& para, const SrcIdnMethod& method) {
	switch (method) {
		case SSE:	return calc_source_SSE(contact_network, nodes_infected);
		case SSEBFS:return calc_source_SSEBFS(contact_network, nodes_infected);
		case SJC:	return calc_source_SJC(contact_network, nodes_infected);
		case JCE:	return calc_source_JCE(contact_network, nodes_infected);
		case RG:	return calc_source_RG(contact_network, nodes_infected, para);
		case DA:	return calc_source_DA(contact_network, nodes_infected);
		case UB:	return calc_source_UB(contact_network, nodes_infected, para.get_ub_r());
		case SLEUTH:return calc_source_NetSleuth(contact_network, nodes_infected);
		default:	std::cerr << "Error for illegal method." << std::endl; exit(-1);
	}
}

/*
NodeSet SourceIdentification::calc_source(const Network& contact_network, const NodeSet& nodes_been_infected, const int& max_source_count, const SrcIdnMethod& method) {
	//TODO
	return NodeSet();
}
*/
std::string SourceIdentification::calc_source_SSE(Network& contact_network, const NodeSet& nodes_been_infected) {
	std::shared_ptr<UndirectedGraph> merged_graph = contact_network.get_merged_graph();
	std::shared_ptr<UndirectedGraph> infection_graph = merged_graph -> get_sub_graph(nodes_been_infected);
	int n = nodes_been_infected.size();
	if (infection_graph->is_tree() == false) {
		NodeVec infected_nodes = Util::unset2vec(nodes_been_infected);
		std::vector<double> candidate_values;
		for (int i = 0; i < n; ++i) {
			std::string root_name = infected_nodes[i];
			std::shared_ptr<BFSTree> tree = std::make_shared<BFSTree>(*infection_graph, root_name);
			std::vector<double> rumor_centralities = SSE_calc_rumor_centralities(tree);
			double rumor_centrality = rumor_centralities[tree->get_node_index(root_name)];
			NodeVec bfs_nodes = infection_graph->get_bfs_nodes(root_name);
			long permitted_permutation_count = SSE_calc_permitted_permutation_count(merged_graph, bfs_nodes);
			candidate_values.push_back(rumor_centrality / permitted_permutation_count);
		}
		int max_index = Util::getMaxIndex(candidate_values);
		return infected_nodes[max_index];
	} else {
		int root_index = Util::gen_rand_int(n);
		std::string root_name = infection_graph->get_node_name(root_index);
		std::shared_ptr<BFSTree> tree = std::make_shared<BFSTree>(*infection_graph, root_name);
		std::vector<double> rumor_centralities = SSE_calc_rumor_centralities(tree);
		int source_index = Util::getMaxIndex(rumor_centralities);
		return tree->get_node_name(source_index);
	}
}

std::string SourceIdentification::calc_source_SSEBFS(Network& contact_network, const NodeSet& nodes_been_infected) {
	std::shared_ptr<UndirectedGraph> infection_graph = contact_network.get_merged_graph() -> get_sub_graph(nodes_been_infected);
	int n = infection_graph->get_node_size();
	std::vector<double> max_rumor_centralities(n, -1);
	for (int root_index = 0; root_index < n; ++root_index) {
		std::string root_name = infection_graph->get_node_name(root_index);
		std::shared_ptr<BFSTree> tree = std::make_shared<BFSTree>(*infection_graph, root_name);
		std::vector<double> rumor_centralities = SSE_calc_rumor_centralities(tree);
		int tree_source_index = Util::getMaxIndex(rumor_centralities);
		int source_index = infection_graph->get_node_index(tree->get_node_name(tree_source_index));
		double max_value = Util::getMax(rumor_centralities);
		if (max_rumor_centralities[source_index] < max_value) {
			max_rumor_centralities[source_index] = max_value;
		}
	}
	int source_index = Util::getMaxIndex(max_rumor_centralities);
	std::string res = infection_graph->get_node_name(source_index);
	return res;
}

std::string SourceIdentification::calc_source_SJC(Network& contact_network, const NodeSet& nodes_being_infected) {
	std::shared_ptr<UndirectedGraph> merged_graph = contact_network.get_merged_graph();
	if (merged_graph->is_tree() == false) {
		std::cerr << "Error because the graph is not a tree." << std::endl;
		exit(-1);
	}
	int n = merged_graph->get_node_size();
	int m = nodes_being_infected.size();
//	int total_received_count = 0; // when the algorithm finishes, this number should be m*n.
	std::vector<std::vector<int>> received_message_time(n, std::vector<int>(m, -1)); // record whether node i received the message of infected node j.
	std::vector<int> received_message_count(n, 0);
	// set the infecting nodes and map
	std::vector<std::string> infecting_nodes;
	std::unordered_map<std::string, int> infecting_map;
	for (NodeSet::const_iterator it = nodes_being_infected.begin(); it != nodes_being_infected.end(); ++it) {
		infecting_map.insert(std::make_pair(*it, infecting_nodes.size()));
		infecting_nodes.push_back(*it);
	}
	// set the initial status
	for (int i = 0; i < m; ++i) {
		int j = merged_graph->get_node_index(infecting_nodes[i]);
		received_message_time[j][i] = 0;
		received_message_count[j] = 1;
//		total_received_count++;
	}
	int t = 0;
	while (!SJC_check_finished(received_message_count, m)) {
		++t;
		for (int i = 0; i < n; ++i) { // broadcast from every node
			for (int j = 0; j < m; ++j) {
				if (received_message_time[i][j] < t) { // broadcast every received message
					std::string from = merged_graph->get_node_name(i);
					NeighborList nl = merged_graph->get_neighbor_list(from);
					for (NeighborList::iterator it = nl.begin(); it != nl.end(); ++it) {
						// broadcast to every neighbor
						std::string to = it->get_name();
						int k = merged_graph->get_node_index(to);
						if (received_message_time[k][j] < 0) {
							// the neighbor has not received the message
							received_message_time[k][j] = t;
							received_message_count[k]++;
//							total_received_count++;
						}
					}
				}
			}
		}
	}
	std::vector<int> received_nodes;
	for (int i = 0; i < n; ++i) {
		if (received_message_count[i] >= m) {
			received_nodes.push_back(i);
		}
	}
	std::string res;
	if (received_nodes.size() > 1) {
		std::vector<int> time_sum(received_nodes.size(), 0);
		int min_value = INT_MAX, min_index = -1;
		for (int i = 0; i < received_nodes.size(); ++i) {
			for (int j = 0; j < m; ++j) {
				time_sum[i] += received_message_time[received_nodes[i]][j];
			}
			if (time_sum[i] < min_value) {
				min_index = received_nodes[i];
				min_value = time_sum[i];
			}
		}
		res = merged_graph->get_node_name(min_index);
	} else {
		res = merged_graph->get_node_name(received_nodes[0]);
	}
	return res;
}

std::string SourceIdentification::calc_source_JCE(Network& contact_network, const NodeSet& nodes_observed_infected) {
	std::shared_ptr<UndirectedGraph> merged_graph = contact_network.get_merged_graph();
	if (merged_graph->is_tree() == false) {
		std::cerr << "Error because the graph is not a tree." << std::endl;
		exit(-1);
	}
	if (nodes_observed_infected.size() == 0) {
		std::cerr << "Error for no nodes observed to construct a spanning tree." << std::endl;
		exit(-1);
	}
	// Generate a temp spanning tree to find root candidates.
	std::string temp_root = *(nodes_observed_infected.begin());
	std::shared_ptr<BFSTree> spanning_tree = std::make_shared<BFSTree>(*merged_graph, temp_root, nodes_observed_infected);
	NodeVec root_candidates;
	for (int i = 0; i < spanning_tree->get_node_size(); ++i) {
		std::string nodename = merged_graph->get_node_name(i);
		if (merged_graph->get_degree(nodename) > 1) {
			root_candidates.push_back(nodename);
		}
	}
	// Find a root from the candidates.
	int n = root_candidates.size();
	if (n == 0) {
		std::cerr << "Error finding no root candidate." << std::endl;
		exit(-1);
	}
	std::string root_name = root_candidates[Util::gen_rand_int(n)];

	std::shared_ptr<BFSTree> tree = std::make_shared<BFSTree>(*spanning_tree, root_name);
	n = tree->get_node_size();
	int root_index = tree->get_node_index(root_name);
	std::vector<JceInfo> jce_info(n, {-1,-1,-1,-1});
	JCE_upward_message_passing(tree, root_index, jce_info);
	std::string res = JCE_downward_message_passing(tree, root_index, jce_info);
	return res;
}

std::string SourceIdentification::calc_source_RG(Network& contact_network, const NodeSet& nodes_observed_infected, const Parameter& para) {
	double explicite_rate = 1;
	std::shared_ptr<UndirectedGraph> merged_graph = contact_network.get_merged_graph();
	int merged_nodes = merged_graph->get_node_size();
	std::vector<double> criteria;
	double infect_rate = DiseaseDynamics::get_infect_prob(para);
	for (int i = 0; i < merged_nodes; ++i) {
		std::string root_name = merged_graph->get_node_name(i);
		std::shared_ptr<BFSTree> tree = std::make_shared<BFSTree>(*merged_graph, root_name, nodes_observed_infected);
		std::shared_ptr<UndirectedGraph> induced_graph = merged_graph->get_sub_graph(tree->get_node_names());
		int Fv_star = RG_for_every_tree(tree, induced_graph);
		int n = tree->get_node_size();
		//criteria.push_back(n * log(infect_rate) + (Fv_star - n) * log(1 - infect_rate) + (n - nodes_observed_infected.size()) * log(1 - explicite_rate));
		criteria.push_back(n * log(infect_rate) + (Fv_star - n) * log(1 - infect_rate));
	}
	int max_index = Util::getMaxIndex(criteria);
	std::string res = merged_graph->get_node_name(max_index);
	return res;
}

std::string SourceIdentification::calc_source_DA(Network& contact_network, const NodeSet& nodes_been_infected) {
	NodeSet sources = calc_source_DA(contact_network, nodes_been_infected, 1);
	return *(sources.begin());
}

NodeSet SourceIdentification::calc_source_DA(Network& contact_network, const NodeSet& nodes_been_infected, const int& source_count) {
	std::shared_ptr<UndirectedGraph> merged_graph = contact_network.get_merged_graph();
	std::shared_ptr<UndirectedGraph> infection_graph = merged_graph->get_sub_graph(nodes_been_infected);
	int n = infection_graph->get_node_size();
	NodeSet res;
	if (n == 1) {
		res.insert(infection_graph->get_node_name(0));
		return res;
	}
	AdjacencyMatrix adjacency_matrix = infection_graph->get_adjacency_matrix();
	std::vector<double> original_eigenvalues = UtilGsl::calcSymmMatrixEigenvalues(adjacency_matrix);
	double original_eigenvalue_max = original_eigenvalues[0];
	std::vector<double> temp_eigenvalues, new_max_eigenvalues(n, 0);
	for (int i = 0; i < n; ++i) {
		AdjacencyMatrix temp_matrix = Util::eliminateFromSquareMatrix(adjacency_matrix, i);
		temp_eigenvalues = UtilGsl::calcSymmMatrixEigenvalues(temp_matrix);
		new_max_eigenvalues[i] = fabs(temp_eigenvalues[0] - original_eigenvalue_max) / original_eigenvalue_max;
	}
	std::vector<int> sorted_indices = Util::getSortedIndices(new_max_eigenvalues, false);
	for (int i = 0; i < source_count; ++i) {
		res.insert(infection_graph->get_node_name(sorted_indices[i]));
	}
	return res;
}

std::string SourceIdentification::calc_source_UB(Network& contact_network, const NodeSet& nodes_been_infected, const double& ub_r) {
	std::shared_ptr<UndirectedGraph> merged_graph = contact_network.get_merged_graph();
	std::shared_ptr<UndirectedGraph> infection_graph = merged_graph->get_sub_graph(nodes_been_infected);
//	std::vector<double> betweennesses = infection_graph->get_betweenness(infection_graph->get_node_names_vec());
	std::vector<double> betweennesses = infection_graph->get_betweenness_all();
	int n = infection_graph->get_node_size();
	std::vector<double> ubs(n, 0);
	for (int i = 0; i < n; ++i) {
		ubs[i] = betweennesses[i] / pow(infection_graph->get_degree(infection_graph->get_node_name(i)), ub_r);
	}
	int max_index = Util::getMaxIndex(ubs);
	return infection_graph->get_node_name(max_index);
}

std::string SourceIdentification::calc_source_DMP_SIR_unknowntime(Network& contact_network, const std::vector<DiseaseStage>& sim_res, const int& max_sim_days, const Parameter& para) {
	if (para.get_net_type() != STATIC) {
		std::cerr << "Error because the network is not STATIC." << std::endl;
		exit(-1);
	}
	if (para.get_disease() != SIR) {
		std::cerr << "Error because disease is not SIR." << std::endl;
		exit(-1);
	}
	if (max_sim_days <= 0) {
		std::cerr << "Error for " << max_sim_days << " simulating days." << std::endl;
		exit(-1);
	}
	std::vector<double> partition_functions;
	std::vector<std::vector<double>> likelihoods;
	partition_functions = DMP_calc_partition_function(contact_network, max_sim_days, sim_res, para, likelihoods);
	int sim_day_index = Util::getMaxIndex(partition_functions);
	int n = contact_network.get_node_size();
	int max_likelihood_index = -1;
	double max_likelihood_value = 0;
	for (int i = 0; i < n; ++i) {
		if (likelihoods[i][sim_day_index] > max_likelihood_value) {
			max_likelihood_value = likelihoods[i][sim_day_index];
			max_likelihood_index = i;
		}
	}
	int seed_index = max_likelihood_index;
	return contact_network.get_node_name(seed_index);
}

std::string SourceIdentification::calc_source_DMP_SIR_knowntime(Network& contact_network, const std::vector<DiseaseStage>& sim_res, const int& days, const Parameter& para) {
	if (para.get_net_type() != STATIC) {
		std::cerr << "Error because the network is not STATIC." << std::endl;
		exit(-1);
	}
	if (para.get_disease() != SIR) {
		std::cerr << "Error because disease is not SIR." << std::endl;
		exit(-1);
	}
	int n = contact_network.get_node_size();
	if (days <= 0) {
		std::cerr << "Error for " << days << " simulating days." << std::endl;
		exit(-1);
	}
	std::vector<double> partition_likelihoods;
	DMP_calc_partition_function_with_same_sim_days(contact_network, days, sim_res, para, partition_likelihoods);
	int seed_index = Util::getMaxIndex(partition_likelihoods);
	return contact_network.get_node_name(seed_index);
}

std::string SourceIdentification::calc_source_MCSM_unknowntime(Network& contact_network, const NodeSet& nodes_been_infected, const int& max_sim_days, const Parameter& para) {
	bool fixed_sim_days = false;
	return calc_source_MCSM_diffsimdays(contact_network, nodes_been_infected, max_sim_days, fixed_sim_days, para);
}

std::string SourceIdentification::calc_source_MCSM_knowntime(Network& contact_network, const NodeSet& nodes_been_infected, const int& sim_days, const Parameter& para) {
	bool fixed_sim_days = true;
	return calc_source_MCSM_diffsimdays(contact_network, nodes_been_infected, sim_days, fixed_sim_days, para);
}

std::string SourceIdentification::calc_source_NetSleuth(Network& contact_network, const NodeSet& nodes_been_infected) {
	std::shared_ptr<UndirectedGraph> merged_graph = contact_network.get_merged_graph();
	std::vector<std::vector<double>> laplacian_matrix = merged_graph->get_laplacian_matrix();
	IndexSet infected_index_set;
	std::for_each(nodes_been_infected.begin(), nodes_been_infected.end(),
			[merged_graph, &infected_index_set] (const std::string& x) { infected_index_set.insert(merged_graph->get_node_index(x)); });
	std::vector<std::vector<double>> infected_laplacian_matrix = Util::getSubMatrix(laplacian_matrix, infected_index_set);
	std::vector<std::vector<double>> eigenvectors = UtilGsl::calcSymmMatrixEigenvectors(infected_laplacian_matrix);
	int seed_index_index = Util::getMaxIndex(eigenvectors.back());
	IndexVec infected_index_vec = Util::unset2vec(infected_index_set);
	std::sort(infected_index_vec.begin(), infected_index_vec.end());
	int seed_index = infected_index_vec[seed_index_index];
	return merged_graph->get_node_name(seed_index);
}

/*
std::string SourceIdentification::calc_source_DMP(Network& contact_network, const NodeSet& nodes_been_infected, const int& max_sim_days, const Parameter& para) {
	if (para.get_net_type() != STATIC) {
		std::cerr << "Error because the network is not STATIC." << std::endl;
		exit(-1);
	}
	if (para.get_disease() != SIR) {
		std::cerr << "Error because disease is not SIR." << std::endl;
		exit(-1);
	}
	int n = contact_network.get_node_size();
	if (max_sim_days <= 0) {
		std::cerr << "Error for " << max_sim_days << " simulating days." << std::endl;
		exit(-1);
	}
	//TODO
}



std::string SourceIdentification::calc_source_DMP(Network& contact_network, const NodeSet& nodes_been_infected, const int& days, const Parameter& para) {
	//TODO
}
*/
long SourceIdentification::SSE_calc_permitted_permutation_count(const std::shared_ptr<UndirectedGraph>& graph, const NodeVec& bfs_nodes) {
	int degree = graph->get_degree(bfs_nodes[0]);
	int sum = degree;
	long res = sum;
	int n = bfs_nodes.size();
	for (int i = 1; i < n; ++i) {
		degree = graph->get_degree(bfs_nodes[i]);
		sum += degree - 2;
		res *= sum;
	}
	if (res == 0) {
		res = 1;
	}
	return res;
}

std::vector<double> SourceIdentification::SSE_calc_rumor_centralities(const std::shared_ptr<BFSTree>& tree) {
	int n = tree->get_node_size();
	std::vector<int> t(n, 0), p(n, 0);
	std::vector<double> res(n, 0);
	int root_index = tree->get_node_index(tree->get_root());
	SSE_calc_t_p(tree, t, p, root_index);
	SSE_calc_r(tree, t, p, res, root_index);
	return res;
}

void SourceIdentification::SSE_calc_t_p(const std::shared_ptr<BFSTree>& tree, std::vector<int>& t, std::vector<int>& p, const int& cur_index) {
	std::string cur_node = tree->get_node_name(cur_index);
	if (t[cur_index] > 0 && p[cur_index] > 0) {
		return;
	} else if (tree->is_leaf(cur_node)) {
		t[cur_index] = 1;
		p[cur_index] = 1;
	} else {
		int sum = 0, prod = 1;
		NodeSet children = tree->get_children(cur_node);
		for (NodeSet::iterator it = children.begin(); it != children.end(); ++it) {
			int child_index = tree->get_node_index(*it);
			SSE_calc_t_p(tree, t, p, child_index);
			sum += t[child_index];
			prod *= p[child_index];
		}
		t[cur_index] = sum + 1;
		p[cur_index] = prod * t[cur_index];
	}
}

void SourceIdentification::SSE_calc_r(const std::shared_ptr<BFSTree>& tree, const std::vector<int>& t, const std::vector<int>& p, std::vector<double>& r, const int& cur_index) {
	std::string cur_node = tree->get_node_name(cur_index);
	int n = tree->get_node_size();
	NodeSet children = tree->get_children(cur_node);
	if (tree->is_root(cur_node)) {
		int prod = 1;
		for (NodeSet::iterator it = children.begin(); it != children.end(); ++it) {
			int child_index = tree->get_node_index(*it);
			prod *= p[child_index];
		}
		r[cur_index] = 1.0 * Util::factorial(n - 1) / prod;
	} else {
		int parent_index = tree->get_node_index(tree->get_parent(cur_node));
		r[cur_index] = r[parent_index] * t[cur_index] / (n - t[cur_index]);
	}
	for (NodeSet::iterator it = children.begin(); it != children.end(); ++it) {
		int child_index = tree->get_node_index(*it);
		SSE_calc_r(tree, t, p, r, child_index);
	}
}


bool SourceIdentification::SJC_check_finished(const std::vector<int>& received_message_count, const int& target) {
	int n = received_message_count.size();
	for (int i = 0; i < n; ++i) {
		if (received_message_count[i] >= target) {
			return true;
		}
	}
	return false;
}

void SourceIdentification::JCE_upward_message_passing(const std::shared_ptr<BFSTree>& tree, const int& cur_index, std::vector<JceInfo>& jce_info) {
	std::string cur_name = tree->get_node_name(cur_index);
	if (tree->is_leaf(cur_name)) {
		jce_info[cur_index].longest_path = 1;
	} else {
		std::vector<int> children_paths;
		NodeSet children = tree->get_children(cur_name);
		std::vector<std::string> children_names;
		for(NodeSet::iterator it = children.begin(); it != children.end(); ++it) {
			children_names.push_back(*it);
			int children_index = tree->get_node_index(*it);
			JCE_upward_message_passing(tree, children_index, jce_info);
			children_paths.push_back(jce_info[children_index].longest_path);
		}
		int first_index = Util::getMaxIndex(children_paths);
		std::string first_name = children_names[first_index];
		jce_info[cur_index].first_child_index = tree->get_node_index(first_name);
		std::sort(children_paths.begin(), children_paths.end());
		int n = children_paths.size();
		int maximum = children_paths[n-1];
		int second_maximum = 0;
		if (n > 1) {
			second_maximum = children_paths[n-2];
		}
		jce_info[cur_index].first_child_path = maximum;
		jce_info[cur_index].second_child_path = second_maximum;
		jce_info[cur_index].longest_path = 1 + maximum;
	}
	return;
}


std::string SourceIdentification::JCE_downward_message_passing(const std::shared_ptr<BFSTree>& tree, const int& cur_index, std::vector<JceInfo>& jce_info) {
	std::string cur_name = tree->get_node_name(cur_index);
	if (tree->is_root(cur_name) == false) {
		std::string parent_name = tree->get_parent(cur_name);
		int parent_index = tree->get_node_index(parent_name);
		jce_info[cur_index].second_child_path = std::max(jce_info[cur_index].second_child_path, jce_info[parent_index].g_value);
	}
	if (jce_info[cur_index].first_child_path - jce_info[cur_index].second_child_path <= 1) {
		return tree->get_node_name(cur_index);
	} else {
		jce_info[cur_index].g_value = jce_info[cur_index].second_child_path + 1;
		return JCE_downward_message_passing(tree, jce_info[cur_index].first_child_index, jce_info);
	}
}

int SourceIdentification::RG_for_every_tree(std::shared_ptr<BFSTree>& tree, std::shared_ptr<UndirectedGraph>& induced_graph) {
	int n = tree->get_node_size();
	std::vector<int> depths = tree->calc_depths();
	std::vector<int> heights = tree->calc_heights();
	std::vector<int> heights_other(n, 0);
	std::string root_name = tree->get_root();
	int root_index = tree->get_node_index(root_name);
	int root_height = heights[root_index];
	for (int d = root_height; d > 0; --d) {
		std::vector<int> x_indices = RG_indices_with_same_depth_and_increasing_height(depths, d, heights);
		for (int x_indices_i = 0; x_indices_i < x_indices.size(); ++x_indices_i) {
			int x_index = x_indices[x_indices_i];
			std::string x_name = tree->get_node_name(x_index);
			NodeSet Y;
			NeighborList nl = induced_graph->get_neighbor_list(x_name);
			for (NeighborList::iterator nl_it = nl.begin(); nl_it != nl.end(); ++nl_it) {
				std::string y_name = nl_it->get_name();
				int y_index = tree->get_node_index(y_name);
			//	if (depths[y_index] <= heights[root_index] - heights[x_index] - 1) {
				if (depths[y_index] <= depths[x_index] - 1) {
					Y.insert(y_name);
					if (y_name.compare(tree->get_parent(x_name)) == 0) {
						heights_other[y_index] = RG_height_with_exception(tree, y_name, x_name, heights);
					} else {
						heights_other[y_index] = heights[y_index];
					} // end if
				} // end if
			} // end for
			// Choose y in Y with largest depth and height_other.
			int y_max_index, y_max_depth = -1, y_max_height_other = -1;
			for (NodeSet::iterator it = Y.begin(); it != Y.end(); ++it) {
				int y_index = tree->get_node_index(*it);
				if (depths[y_index] > y_max_depth) {
					y_max_index = y_index;
					y_max_depth = depths[y_index];
					y_max_height_other = heights_other[y_index];
				} else if (depths[y_index] == y_max_depth 
						&& heights_other[y_index] > y_max_height_other) {
					y_max_index = y_index;
					y_max_height_other = heights_other[y_index];
				} // end if
			} // end for
			int y_index = y_max_index;
			// Rewire x from its parent to y.
			std::string y_name = tree->get_node_name(y_index);
			tree->rewire(x_name, tree->get_parent(x_name), y_name);
			// Modify depth of all nodes in the subtree(x).
			int delta_x = depths[y_index] + 1 - depths[x_index];
			if (delta_x != 0) {
				NodeSet subtree_nodes = tree->get_subtree_nodes(x_name);
				for (NodeSet::iterator it = subtree_nodes.begin(); it != subtree_nodes.end(); ++it) {
					int sub_index = tree->get_node_index(*it);
					depths[sub_index] += delta_x;
				} // end for
			} // end if
			// Modify y to its ancestors.
			int delta_y = std::max(0, heights[x_index] + 1 - heights[y_index]);
			while (delta_y != 0) {
				heights[y_index] += delta_y;
				std::string y_parent_name = tree->get_parent(y_name);
				int y_parent_index = tree->get_node_index(y_parent_name);
				delta_y = std::max(0, heights[y_index] + 1 - heights[y_parent_index]);
				y_name = y_parent_name;
				y_index = y_parent_index;
			} // end while
		} // end for
	} // end for
	int res = 2 * heights[root_index];
	NodeSet all_nodes = tree->get_node_names();
	for (NodeSet::iterator it = all_nodes.begin(); it != all_nodes.end(); ++it) {
		res += (tree->get_degree(*it) - 2) * heights[tree->get_node_index(*it)];
	}
	return res;
}


std::vector<int> SourceIdentification::RG_indices_with_same_depth_and_increasing_height(const std::vector<int>& depths, const int& d, const std::vector<int>& heights) {
	struct mystruct {
		int index;
		int height;
	};
	struct myclass {
		bool operator() (const mystruct& i, const mystruct& j) {
			return (i.height < j.height);
		}
	} myobject;

	std::vector<mystruct> heights_with_indices;
	for (int i = 0; i < depths.size(); ++i) {
		if (depths[i] == d) {
			heights_with_indices.push_back({i, heights[i]});
		}
	}
	std::sort(heights_with_indices.begin(), heights_with_indices.end(), myobject);
	std::vector<int> res;
	for (int i = 0; i < heights_with_indices.size(); ++i) {
		res.push_back(heights_with_indices[i].index);
	}
	return res;
}

int SourceIdentification::RG_height_with_exception(const std::shared_ptr<BFSTree>& tree, const std::string& y_name, const std::string& x_name, const std::vector<int>& heights) {
	NodeSet children = tree->get_children(y_name);
	std::vector<int> candidate_heights;
	for (NodeSet::iterator it = children.begin(); it != children.end(); ++it) {
		if (x_name.compare(*it) != 0) {
			int index = tree->get_node_index(*it);
			candidate_heights.push_back(heights[index] + 1);
		}
	}
	if (candidate_heights.size() == 0) {
		return 0;
	} else {
		return Util::getMax(candidate_heights);
	}
}

int SourceIdentification::RG_height_with_exception(const std::shared_ptr<BFSTree>& tree, const std::string& y_name, const std::string& x_name) {
	if (y_name.compare(x_name) == 0) {
		return -1;
	} else if (tree->is_leaf(y_name)) {
		return 0;
	} else {
		NodeSet children = tree->get_children(y_name);
		int children_height_max = -1;
		for (NodeSet::iterator it = children.begin(); it != children.end(); ++it) {
			int cur_height = RG_height_with_exception(tree, *it, x_name);
			if (cur_height > children_height_max) {
				children_height_max = cur_height;
			}
		}
		return 1 + children_height_max;
	}
}

std::vector<double> SourceIdentification::DMP_calc_partition_function(Network& contact_network, const int& max_sim_days, const std::vector<DiseaseStage>& res_stages, const Parameter& para, std::vector<std::vector<double>>& likelihoods) {
	std::vector<double> res(max_sim_days, 0);
	int n = contact_network.get_node_size();
	likelihoods.clear();
	likelihoods.resize(n, std::vector<double>(max_sim_days, 0)); // row: seeds; column: sim_days
	for (int i = 0; i < n; ++i) {
		std::string seed = contact_network.get_node_name(i);
		likelihoods[i] = DMP_calc_likelihoods_with_same_seed(contact_network, seed, max_sim_days, res_stages, para);
	}
	for (int i = 0; i < max_sim_days; ++i) {
		for (int j = 0; j < n; ++j) {
			res[i] += likelihoods[j][i];
		}
	}
	return res;
}

double SourceIdentification::DMP_calc_partition_function_with_same_sim_days(Network& contact_network, const int& sim_days, const std::vector<DiseaseStage>& res_stages, const Parameter& para, std::vector<double>& likelihoods) {
	likelihoods.clear();
	int n = contact_network.get_node_size();
	likelihoods.resize(n, 0);
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		double likelihood = DMP_calc_likelihood_with_fixed_seed_and_simdays(contact_network, contact_network.get_node_name(i), sim_days, res_stages, para);
		sum += likelihood;
		likelihoods[i] = likelihood;
	}
	return sum;
}

std::vector<double> SourceIdentification::DMP_calc_likelihoods_with_same_seed(Network& contact_network, const std::string& seed_node, const int& sim_days, const std::vector<DiseaseStage>& res_stages, const Parameter& para) {
	std::shared_ptr<DMPNeiMessage> init_msg_ptr = std::make_shared<DMPNeiMessage>(DMP_init_message_step_stat_net(contact_network, seed_node));
	std::shared_ptr<DMPNeiMessage> new_msg_ptr = init_msg_ptr, old_msg_ptr;
	DMPMsgVecPtr pss0 = DMP_init_pss0_step_stat_net(contact_network, seed_node);
	std::vector<double> res(sim_days, 0);
	for (int t = 0; t < sim_days; ++t) {
		old_msg_ptr = new_msg_ptr;
		new_msg_ptr = std::make_shared<DMPNeiMessage>(DMP_pass_message_step_stat_net(contact_network, *old_msg_ptr, pss0, para));
		res[t] = DMP_calc_likelihood(res_stages, *new_msg_ptr);
	}
	return res;
}

double SourceIdentification::DMP_calc_likelihood_with_fixed_seed_and_simdays(Network& contact_network, const std::string& seed_node, const int& sim_days, const std::vector<DiseaseStage>& res_stages, const Parameter& para) {
	std::shared_ptr<DMPNeiMessage> init_msg_ptr = std::make_shared<DMPNeiMessage>(DMP_init_message_step_stat_net(contact_network, seed_node));
	std::shared_ptr<DMPNeiMessage> new_msg_ptr = init_msg_ptr, old_msg_ptr;
	DMPMsgVecPtr pss0 = DMP_init_pss0_step_stat_net(contact_network, seed_node);
	for (int t = 0; t < sim_days; ++t) {
		old_msg_ptr = new_msg_ptr;
		new_msg_ptr = std::make_shared<DMPNeiMessage>(DMP_pass_message_step_stat_net(contact_network, *old_msg_ptr, pss0, para));
	}
	return DMP_calc_likelihood(res_stages, *new_msg_ptr);
}

/*
void SourceIdentification::DMP_pass_message_day(Network& contact_network, const DMPMessage& pre_message, const int& day, const DMPOneMessage& pss0, DMPMessage& next_message) {
	int parts = contact_network.get_parts();
	for (int i = 0; i < parts; ++i) {
		//TODO
	}
}
*/

DMPNeiMessage SourceIdentification::DMP_init_message_step_stat_net(Network& contact_network, const std::string& seed_node) {
	int n = contact_network.get_node_size();
	int seed_index = contact_network.get_node_index(seed_node);
	if (seed_index < 0) {
		std::cerr << "Error for node \"" << seed_node << "\" not exist." << std::endl;
		exit(-1);
	}
	DMPNeiMessage res = DMP_init_nei_msg(n);
	for (int i = 0; i < n; ++i) {
		std::string node = contact_network.get_node_name(i);
		NeighborList nl = contact_network.get_neighbor_list(node, 0, 0);
		for (NeighborList::const_iterator it = nl.begin(); it != nl.end(); ++it) {
			std::string j_node = it->get_name();
			DMPMsg2Nei& temp = DMP_set_msg2nei(res, i, j_node);
			temp.to_nei_theta = 1;
			temp.to_nei_phi = (i == seed_index ? 1 : 0);
			temp.to_nei_psd = 1 - temp.to_nei_phi;
		}
		res.prs->at(i) = 0;
		res.pis->at(i) = (i == seed_index ? 1 : 0);
		res.pss->at(i) = 1 - res.pis->at(i);
	}
	return res;
}

DMPMsgVecPtr SourceIdentification::DMP_init_pss0_step_stat_net(Network& contact_network, const std::string& seed_node) {
	int n = contact_network.get_node_size();
	int seed_index = contact_network.get_node_index(seed_node);
	if (seed_index < 0) {
		std::cerr << "Error for node \"" << seed_node << "\" not exist." << std::endl;
		exit(-1);
	}
	DMPMsgVecPtr res = std::make_shared<DMPMsgVec>(n, 1);
	res->at(seed_index) = 0;
	return res;
}

DMPNeiMessage SourceIdentification::DMP_pass_message_step_stat_net(Network& contact_network, const DMPNeiMessage& pre_messages, const DMPMsgVecPtr& pss0, const Parameter& para) {
	int n = contact_network.get_node_size();
	DMPNeiMessage res = DMP_init_nei_msg(n);
	DMP_calc_new_thetas_stat_net(contact_network, pre_messages, para, res);
	DMP_calc_new_psds_stat_net(contact_network, pss0, res);
	DMP_calc_new_phis_stat_net(contact_network, pre_messages, para, res);
	DMP_calc_new_prs_stat_net(contact_network, pre_messages, para, res);
	DMP_calc_new_pss_stat_net(contact_network, pss0, res);
	DMP_calc_new_pis_stat_net(contact_network, res);
	return res;
}

double SourceIdentification::DMP_calc_likelihood(const std::vector<DiseaseStage>& stages, const DMPNeiMessage& messages) {
	double res = 1.0, prob;
	for (int i = 0; i < stages.size(); ++i) {
		switch (stages[i]) {
			case SUSCEPTIBLE:	prob = messages.pss->at(i); break;
			case INFECTIOUS:	prob = messages.pis->at(i); break;
			case RECOVERED:		prob = messages.prs->at(i); break;
			default: {
						 std::cerr << "Error finding undefined disease stage." << std::endl;
						 exit(-1);
					 }
		}
		res *= prob;
	}
	return res;
}

void SourceIdentification::DMP_init_mat(const int& n, DMPMsgMatPtr& mp) {
	mp = std::make_shared<DMPMsgMat>(n, DMPMsgVec(n, 1));
}

void SourceIdentification::DMP_init_vec(const int& n, DMPMsgVecPtr& vp) {
	vp = std::make_shared<DMPMsgVec>(n, 1);
}

DMPMatMessage SourceIdentification::DMP_init_mat_msg(const int& n, const DMPMatMessage& pre_msg) {
	DMPMatMessage res;
	DMP_init_mat(n, res.thetas);
	DMP_init_mat(n, res.phis);
	res.pre_psds = pre_msg.psds;
	DMP_init_mat(n, res.psds);
	DMP_init_vec(n, res.pss);
	DMP_init_vec(n, res.prs);
	DMP_init_vec(n, res.pis);
	return res;
}

DMPNeiMessage SourceIdentification::DMP_init_nei_msg(const int& n) {
	DMPNeiMessage res;
	res.neighbors = std::make_shared<DMPMsgNeiListVec>(n, DMPMsgNeiList());
	res.pss = std::make_shared<DMPMsgVec>(n, 0);
	res.prs = std::make_shared<DMPMsgVec>(n, 0);
	res.pis = std::make_shared<DMPMsgVec>(n, 0);
	return res;
}

const DMPMsg2Nei& SourceIdentification::DMP_get_msg2nei(const DMPNeiMessage& messages, const int& from_index, const std::string& to_node) {
	DMPMsgNeiList& msg_list = messages.neighbors->at(from_index);
	for (DMPMsgNeiList::const_iterator it = msg_list.begin(); it != msg_list.end(); ++it) {
		if (it->nei_name.compare(to_node) == 0) {
			return *it;
		}
	}
	std::cerr << "Error for finding old message to neighbor \"" << to_node << "\" failed." << std::endl;
	exit(-1);
}

DMPMsg2Nei& SourceIdentification::DMP_set_msg2nei(DMPNeiMessage& messages, const int& from_index, const std::string& to_node) {
	DMPMsgNeiList& msg_list = messages.neighbors->at(from_index);
	for (DMPMsgNeiList::iterator it = msg_list.begin(); it != msg_list.end(); ++it) {
		if (it->nei_name.compare(to_node) == 0) {
			return *it;
		}
	}
	msg_list.push_back({to_node, -1, -1, -1}); // We set default new values to be -1.
	return msg_list.back();
}

void SourceIdentification::DMP_calc_new_thetas_stat_net(Network& contact_network, const DMPNeiMessage& old_messages, const Parameter& para, DMPNeiMessage& new_messages) {
	int n = contact_network.get_node_size();
	DMPMsgNeiPtr new_msg_neighbors = new_messages.neighbors;
	for (int i = 0; i < n; ++i) {
		std::string cur_node = contact_network.get_node_name(i);
		NeighborList neighbor_list = contact_network.get_neighbor_list(cur_node, 0, 0);
		DMPMsgNeiList& msg_list = new_msg_neighbors->at(i);
		for (NeighborList::const_iterator it = neighbor_list.begin(); it != neighbor_list.end(); ++it) {
			std::string nei_name = it->get_name();
			double weight = it->get_weight();
			double lambda = DiseaseDynamics::get_infect_prob(weight, para);
			// cur_node --> nei_name
			DMPMsgNeiList::iterator jt;
			double new_theta = DMP_calc_new_theta_stat_net(old_messages, i, nei_name, lambda);
			DMP_set_msg2nei(new_messages, i, nei_name).to_nei_theta = new_theta;
			/*
			for (jt = msg_list.begin(); jt != msg_list.end(); ++jt) {
				if (jt->nei_name.compare(nei_name) == 0) {
					jt->to_nei_theta = new_theta;
					break;
				}
			}
			if (jt == msg_list.end()) {
				msg_list.push_back({nei_name, new_theta, -1, -1}); // We take default new values to be -1.
			}
			*/
		}
	}
}


double SourceIdentification::DMP_calc_new_theta_stat_net(const DMPNeiMessage& old_messages, const int& from_index, const std::string& to_node, const double& lambda) {
	const DMPMsg2Nei& msg2nei = DMP_get_msg2nei(old_messages, from_index, to_node);
	return msg2nei.to_nei_theta - lambda * msg2nei.to_nei_phi;
}

double SourceIdentification::DMP_get_theta_stat_net(const DMPNeiMessage& messages, const int& from_index, const std::string& to_node) {
	return DMP_get_msg2nei(messages, from_index, to_node).to_nei_theta;
}

void SourceIdentification::DMP_calc_new_psds_stat_net(Network& contact_network, const DMPMsgVecPtr& pss0, DMPNeiMessage& new_messages) {
	int n = contact_network.get_node_size();
	for (int i = 0; i < n; ++i) {
		std::string i_node = contact_network.get_node_name(i);
		NeighborList neighbor_list = contact_network.get_neighbor_list(i_node, 0, 0);
		std::unordered_map<int, double> k2theta_map;
		for (NeighborList::const_iterator kt = neighbor_list.begin(); kt != neighbor_list.end(); ++kt) {
			std::string k_node = kt->get_name();
			int k_index = contact_network.get_node_index(k_node);
			double theta_k2i = DMP_get_theta_stat_net(new_messages, k_index, i_node);
			k2theta_map.insert(std::make_pair(k_index, theta_k2i));
		}
		for (NeighborList::const_iterator jt = neighbor_list.begin(); jt != neighbor_list.end(); ++jt) {
			std::string j_node = jt->get_name();
			int j_index = contact_network.get_node_index(j_node);
			double new_psd = DMP_calc_new_psd_stat_net(j_index, pss0->at(i), k2theta_map);
			DMP_set_msg2nei(new_messages, i, j_node).to_nei_psd = new_psd;
		}
	}
}

double SourceIdentification::DMP_calc_new_psd_stat_net(const int& j_index, const double& pss0, const std::unordered_map<int, double>& k2theta_map) {
	double res = pss0;
	for (std::unordered_map<int, double>::const_iterator it = k2theta_map.begin(); it != k2theta_map.end(); ++it) {
		if (it->first != j_index) {
			res *= it->second;
		}
	}
	return res;
}

void SourceIdentification::DMP_calc_new_phis_stat_net(Network& contact_network, const DMPNeiMessage& old_messages, const Parameter& para, DMPNeiMessage& new_messages) {
	int n = contact_network.get_node_size();
	for (int k = 0; k < n; ++k) {
		std::string k_node = contact_network.get_node_name(k);
		NeighborList neighbor_list = contact_network.get_neighbor_list(k_node, 0, 0);
		for (NeighborList::const_iterator it = neighbor_list.begin(); it != neighbor_list.end(); ++it) {
			std::string i_node = it->get_name();
			double weight = it->get_weight();
			double lambda = DiseaseDynamics::get_infect_prob(weight, para);
			double mu = DiseaseDynamics::get_recover_prob(para);
			double new_phi = DMP_calc_new_phi_stat_net(old_messages, new_messages, k, i_node, lambda, mu);
			DMP_set_msg2nei(new_messages, k, i_node).to_nei_phi = new_phi;
		}
	}
}

double SourceIdentification::DMP_calc_new_phi_stat_net(const DMPNeiMessage& old_messages, const DMPNeiMessage& new_messages, const int& from_index, const std::string& to_node, const double& lambda, const double& mu) {
	const DMPMsg2Nei& old_msg2nei = DMP_get_msg2nei(old_messages, from_index, to_node);
	const DMPMsg2Nei& new_msg2nei = DMP_get_msg2nei(new_messages, from_index, to_node);
	return (1 - lambda) * (1 - mu) * old_msg2nei.to_nei_phi - (new_msg2nei.to_nei_psd - old_msg2nei.to_nei_psd);
}

void SourceIdentification::DMP_calc_new_prs_stat_net(const Network& contact_network, const DMPNeiMessage& old_messages, const Parameter& para, DMPNeiMessage& new_messages) {
	int n = contact_network.get_node_size();
	double mu = DiseaseDynamics::get_recover_prob(para);
	for (int i = 0; i < n; ++i) {
		new_messages.prs->at(i) = DMP_calc_new_pr_stat_net(old_messages.prs->at(i), old_messages.pis->at(i), mu);
	}
}

double SourceIdentification::DMP_calc_new_pr_stat_net(const double& old_pr, const double& old_pi, const double& mu) {
	return old_pr + mu * old_pi;
}

void SourceIdentification::DMP_calc_new_pss_stat_net(Network& contact_network, const DMPMsgVecPtr& pss0, DMPNeiMessage& new_messages) {
	int n = contact_network.get_node_size();
	for (int i = 0; i < n; ++i) {
		new_messages.pss->at(i) = DMP_calc_new_ps_stat_net(contact_network, pss0->at(i), new_messages, i);
	}
}

double SourceIdentification::DMP_calc_new_ps_stat_net(Network& contact_network, const double& pss0, const DMPNeiMessage& new_messages, const int& index) {
	std::string node = contact_network.get_node_name(index);
	NeighborList nl = contact_network.get_neighbor_list(node, 0, 0);
	if (nl.size() == 0) {
		return pss0;
	} else {
		std::string j_node = nl.begin()->get_name();
		int j_index = contact_network.get_node_index(j_node);
		return DMP_get_msg2nei(new_messages, index, j_node).to_nei_psd * DMP_get_theta_stat_net(new_messages, j_index, node);
	}
}

void SourceIdentification::DMP_calc_new_pis_stat_net(const Network& contact_network, DMPNeiMessage& new_messages) {
	int n = contact_network.get_node_size();
	for (int i = 0; i < n; ++i) {
		new_messages.pis->at(i) = 1 - new_messages.pss->at(i) - new_messages.prs->at(i);
	}
}

std::string SourceIdentification::calc_source_MCSM_diffsimdays(Network& contact_network, const NodeSet& nodes_been_infected, const int& max_sim_days, const bool& fixed_sim_days, const Parameter& para) {
	NodeVec potential_seeds = Util::unset2vec(nodes_been_infected);
	int potential_seeds_count = potential_seeds.size();
	int basic_repeat_times = potential_seeds_count * log(potential_seeds_count);
/*	if (!fixed_sim_days) {
		basic_repeat_times *= potential_seeds_count;
	}
	*/
	std::vector<std::vector<double>> similarities = MCSM_calc_similarities_from_simulations(contact_network, potential_seeds, basic_repeat_times, max_sim_days, fixed_sim_days, para);
	int a_reci = 1024 * 32;
	int seed_index_map;
	while (a_reci > 1) {
		std::vector<double> likelihoods1, likelihoods2;
		MCSM_calc_likelihoods(similarities, a_reci, likelihoods1, likelihoods2);
		seed_index_map = Util::getMaxIndex(likelihoods2);
		double posterior2 = MCSM_calc_posterior(likelihoods2, seed_index_map);
		double posterior1 = MCSM_calc_posterior(likelihoods1, seed_index_map);
		if (fabs(posterior2 - posterior1) <= 0.05) {
			break;
		}
		a_reci /= 2;
	}
	return potential_seeds[seed_index_map];
}


std::vector<std::vector<double>> SourceIdentification::MCSM_calc_similarities_from_simulations(Network& contact_network, const NodeVec& potential_seeds, const int& basic_repeat_times, const int& sim_days, const bool& fixed_sim_days, const Parameter& para) {
	int n = contact_network.get_node_size();
	int seeds_count = potential_seeds.size();
	std::vector<std::vector<double>> res(seeds_count, std::vector<double>(2 * basic_repeat_times, 0));
	NodeSet standard_infected = Util::vec2unset(potential_seeds);
	Parameter* para_ptr = new Parameter(para);
	for (int i = 0; i < seeds_count; ++i) {
		std::string seed = potential_seeds[i];
		int seed_index = contact_network.get_node_index(seed);
		for (int j = 0; j < 2 * basic_repeat_times; ++j) {
			if (fixed_sim_days == false) {
				para_ptr->set_max_sim_days(1 + Util::gen_rand_int(sim_days));
			}
			Simulator simulator(*para_ptr, n, {seed_index});
			std::vector<DiseaseStage> stages = simulator.get_sim_res(contact_network, *para_ptr);
			NodeSet infected_set = simulator.get_nodeset_been_infected_from_sim_res(stages, para_ptr->get_disease(), contact_network);
			res[i][j] = Util::calcJaccardSimilarity(infected_set, standard_infected);
		}
	}
	return res;
}

void SourceIdentification::MCSM_calc_likelihoods(const std::vector<std::vector<double>>& similarities, const int& a_reci, std::vector<double>& half_res, std::vector<double>& total_res) {
	int seeds, repeat_times;
	Util::getMatrixSize(similarities, seeds, repeat_times);
	half_res.clear();
	total_res.clear();
	half_res.resize(seeds, 0);
	total_res.resize(seeds, 0);
	for (int i = 0; i < seeds; ++i) {
		double sum = 0;
		int half_repeat_times = repeat_times / 2;
		int j = 0;
		for (; j < half_repeat_times; ++j) {
			double temp = a_reci * (1-similarities[i][j]);
			sum += exp(- temp * temp);
		}
		half_res[i] = sum / half_repeat_times;
		for (; j < repeat_times; ++j) {
			double temp = a_reci * (1-similarities[i][j]);
			sum += exp(- temp * temp);
		}
		total_res[i] = sum / repeat_times;
	}
}

double SourceIdentification::MCSM_calc_posterior(const std::vector<double>& likelihoods, const int& seed_index) {
	if (seed_index < 0 || seed_index >= likelihoods.size()) {
		std::cerr << "Error for illegal seed_index = " << seed_index << std::endl;
		exit(-1);
	}
	double sum = 0;
	for (int i = 0; i < likelihoods.size(); ++i) {
		sum += likelihoods[i];
	}
	return likelihoods[seed_index] / sum;
}


