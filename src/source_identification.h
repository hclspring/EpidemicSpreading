#ifndef NETWORKPROJECT_SOURCEIDENTIFICATION_H_
#define NETWORKPROJECT_SOURCEIDENTIFICATION_H_

#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>
#include <list>
#include <memory>
#include "constant.h"

class Network;
class UndirectedGraph;
class BFSTree;
class Parameter;
class NeighborInfo;
class DiseaseDynamics;

struct JceInfo {
public:
	int first_child_path;
	int second_child_path;
	int first_child_index;
	int longest_path;
	int g_value;
};

typedef std::list<NeighborInfo> NeighborList;
typedef std::unordered_set<std::string> NodeSet;
typedef std::vector<std::string> NodeVec;
typedef std::vector<JceInfo> JceVec;
typedef std::vector<double> AdjacencyRow;
typedef std::vector<AdjacencyRow> AdjacencyMatrix;
typedef std::vector<double> DMPMsgVec;
typedef std::shared_ptr<DMPMsgVec> DMPMsgVecPtr;
typedef std::vector<std::vector<double>> DMPMsgMat;
typedef std::shared_ptr<DMPMsgMat> DMPMsgMatPtr;
struct DMPMsg2Nei {
	std::string nei_name;
	double to_nei_theta;
	double to_nei_phi;
	double to_nei_psd;
};
typedef std::list<DMPMsg2Nei> DMPMsgNeiList;
typedef std::vector<DMPMsgNeiList> DMPMsgNeiListVec;
typedef std::shared_ptr<DMPMsgNeiListVec> DMPMsgNeiPtr;

struct DMPNeiMessage {
	DMPMsgNeiPtr neighbors;
	DMPMsgVecPtr pss;
	DMPMsgVecPtr prs;
	DMPMsgVecPtr pis;
};

struct DMPMatMessage {
	DMPMsgMatPtr thetas;
	DMPMsgMatPtr phis;
	DMPMsgMatPtr pre_psds;
	DMPMsgMatPtr psds;
	DMPMsgVecPtr pss;
	DMPMsgVecPtr prs;
	DMPMsgVecPtr pis;
};

class SourceIdentification {
private:
//	UndirectedGraph* _infection_graph;

public:
//	SourceIdentification(const Network& contact_network, const NodeSet& infected_nodes);
//	~SourceIdentification();

	// calculate identification effect on one simulation
	static double calc_error_distance(Network& contact_network, const std::string& real_source, const std::string& inferred_source);
	static double calc_error_distance(Network& contact_network, const NodeSet& real_source, const NodeSet& inferred_source, const double& diff_count_penalty);
	// calculate identification effect on repeated simulations
	static double calc_detection_rate(const std::string& real_source, const std::string& inferred_source);
//	static double calc_detection_rate(const NodeSet& real_source, const NodeSet& inferred_source);

public:
	// calculate sources
	// Total: 
	//		When known_sim_days = false, max_sim_days in the identification method will be (10 * para.get_max_sim_days());
	//		For fair, JCE & RG will know all people's stages.
	static std::string calc_source(Network& contact_network, const std::vector<DiseaseStage>& sim_res, const Parameter& para, const SrcIdnMethod& method, const bool& known_sim_days);
	static std::string calc_source_with_infection_graph(Network& contact_network, const NodeSet& nodes_been_infected, const Parameter& para, const SrcIdnMethod& method);
	//static NodeSet calc_source(const Network& contact_network, const NodeSet& nodes_been_infected, const int& max_source_count, const SrcIdnMethod& method);

	// SSE: tree
	static std::string calc_source_SSE(Network& contact_network, const NodeSet& nodes_been_infected);
	// SSEBFS: 
	static std::string calc_source_SSEBFS(Network& contact_network, const NodeSet& nodes_been_infected);
	// SJC: tree
	static std::string calc_source_SJC(Network& contact_network, const NodeSet& nodes_being_infected);
	// JCE: tree
	static std::string calc_source_JCE(Network& contact_network, const NodeSet& nodes_observed_infected);
	// RG:
	static std::string calc_source_RG(Network& contact_network, const NodeSet& nodes_observed_infected, const Parameter& para);
	// DA:
	static std::string calc_source_DA(Network& contact_network, const NodeSet& nodes_been_infected);
	static NodeSet calc_source_DA(Network& contact_network, const NodeSet& nodes_been_infected, const int& source_count);
	// UB:
	static std::string calc_source_UB(Network& contact_network, const NodeSet& nodes_been_infected, const double& ub_r);
	// DMP: SIR, STATIC
	static std::string calc_source_DMP_SIR_unknowntime(Network& contact_network, const std::vector<DiseaseStage>& sim_res, const int& max_sim_days, const Parameter& para);
	static std::string calc_source_DMP_SIR_knowntime(Network& contact_network, const std::vector<DiseaseStage>& sim_res, const int& days, const Parameter& para);
	// DMP: STATIC
//	static std::string calc_source_DMP(Network& contact_network, const NodeSet& nodes_been_infected, const int& max_sim_days, const Parameter& para);
//	static std::string calc_source_DMP(Network& contact_network, const NodeSet& nodes_been_infected, const int& days, const Parameter& para);
	// MCSM
	static std::string calc_source_MCSM_unknowntime(Network& contact_network, const NodeSet& nodes_been_infected, const int& max_sim_days, const Parameter& para);
	static std::string calc_source_MCSM_knowntime(Network& contact_network, const NodeSet& nodes_been_infected, const int& sim_days, const Parameter& para);
	// NetSleuth
	static std::string calc_source_NetSleuth(Network& contact_network, const NodeSet& nodes_been_infected);

private:
	// SSE
	static long SSE_calc_permitted_permutation_count(const std::shared_ptr<UndirectedGraph>& graph, const NodeVec& bfs_nodes);
	static std::vector<double> SSE_calc_rumor_centralities(const std::shared_ptr<BFSTree>& tree);
	static void SSE_calc_t_p(const std::shared_ptr<BFSTree>& tree, std::vector<int>& t, std::vector<int>& p, const int& cur_index);
	static void SSE_calc_r(const std::shared_ptr<BFSTree>& tree, const std::vector<int>& t, const std::vector<int>& p, std::vector<double>& r, const int& cur_index);
	// SJC
	static bool SJC_check_finished(const std::vector<int>& received_message_count, const int& target);
	// JCE
	static void JCE_upward_message_passing(const std::shared_ptr<BFSTree>& tree, const int& root_index, std::vector<JceInfo>& jce_info);
	static std::string JCE_downward_message_passing(const std::shared_ptr<BFSTree>& tree, const int& root_index, std::vector<JceInfo>& jce_info);
	// RG
	static int RG_for_every_tree(std::shared_ptr<BFSTree>& tree, std::shared_ptr<UndirectedGraph>& induced_graph);
	static std::vector<int> RG_indices_with_same_depth_and_increasing_height(const std::vector<int>& depths, const int& d, const std::vector<int>& heights);
	static int RG_height_with_exception(const std::shared_ptr<BFSTree>& tree, const std::string& y_name, const std::string& x_name, const std::vector<int>& heights);
	static int RG_height_with_exception(const std::shared_ptr<BFSTree>& tree, const std::string& y_name, const std::string& x_name);
	// DMP
	// 2 main functions:
	// DMP_calc_partition_function: used when simulation time is unknown.
	//			return a vector of partition, i-th value denotes partition function after simulating (i+1) days.
	//			parameter "likelihoods" is a matrix, with rows denoting diff seeds, and columns denotes diff simulating days.
	static std::vector<double> DMP_calc_partition_function(Network& contact_network, const int& max_sim_days, const std::vector<DiseaseStage>& res_stages, const Parameter& para, std::vector<std::vector<double>>& likelihoods);
	// DMP_calc_partition_function_with_same_sim_days: used when simulation time is known.
	//			parameter "likelihoods" is a vector, i-th value denotes likelihood when i-th node is the seed.
	static double DMP_calc_partition_function_with_same_sim_days(Network& contact_network, const int& sim_days, const std::vector<DiseaseStage>& res_stages, const Parameter& para, std::vector<double>& likelihoods);
	// DMP likelihoods:
	static std::vector<double> DMP_calc_likelihoods_with_same_seed(Network& contact_network, const std::string& seed_node, const int& max_sim_days, const std::vector<DiseaseStage>& res_stages, const Parameter& para); // This function returns a vector with (max_sim_days) likelihoods.
	static double DMP_calc_likelihood_with_fixed_seed_and_simdays(Network& contact_network, const std::string& seed_node, const int& sim_days, const std::vector<DiseaseStage>& res_stages, const Parameter& para);
	// DMP message initialization and passing
//	static void DMP_pass_message_day(Network& contact_network, const DMPMessage& pre_message, const int& day, const DMPOneMessage& pss0, DMPMessage& next_message);
	static DMPNeiMessage DMP_init_message_step_stat_net(Network& contact_network, const std::string& seed_node);
	static DMPMsgVecPtr DMP_init_pss0_step_stat_net(Network& contact_network, const std::string& seed_node);
	static DMPNeiMessage DMP_pass_message_step_stat_net(Network& contact_network, const DMPNeiMessage& pre_messages, const DMPMsgVecPtr& pss0, const Parameter& para);
	static double DMP_calc_likelihood(const std::vector<DiseaseStage>& stages, const DMPNeiMessage& messages);
	// DMP initialization in every step
	static void DMP_init_mat(const int& n, DMPMsgMatPtr& mp);
	static void DMP_init_vec(const int& n, DMPMsgVecPtr& vp);
	static DMPMatMessage DMP_init_mat_msg(const int& n, const DMPMatMessage& pre_msg);
	static DMPNeiMessage DMP_init_nei_msg(const int& n);
	// DMP get/set neighbor message
	static const DMPMsg2Nei& DMP_get_msg2nei(const DMPNeiMessage& messages, const int& from_index, const std::string& to_node);
	static DMPMsg2Nei& DMP_set_msg2nei(DMPNeiMessage& messages, const int& from_index, const std::string& to_node);
	// DMP calc theta in every step
	static void DMP_calc_new_thetas_stat_net(Network& contact_network, const DMPNeiMessage& old_messages, const Parameter& para, DMPNeiMessage& new_messages);
	static double DMP_calc_new_theta_stat_net(const DMPNeiMessage& old_messages, const int& from_index, const std::string& to_node, const double& lambda);
	static double DMP_get_theta_stat_net(const DMPNeiMessage& messages, const int& from_index, const std::string& to_node);
	// DMP calc psd in every step
	static void DMP_calc_new_psds_stat_net(Network& contact_network, const DMPMsgVecPtr& pss0, DMPNeiMessage& new_messages);
	static double DMP_calc_new_psd_stat_net(const int& j_index, const double& pss0, const std::unordered_map<int, double>& k2theta_map);
	// DMP calc phi in every step
	static void DMP_calc_new_phis_stat_net(Network& contact_network, const DMPNeiMessage& old_messages, const Parameter& para, DMPNeiMessage& new_messages);
	static double DMP_calc_new_phi_stat_net(const DMPNeiMessage& old_messages, const DMPNeiMessage& new_messages, const int& from_index, const std::string& to_node, const double& lambda, const double& mu);
	// DMP calc prs in every step
	static void DMP_calc_new_prs_stat_net(const Network& contact_network, const DMPNeiMessage& old_messages, const Parameter& para, DMPNeiMessage& new_messages);
	static double DMP_calc_new_pr_stat_net(const double& old_pr, const double& old_pi, const double& mu);
	// DMP calc pss in every step
	static void DMP_calc_new_pss_stat_net(Network& contact_network, const DMPMsgVecPtr& pss0, DMPNeiMessage& new_messages);
	static double DMP_calc_new_ps_stat_net(Network& contact_network, const double& pss0, const DMPNeiMessage& new_messages, const int& index);
	// DMP calc psi in every step
	static void DMP_calc_new_pis_stat_net(const Network& contact_network, DMPNeiMessage& new_messages);
	// MCSM
	static std::string calc_source_MCSM_diffsimdays(Network& contact_network, const NodeSet& nodes_been_infected, const int& max_sim_days, const bool& fixed_sim_days, const Parameter& para);
	static std::vector<std::vector<double>> MCSM_calc_similarities_from_simulations(Network& contact_network, const NodeVec& potential_seeds, const int& basic_repeat_times, const int& sim_days, const bool& fixed_sim_days, const Parameter& para);
	static void MCSM_calc_likelihoods(const std::vector<std::vector<double>>& similarities, const int& a_reci, std::vector<double>& half_likelihoods, std::vector<double>& total_likelihoods);
	static double MCSM_calc_posterior(const std::vector<double>& likelihoods, const int& seed_index);
};

#endif // NETWORKPROJECT_SOURCEIDENTIFICATION_H_
