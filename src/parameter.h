#ifndef NETWORKPROJECT_PARAMETER_H_
#define NETWORKPROJECT_PARAMETER_H_

#include <unistd.h>
#include <getopt.h>
#include <string>
#include <cstring>
#include <iostream>

#include "constant.h"

class UtilConstant;

class Parameter {
private:
	bool _help {false};
	NetType _net_type {STATIC};
	std::string _net_inroot {""};
	std::string _net_injson {""}; // only useful when _net_type != STATIC
	std::string _net_volunteers {""}; // when not empty, there will be only these nodes in the network(s)
	std::string _out_dir {""}; // it can also be used to denote the output file
	int _part_str_length {3};
	DiseaseModel _disease {SI};
	double _infect_rate {0.0}; // S->I or S->E
	double _infect_rate_seconds {1.0};
	double _infectious_rate {0.0}; // E->I
	double _infectious_rate_seconds {1.0};
	double _recover_rate {0.0}; // I->R or I->S
	double _recover_rate_seconds {1.0};
	
	double _seconds_per_weight {1.0};
	double _seconds_per_step {1.0};
	int _source_count {1};
	int _max_sim_days {-1};
	double _snapshot_coverage {2}; // greater than 1 means that all nodes are covered
	int _repeat_times {1};
	SrcIdnMethod _source_identification_method {SSEBFS};
	double _ub_r {0.85};
	bool _source_identification_knowntime {false};

	int _start_part {92};
	int _end_part {110};
	int _last_parts_threshold {4};

	// about calculating features
	bool _calc_edges {false};

	// about large component evolve
	std::string _net_friendship {""};
	std::string _fd_func {""};

	double _evolve_para_alpha {0.5};
	double _evolve_para_a {0.5};
	double _evolve_para_b {0.5};

	std::string _if {""};
	std::string _ie {""};

public:
	Parameter();
	Parameter(const Parameter& para);
	Parameter(int argc, char*const* argv, const char* shortopts, const struct option* longopts);
	void check_notnull();

	bool get_help() const { return _help; }
	NetType get_net_type() const { return _net_type; }
	std::string get_net_inroot() const { return _net_inroot; }
	std::string get_net_injson() const { return _net_injson; }
	std::string get_net_volunteers() const { return _net_volunteers; }
	std::string get_out_dir() const { return _out_dir; }
	int get_part_str_length() const { return _part_str_length; }
	DiseaseModel get_disease() const { return _disease; }
	double get_infect_rate() const { return _infect_rate; }
	double get_infect_rate_seconds() const { return _infect_rate_seconds; }
	double get_infectious_rate() const { return _infectious_rate; }
	double get_infectious_rate_seconds() const { return _infectious_rate_seconds; }
	double get_recover_rate() const { return _recover_rate; }
	double get_recover_rate_seconds() const { return _recover_rate_seconds; }
	double get_seconds_per_weight() const { return _seconds_per_weight; }
	double get_seconds_per_step() const { return _seconds_per_step; }
	int get_source_count() const { return _source_count; }
	int get_max_sim_days() const { return _max_sim_days; }
	double get_snapshot_coverage() const { return _snapshot_coverage; }
	int get_repeat_times() const { return _repeat_times; }
	SrcIdnMethod get_source_identification_method() const { return _source_identification_method; }
	double get_ub_r() const { return _ub_r; }
	bool get_source_identification_knowntime() const { return _source_identification_knowntime; }
	int get_start_part() const { return _start_part; }
	int get_end_part() const { return _end_part; }
	int get_last_parts_threshold() const { return _last_parts_threshold; }
	bool get_calc_edges() const { return _calc_edges; }
	std::string get_net_friendship() const { return _net_friendship; }
	std::string get_fd_func() const { return _fd_func; }
	double get_evolve_para_alpha () const { return _evolve_para_alpha; }
	double get_evolve_para_a () const { return _evolve_para_a; }
	double get_evolve_para_b () const { return _evolve_para_b; }
	std::string get_if() const { return _if; }
	std::string get_ie() const { return _ie; }

public:
	/* 
	 * FOR all set functions:
	 *     IF opt_val can be parsed correctly:
	 *	       return 0; 
	 *	   ELSE:
	 *	       return -1;
	 */
	//int set_para(const OptionKey& opt_key, const std::string& opt_val);
	int set_para(OptionKey opt_key, const std::string& opt_val);
	//int set_para(const OptionKey& opt_key);
	int set_para(OptionKey opt_key);

	int set_help();
	int set_net_type(const std::string& val);
	int set_net_inroot(const std::string& val);
	int set_net_injson(const std::string& val);
	int set_net_volunteers(const std::string& val);
	int set_out_dir(const std::string& val);
	int set_part_str_length(const std::string& val);
	int set_disease(const std::string& val);
	int set_infect_rate(const std::string& val);
	int set_infect_rate_seconds(const std::string& val);
	int set_infectious_rate(const std::string& val);
	int set_infectious_rate_seconds(const std::string& val);
	int set_recover_rate(const std::string& val);
	int set_recover_rate_seconds(const std::string& val);
	int set_seconds_per_weight(const std::string& val);
	int set_seconds_per_step(const std::string& val);
	int set_source_count(const std::string& val);
	int set_max_sim_days(const std::string& val);
	int set_snapshot_coverage(const std::string& val);
	int set_repeat_times(const std::string& val);
	int set_source_identification_method(const std::string& val);
	int set_ub_r(const std::string& val);
	int set_source_identification_knowntime(const std::string& val);
	int set_start_part(const std::string& val);
	int set_end_part(const std::string& val);
	int set_last_parts_threshold(const std::string& val);
	int set_calc_edges();
	int set_net_friendship(const std::string& val);
	int set_fd_func(const std::string& val);
	int set_evolve_para_alpha(const std::string& val);
	int set_evolve_para_a(const std::string& val);
	int set_evolve_para_b(const std::string& val);
	int set_if(const std::string& val);
	int set_ie(const std::string& val);

public:
	int set_net_type(const NetType& net_type);
	int set_part_str_length(const int& part_str_length);
	int set_disease(const DiseaseModel& disease);
	int set_infect_rate(const double& infect_rate);
	int set_infect_rate_seconds(const double& infect_rate_seconds);
	int set_infectious_rate(const double& infectious_rate);
	int set_infectious_rate_seconds(const double& infectious_rate_seconds);
	int set_recover_rate(const double& recover_rate);
	int set_recover_rate_seconds(const double& recover_rate_seconds);
	int set_seconds_per_weight(const double& seconds_per_weight);
	int set_seconds_per_step(const double& seconds_per_step);
	int set_source_count(const int& source_count);
	int set_max_sim_days(const int& max_sim_days);
	int set_snapshot_coverage(const double& snapshot_coverage);
	int set_repeat_times(const int& repeat_times);
	int set_source_identification_method(const SrcIdnMethod& method);
	int set_ub_r(const double& ub_r);
	int set_source_identification_knowntime(const bool& knowntime);
	int set_start_part(const int& start_part);
	int set_end_part(const int& end_part);
	int set_last_parts_threshold(const int& last_parts_threshold);
	int set_evolve_para_alpha(const double& evolve_para_alpha);
	int set_evolve_para_a(const double& evolve_para_a);
	int set_evolve_para_b(const double& evolve_para_b);
};


#endif // NETWORKPROJECT_PARAMETER_H_
