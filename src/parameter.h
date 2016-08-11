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
	NetType _net_type {STATIC};
	std::string _net_inroot {""};
	std::string _net_injson {""}; // only useful when _net_type != STATIC
	std::string _net_volunteers {""}; // when not empty, there will be only these nodes in the network(s)
	std::string _out_dir {""}; // it can also be used to denote the output file
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

public:
	Parameter();
	Parameter(const Parameter& para);
	Parameter(int argc, char*const* argv, const char* shortopts, const struct option* longopts);
	void check_notnull();

	NetType get_net_type() const { return _net_type; }
	std::string get_net_inroot() const { return _net_inroot; }
	std::string get_net_injson() const { return _net_injson; }
	std::string get_net_volunteers() const { return _net_volunteers; }
	std::string get_out_dir() const { return _out_dir; }
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

private:
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
	
	int set_net_type(const std::string& val);
	int set_net_inroot(const std::string& val);
	int set_net_injson(const std::string& val);
	int set_net_volunteers(const std::string& val);
	int set_out_dir(const std::string& val);
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

public:
	int set_net_type(const NetType& net_type);
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
};


#endif // NETWORKPROJECT_PARAMETER_H_
