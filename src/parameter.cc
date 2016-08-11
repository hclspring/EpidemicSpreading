#include "parameter.h"

#include <string>
#include <boost/lexical_cast.hpp>
#include "util.h"
#include "util_constant.h"
#include "constant.h"

Parameter::Parameter() {}

Parameter::Parameter(const Parameter& para) {
	this->_net_type = para._net_type;
	this->_net_inroot = para._net_inroot;
	this->_net_injson = para._net_injson;
	this->_net_volunteers = para._net_volunteers;
	this->_out_dir = para._out_dir;
	this->_disease = para._disease;
	this->_infect_rate = para._infect_rate;
	this->_infect_rate_seconds = para._infect_rate_seconds;
	this->_infectious_rate = para._infectious_rate;
	this->_infectious_rate_seconds = para._infectious_rate_seconds;
	this->_recover_rate = para._recover_rate;
	this->_recover_rate_seconds = para._recover_rate_seconds;
	this->_seconds_per_weight = para._seconds_per_weight;
	this->_seconds_per_step = para._seconds_per_step;
	this->_source_count = para._source_count;
	this->_max_sim_days = para._max_sim_days;
	this->_snapshot_coverage = para._snapshot_coverage;
	this->_repeat_times = para._repeat_times;
	this->_source_identification_method = para._source_identification_method;
	this->_ub_r = para._ub_r;
	this->_source_identification_knowntime = para._source_identification_knowntime;
}

Parameter::Parameter(int argc, char*const* argv, 
		const char* shortopts, const struct option* longopts) {
	int opt_key;
	while ((opt_key = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1) {
		if (optarg != NULL) {
			set_para(UtilConstant::int2OptionKey(opt_key), std::string(optarg));
		} else {
			set_para(UtilConstant::int2OptionKey(opt_key));
		}
	}
}

int Parameter::set_para(OptionKey opt_key, const std::string& opt_val) {
	switch (opt_key) {
		case OPT_NET_TYPE:					return set_net_type(opt_val);
		case OPT_NET_INROOT:				return set_net_inroot(opt_val);
		case OPT_NET_INJSON:				return set_net_injson(opt_val);
		case OPT_NET_VOLUNTEERS:			return set_net_volunteers(opt_val);
		case OPT_OUT_DIR:					return set_out_dir(opt_val);
		case OPT_DISEASE:					return set_disease(opt_val);
		case OPT_INFECT_RATE:				return set_infect_rate(opt_val);
		case OPT_INFECT_RATE_SECONDS:		return set_infect_rate_seconds(opt_val);
		case OPT_INFECTIOUS_RATE:			return set_infectious_rate(opt_val);
		case OPT_INFECTIOUS_RATE_SECONDS:	return set_infectious_rate(opt_val);
		case OPT_RECOVER_RATE:				return set_recover_rate(opt_val);
		case OPT_RECOVER_RATE_SECONDS:		return set_recover_rate_seconds(opt_val);
		case OPT_SECONDS_PER_WEIGHT:		return set_seconds_per_weight(opt_val);
		case OPT_SECONDS_PER_STEP:			return set_seconds_per_step(opt_val);
		case OPT_SOURCE_COUNT:				return set_source_count(opt_val);
		case OPT_MAX_SIM_DAYS:				return set_max_sim_days(opt_val);
		case OPT_SNAPSHOT_COVERAGE:			return set_snapshot_coverage(opt_val);
		case OPT_REPEAT_TIMES:				return set_repeat_times(opt_val);
		case OPT_SRC_IDN_METHOD:			return set_source_identification_method(opt_val);
		case OPT_UB_R:						return set_ub_r(opt_val);
		case OPT_SRC_IDN_KNOWNTIME:			return set_source_identification_knowntime(opt_val);
		default:							return -1;
	}
}

int Parameter::set_para(OptionKey opt_key) {
	switch (opt_key) {
		default:	return -1;
	}
}

int Parameter::set_net_type(const std::string& val) {
	_net_type = UtilConstant::toNetType(val);
	return 0;
}
int Parameter::set_net_inroot(const std::string& val) {
	_net_inroot = val;
	return 0;
}
int Parameter::set_net_injson(const std::string& val) {
	_net_injson = val;
	return 0;
}
int Parameter::set_net_volunteers(const std::string& val) {
	_net_volunteers = val;
	return 0;
}
int Parameter::set_out_dir(const std::string& val) {
	_out_dir = val;
	return 0;
}
int Parameter::set_disease(const std::string& val) {
	_disease = UtilConstant::toDiseaseModel(val);
	return 0;
}
int Parameter::set_infect_rate(const std::string& val) {
	_infect_rate = boost::lexical_cast<double>(val);
	return 0;
}
int Parameter::set_infect_rate_seconds(const std::string& val) {
	_infect_rate_seconds = boost::lexical_cast<double>(val);
	return 0;
}
int Parameter::set_infectious_rate(const std::string& val) {
	_infectious_rate = boost::lexical_cast<double>(val);
	return 0;
}
int Parameter::set_infectious_rate_seconds(const std::string& val) {
	_infectious_rate_seconds = boost::lexical_cast<double>(val);
	return 0;
}
int Parameter::set_recover_rate(const std::string& val) {
	_recover_rate = boost::lexical_cast<double>(val);
	return 0;
}
int Parameter::set_recover_rate_seconds(const std::string& val) {
	_recover_rate_seconds = boost::lexical_cast<double>(val);
	return 0;
}
int Parameter::set_seconds_per_weight(const std::string& val) {
	_seconds_per_weight = boost::lexical_cast<double>(val);
	return 0;
}
int Parameter::set_seconds_per_step(const std::string& val) {
	_seconds_per_step = boost::lexical_cast<double>(val);
	return 0;
}
int Parameter::set_source_count(const std::string& val) {
	_source_count = boost::lexical_cast<int>(val);
	return 0;
}
int Parameter::set_max_sim_days(const std::string& val) {
	_max_sim_days = boost::lexical_cast<int>(val);
	return 0;
}
int Parameter::set_snapshot_coverage(const std::string& val) {
	_snapshot_coverage = boost::lexical_cast<double>(val);
	return 0;
}
int Parameter::set_repeat_times(const std::string& val) {
	_repeat_times = boost::lexical_cast<int>(val);
	return 0;
}
int Parameter::set_source_identification_method(const std::string& val) {
	_source_identification_method = UtilConstant::toSrcIdnMethod(val);
	return 0;
}
int Parameter::set_ub_r(const std::string& val) {
	_ub_r = boost::lexical_cast<int>(val);
	return 0;
}
int Parameter::set_source_identification_knowntime(const std::string& val) {
	_source_identification_knowntime = UtilConstant::toBool(val);
	return 0;
}

int Parameter::set_net_type(const NetType& net_type) {
	this->_net_type = net_type;
	return 0;
}
int Parameter::set_disease(const DiseaseModel& disease) {
	this->_disease = disease;
	return 0;
}
int Parameter::set_infect_rate(const double& infect_rate) {
	this->_infect_rate = infect_rate;
	return 0;
}
int Parameter::set_infect_rate_seconds(const double& infect_rate_seconds) {
	this->_infect_rate_seconds = infect_rate_seconds;
	return 0;
}
int Parameter::set_infectious_rate(const double& infectious_rate) {
	this->_infectious_rate = infectious_rate;
	return 0;
}
int Parameter::set_infectious_rate_seconds(const double& infectious_rate_seconds) {
	this->_infectious_rate_seconds = infectious_rate_seconds;
	return 0;
}
int Parameter::set_recover_rate(const double& recover_rate) {
	this->_recover_rate = recover_rate;
	return 0;
}
int Parameter::set_recover_rate_seconds(const double& recover_rate_seconds) {
	this->_recover_rate_seconds = recover_rate_seconds;
	return 0;
}
int Parameter::set_seconds_per_weight(const double& seconds_per_weight) {
	this->_seconds_per_weight = seconds_per_weight;
	return 0;
}
int Parameter::set_seconds_per_step(const double& seconds_per_step) {
	this->_seconds_per_step = seconds_per_step;
	return 0;
}
int Parameter::set_source_count(const int& source_count) {
	this->_source_count = source_count;
	return 0;
}
int Parameter::set_max_sim_days(const int& max_sim_days) {
	this->_max_sim_days = max_sim_days;
	return 0;
}
int Parameter::set_snapshot_coverage(const double& snapshot_coverage) {
	this->_snapshot_coverage = snapshot_coverage;
	return 0;
}
int Parameter::set_repeat_times(const int& repeat_times) {
	this->_repeat_times = repeat_times;
	return 0;
}
int Parameter::set_source_identification_method(const SrcIdnMethod& method) {
	this->_source_identification_method = method;
	return 0;
}
int Parameter::set_ub_r(const double& ub_r) {
	this->_ub_r = ub_r;
	return 0;
}
int Parameter::set_source_identification_knowntime(const bool& knowntime) {
	this->_source_identification_knowntime = knowntime;
	return 0;
}


