#include "util_constant.h"

#include <iostream>
#include <boost/algorithm/string/case_conv.hpp>
#include "util.h"

NetType UtilConstant::toNetType(const std::string& val) {
	std::string lowercase = boost::to_lower_copy(val);
	if (lowercase.compare("static") == 0) {
		return STATIC;
	} else if (lowercase.compare("dynamic") == 0) {
		return DYNAMIC;
	} else if (lowercase.compare("pearl") == 0) {
		return PEARL;
	} else {
		std::cerr << "Error parsing string: " << val << std::endl;
		exit(-1);
	}
}

DiseaseModel UtilConstant::toDiseaseModel(const std::string& val) {
	std::string lowercase = boost::to_lower_copy(val);
	if (lowercase.compare("si") == 0) {
		return SI;
	} else if (lowercase.compare("sis") == 0) {
		return SIS;
	} else if (lowercase.compare("sir") == 0) {
		return SIR;
	} else if (lowercase.compare("seir") == 0) {
		return SEIR;
	} else if (lowercase.compare("seir2019") == 0) {
		return SEIR2019;
	} else {
		std::cerr << "Error parsing string: " << val << std::endl;
		exit(-1);
	}
}

SrcIdnMethod UtilConstant::toSrcIdnMethod(const std::string& val) {
	std::string lowercase = boost::to_lower_copy(val);
	if (lowercase.compare("sse") == 0) {
		return SSE;
	} else if (lowercase.compare("ssebfs") == 0) {
		return SSEBFS;
	} else if (lowercase.compare("tse") == 0) {
		return TSE;
	} else if (lowercase.compare("msep") == 0) {
		return MSEP;
	} else if (lowercase.compare("msepbfs") == 0) {
		return MSEPBFS;
	} else if (lowercase.compare("urcc") == 0) {
		return URCC;
	} else if (lowercase.compare("sjc") == 0) {
		return SJC;
	} else if (lowercase.compare("jce") == 0) {
		return JCE;
	} else if (lowercase.compare("mjc") == 0) {
		return MJC;
	} else if (lowercase.compare("rg") == 0) {
		return RG;
	} else if (lowercase.compare("da") == 0) {
		return DA;
	} else if (lowercase.compare("ub") == 0) {
		return UB;
	} else if (lowercase.compare("aub") == 0) {
		return AUB;
	} else if (lowercase.compare("dmp") == 0) {
		return DMP;
	} else if (lowercase.compare("bp") == 0) {
		return BP;
	} else if (lowercase.compare("mcsm") == 0) {
		return MCSM;
	} else if (lowercase.compare("sleuth") == 0 || lowercase.compare("netsleuth") == 0) {
		return SLEUTH;
	} else {
		std::cerr << "Error parsing string: " << val << std::endl;
		exit(-1);
	}
}

bool UtilConstant::toBool(const std::string& val) {
	std::string lowercase = boost::to_lower_copy(val);
	if (lowercase.compare("true") == 0) {
		return true;
	} else {
		return false;
	}
}

OptionKey UtilConstant::int2OptionKey(const int& val) {
	std::vector<OptionKey> vec = {
		OPT_HELP, // 'h'
		OPT_NET_TYPE,
		OPT_NET_INROOT,
		OPT_NET_INJSON,
		OPT_NET_VOLUNTEERS,
		OPT_OUT_DIR,
		OPT_PART_STR_LENGTH,
		OPT_DISEASE,
		OPT_INFECT_RATE,
		OPT_INFECT_RATE_SECONDS,
		OPT_INFECTIOUS_RATE,
		OPT_INFECTIOUS_RATE_SECONDS,
		OPT_RECOVER_RATE,
		OPT_RECOVER_RATE_SECONDS,
		OPT_SECONDS_PER_WEIGHT,
		OPT_SECONDS_PER_STEP,
		OPT_SOURCE_COUNT,
		OPT_SNAPSHOT_COVERAGE,
		OPT_MAX_SIM_DAYS,
		OPT_REPEAT_TIMES,
		OPT_SRC_IDN_METHOD,
		OPT_UB_R,
		OPT_SRC_IDN_KNOWNTIME,
		OPT_START_PART,
		OPT_END_PART,
		OPT_LAST_PARTS_THRESHOLD,
		OPT_CALC_EDGES,
		OPT_NET_FRIENDSHIP,
		OPT_FD_FUNC,
		OPT_EVOLVE_PARA_ALPHA,
		OPT_EVOLVE_PARA_A,
		OPT_EVOLVE_PARA_B,
		OPT_IF,
		OPT_IE,
		OPT_MERGE_PARTS
	};
	if (val >= 'h' && val < 'h' + vec.size()) {
		return vec[val - 'h'];
	} else {
		std::cerr << "Error: unknown OptionKey." << std::endl;
		exit(-1);
	}
}

DiseaseStage UtilConstant::int2DiseaseStage(const int& val) {
	switch (val) {
		case 0: return SUSCEPTIBLE;
		case 1: return EXPOSED;
		case 2: return INFECTIOUS;
		case 3: return RECOVERED;
		default: std::cerr << "Error: unknown DiseaseStage." << std::endl; exit(-1);
	}
}

bool UtilConstant::canBeInfected(const DiseaseStage& stage) {
	return stage == SUSCEPTIBLE;
}

bool UtilConstant::isInfected(const DiseaseModel& disease, const DiseaseStage& stage) {
	switch (stage) {
		case INFECTIOUS: return true;
		case EXPOSED:
			 switch (disease) {
				 case SI:
				 case SIR:
				 case SIS:
					 return false;
				 case SEIR:
				 case SEIR2019:
					 return true;
				 default:
					 std::cerr << "Error: unknown DiseaseModel." << std::endl;
					 exit(-1);
			 }
		default: return false;
	}
}

bool UtilConstant::isInfectious(const DiseaseModel& disease, const DiseaseStage& stage) {
	if (stage == INFECTIOUS) return true;
	else if (stage == EXPOSED && disease == SEIR2019) return true;
	else return false;
}

bool UtilConstant::hasBeenInfected(const DiseaseModel& disease, const DiseaseStage& stage) {
	switch (disease) {
		case SI: if (stage == INFECTIOUS) return true; else return false;
		case SIR: if (stage == INFECTIOUS || stage == RECOVERED) return true; else return false;
		case SIS: if (stage == INFECTIOUS) return true; else return false;
		case SEIR:
		case SEIR2019: if (stage == EXPOSED || stage == INFECTIOUS || stage == RECOVERED) {
					   return true;
				   } else {
					   return false;
				   }
		default:
			std::cerr << "Error: unknown DiseaseModel." << std::endl;
			exit(-1);
	}
}

std::vector<DiseaseStage> UtilConstant::getUnstableStages(const DiseaseModel& disease) {
	std::vector<DiseaseStage> res;
	switch (disease) {
		case SI:
		case SIS:
		case SIR:
			res.push_back(INFECTIOUS); break;
		case SEIR:
		case SEIR2019:
			res.push_back(INFECTIOUS); res.push_back(EXPOSED); break;
		default:
			std::cerr << "Error: unknown DiseaseModel." << std::endl;
			exit(-1);
	}
	return res;
}

bool UtilConstant::isStable(
		const DiseaseModel& disease, 
		const std::vector<DiseaseStage>& stages, 
		const std::vector<DiseaseStage>& unstable_stages) {
	for (int i = 0; i < stages.size(); ++i) {
		if (Util::contains(unstable_stages, stages[i])) {
			return false;
		}
	}
	return true;
}

DiseaseStage UtilConstant::getInitialInfectedStage(const DiseaseModel& disease) {
	switch (disease) {
		case SI:
		case SIS:
		case SIR:
			return INFECTIOUS;
		case SEIR:
		case SEIR2019:
			return EXPOSED;
		default:
			std::cerr << "Error: unknown disease." << std::endl;
			exit(-1);
	}
}
