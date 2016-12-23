#ifndef NETWORKPROJECT_CONSTANT_H_
#define NETWORKPROJECT_CONSTANT_H_

#include <cstddef>

enum NetType {
	STATIC,
	DYNAMIC,
	PEARL,
};

enum DiseaseModel {
	SI,
	SIS,
	SIR,
	SEIR,
};

const size_t OPTION_KEY_MAXSIZE = 100;
/*
 * NOTE: when adding OptionKey values, you must add lines in the function UtilConstant::int2OptionKey(const int&) in file "util_constant.cc".
 */
enum OptionKey {
	OPT_HELP = 'h',
//	OPT_NET_DIM = 0, // dimension of the input network(s), possible values are: 0 (single), 1 (dynamic), 2 (like pearl)
	OPT_NET_TYPE, // type of the input network, possible values: STATIC, DYNAMIC, PEARL
	OPT_NET_INROOT, // input of the network: filename for common-type, and root-directory for pearl-type
	OPT_NET_INJSON, // json file of the input network(s), which records directories and files of the networks in detail
	OPT_NET_VOLUNTEERS, // a file contains all nodes of the network(s). when this option exists, the network(s) contain(s) exactly these volunteers as nodes.
	OPT_OUT_DIR, // output of the results: dir or filename
	OPT_PART_STR_LENGTH, //  length of the part string (usually be 3 for ordinary pearl networks)
	OPT_DISEASE, // disease model of the simulation, possible values: SI, SIR, SIS
	OPT_INFECT_RATE, // beta parameter of the transmission (S->I), double
	OPT_INFECT_RATE_SECONDS, // time duration for beta, double
	OPT_INFECTIOUS_RATE, // delta parameter of the transmission (E->I), double
	OPT_INFECTIOUS_RATE_SECONDS, // time duration for delta, double
	OPT_RECOVER_RATE, // gamma parameter of the transmission (I->R), double
	OPT_RECOVER_RATE_SECONDS, // time duration for gamma, double
	OPT_SECONDS_PER_WEIGHT, // how many seconds does one weight in the network mean
	OPT_SECONDS_PER_STEP, // how many seconds does one step (or part) in the network mean
	OPT_SOURCE_COUNT, // number of infection sources, int
	OPT_SNAPSHOT_COVERAGE, // coverage of snapshot, double; note that when it is greater than 1, it means the snapshot contains all nodes.
	OPT_MAX_SIM_DAYS, // number of maximal simulating days, -1 for infty
	OPT_REPEAT_TIMES,
	OPT_SRC_IDN_METHOD,
	OPT_UB_R, // parameter r of unbiased betweenness
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
};


const size_t DISEASE_STAGE_SIZE = 4;
enum DiseaseStage {
	SUSCEPTIBLE = 0,
	EXPOSED,
	INFECTIOUS,
	RECOVERED,
};

enum SrcIdnMethod {
	SSE = 0,
	SSEBFS,
	TSE,
	MSEP,
	MSEPBFS,
	URCC,
	SJC,
	JCE,
	MJC,
	RG,
	DA,
	UB,
	AUB,
	DMP,
	BP,
	MCSM,
	SLEUTH,
};

#endif // NETWORKPROJECT_CONSTANT_H_
