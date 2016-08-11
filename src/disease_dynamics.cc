#include "disease_dynamics.h"

#include <cstdlib>
#include <climits>
#include "parameter.h"
#include "util.h"

DiseaseStage DiseaseDynamics::get_next_status_byself(const Parameter& para, const DiseaseStage& cur_stage, const double& seconds) {
	switch (para.get_disease()) {
		case SI: return si_get_next_status_byself(para, cur_stage, seconds);
		case SIS: return sis_get_next_status_byself(para, cur_stage, seconds);
		case SIR: return sir_get_next_status_byself(para, cur_stage, seconds);
		default: std::cerr << "Error with undefined disease model." << std::endl; exit(-1);
	}
	return cur_stage;
}

DiseaseStage DiseaseDynamics::get_next_status_bycontact(const Parameter& para, const DiseaseStage& cur_stage, const std::list<double>& seconds) {
	switch (para.get_disease()) {
		case SI: return si_get_next_status_bycontact(para, cur_stage, seconds);
		case SIS: return sis_get_next_status_bycontact(para, cur_stage, seconds);
		case SIR: return sir_get_next_status_bycontact(para, cur_stage, seconds);
		default: std::cerr << "Error with undefined disease model." << std::endl; exit(-1);
	}
	return cur_stage;
}

DiseaseStage DiseaseDynamics::si_get_next_status_byself(const Parameter& para, const DiseaseStage& cur_stage, const double& seconds) {
	return cur_stage;
}

DiseaseStage DiseaseDynamics::si_get_next_status_bycontact(const Parameter& para, const DiseaseStage& cur_stage, const std::list<double>& seconds) {
	switch (cur_stage) {
		case SUSCEPTIBLE: 
			if (Util::gen_rand_double() < get_prob(para.get_infect_rate(), para.get_infect_rate_seconds(), seconds)) {
				return INFECTIOUS;
			} else {
				return SUSCEPTIBLE;
			}
		default: return cur_stage;
	}
}

DiseaseStage DiseaseDynamics::sis_get_next_status_byself(const Parameter& para, const DiseaseStage& cur_stage, const double& seconds) {
	switch (cur_stage) {
		case INFECTIOUS:
			if (Util::gen_rand_double() < get_prob(para.get_recover_rate(), para.get_recover_rate_seconds(), seconds)) {
				return SUSCEPTIBLE;
			} else {
				return INFECTIOUS;
			}
		default: return cur_stage;
	}
}

DiseaseStage DiseaseDynamics::sis_get_next_status_bycontact(const Parameter& para, const DiseaseStage& cur_stage, const std::list<double>& seconds) {
	switch (cur_stage) {
		case SUSCEPTIBLE:
			if (Util::gen_rand_double() < get_prob(para.get_infect_rate(), para.get_infect_rate_seconds(), seconds)) {
				return INFECTIOUS;
			} else {
				return SUSCEPTIBLE;
			}
		default: return cur_stage;
	}
}

DiseaseStage DiseaseDynamics::sir_get_next_status_byself(const Parameter& para, const DiseaseStage& cur_stage, const double& seconds) {
	switch (cur_stage) {
		case INFECTIOUS:
			if (Util::gen_rand_double() < get_prob(para.get_recover_rate(), para.get_recover_rate_seconds(), seconds)) {
				return RECOVERED;
			} else {
				return INFECTIOUS;
			}
		default: return cur_stage;
	}
}

DiseaseStage DiseaseDynamics::sir_get_next_status_bycontact(const Parameter& para, const DiseaseStage& cur_stage, const std::list<double>& seconds) {
	switch (cur_stage) {
		case SUSCEPTIBLE:
			if (Util::gen_rand_double() < get_prob(para.get_infect_rate(), para.get_infect_rate_seconds(), seconds)) {
				return INFECTIOUS;
			} else {
				return SUSCEPTIBLE;
			}
		default: return cur_stage;
	}
}


double DiseaseDynamics::get_prob(const double& para_prob, const double& para_sec, const double& seconds) {
	return 1 - pow(1 - para_prob, seconds / para_sec);
}
	
double DiseaseDynamics::get_prob(const double& para_prob, const double& para_sec, const std::list<double>& seconds) {
	double res = 1.0;
	std::for_each(seconds.begin(), seconds.end(),
			[&para_prob, &para_sec, &res] (const double& x) {
				res = res * (1 - get_prob(para_prob, para_sec, x));
			}
			);
	return 1 - res;
}

double DiseaseDynamics::get_infect_prob(const double& weight, const Parameter& para) {
	return get_prob(para.get_infect_rate(), para.get_infect_rate_seconds(), para.get_seconds_per_weight() * weight);
}

double DiseaseDynamics::get_infect_prob(const Parameter& para) {
	return get_prob(para.get_infect_rate(), para.get_infect_rate_seconds(), para.get_seconds_per_step());
}

double DiseaseDynamics::get_recover_prob(const Parameter& para) {
	return get_prob(para.get_recover_rate(), para.get_recover_rate_seconds(), para.get_seconds_per_step());
}

