#include <cmath>
#include <vector>
#include <list>

#include "constant.h"

class Parameter;

class DiseaseDynamics {
public:
//	static DiseaseStage get_next_status(const Parameter& para, const DiseaseStage& cur_stage, const std::vector<double>& infectious_neighbor_seconds);



public:
	static DiseaseStage get_next_status_byself(const Parameter& para, const DiseaseStage& cur_stage, const double& seconds);
	static DiseaseStage get_next_status_bycontact(const Parameter& para, const DiseaseStage& cur_stage, const std::list<double>& seconds);

private:
	static DiseaseStage si_get_next_status_byself(const Parameter& para, const DiseaseStage& cur_stage, const double& seconds);
	static DiseaseStage si_get_next_status_bycontact(const Parameter& para, const DiseaseStage& cur_stage, const std::list<double>& seconds);
	static DiseaseStage sis_get_next_status_byself(const Parameter& para, const DiseaseStage& cur_stage, const double& seconds);
	static DiseaseStage sis_get_next_status_bycontact(const Parameter& para, const DiseaseStage& cur_stage, const std::list<double>& seconds);
	static DiseaseStage sir_get_next_status_byself(const Parameter& para, const DiseaseStage& cur_stage, const double& seconds);
	static DiseaseStage sir_get_next_status_bycontact(const Parameter& para, const DiseaseStage& cur_stage, const std::list<double>& seconds);
	static DiseaseStage seir_get_next_status_byself(const Parameter& para, const DiseaseStage& cur_stage, const double& seconds); // original SEIR model (uninfectious in E state)
	static DiseaseStage seir_get_next_status_bycontact(const Parameter& para, const DiseaseStage& cur_stage, const std::list<double>& seconds); // original SEIR model (uninfectious in E state)
	static DiseaseStage seir2019_get_next_status_byself(const Parameter& para, const DiseaseStage& cur_stage, const double& seconds); // new SEIR model adapted to 2019-nCov (infectious in E state)
	static DiseaseStage seir2019_get_next_status_bycontact(const Parameter& para, const DiseaseStage& cur_stage, const std::list<double>& seconds); // new SEIR model adapted to 2019-nCov (infectious in E state)


	static double get_prob(const double& para_prob, const double& para_sec, const double& seconds);
	static double get_prob(const double& para_prob, const double& para_sec, const std::list<double>& seconds);

public:
	static double get_infect_prob(const double& weight, const Parameter& para);
	static double get_infect_prob(const Parameter& para);
	static double get_recover_prob(const Parameter& para);
};
