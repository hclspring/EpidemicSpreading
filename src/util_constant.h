#ifndef NETWORKPROJECT_UTILCONSTANT_H_
#define NETWORKPROJECT_UTILCONSTANT_H_

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>

#include "constant.h"

class UtilConstant {
public:
	static NetType toNetType(const std::string&);
	static DiseaseModel toDiseaseModel(const std::string&);
	static SrcIdnMethod toSrcIdnMethod(const std::string&);
	static bool toBool(const std::string&);
	static OptionKey int2OptionKey(const int&);
	static DiseaseStage int2DiseaseStage(const int&);
//	bool canBeInfected(const DiseaseModel& disease, const DiseaseStage& stage);
	static bool canBeInfected(const DiseaseStage& stage);
	static bool isInfected(const DiseaseModel& disease, const DiseaseStage& stage);
	static bool isInfectious(const DiseaseModel& disease, const DiseaseStage& stage);
	static bool hasBeenInfected(const DiseaseModel& disease, const DiseaseStage& stage);
	static std::vector<DiseaseStage> getUnstableStages(const DiseaseModel& disease);
	static bool isStable(const DiseaseModel& disease, const std::vector<DiseaseStage>& stages, const std::vector<DiseaseStage>& unstable_stages);
	static DiseaseStage getInitialInfectedStage(const DiseaseModel& disease);
};

#endif // NETWORKPROJECT_UTILCONSTANT_H_
