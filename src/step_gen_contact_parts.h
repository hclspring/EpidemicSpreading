#ifndef NETWORKPROJECT_STEP_GEN_CONTACT_PARTS_H_
#define NETWORKPROJECT_STEP_GEN_CONTACT_PARTS_H_

#include <unistd.h>
#include <getopt.h>
#include <map>
#include <iostream>
#include <string>

#include "runner.h"

class RunnerManager;
class Parameter;

class StepGenContactParts: public Runner {
private:
	static StepGenContactParts _step_gen_contact_parts;

	StepGenContactParts();

public:
	virtual void help();
	virtual int run(const Parameter& para);

private:
//	void write_result(const std::string& res_dir, const std::vector<SourceIdentificationResV2>& res_structs);
};


#endif // NETWORKPROJECT_STEP_GEN_CONTACT_PARTS_H_
