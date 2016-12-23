#ifndef NETWORKPROJECT_STEP_TEST_H_
#define NETWORKPROJECT_STEP_TEST_H_

#include <unistd.h>
#include <getopt.h>
#include <map>
#include <iostream>
#include <string>

#include "runner.h"

class RunnerManager;
class Parameter;

class StepTest: public Runner {
private:
	static StepTest _step_test;

	StepTest();

public:
	virtual void help();
	virtual int run(const Parameter& para);
};


#endif // NETWORKPROJECT_STEP_TEST_H_
