#include <unistd.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <igraph.h>

#include "runner_manager.h"
#include "runner.h"
//#include "step_test.h"
#include "constant.h"
#include "parameter.h"

int main(int argc, char** argv) {
	igraph_i_set_attribute_table(&igraph_cattribute_table);
	srand(time(NULL));

	if (argc < 2) {
		return RunnerManager::help();
	}

	RunnerPtr runner = RunnerManager::instance()->get_runner(std::string(argv[1]));
	if (runner == NULL) {
		return RunnerManager::help();
	}

	Parameter para(argc - 1, argv + 1, runner->get_short_options().c_str(), runner->get_long_options());
	if (para.get_help()) {
		runner->help();
	} else {
		runner->run(para);
	}

	return 0;
}
