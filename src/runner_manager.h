#ifndef NETWORKPROJECT_RUNNER_MANAGER_H_
#define NETWORKPROJECT_RUNNER_MANAGER_H_

#include <iostream>
#include <string>
#include <map>

class Runner;
typedef Runner* RunnerPtr;

class RunnerManager {
private:
	std::map<std::string, RunnerPtr> _runners;

public:
	// 根据命令返回适用的runner
	RunnerPtr get_runner(const std::string& name) const {
		return _runners.find(name) == _runners.end() ? NULL : _runners.find(name)->second;
	}

	void install(const std::string& name, RunnerPtr p) {
		_runners.insert(std::make_pair(name, p));
	}

	static RunnerManager* instance() {
		static RunnerManager manager;
		return &manager;
	}

	static int help() {
		std::cout << "Helps" << std::endl;
		return 1;
	}

};

#endif // NETWORKPROJECT_RUNNER_MANAGER_H_
