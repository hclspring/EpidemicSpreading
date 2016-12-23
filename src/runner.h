#ifndef NETWORKPROJECT_RUNNER_H_
#define NETWORKPROJECT_RUNNER_H_

#include <iostream>
#include <string>
#include <map>
#include <vector>

class Parameter;

class Runner {
protected:
	std::string _short_options;
	struct option* _long_options;

public:

	Runner() {}
	virtual ~Runner() {}

	virtual int run(const Parameter& para) = 0;
	virtual void help() = 0;

	void set_short_options(const std::string& so) {
		_short_options = so;
	}
	void set_long_options(struct option* lo) {
		_long_options = lo;
	}
	const std::string& get_short_options() const {
		return _short_options;
	}
	const struct option* get_long_options() const {
		return _long_options;
	}
};


#endif //NETWORKPROJECT_RUNNER_H_
