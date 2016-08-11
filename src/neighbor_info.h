#ifndef NETWORKPROJECT_NEIGHBORINFO_H_
#define NETWORKPROJECT_NEIGHBORINFO_H_

#include <string>

class NeighborInfo {
private:
	std::string _name;
	double _weight;

public:
	// constructor
	NeighborInfo(const std::string& tn, const double& w):
		_name(tn), _weight(w) {}
	NeighborInfo() : NeighborInfo("", 0) {}
	~NeighborInfo() {}

	// set functions
	void set_name(const std::string& tn) { _name = tn; }
	void set_weight(const double& w) { _weight = w; }

	// get functions
	std::string get_name() const { return _name; }
	double get_weight() const { return _weight; }

};

#endif // NETWORKPROJECT_NEIGHBORINFO_H_
