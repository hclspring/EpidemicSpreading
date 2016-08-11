#ifndef _NETWORKPROJECT_NODEDATA_H_
#define _NETWORKPROJECT_NODEDATA_H_

#include <string>


class NodeData {
private:
	std::string _label;

public:
	NodeData(): NodeData("") {}
	NodeData(const std::string& label): _label(label) {}
	~NodeData() { _label.clear(); }

	void set_label(const std::string& label) { _label = label; }
	std::string get_label() { return _label; }
};

#endif // _NETWORKPROJECT_NODEDATA_H_
