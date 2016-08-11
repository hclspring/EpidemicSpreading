#ifndef NETWORKPROJECT_IGRAPH_H_
#define NETWORKPROJECT_IGRAPH_H_

#include <vector>
#include <list>
#include <unordered_set>

#include <igraph.h>

class UtilIgraph {
	static igraph_vector_t getIgraphVector(const std::vector<double>& nums);
};


#endif // NETWORKPROJECT_IGRAPH_H_
