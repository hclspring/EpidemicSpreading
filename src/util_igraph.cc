#include "util_igraph.h"


igraph_vector_t UtilIgraph::getIgraphVector(const std::vector<double>& nums) {
	igraph_vector_t res;
	igraph_vector_init(&res, nums.size());
	for (int i = 0; i < nums.size(); ++i) {
		VECTOR(res)[i] = nums[i];
	}
	return res;
}



