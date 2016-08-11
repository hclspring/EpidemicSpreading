#ifndef NETWORKPROJECT_UTILGSL_H_
#define NETWORKPROJECT_UTILGSL_H_

#include <vector>
#include <iostream>

#include "util.h"

class UtilGsl {
public:
	static std::vector<double> calcSymmMatrixEigenvalues(const std::vector<std::vector<double>>& symm_matrix);
	static std::vector<std::vector<double>> calcSymmMatrixEigenvectors(const std::vector<std::vector<double>>& symm_matrix);

private:
	template<typename T>
	static double* getMatrixArray(const std::vector<std::vector<T>>& matrix) {
		if (Util::isMatrix(matrix) == false) {
			std::cerr << "Error for finding the input not being a matrix." << std::endl;
			exit(-1);
		}
		int rows, columns, total_size;
		total_size = Util::getMatrixSize(matrix, rows, columns);
		double* res = new double[total_size];
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				res[i * columns + j] = static_cast<double>(matrix[i][j]);
			}
		}
		return res;
	}
};

#endif // NETWORKPROJECT_UTILGSL_H_
