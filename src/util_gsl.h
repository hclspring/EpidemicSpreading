#ifndef NETWORKPROJECT_UTILGSL_H_
#define NETWORKPROJECT_UTILGSL_H_

#include <vector>
#include <iostream>
#include <functional>

#include <gsl/gsl_multimin.h>

#include "util.h"

class UtilGsl {
public:
	typedef double FunctionBeforeNM(const std::vector<double>&, void*);
	typedef FunctionBeforeNM* FuncBeforePtr;
	typedef double FunctionInNM(const gsl_vector*, void*);
	typedef FunctionInNM* FuncInPtr;
	// About calculate eigenvalues and eigenvectors
	static std::vector<double> calcSymmMatrixEigenvalues(const std::vector<std::vector<double>>& symm_matrix);
	static std::vector<std::vector<double>> calcSymmMatrixEigenvectors(const std::vector<std::vector<double>>& symm_matrix);

	// About Nelder-Mead (Simplex) algorithm
	//static double getMinWithNelderMead(const std::vector<double>& init_solution, double (*f)(const std::vector<double>&), std::vector<double>& best_solution);
	static double getMinWithNelderMead(
		const std::vector<double>& init_solution, 
		const std::vector<double>& step_size, 
		const int& iter_max,
		const double& simplex_size_max,
//		std::function<double(const std::vector<double>&)> f,
//		double (*f)(const std::vector<double>&), // 此行改成下一行
		void** f_para, // f_para points to an array of pointers, where the first pointer points to a function (double f(const vector<double>&, void*)), and the second pointer points to the second parameter of the function.
		const bool& iteration_info_wanted,
		std::vector<double>& best_solution
);

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

	static gsl_vector* get_gsl_vector(const std::vector<double>& x);
	static std::vector<double> get_std_vector(const gsl_vector*);

	static double gsl_nelder_mead_eval(const gsl_vector* v, void* para);

//	static FuncInPtr getNMfunc(FuncBeforePtr, const std::vector<double>&);

};

#endif // NETWORKPROJECT_UTILGSL_H_
