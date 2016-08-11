#include "util_gsl.h"

#include <vector>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

#include "util.h"

std::vector<double> UtilGsl::calcSymmMatrixEigenvalues(const std::vector<std::vector<double>>& symm_matrix) {
	if (Util::isSymmetricMatrix(symm_matrix) == false) {
		std::cerr << "Error for finding the matrix non-symmetric." << std::endl;
		exit(-1);
	}
	std::vector<double> eigenvalues;
	int n;
	Util::getMatrixSize(symm_matrix, n, n);
	double* matrix_data = getMatrixArray(symm_matrix);
	gsl_matrix_view m = gsl_matrix_view_array(matrix_data, n, n);
	gsl_vector* eval = gsl_vector_alloc(n);
	gsl_matrix* evec = gsl_matrix_alloc(n, n);
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(n);
	gsl_eigen_symmv(&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

	{
		int i;
		for (i = 0; i < n; ++i) {
			double eval_i = gsl_vector_get(eval, i);
			eigenvalues.push_back(eval_i);
			gsl_vector_view evec_i = gsl_matrix_column(evec, i);
		}
	}

	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	free(matrix_data);
	return eigenvalues;
}

std::vector<std::vector<double>> UtilGsl::calcSymmMatrixEigenvectors(const std::vector<std::vector<double>>& symm_matrix) {
	if (Util::isSymmetricMatrix(symm_matrix) == false) {
		std::cerr << "Error for finding the matrix non-symmetric." << std::endl;
		exit(-1);
	}
	std::vector<double> eigenvalues;
	std::vector<std::vector<double>> eigenvectors;
	int n;
	Util::getMatrixSize(symm_matrix, n, n);
	double* matrix_data = getMatrixArray(symm_matrix);
	gsl_matrix_view m = gsl_matrix_view_array(matrix_data, n, n);
	gsl_vector* eval = gsl_vector_alloc(n);
	gsl_matrix* evec = gsl_matrix_alloc(n, n);
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(n);
	gsl_eigen_symmv(&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

	{
		eigenvectors.resize(n, std::vector<double>(n, 0));
		int i;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				eigenvectors[i][j] = gsl_matrix_get(evec, j, i);
			}
		}
	}

	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	free(matrix_data);
	return eigenvectors;
}


