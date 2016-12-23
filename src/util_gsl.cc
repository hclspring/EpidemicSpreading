#include "util_gsl.h"

#include <vector>
#include <iostream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
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


double UtilGsl::getMinWithNelderMead(
		const std::vector<double>& init_solution, 
		const std::vector<double>& step, 
		const int& iter_max,
		const double& simplex_size_max,
//		double (*out_f)(const std::vector<double>&),
		void** f_para,
		const bool& iteration_info_wanted,
		std::vector<double>& best_solution
) {
	const gsl_multimin_fminimizer_type *algo = gsl_multimin_fminimizer_nmsimplex2;
	Util::checkTrue(init_solution.size() == step.size(), "Error: init_solution.size() != step.size().");
	int n = init_solution.size();

	// starting point
	gsl_vector *x = get_gsl_vector(init_solution);
	gsl_vector *step_size = get_gsl_vector(step);

	gsl_multimin_function min_func;
	min_func.n = n;
	min_func.f = &gsl_nelder_mead_eval;
	min_func.params = (void*) f_para;

	using namespace std;

	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(algo, n);
	gsl_multimin_fminimizer_set(s, &min_func, x, step_size);
	size_t iter = 0;
	int status = GSL_CONTINUE;
	while (status == GSL_CONTINUE && iter < iter_max) {
		status = gsl_multimin_fminimizer_iterate(s);

		if (status) break;

		double simplex_size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(simplex_size, simplex_size_max);
	
		if (iteration_info_wanted) {
			cout << "Iteration\t" << iter << ":\t";
			cout << "x = ";
			for (int i = 0; i < 3; ++i) {
				cout << gsl_vector_get(s->x, i) << ", ";
			}
			cout << "y = " << s->fval << "; simplex size = " << simplex_size << endl;
		}

//		cout << simplex_size << "   " << status << endl;
		++iter;
	}

	gsl_vector_free(x);
	gsl_vector_free(step_size);

	best_solution = get_std_vector(s->x);
	double res = s->fval;
	gsl_multimin_fminimizer_free(s);

	return res;
}


gsl_vector* UtilGsl::get_gsl_vector(const std::vector<double>& x) {
	gsl_vector* res = gsl_vector_alloc(x.size());
	for (size_t i = 0; i < x.size(); ++i) {
		gsl_vector_set(res, i, x[i]);
	}
	return res;
}

std::vector<double> UtilGsl::get_std_vector(const gsl_vector* x) {
	std::vector<double> res(x->size, 0);
	for (size_t i = 0; i < x->size; ++i) {
		res[i] = gsl_vector_get(x, i);
	}
	return res;
}


/*
UtilGsl::ConvertInNM UtilGsl::getNMfunc(UtilGsl::ConvertBeforeNM, const std::vector<double>& v);
//double UtilGsl::(*getNMfunc(UtilGsl::ConvertBeforeNM)) (const gsl_vector* v, void* para) {
	return ConvertBeforeNM(get_std_vector(v));
}
*/

double UtilGsl::gsl_nelder_mead_eval(const gsl_vector* v, void* para) {
	void** new_para = (void**) para;
	FuncBeforePtr f = (FuncBeforePtr) new_para[0];
	return (*f)(get_std_vector(v), new_para[1]);
}

/*
UtilGsl::FuncInPtr UtilGsl::getNMfunc(UtilGsl::FuncBeforePtr f, const std::vector<double>& x) {
	return f();
}
*/
