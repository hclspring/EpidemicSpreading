#ifndef NETWORKPROJECT_UTIL_H_
#define NETWORKPROJECT_UTIL_H_

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <climits>
#include <memory>
#include <cmath>

#include "unistd.h"

class Util {
public:	
	// Check condition and give output if the condition is not satisfied
	static void checkTrue(bool condition_to_satisfy, const std::string& error_info);
	static void checkFalse(bool condition_not_to_satisfy, const std::string& error_info);

	template<typename T>
	static void checkPointer(T* p, const std::string& p_name) {
		checkTrue(p, "Error using null pointer " + p_name);
	}

	template<typename T>
	static void checkNotNull(T* p, const std::string& p_name) {
		checkPointer(p, p_name);
	}

	// Deal with strings
	static int replace_first(std::string& original, const std::string& source, const std::string& target);
	static void replace_all(std::string& original, const std::string& source, const std::string& target);
	static std::vector<std::string> addSamePrefix(const std::string&, const std::vector<std::string>&);
	static std::vector<int> parseIntegers(const std::string&, const char& split);
//	static std::vector<int> expandInterval(const int& begin, const int& end);
	static std::vector<int> expandInterval(const std::vector<std::string>&, const char& split);
	static std::string getPartString(const int& val, const int& min_length);

	// Deal with files
	static std::vector<std::string> readLines(const std::string&);
	static void writeLines(const std::vector<std::string>& lines, std::ostream& os);
	static void writeLines(const std::vector<std::string>& lines, const std::string& filename);

	template <typename T>
	static void writeVector(const std::vector<T>& vec, const std::string& filename) {
		std::ofstream ofs(filename.c_str());
		checkFalse(ofs.fail(), "Error writing file: " + filename);
		for (int i = 0; i < vec.size(); ++i) {
			ofs << vec[i] << std::endl;
		}
		ofs.close();
	}

	// Vector --> Line(s)
	template <typename T>
	static std::string getLine(const std::vector<T>& vec, const std::string& split) {
		std::string res;
		if (vec.size() > 0) {
			res.append(std::to_string(vec[0]));
		}
		for (int i = 1; i < vec.size(); ++i) {
			res.append(split);
			res.append(std::to_string(vec[i]));
		}
		return res;
	}
	
	template <typename T>
	static std::vector<std::string> getLines(const std::vector<std::vector<T>>& vec_vec, const std::string& split) {
		std::vector<std::string> res;
		for (int i = 0; i < vec_vec.size(); ++i) {
			res.push_back(getLine(vec_vec[i], split));
		}
		return res;
	}


	// Deal with file system
	static std::vector<std::string> getAllDirPaths(const std::string& root_dir, const int& layer);
	static std::vector<std::string> getAllSubDirPaths(const std::string& root_dir, const int& layer);
	static bool existDir(const std::string& dir);
	static std::vector<std::string> getAllSubs(const std::string& dir);
	static std::vector<std::string> getAllSortedSubs(const std::string& dir);
	static std::vector<std::string> getAllSubDirs(const std::string& dir);
	static std::vector<std::string> getAllSortedSubDirs(const std::string& dir);
	static std::vector<std::string> getAllSubFiles(const std::string& dir);
	static std::vector<std::string> getAllSortedSubFiles(const std::string& dir);

	// Deal with vectors and unordered_sets
	template<typename Type> static std::unordered_set<Type> vec2unset(
			const std::vector<Type>& input)
	{
		std::unordered_set<Type> res;
		std::for_each(
			input.begin(),
			input.end(),
			[&res] (const Type& x) { res.insert(x); }
		);
		return res;
	}

	template<typename Type> static std::vector<Type> unset2vec(
			const std::unordered_set<Type>& input)
	{
		std::vector<Type> res;
		std::for_each(
			input.begin(),
			input.end(),
			[&res] (const Type& x) { res.push_back(x); }
		);
		return res;
	}

	template<typename Type> 
	static std::vector<Type> getVector(const std::unordered_set<Type>& input) {
		std::vector<Type> res;
		std::for_each(input.begin(), input.end(),
				[&res] (const Type& x) { res.push_back(x); } );
	}

	// Judge whether a vector contains a specified value.
	template<typename Type>
	static bool contains(const std::vector<Type>& vec, const Type& x) {
		for (int i = 0; i < vec.size(); ++i) {
			if (vec[i] == x) {
				return true;
			}
		}
		return false;
	}

	// Given a vector with different values, return a map with indices of the same value are merged together into a vector.
	template <typename T>
	static std::unordered_map<T, std::shared_ptr<std::vector<int>>>
	get_vecptrmap(const std::vector<T>& vec) {
		std::unordered_map<T, std::shared_ptr<std::vector<int>>> res;
		for (int i = 0; i < vec.size(); ++i) {
			if (res.find(vec[i]) == res.end()) {
				std::vector<int> a(1, i);
				res.insert(std::make_pair(vec[i], std::make_shared<std::vector<int>>(a)));
			} else {
				res[vec[i]]->push_back(i);
			}
		}
		return res;
	}

	// Deal with set calculation
	template<typename Type> static std::unordered_set<Type> getIntersection(
			const std::unordered_set<Type>& setA,
			const std::unordered_set<Type>& setB)
	{
		//iterate the smaller set
		if (setA.size() <= setB.size()) {
			std::unordered_set<Type> intersectionSet;
			std::for_each(
				setA.begin(),
				setA.end(),
				[&setB, &intersectionSet] (const Type& a) {
					if (setB.count(a) > 0) { intersectionSet.insert(a);	}
				}
			);
			return intersectionSet;
		} else {
			return getIntersection(setB, setA);
		}
	}
	
	template<typename Type>	static std::unordered_set<Type> getUnion(
			const std::unordered_set<Type>& setA,
			const std::unordered_set<Type>& setB)
	{
		//iterate the smaller set
		if (setA.size() <= setB.size()) {
			std::unordered_set<Type> unionSet(setB);
			std::for_each(
				setA.begin(),
				setA.end(),
				[&unionSet] (const Type& a) { unionSet.insert(a); }
			);
			return unionSet;
		} else {
			return getUnion(setB, setA);
		}
	}

	template<typename Type> static double calcJaccardSimilarity(
			const std::unordered_set<Type>& setA,
			const std::unordered_set<Type>& setB)
	{
		return 1.0 * getIntersection(setA, setB).size() / getUnion(setA, setB).size();
	}

	template<typename Type>	static std::unordered_set<Type> getDiff(
			const std::unordered_set<Type>& setA,
			const std::unordered_set<Type>& setB)
	{
		std::unordered_set<Type> diffSet;
		std::for_each(
			setA.begin(),
			setA.end(),
			[&setB, &diffSet] (const Type& a) {
				if(setB.count(a) == 0) { diffSet.insert(a); }
			}
		);
		return diffSet;
	}

	// Deal with min/max of a vector
	template<typename Type>
	static Type getMin(const std::vector<Type>& input) {
		if (input.size() <= 0) {
			std::cerr << "Error with vector size 0." << std::endl;
			exit(-1);
		}
		Type res = input[0];
		for (int i = 1; i < input.size(); ++i) {
			res = std::min(res, input[i]);
		}
		return res;
	}
	
	template<typename Type>
	static Type getMax(const std::vector<Type>& input) {
		if (input.size() <= 0) {
			std::cerr << "Error with vector size 0." << std::endl;
			exit(-1);
		}
		Type res = input[0];
		for (int i = 1; i < input.size(); ++i) {
			res = std::max(res, input[i]);
		}
		return res;
	}

	template<typename Type>
	static int getMinIndex(const std::vector<Type>& input) {
		if (input.size() <= 0) {
			std::cerr << "Error with vector size 0." << std::endl;
			exit(-1);
		}
		Type min_value = input[0];
		int min_index = 0;
		for (int i = 1; i < input.size(); ++i) {
			if (input[i] < min_value) {
				min_index = i;
				min_value = input[i];
			}
		}
		return min_index;
	}
	
	template<typename Type>
	static int getMaxIndex(const std::vector<Type>& input) {
		if (input.size() <= 0) {
			std::cerr << "Error with vector size 0." << std::endl;
			exit(-1);
		}
		Type max_value = input[0];
		int max_index = 0;
		for (int i = 1; i < input.size(); ++i) {
			if (input[i] > max_value) {
				max_index = i;
				max_value = input[i];
			}
		}
		return max_index;
	}

	// Deal with random functions
	static double gen_rand_double();
	static int gen_rand_int(const int& range);
	static std::vector<double> gen_rand_vector_double(const int& size, const double& lower_bound, const double& upper_bound);
	static std::unordered_set<int> gen_rand_indices(const int& range, const int& count);

	// Deal with maths
	static double factorial(int n);
	static double choose(const int& n, const int& k);

	// Deal with a matrix.
	template<typename Type>
	static bool isMatrix(const std::vector<std::vector<Type>>& input) {
		int n = input.size();
		if (n == 0) {
			return false;
		}
		int m = input[0].size();
		for (int i = 1; i < n; ++i) {
			if (input[i].size() != m) {
				return false;
			}
		}
		return true;
	}

	template<typename Type>
	static bool isSquareMatrix(const std::vector<std::vector<Type>>& input) {
		if (!isMatrix(input)) {
			return false;
		}
		int n = input.size();
		if (input[0].size() != n) {
			return false;
		} else {
			return true;
		}
	}
	
	template<typename Type>
	static bool isSymmetricMatrix(const std::vector<std::vector<Type>>& input) {
		if (!isSquareMatrix(input)) {
			return false;
		}
		int n = input.size();
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				if (input[i][j] != input[j][i]) {
					return false;
				}
			}
		}
		return true;
	}

	template<typename T>
	static int getMatrixSize(const std::vector<std::vector<T>>& matrix, int& rows, int& columns) {
		if (isMatrix(matrix) == false) {
			std::cerr << "Error for finding the input not being a matrix." << std::endl;
			exit(-1);
		}
		rows = matrix.size();
		columns = matrix[0].size();
		return rows * columns;
	}

	template<typename Type>
	static std::vector<std::vector<Type>> eliminateFromSquareMatrix(const std::vector<std::vector<Type>>& input, const int& index_to_eliminate) {
		if (!isSquareMatrix(input)) {
			std::cerr << "Error for finding the matrix non-square." << std::endl;
			exit(-1);
		}
		int n = input.size();
		if (index_to_eliminate < 0 || index_to_eliminate >= n) {
			std::cerr << "Error for illegal index to eliminate: " << index_to_eliminate << std::endl;
			exit(-1);
		}
		std::vector<std::vector<Type>> res(n - 1, std::vector<Type>(n - 1));
		for (int i = 0; i < index_to_eliminate; ++i) {
			for (int j = 0; j < index_to_eliminate; ++j) {
				res[i][j] = input[i][j];
			}
			for (int j = index_to_eliminate + 1; j < n; ++j) {
				res[i][j-1] = input[i][j];
			}
		}
		for (int i = index_to_eliminate + 1; i < n; ++i) {
			for (int j = 0; j < index_to_eliminate; ++j) {
				res[i-1][j] = input[i][j];
			}
			for (int j = index_to_eliminate + 1; j < n; ++j) {
				res[i-1][j-1] = input[i][j];
			}
		}
		return res;
	}

	template<typename Type>
	static std::vector<std::vector<Type>> getSubMatrix(
			const std::vector<std::vector<Type>>& input, 
			const std::unordered_set<int>& indices) {
		if (!isMatrix(input)) {
			std::cerr << "Error for finding input not an matrix." << std::endl;
			exit(-1);
		}
		std::vector<std::vector<Type>> res;
		for (int i = 0; i < input.size(); ++i) {
			if (indices.find(i) != indices.end()) {
				std::vector<Type> temp;
				for (int j = 0; j < input[i].size(); ++j) {
					if (indices.find(j) != indices.end()) {
						temp.push_back(input[i][j]);
					}
				}
				res.push_back(temp);
			}
		}
		return res;
	}

	template<typename T>
	static std::vector<std::vector<T>> getDiffMatrix(
			const std::vector<std::vector<T>>& x1,
			const std::vector<std::vector<T>>& x2) {
		int r1, c1, r2, c2;
		getMatrixSize(x1, r1, c1);
		getMatrixSize(x2, r2, c2);
		checkTrue(r1 == r2 && c1 == c2, "Error: Cannot calculate difference of two matrices with different sizes.");
		std::vector<std::vector<T>> res(r1, std::vector<T>(c1));
		for (size_t i = 0; i < r1; ++i) {
			for (size_t j = 0; j < c1; ++j) {
				res[i][j] = x1[i][j] - x2[i][j];
			}
		}
		return res;
	}

	template<typename T>
	static double getFrobeniusNorm(const std::vector<std::vector<T>>& m) {
		checkTrue(isMatrix(m), "Error: input is not a matrix.");
		double sum = 0;
		for (size_t i = 0; i < m.size(); ++i) {
			for (size_t j = 0; j < m[i].size(); ++j) {
				sum += m[i][j] * m[i][j];
			}
		}
		return sqrt(sum);
	}

	template <typename T>
	static std::vector<std::vector<T>> getMin(const std::vector<std::vector<T>>& m1,
			const std::vector<std::vector<T>>& m2) {
		int r1, c1, r2, c2;
		getMatrixSize(m1, r1, c1);
		getMatrixSize(m2, r2, c2);
		checkTrue(r1 == r2 && c1 == c2, "Error: Cannot calculate minimum of two matrices with different sizes.");
		std::vector<std::vector<T>> res(r1, std::vector<T>(c1));
		for (size_t i = 0; i < r1; ++i) {
			for (size_t j = 0; j < c1; ++j) {
				res[i][j] = std::min(m1[i][j], m2[i][j]);
			}
		}
		return res;
	}

	template <typename T>
	static int count(const std::vector<std::vector<T>>& m, const T& x) {
		int res = 0;
		for (int i = 0; i < m.size(); ++i) {
			for (int j = 0; j < m[i].size(); ++j) {
				if (m[i][j] == x) {
					++res;
				}
			}
		}
		return res;
	}

	template <typename T>
	static int count(const std::vector<std::vector<T>>& m, const T& x, const T& error) {
		checkTrue(error > 0, "Error: error must be a positive small number.");
		int res = 0;
		for (int i = 0; i < m.size(); ++i) {
			for (int j = 0; j < m[i].size(); ++j) {
				if (m[i][j] - x > -error && m[i][j] -x < error) {
					++res;
				}
			}
		}
		return res;
	}

	// Get system usage.

	static std::string getMemoryPeakUsage();

	// Some calculations of a vector.
	static double getMean(const std::vector<double>& input);
	static double getDeviation(const std::vector<double>& input, const double& mean);

	template<typename Type>
	static Type getSum(const std::vector<Type>& input) {
		Type sum = 0;
		for (int i = 0; i < input.size(); ++i) {
			sum += input[i];
		}
		return sum;
	}
	
	template<typename Type>
	static std::vector<int> getSortedIndices(const std::vector<Type>& input, const bool& ascending) {
		struct myclass {
			int index;
			Type value;
		};
		int n = input.size();
		std::vector<myclass> myobjects(n);
		for (int i = 0; i < n; ++i) {
			myobjects[i].index = i;
			myobjects[i].value = input[i];
		}
		struct MyCompare {
			bool _ascending;
			MyCompare(const bool& x) {
				_ascending = x;
			}
			bool operator() (const myclass& i, const myclass& j) const {
				if (_ascending) {
					return i.value < j.value;
				} else {
					return i.value > j.value;
				}
			}
		};
		MyCompare mycompare(ascending);
		std::sort(myobjects.begin(), myobjects.end(), mycompare);
		std::vector<int> res(n);
		for(int i = 0; i < n; ++i) {
			res[i] = myobjects[i].index;
		}
		return res;
	}



private:
	// get all sub directories/files under "dir"
	// "type" is used to denote whether directories or files are wanted. 4 for directories, 8 for files, and 12 for both.
	static std::vector<std::string> getAllSubs(const std::string& dir, const int& type);
};

#endif // NETWORKPROJECT_UTIL_H_
