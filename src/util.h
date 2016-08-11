#ifndef NETWORKPROJECT_UTIL_H_
#define NETWORKPROJECT_UTIL_H_

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <climits>
#include "unistd.h"

class Util {
public:
	static std::vector<std::string> addSamePrefix(const std::string&, const std::vector<std::string>&);
	static std::vector<std::string> readLines(const std::string&);

	static void checkNotNull(void* p, const std::string& p_name);

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

	template<typename Type> 
	static std::vector<Type> getVector(const std::unordered_set<Type>& input) {
		std::vector<Type> res;
		std::for_each(input.begin(), input.end(),
				[&res] (const Type& x) { res.push_back(x); } );
	}

	template<typename Type>
	static bool contains(const std::vector<Type>& vec, const Type& x) {
		for (int i = 0; i < vec.size(); ++i) {
			if (vec[i] == x) {
				return true;
			}
		}
		return false;
	}
	
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


	static double gen_rand_double();
	static int gen_rand_int(const int& range);
	static std::unordered_set<int> gen_rand_indices(const int& range, const int& count);
	static int factorial(int n);

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

	static std::string getMemoryPeakUsage();
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
};

#endif // NETWORKPROJECT_UTIL_H_
