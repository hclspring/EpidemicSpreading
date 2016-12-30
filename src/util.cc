#include "util.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <dirent.h>

void Util::checkTrue(bool condition_to_satisfy, const std::string& error_info) {
	if (!condition_to_satisfy) {
		std::cerr << error_info << std::endl;
		exit(-1);
	}
}

void Util::checkFalse(bool condition_not_to_satisfy, const std::string& error_info) {
	checkTrue(!condition_not_to_satisfy, error_info);
}

int Util::replace_first(std::string& original, const std::string& source, const std::string& target) {
	int pos = 0;
	if ((pos = original.find(source.c_str())) != std::string::npos) {
		original.replace(pos, source.size(), target);
		return pos;
	} else {
		return -1;
	}
}
	
void Util::replace_all(std::string& original, const std::string& source, const std::string& target) {
	int pos = -1;
	do {
		pos = replace_first(original, source, target);
	} while (pos >= 0);
}

std::vector<std::string> Util::addSamePrefix(const std::string& prefix, const std::vector<std::string>& strs) {
	std::vector<std::string> res(strs.size());
	for (size_t i = 0; i < strs.size(); ++i) {
		res[i] = prefix;
		res[i].append(strs[i]);
	}
	return res;
}

std::vector<int> Util::parseIntegers(const std::string& str, const char& split) {
	std::vector<int> res;
	int index_split = 0;
	int next_split = str.find_first_of(split, index_split);
	while (next_split != std::string::npos) {
		res.push_back(std::stoi(str.substr(index_split, next_split - index_split)));
		index_split = next_split + 1;
		next_split = str.find_first_of(split, index_split);
	}
	res.push_back(std::stoi(str.substr(index_split)));
	return res;
}

/*
std::vector<int> Util::expandInterval(const int& begin, const int& end) {
	Util::checkTrue(begin < end
	if (begin > end) {
		return expandInterval(end, begin);
	} else {
		std::vector<int> res;
		for (int i = begin; i <= end; ++i) {
			res.push_back(i);
		}
		return res;
	}
}
*/
std::vector<int> Util::expandInterval(const std::vector<std::string>& strs, const char& split) {
	std::vector<int> res;
	for (int i = 0; i < strs.size(); ++i) {
		std::vector<int> temp = parseIntegers(strs[i], split);
		if (temp.size() == 1) {
			res.push_back(std::stoi(strs[i]));
		} else if (temp.size() == 2) {
			Util::checkTrue(temp[0] <= temp[1], "Error: first integer > second integer in the string \"" + strs[i] + "\".");
			for (int j = temp[0]; j <= temp[1]; ++j) {
				res.push_back(j);
			}
		} else {
			Util::checkTrue(false, "Error: parsing string \"" + strs[i] + "\".");
		}
	}
	return res;
}

std::string Util::getPartString(const int& val, const int& min_length){
	checkTrue(val >= 0, "Error because of trying to get a part string from a negative integer.");
	std::string res = std::to_string(val);
	if (res.length() < min_length) {
		std::string res1;
		for (int i = 0; i < min_length - res.length(); ++i) {
			res1.push_back('0');
		}
		res1.append(res);
		return res1;
	} else {
		return res;
	}
}

std::vector<std::string> Util::readLines(const std::string& infile) {
	std::vector<std::string> res;
	std::ifstream ifs(infile.c_str());
	checkFalse(ifs.fail(), "Error reading file: " + infile);
	std::string line;
	while (getline(ifs, line, '\n')) {
		res.push_back(line);
	}
	ifs.close();
	return res;
}

void Util::writeLines(const std::vector<std::string>& lines, std::ostream& os) {
	for (int i = 0; i < lines.size(); ++i) {
		os << lines[i] << std::endl;
	}
}

void Util::writeLines(const std::vector<std::string>& lines, const std::string& filename) {
	std::ofstream ofs(filename.c_str());
	checkFalse(ofs.fail(), "Error writing file: " + filename);
	writeLines(lines, ofs);
	ofs.close();
}

std::vector<std::string> Util::getAllDirPaths(const std::string& root_dir, const int& layer) {
	if (layer <= 0) {
		return std::vector<std::string> { root_dir + "/" };
	} else {
		std::vector<std::string> subdirs = getAllSortedSubDirs(root_dir);
		std::vector<std::string> res, temp_res;
		for (int i = 0; i < subdirs.size(); ++i) {
			temp_res = getAllDirPaths(root_dir + "/" + subdirs[i], layer - 1);
			for (int j = 0; j < temp_res.size(); ++j) {
				res.push_back(temp_res[j]);
			}
		}
		return res;
	}
}

std::vector<std::string> Util::getAllSubDirPaths(const std::string& root_dir, const int& layer) {
	if (layer <= 0) {
		return std::vector<std::string> { "" };
	} else {
		std::vector<std::string> subdirs = getAllSortedSubDirs(root_dir);
		std::vector<std::string> res, temp_res;
		for (int i = 0; i < subdirs.size(); ++i) {
			temp_res = getAllSubDirPaths(root_dir + "/" + subdirs[i], layer - 1);
			for (int j = 0; j < temp_res.size(); ++j) {
				res.push_back(subdirs[i] + "/" + temp_res[j]);
			}
		}
		return res;
	}
}


bool Util::existDir(const std::string& dir) {
	return opendir(dir.c_str());
}

std::vector<std::string> Util::getAllSubs(const std::string& dir) {
	return getAllSubs(dir, 12);
}

std::vector<std::string> Util::getAllSortedSubs(const std::string& dir) {
	std::vector<std::string> res = getAllSubs(dir);
	std::sort(res.begin(), res.end());
	return res;
}

std::vector<std::string> Util::getAllSubDirs(const std::string& dir) {
	return getAllSubs(dir, 4);
}

std::vector<std::string> Util::getAllSortedSubDirs(const std::string& dir) {
	std::vector<std::string> res = getAllSubDirs(dir);
	std::sort(res.begin(), res.end());
	return res;
}

std::vector<std::string> Util::getAllSubFiles(const std::string& dir) {
	return getAllSubs(dir, 8);
}
	
std::vector<std::string> Util::getAllSortedSubFiles(const std::string& dir) {
	std::vector<std::string> res = getAllSubFiles(dir);
	std::sort(res.begin(), res.end());
	return res;
}

double Util::gen_rand_double() {
	return (1.0 * rand()) / RAND_MAX;
}

int Util::gen_rand_int(const int& range) {
	if (range <= 0) {
		std::cerr << "Error generating random integer in range [0, " << range << ")." << std::endl;
		exit(-1);
	}
	return rand() % range;
}

std::vector<double> Util::gen_rand_vector_double(const int& size, const double& lower_bound, const double& upper_bound) {
	std::vector<double> res(size);
	for (int i = 0; i < size; ++i) {
		res[i] = (upper_bound - lower_bound) * gen_rand_double() + lower_bound;
	}
	return res;
}

std::unordered_set<int> Util::gen_rand_indices(const int& range, const int& count) {
	if (range <= 0 || count <= 0 || count > range) {
		std::cerr << "Error generating " << count << " random indices in range [0, " << range << ")." << std::endl;
		exit(-1);
	}
	std::unordered_set<int> res;
	if (count < sqrt(range)) {
		int cur;
		for (int i = 0; i < count; ++i) {
			do {
				cur = gen_rand_int(range);
			} while (res.count(cur) > 0);
			res.insert(cur);
		}
	} else if (count < range - sqrt(range)) {
		std::vector<int> temp(range, 0);
		for (int i = 0; i < range; ++i) {
			temp[i] = i;
		}
		random_shuffle(temp.begin(), temp.end());
		for (int i = 0; i < count; ++i) {
			res.insert(temp[i]);
		}
	} else {
		std::unordered_set<int> other_indices = gen_rand_indices(range, range - count);
		for (int i = 0; i < range; ++i) {
			if (other_indices.count(i) == 0) {
				res.insert(i);
			}
		}
	}
	return res;
}

double Util::factorial(int n) {
	if (n < 0) {
		std::cerr << "Error calculation the factorial of " << n << "." << std::endl;
		exit(-1);
	} else if (n == 0 || n == 1) {
		return 1;
	} else {
		return n * factorial(n - 1);
	}
}

double Util::choose(const int& n, const int& k) {
	checkTrue(k >= 0 && n >= k, "Error: try to calculate illegal choose number.");
	if (k > n / 2) {
		return choose(n, n - k);
	} else if (n == 0) {
		return 1;
	} else {
		std::vector<double> a(k + 1, 1);
		for (int nn = 2; nn <= n; ++nn) {
			for (int kk = std::min(nn - 1, k); kk >= 1; --kk) {
				a[kk] += a[kk-1];
			}
//			std::cout << a[k] << std::endl;
		}
		return a[k];
	}
}



std::string Util::getMemoryPeakUsage() {
	pid_t pid = getpid();
	std::stringstream sstr;
	sstr << pid;
	std::string pid_str;
	sstr >> pid_str;
	std::string status_file = "/proc/" + pid_str + "/status";
	std::ifstream ifs(status_file.c_str());
	std::string result, line, name;	
	while (getline(ifs, line, '\n')) {
		std::stringstream sstr;
		sstr << line;
		sstr >> name;
		if (name.compare("VmPeak:") == 0) {
			sstr >> line;
			result = line; // 数值
			sstr >> line;
			result.append(" ").append(line); // 单位
			ifs.close();
			return result;
		}
	}

	std::cerr << "Error to find VmPeak line." << std::endl;
	ifs.close();
	return "error";
}

double Util::getMean(const std::vector<double>& input) {
	int n = input.size();
	if (n <= 0) {
		std::cerr << "Error for no numbers in the vector." << std::endl;
		exit(-1);
	}
	double sum = 0;
	std::for_each(input.begin(), input.end(),
			[&sum] (const double& x) { sum += x; });
	return sum / n;
}

double Util::getDeviation(const std::vector<double>& input, const double& mean) {
	int n = input.size();
	if (n <= 0) {
		std::cerr << "Error for no numbers in the vector." << std::endl;
		exit(-1);
	}
	double sum = 0;
	std::for_each(input.begin(), input.end(),
			[&sum, mean] (const double& x) { sum += (x - mean) * (x - mean); });
	return sqrt(sum / n);
}

std::vector<std::string> Util::getAllSubs(const std::string& root_dir, const int& type) {
	struct dirent *ent = NULL;
	DIR *pDir = opendir(root_dir.c_str());
	checkNotNull(pDir, "Error for the directory " + root_dir + " does not exist.");
	std::vector<std::string> res;
	while (ent = readdir(pDir)) {
		if (ent->d_type & type) {
			std::string sub(ent->d_name);
			if (sub.compare(".") != 0 && sub.compare("..") != 0) {
				res.push_back(sub);
			}
		}
	}
	return res;
}


