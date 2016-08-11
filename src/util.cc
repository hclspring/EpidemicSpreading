#include "util.h"

#include "sstream"

std::vector<std::string> Util::addSamePrefix(const std::string& prefix, const std::vector<std::string>& strs) {
	std::vector<std::string> res(strs.size());
	for (size_t i = 0; i < strs.size(); ++i) {
		res[i] = prefix;
		res[i].append(strs[i]);
	}
	return res;
}

std::vector<std::string> Util::readLines(const std::string& infile) {
	std::vector<std::string> res;
	std::ifstream ifs(infile.c_str());
	if (ifs.fail()) {
		std::cerr << "Error reading file: " << infile << std::endl;
		exit(-1);
	}
	std::string line;
	while (getline(ifs, line, '\n')) {
		res.push_back(line);
	}
	ifs.close();
	return res;
}

void Util::checkNotNull(void* p, const std::string& p_name) {
	if (p == NULL) {
		std::cerr << "Error using null pointer " << p_name << std::endl;
		exit(-1);
	}
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

int Util::factorial(int n) {
	if (n < 0) {
		std::cerr << "Error calculation the factorial of " << n << "." << std::endl;
		exit(-1);
	} else if (n == 0 || n == 1) {
		return 1;
	} else {
		return n * factorial(n - 1);
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

