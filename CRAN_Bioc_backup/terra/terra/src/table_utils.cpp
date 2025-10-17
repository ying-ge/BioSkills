#include <map>
#include <vector>
#include <cmath>
#include <algorithm>


std::map<double, size_t> table(std::vector<double> &v) {
	std::map<double, size_t> count;
	for_each( v.begin(), v.end(), [&count]( double val ){
			if(!std::isnan(val)) count[val]++;
		}
	);
	return count;
}


std::map<double, size_t> combine_tables(std::map<double, size_t> &x, std::map<double, size_t> &y) {
	for(auto p : y) {
		x[p.first] += p.second;
	}
	return(x);
}


std::vector<double> table2vector(std::map<double, size_t> &x) {
	std::vector<std::vector<double>> out(2);
	for( auto p : x ) {
		out[0].push_back(p.first);
		out[1].push_back(p.second);
	}
	out[0].insert(out[0].end(), out[1].begin(), out[1].end());
	return out[0];
}

std::vector<std::vector<double>> table2vector2(std::map<double, size_t> &x) {
	std::vector<std::vector<double>> out(2);
	for( auto p : x ) {
		out[0].push_back(p.first);
		out[1].push_back(p.second);
	}
	return out;
}
