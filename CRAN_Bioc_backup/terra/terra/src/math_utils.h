// Copyright (c) 2018-2023  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.
#include <cmath>
#include <limits>
#include <vector>
#include <random>



double modal_value(std::vector<double> values, unsigned ties, bool narm, std::default_random_engine rgen, std::uniform_real_distribution<double> dist);

void na_omit(std::vector<double> &x);

bool is_equal(double a, double b, double tolerance=10.0);
bool about_equal(double a, double b, double tolerance);
bool is_equal_relative(double a, double b, double tolerance);
bool is_equal_range(double x, double y, double range, double tolerance);
void vector_minmax(std::vector<double> v, double &min, int &imin, double &max, int &imax);
double roundn(double x, int n);
double signif(double x, unsigned n);


template <typename Iterator>
void minmax(Iterator start, Iterator end, double &vmin, double &vmax, double exclude) {
    vmin = std::numeric_limits<double>::max();
    vmax = std::numeric_limits<double>::lowest();
    bool none = true;
	if (std::isnan(exclude)) {
		for (Iterator v = start; v !=end; ++v) {
			if (!std::isnan(*v)) {
				if (*v > vmax) {
					vmax = *v;
					none = false;
				}
				if (*v < vmin) {
					vmin = *v;
				}
			}
		}
	} else {
		for (Iterator v = start; v !=end; ++v) {
			if ((!std::isnan(*v)) && (exclude != *v)) {
				if (*v > vmax) {
					vmax = *v;
					none = false;
				}
				if (*v < vmin) {
					vmin = *v;
				}
			}
		}
	}
    if (none) {
        vmin = NAN;
        vmax = NAN;
    }
}




template <typename T> 
void sort_unique_2d(std::vector<T> &x, std::vector<T> &y) {

	std::vector<std::vector<T>> v(x.size());
	for (size_t i=0; i<v.size(); i++) {
		v[i] = {x[i], y[i]};
	}

    std::sort(v.begin(), v.end());
	v.erase(std::unique(v.begin(), v.end()), v.end());
	x.resize(v.size());
	y.resize(v.size());
	for (size_t i=0; i<x.size(); i++) {
		x[i] = v[i][0];
		y[i] = v[i][1];
	}
}
