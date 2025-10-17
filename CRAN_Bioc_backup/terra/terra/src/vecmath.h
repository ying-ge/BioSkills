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

#ifndef VECMATH_GUARD
#define VECMATH_GUARD

#include <functional>
#include <iterator>
#include <string>
#include <type_traits>
#include <vector>
#include "NA.h"
#include <math.h>

#include <map>
#include <algorithm>
#include <unordered_map>


bool haveFun(std::string fun);
std::function<double(std::vector<double>&, bool)> getFun(std::string fun);
bool bany(const std::vector<bool>& v);
bool ball(const std::vector<bool>& v);


template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size();
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}


static inline double interpolate(double x, double y1, double y2, double x1, double x2) {
	double denom = (x2-x1);
	return y1 + (x-x1) * (y2-y1)/denom;
}



static inline std::vector<double> vquantile(std::vector<double> v, const std::vector<double>& probs, bool narm) {
	size_t n = v.size();
    if (n==0) {
        return std::vector<double>(probs.size(), NAN);
    }
    if (n == 1) {
        return std::vector<double>(probs.size(), v[0]);
    }

	//na_omit(v);
	v.erase(std::remove_if(std::begin(v), std::end(v),
        [](const double& value) { return std::isnan(value); }),
        std::end(v));

	if (((!narm) && (v.size() < n)) || v.empty()) {
		return std::vector<double>(probs.size(), NAN);
	}

	n = v.size();
    std::sort(v.begin(), v.end());

	size_t pn = probs.size();
	std::vector<double> q(pn);

    for (size_t i = 0; i < pn; ++i) {
		double x = probs[i] * (double)(n-1);
		size_t x1 = (size_t)std::floor(x);
		size_t x2 = (size_t)std::ceil(x);
		if (x1 == x2) {
			q[i] = v[x1];
		} else {
			q[i] = interpolate(x, v[x1], v[x2], (double)x1, (double)x2);
		}
    }
    return q;
}



template <typename T>
std::vector<T> vunique(std::vector<T> d) {
	std::sort(d.begin(), d.end());
	d.erase(std::unique(d.begin(), d.end()), d.end());
	return d;
}

template <typename T>
std::vector<std::string> vtostring(std::vector<T>& v) {
	std::vector<std::string> s;
	std::transform(std::begin(v),
           std::end(v), std::back_inserter(s),
           [](double d) { return std::to_string(d); } 
        );
	return s;
}



template <typename T>
T vmedian(std::vector<T>& v, bool narm) {
	size_t n = v.size();
	std::vector<T> vv;
	vv.reserve(n);
	for (size_t i=0; i<n; i++) {
        if (!is_NA(v[i])) {
            vv.push_back(v[i]);
        } else if (!narm) {
            return NA<T>::value;
        }
	}
	n = vv.size();
	if (n == 0) {
		return(NA<T>::value);
	}
	if (n == 1) {
		return(vv[0]);
	}

	size_t n2 = n / 2;
	if (n % 2) {
		std::nth_element(vv.begin(), vv.begin()+n2, vv.end());
		return vv[n2];
	} else {
		std::sort(vv.begin(), vv.end());
		return (vv[n2] + vv[n2-1]) / 2;
	}
	
}



template <typename T>
T vsum(const std::vector<T>& v, bool narm) {
	T x = v[0];
	if (narm) {		
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(x)) {
				x = v[i];
			} else if (!is_NA(v[i])) {
				x += v[i];
			}
		}
	} else {
		if (is_NA(x)) {
			return(x);
		}
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(v[i])) {
				x = NA<T>::value;
				return(x);
			} else {
				x += v[i];
			}
		}
	}
	return x;
}

template <typename T>
T vsum2(const std::vector<T>& v, bool narm) {
	T x = v[0] * v[0];
	if (narm) {		
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(x)) {
				x = v[i] * v[i];
			} else if (!is_NA(v[i])) {
				x += v[i] * v[i];
			}
		}
	} else {
		if (is_NA(v[0])) {
			return(v[0]);
		}
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(v[i])) {
				x = NA<T>::value;
				return(x);
			} else {
				x += v[i] * v[i];
			}
		}
	}
	return x;
}


template <typename T>
T vprod(const std::vector<T>& v, bool narm) {
	T x = v[0];
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(x)) {
				x = v[i];
			} else if (!is_NA(v[i])) {
				x *= v[i];
			}
		}
	} else {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(x)) {
				if (is_NA(v[i])) {
					x = NA<T>::value;
					return(x);
				} else {
					x *= v[i];
				}
			}
		}
	}
	return x;
}



template <typename T>
double vmean(const std::vector<T>& v, bool narm) {
	double x = 0;
	unsigned d = 0;
	if (narm) {
		for (size_t i=0; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				x += v[i];
				d++;
			}
		}
	} else {
		for (size_t i=0; i<v.size(); i++) {
			if (is_NA(v[i])) {
				return(NAN);
			} else {
				x += v[i];
				d++;
			}
		}
	}
	if (d > 0) {
		x /= (double) d;
	} else {
		x = NAN;
	}
	return x;
}

template <typename T>
double vsd(const std::vector<T>& v, bool narm) {
	double m = vmean(v, narm);
	if (std::isnan(m)) return m;
	double x = 0;
	size_t n = 0;
	for (size_t i=0; i<v.size(); i++) {
		if (!is_NA(v[i])) {
			double d = (v[i] - m);
			x += (d * d);
			n++;
		}
	}
	n--;
	//if (n==0) return NAN;
	x = sqrt(x / n);
	return x;
}



template <typename T>
double vsdpop(const std::vector<T>& v, bool narm) {
	double m = vmean(v, narm);
	if (std::isnan(m)) return m;
	//double x = v[0];
	double x = 0;
	size_t n = 0;
	for (size_t i=0; i<v.size(); i++) {
		if (!is_NA(v[i])) {
			double d = (v[i] - m);
			x += d * d;
			n++;
		}
	}
	x = sqrt(x / n);
	return x;
}



template <typename T>
T vmin(const std::vector<T>& v, bool narm) {
	T x = v[0];
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (is_NA(x)) {
					x = v[i];
				} else {
					x = std::min(x, v[i]);
				}
			}
		}
	} else {
		if (is_NA(x)) return x;
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(v[i])) {
				return NA<T>::value;
			} else {
				x = std::min(x, v[i]);
			}
		}
	}
	return x;
}


template <typename T>
T vfirst(const std::vector<T>& v, bool narm) {
	if (narm) {
		for (size_t i=0; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				return v[i];
			}
		}
	} 
	return v[0];
}


template <typename T>
T vmax(const std::vector<T>& v, bool narm) {
	T x = v[0];
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (is_NA(x)) {
					x = v[i];
				} else {
					x = std::max(x, v[i]);
				}
			}
		}
	} else {
		if (is_NA(x)) return x;
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(v[i])) {
				return NA<T>::value;
			} else {
				x = std::max(x, v[i]);
			}
		}
	}
	return x;
}


template <typename T>
double vwhich(const std::vector<T>& v, bool narm) {
	double out;
	for (size_t i=0; i<v.size(); i++) {
		if ((!is_NA(v[i])) && v[i] != 0) {
			out = i+1;
			return out;
		}
	}
	out = NAN;
	return out;
}



template <typename T>
T vwhichmin(const std::vector<T>& v, bool narm) {
	T x = v[0];
	T out;
	if (is_NA(x)) {
		out = NA<T>::value;
	} else {
		out = 0;		
	}
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (is_NA(out)) {
					x = v[i];
					out = i;
				} else if (v[i] < x) {
					x = v[i];
					out = i;
				}
			}
		}
	} else {
		if (is_NA(x)) { return out; }
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(v[i])) {
				return NA<T>::value;
			} else {
				if (v[i] < x) {
					x = v[i];
					out = i;
				}
			}
		}
	}
	if (is_NA(out)) {
		return out;
	} else {
		return (out + 1);  // +1 for R
	}	
}


template <typename T>
T vwhichmax(const std::vector<T>& v, bool narm) {
	T x = v[0];
	T out;
	if (is_NA(x)) {
		out = NA<T>::value;
	} else {
		out = 0;		
	}
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (is_NA(out)) {
					x = v[i];
					out = i;
				} else if (v[i] > x) {
					x = v[i];
					out = i;
				}
			}
		}
	} else {
		if (is_NA(x)) { return out; }
		for (size_t i=0; i<v.size(); i++) {
			if (is_NA(v[i])) {
				return NA<T>::value;
			} else {
				if (v[i] > x) {
					x = v[i];
					out = i;
				}
			}
		}
	}
	if (is_NA(out)) {
		return out;
	} else {
		return (out + 1);  // +1 for R
	}	
}


// problematic; should be ok for int and float but
// won't work with bool values (nodata == 0)
template <typename T>
T vall(const std::vector<T>& v, bool narm) {
	T x;
	if (narm) {
		x = NA<T>::value;
        for (size_t i=0; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (v[i] == 0) {
					x = 0;
					break;
				} else {
					x = 1;					
				}
			}
        }
		//x = x < 0 ? NA<T>::value : x;
    } else {
		x = 1;
        for (size_t i=0; i<v.size(); i++) {
            if (is_NA(v[i]) || (v[i] == 0)) {
                x = v[i];
                break;
			}
		}
	}
	return x;
}



template <typename T>
T vany(const std::vector<T>& v, bool narm) {
	T x = NA<T>::value;
	x = 0;
	if (narm) {
		for (size_t i=0; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (v[i] != 0) {
					x = 1;
					break;
				} 
			}
		}
	} else {
		for (size_t i=0; i<v.size(); i++) {
			if (is_NA(v[i])) {
				x = NA<T>::value;
				break;				
			} else if (v[i] != 0) {
				x = 1;
			}
		}
	}
	return x;
}



template <typename T>
std::vector<T> vrange(const std::vector<T>& v, bool narm) {
	
	std::vector<T> x = { v[0], v[0] };

	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (is_NA(x[0])) {
					x[0] = v[i];
					x[1] = v[i];
				} else {
					x[0] = std::min(x[0], v[i]);
					x[1] = std::max(x[1], v[i]);
				}
			}
		}
	} else {
		if (is_NA(x[0])) { return(x); }
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(v[i])) {
				x[0] = NA<T>::value;
				x[1] = NA<T>::value;
				return(x);
			} else {
				x[0] = std::min(x[0], v[i]);
				x[1] = std::max(x[1], v[i]);
			}
		}
	}
	return x;
}



template <typename T>
T vmodal_old(std::vector<T>& v, bool narm) {

	size_t n = v.size();
    std::vector<unsigned> counts(n, 0);

	std::sort(v.begin(), v.end());

    for (size_t i = 0; i < n; ++i) {
        counts[i] = 0;
        size_t j = 0;
        while ((j < i) && (v[i] != v[j])) {
            ++j;
        }
        ++(counts[j]);
    }
	
    size_t maxCount = 0;
	for (size_t i = 1; i < n; ++i) {
		if (counts[i] > counts[maxCount]) {
			maxCount = i;
		}
	}
	
    return v[maxCount];
}


template <typename T>
T vmodal(std::vector<T>& v, bool narm) {

	if (narm) {
		std::map<double, size_t> count;
		for_each( v.begin(), v.end(), [&count]( double val ){
				if(!std::isnan(val)) count[val]++;
			}
		);
		if (count.size() == 0) return NAN;
		
		std::map<double, size_t>::iterator mode =	
			std::max_element(count.begin(), count.end(),[] (const std::pair<double, size_t>& a, 
			const std::pair<double, size_t>& b)->bool{ return a.second < b.second; } );
			
		return mode->first;

	}  else {
	
		std::map<double, size_t> count;
		for(size_t i=0; i<v.size(); i++) {
			if (std::isnan(v[i])) {
				return NAN;
			} else {
				count[v[i]]++;
			}
		}

		std::map<double, size_t>::iterator mode =	
			std::max_element(count.begin(), count.end(),[] (const std::pair<double, size_t>& a, 
			const std::pair<double, size_t>& b)->bool{ return a.second < b.second; } );
			
		return mode->first;
	}
	}





template <typename T>
std::vector<bool> visna(const std::vector<T>& v) {
	std::vector<bool> x(v.size(), false);
	for (size_t i=0; i<v.size(); i++) {
		if (is_NA(v[i])) {
			x[i] = true;
		}
	}
	return x;
}


template <typename T>
std::vector<bool> visnotna(const std::vector<T>& v) {
	std::vector<bool> x(v.size(), true);
	for (size_t i=0; i<v.size(); i++) {
		if (is_NA(v[i])) {
			x[i] = false;
		}
	}
	return x;
}


template <typename T>
bool vany_notna(const std::vector<T>& v) {
	for (size_t i=0; i<v.size(); i++) {
		if (!is_NA(v[i])) {
			return true;
		}
	}
	return false;
}

template <typename T>
bool vany_na(const std::vector<T>& v) {
	for (size_t i=0; i<v.size(); i++) {
		if (is_NA(v[i])) {
			return true;
		}
	}
	return false;
}


template <typename T>
void cumsum(std::vector<T>& v, bool narm) {
    if (narm) {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i])) {
                v[i] = v[i-1];
            } else if (!is_NA(v[i-1])){
                v[i] += v[i-1];
            }
        }
    } else {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i]) || is_NA(v[i-1])) {
                v[i] = NA<T>::value;
            } else {
                v[i] += v[i-1];
            }
        }
    }
}

template <typename T>
void cumprod(std::vector<T>& v, bool narm) {
    if (narm) {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i])) {
                v[i] = v[i-1];
            } else if (!is_NA(v[i-1])){
                v[i] *= v[i-1];
            }
        }
    } else {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i]) || is_NA(v[i-1])) {
                v[i] = NA<T>::value;
            } else {
                v[i] *= v[i-1];
            }
        }
    }
}


template <typename T>
void cummax(std::vector<T>& v, bool narm) {
    if (narm) {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i])) {
                v[i] = v[i-1];
            } else if (!is_NA(v[i-1])){
                v[i] = std::max(v[i], v[i-1]);
            }
        }
    } else {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i]) || is_NA(v[i-1])) {
                v[i] = NA<T>::value;
            } else {
                v[i] = std::max(v[i], v[i-1]);
            }
        }
    }
}


template <typename T>
void cummin(std::vector<T>& v, bool narm) {
    if (narm) {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i])) {
                v[i] = v[i-1];
            } else if (!is_NA(v[i-1])){
                v[i] = std::min(v[i], v[i-1]);
            }
        }
    } else {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i]) || is_NA(v[i-1])) {
                v[i] = NA<T>::value;
            } else {
                v[i] = std::min(v[i], v[i-1]);
            }
        }
    }
}

/*
#include <numeric>

template <typename T>
std::vector<size_t> order(const std::vector<T> &v) {
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);
	stable_sort(idx.begin(), idx.end(),
		[&v](size_t i, size_t j) {return v[i] < v[j];});
	return idx;
}
*/


template <typename T>
double expH(std::vector<T> d, bool narm) {
  
	std::unordered_map<T, unsigned int> counts;
	double s = 0;
	if (narm) {
		for (int v : d) {
			if (!is_NA(v)) {
				counts[v]++;
				s++;
			}
		}
	} else {
		for (int v : d) {
			if (is_NA(v)) {
				return NAN;
			} else {
				counts[v]++;
				s++;
			}
		}
	}
	if (s == 0) {
		return NAN;
	} 

	double sump = 0;
	for (auto const& pair : counts) { 
		double p = pair.second / s;
		sump += p * log(p);
	}
  
	return exp(-sump);
}

#endif

