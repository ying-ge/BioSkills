// Copyright (c) 2018-2025  Robert J. Hijmans
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

#include "spatRaster.h"
#include "distance.h"
#include <limits>
#include <cmath>
#include "geodesic.h"
#include "recycle.h"
#include "math_utils.h"
#include "vecmath.h"
#include "file_utils.h"
#include "string_utils.h"
#include "crs.h"
#include "sort.h"
#include "geosphere.h"

/*
inline void shortDistPoints(std::vector<double> &d, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &px, const std::vector<double> &py, const bool& lonlat, const std::string& method, const double &lindist) {
	if (lonlat) {
		distanceToNearest_lonlat(d, x, y, px, py, lindist, method);
	} else {
		distanceToNearest_plane(d, x, y, px, py, lindist);
	}
}
*/

inline void shortDirectPoints(std::vector<double> &d, std::vector<double> &x, std::vector<double> &y, std::vector<double> &px, std::vector<double> &py, const bool& lonlat, bool &from, bool &degrees, const std::string &method) {
	if (lonlat) {
		directionToNearest_lonlat(d, x, y, px, py, degrees, from, method);
	} else {
		directionToNearest_plane(d, x, y, px, py, degrees, from);
	}
}


/*
bool get_m(double &m, SpatSRS srs, bool lonlat, std::string unit) {
	m = 1;
	if (!lonlat) {
		m = srs.to_meter();
		m = std::isnan(m) ? 1 : m;
	}
	std::vector<std::string> ss {"m", "km"};
	if (std::find(ss.begin(), ss.end(), unit) == ss.end()) {
		return false;
	}
	if (unit == "km")	{
		m /= 1000;
	}
	return true;
}
*/

std::vector<double> dist_bounds(const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& rx, const double& ry, size_t& first, size_t& last, const bool& lonlat, const std::string &method) {

	std::vector<double> d(rx.size(), std::numeric_limits<double>::max());
	size_t oldfirst = first;
	first = vx.size();
	last = 0;


	if (lonlat) {
		std::function<double(double, double, double, double)> dfun;
		if (method == "haversine") {
			dfun = distance_hav;
		} else if (method == "cosine") {
			dfun = distance_cos;			
		} else {
			dfun = distance_geo;
		}

		
		for (size_t i=0; i<rx.size(); i++) {
			size_t thisone = 0;
			for (size_t j=oldfirst; j<vx.size(); j++) {
				double dd = dfun(rx[i], ry, vx[j], vy[j]);
				if (dd < d[i]) {
					d[i] = dd;
					thisone = j;
				}
			}
			first = std::min(thisone, first);
			last  = std::max(thisone, last);
		}
	} else {
		for (size_t i=0; i<rx.size(); i++) {
			size_t thisone = 0;
			for (size_t j=oldfirst; j<vx.size(); j++) {
				double dd = distance_plane(rx[i], ry, vx[j], vy[j]);
				if (dd < d[i]) {
					d[i] = dd;
					thisone = j;
				}
			}
			first = std::min(thisone, first);
			last  = std::max(thisone, last);
		}
	}
	last += 1;
	return d;
}

void dist_only(std::vector<double> &d, const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& rx, const std::vector<double>& ry, const size_t& first, const size_t& last, const bool& lonlat, const std::vector<double>& dlast, bool skip, const std::vector<double>& v, const std::string& method, bool setNA) {

	size_t rxs = rx.size();
	d.reserve(rxs + dlast.size());

	double inf = std::numeric_limits<double>::infinity();

	if (lonlat) {

		if (method == "geo") {

			double dd, azi1, azi2;
			struct geod_geodesic g;
			// get a and f from crs?
			double a = 6378137.0;
			double f = 1/298.257223563;
			geod_init(&g, a, f);

			if (skip) {
				for (size_t i=0; i<rxs; i++) {
					if (std::isnan(v[i])) {
						d.push_back(inf);
						for (size_t j=first; j<last; j++) {
							geod_inverse(&g, ry[i], rx[i], vy[j], vx[j], &dd, &azi1, &azi2);
							if (dd < d[i]) {
								d[i] = dd;
							}
						}
					} else {
						d.push_back(0);
					}
				}
			} else { // no skip

				for (size_t i=0; i<rxs; i++) {
					d.push_back(inf);
					for (size_t j=first; j<last; j++) {
						geod_inverse(&g, ry[i], rx[i], vy[j], vx[j], &dd, &azi1, &azi2);
						if (dd < d[i]) {
							d[i] = dd;
						}
					}
				} 
			}

		} else { // not geo
	
			std::function<double(double, double, double, double)> dfun;
			if (method == "haversine") {
				dfun = distance_hav;
			} else if (method == "cosine") {
				dfun = distance_cos;
			} 
			if (skip) {
				
				for (size_t i=0; i<rxs; i++) {
					if (std::isnan(v[i])) {
						d.push_back(inf);
						for (size_t j=first; j<last; j++) {
							double dd = dfun(rx[i], ry[i], vx[j], vy[j]);
							if (dd < d[i]) {
								d[i] = dd;
							}
						}
					} else {
						d.push_back(0);
					}
				}		
			} else {	
				for (size_t i=0; i<rxs; i++) {
					d.push_back(inf);
					for (size_t j=first; j<last; j++) {
						//double dd = distHaversine(rx[i], ry[i], vx[j], vy[j]);
						double dd = dfun(rx[i], ry[i], vx[j], vy[j]);

						if (dd < d[i]) {
							d[i] = dd;
						}
					}
				}
			} 
		}
		
	} else { // not lonlat
		if (skip) {
			for (size_t i=0; i<rxs; i++) {
				if (std::isnan(v[i])) {
					d.push_back(inf);
					for (size_t j=first; j<last; j++) {
						double dd = distance_plane(rx[i], ry[i], vx[j], vy[j]);
						if (dd < d[i]) {
							d[i] = dd;
						}
					}
				} else {
					d.push_back(0);
				}
			}
		} else {
			for (size_t i=0; i<rxs; i++) {
				d.push_back(inf);
				for (size_t j=first; j<last; j++) {
					double dd = distance_plane(rx[i], ry[i], vx[j], vy[j]);
					if (dd < d[i]) {
						d[i] = dd;
					}
				}
			}
		}
	}

	d.insert(d.end(), dlast.begin(), dlast.end());
	if (skip) {
		for (size_t i=rxs; i< v.size(); i++) {
			if (!std::isnan(v[i])) {
				d[i] = 0;
			}
		}
		if (setNA) {
			double mxval = std::numeric_limits<double>::max();
			for (size_t i=0; i< v.size(); i++) {
				if (v[i] == mxval) {
					d[i] = NAN;
				}
			}
		}
	}
}


SpatRaster SpatRaster::distance_crds(std::vector<double>& x, std::vector<double>& y, const std::string& method, bool skip, bool setNA, std::string unit, double max_dist, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (x.empty()) {
		out.setError("no locations to compute distance from");
		return(out);
	}
	const double toRad = 0.0174532925199433;
	std::vector<std::size_t> pm = sort_order_d(y);
	permute(x, pm);
	permute(y, pm);

	bool lonlat = is_lonlat(); 

	double m=1;
	if (!source[0].srs.m_dist(m, lonlat, unit)) {
		out.setError("invalid unit");
		return(out);
	}

	size_t nc = ncol();
	if (nrow() > 1000) {
		opt.steps = std::max(opt.steps, (size_t) 4);
		opt.progress = opt.progress * 2;
	}
	
 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	std::vector<double> cells;
	std::vector<double> dlast;

	std::vector<int64_t> cols;
	cols.resize(ncol());
	std::iota(cols.begin(), cols.end(), 0);
	std::vector<double> tox = xFromCol(cols);

	if (lonlat && (method != "geo")) {
		for (double &d : x) d *= toRad;
		for (double &d : y) d *= toRad;
		for (double &d : tox) d *= toRad;
	}

	size_t oldfirst = 0;
	size_t first = 0;
	size_t last  = x.size();

	std::vector<double> v;
	if (skip) {
		if (!readStart()) {
			out.setError(getError());
			return(out);
		}
		for (size_t i = 0; i < out.bs.n; i++) {
			cells.resize((out.bs.nrows[i] -1) * nc) ;
			std::iota(cells.begin(), cells.end(), out.bs.row[i] * nc);
			std::vector<std::vector<double>> rxy = xyFromCell(cells);
			double toy = yFromRow(out.bs.row[i] + out.bs.nrows[i] - 1);
			if (lonlat && (method != "geo")) {
				toy *= toRad;
				for (double &d : rxy[0]) d *= toRad;
				for (double &d : rxy[1]) d *= toRad;
			}
			readBlock(v, out.bs, i);
			dlast = dist_bounds(x, y, tox, toy, first, last, lonlat, method);
			std::vector<double> d;
			dist_only(d, x, y, rxy[0], rxy[1], oldfirst, last, lonlat, dlast, true, v, method, setNA);
			oldfirst = first;
			if (m != 1) {
				for (double &v : d) v *= m;
			}
			if (max_dist > 0) {
				for (double &v : d) v = v > max_dist ? NAN : v;
			}
			if (!out.writeBlock(d, i)) return out;
		}
		readStop();
	} else {
		for (size_t i = 0; i < out.bs.n; i++) {
			double toy = yFromRow(out.bs.row[i] + out.bs.nrows[i] - 1);
			cells.resize((out.bs.nrows[i] -1) * nc) ;
			std::iota(cells.begin(), cells.end(), out.bs.row[i] * nc);
			std::vector<std::vector<double>> rxy = xyFromCell(cells);
			if (lonlat && (method != "geo")) {
				toy *= toRad;
				for (double &d : rxy[0]) d *= toRad;
				for (double &d : rxy[1]) d *= toRad;
			}
			dlast = dist_bounds(x, y, tox, toy, first, last, lonlat, method);
			std::vector<double> d;
			dist_only(d, x, y, rxy[0], rxy[1], oldfirst, last, lonlat, dlast, false, v, method, setNA);
			oldfirst = first;
			if (m != 1) {
				for (double &v : d) v *= m;
			}
			if (max_dist > 0) {
				for (double &v : d) v = v > max_dist ? NAN : v;
			}
			if (!out.writeBlock(d, i)) return out;
		}
	}
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::distance_vector(SpatVector p, bool rasterize, std::string unit, const std::string& method, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (source[0].srs.wkt.empty()) {
		out.setError("CRS not defined");
		return(out);
	}
	if (!source[0].srs.is_same(p.srs, false)) {
		out.setError("CRS do not match");
		return(out);
	}
	if (p.empty()) {
		out.setError("no locations to compute distance from");
		return(out);
	}
	bool lonlat = is_lonlat(); 

	if ((unit != "m") && (unit != "km")) {
		out.setError("invalid unit. Must be 'm' or 'km'");
		return(out);
	}
	if ((method != "geo") && (method != "cosine") && (method != "haversine")) {
		out.setError("invalid method. Must be 'geo', 'haversine' or 'cosine'");
		return(out);
	}

	if (rasterize) {
//SpatRaster SpatRaster::distance_rasterize(SpatVector p, double target, double exclude, std::string unit, const std::string& method, SpatOptions &opt) {
//		double target = NAN;

		double exclude = NAN;

		SpatRaster x;
		SpatOptions ops(opt);
		std::string gtype = p.type();
		bool poly = gtype == "polygons";

		x = out.rasterize(p, "", {1}, NAN, false, "", false, false, false, ops);

		if (!lonlat) {
			return x.distance(NAN, 0, false, unit, false, method, false, -1, opt);
		}

		if (poly) {
			x  = x.edges(false, "inner", 8, 0, ops);
			SpatRaster xp = x.replaceValues({0}, {exclude}, 1, false, NAN, false, ops);
			p  = xp.as_points(false, true, false, opt);
		} else {
	//		x = x.edges(false, "inner", 8, NAN, ops);
			p = x.as_points(false, true, false, opt);
		}

		std::vector<std::vector<double>> pxy = p.coordinates();

		if (pxy.empty()) {
			out.setError("no locations to compute from");
			return(out);
		}


		bool setNA = false;
		out = x.distance_crds(pxy[0], pxy[1], method, poly, setNA, unit, -1, opt);

	} else {

		if ((p.type() == "polygons") || (p.type() == "lines")) {
			
			if (p.nrow() > 1) {
				p = p.aggregate(true);
			}
			
			std::vector<double> cells;
			size_t nc = ncol();

			if (!readStart()) {
				out.setError(getError());
				return(out);
			}
			if (!out.writeStart(opt, filenames())) {
				readStop();
				return out;
			}
			for (size_t i = 0; i < out.bs.n; i++) {
				cells.resize(out.bs.nrows[i] * nc) ;
				std::iota(cells.begin(), cells.end(), out.bs.row[i] * nc);
				std::vector<std::vector<double>> rxy = xyFromCell(cells);

				SpatVector pnts;
				pnts.srs = source[0].srs;
				pnts.setPointsGeometry(rxy[0], rxy[1]);
				std::vector<double> d = pnts.distance(p, false, unit, method, false, opt);
					
				if (!out.writeBlock(d, i)) return out;
			}
			readStop();
			out.writeStop();
			
		} else {
		//p = p.aggregate(false);
			if ((method != "geo") && (method != "cosine") && (method != "haversine")) {
				out.setError("invalid method. Must be 'geo', 'cosine' or 'haversine'");
				return(out);
			}

			std::vector<std::vector<double>> pxy = p.coordinates();
			SpatOptions ops(opt);
			bool setNA = false;
			out = distance_crds(pxy[0], pxy[1], method, false, setNA, unit, -1, opt);
		}
	}
	return out;

}


SpatRaster SpatRaster::direction_rasterize(SpatVector p, bool from, bool degrees, double target, double exclude, const std::string &method, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (source[0].srs.wkt.empty()) {
		out.setError("CRS not defined");
		return(out);
	}
	if (!source[0].srs.is_same(p.srs, false)) {
		out.setError("CRS do not match");
		return(out);
	}
	bool lonlat = is_lonlat(); 

	SpatRaster x;
	SpatOptions ops(opt);
	std::string gtype = p.type();
	bool poly = gtype == "polygons";

	x = out.rasterize(p, "", {1}, NAN, false, "", false, false, false, ops);


	if (poly) {
		x  = x.edges(false, "inner", 8, 0, ops);
		SpatRaster xp = x.replaceValues({0}, {exclude}, 1, false, NAN, false, ops);
		p  = xp.as_points(false, true, false, opt);
	} else {
//		x = x.edges(false, "inner", 8, NAN, ops);
		p = x.as_points(false, true, false, opt);
	}

	std::vector<std::vector<double>> pxy = p.coordinates();

	if (pxy.empty()) {
		out.setError("no locations to compute from");
		return(out);
	}
	size_t nc = ncol();
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v;
		std::vector<double> cells(out.bs.nrows[i] * nc) ;
		std::vector<double> vals;
		vals.resize(out.bs.nrows[i] * nc, NAN);

		std::iota(cells.begin(), cells.end(), out.bs.row[i] * nc);

		x.readBlock(v, out.bs, i);
		if (std::isnan(target)) {
			for (size_t j=0; j<v.size(); j++) {
				if (!std::isnan(v[j])) {
					cells[j] = NAN;
				}
			}
		} else {
			for (size_t j=0; j<v.size(); j++) {
				if (v[j] != target) {
					cells[j] = NAN;
					if (std::isnan(v[j])) {
						vals[j] = NAN;
					}
				}
			}
		}
		std::vector<std::vector<double>> xy = xyFromCell(cells);
		shortDirectPoints(vals, xy[0], xy[1], pxy[0], pxy[1], lonlat, from, degrees, method);
		if (!out.writeBlock(vals, i)) return out;
	}

	out.writeStop();
	readStop();
	return(out);
}


/*
SpatRaster SpatRaster::distance_vector(SpatVector p, std::string unit, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (source[0].srs.wkt == "") {
		out.setError("CRS not defined");
		return(out);
	}
	if (!source[0].srs.is_same(p.srs, false) ) {
		out.setError("CRS does not match");
		return(out);
	}

	bool lonlat = is_lonlat(); 
	double m=1;
//	if (!get_m(m, source[0].srs, lonlat, unit)) {
	if (!source[0].srs.m_dist(m, lonlat, unit)) {
		
		out.setError("invalid unit");
		return(out);
	}

	if (p.size() == 0) {
		out.setError("no locations to compute distance from");
		return(out);
	}
	p = p.aggregate(false);

//	bool lonlat = is_lonlat(); // m == 0
	size_t nc = ncol();

 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	std::vector<double> cells;

	for (size_t i = 0; i < out.bs.n; i++) {
		double s = out.bs.row[i] * nc;
		cells.resize(out.bs.nrows[i] * nc) ;
		std::iota(cells.begin(), cells.end(), s);
		std::vector<std::vector<double>> xy = xyFromCell(cells);
		SpatVector pv(xy[0], xy[1], points, "");
		pv.srs = source[0].srs;
		std::vector<double> d = p.distance(pv, false, unit);
		if (p.hasError()) {
			out.setError(p.getError());
			out.writeStop();
			return(out);
		}
		if (m != 1) {
			for (double &v : d) v *= m;
		}
		if (!out.writeBlock(d, i)) return out;
	}
	out.writeStop();
	return(out);
}

*/

SpatRaster SpatRaster::distance(double target, double exclude, bool keepNA, std::string unit, bool remove_zero, const std::string method, bool values, double threshold, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return out;
	}

	SpatOptions ops(opt);
	size_t nl = nlyr();
	if (nl > 1) {
		std::vector<std::string> nms = getNames();
		if (ops.names.size() == nms.size()) {
			nms = opt.names;
		}
		out.source.resize(nl);
		for (size_t i=0; i<nl; i++) {
			std::vector<size_t> lyr = {i};
			SpatRaster r = subset(lyr, ops);
			ops.names = {nms[i]};
			r = r.distance(target, exclude, keepNA, unit, remove_zero, method, values, threshold, ops);
			out.source[i] = r.source[0];
		}
		if (!opt.get_filename().empty()) {
			out = out.writeRaster(opt);
		}
		return out;
	}
	if (!(values || is_lonlat())) { // && std::isnan(target) && std::isnan(exclude)) {
		return proximity(target, exclude, keepNA, unit, false, threshold, remove_zero, opt); 
	}

	bool setNA = false;
	std::vector<std::vector<double>> p;
	if (!std::isnan(exclude)) {
		SpatRaster x;
		if (std::isnan(target)) {
			x = replaceValues({exclude}, {target}, 1, false, NAN, false, ops);
			x = x.edges(false, "inner", 8, 1, ops);
			p = x.as_points_value(1, ops);
			if (p.empty()) {
				return out.init({0}, opt);
			}
			if (values) {
				std::vector<std::vector<double>> vv = extractXY(p[0], p[1], "", false, opt);
				return distance_crds_vals(p[0], p[1], vv[0], method, true, setNA, unit, threshold, opt);				
			} else {
				return distance_crds(p[0], p[1], method, true, setNA, unit, threshold, opt);
			}
		} else {
			x = replaceValues({exclude, target}, {NAN, NAN}, 1, false, NAN, false, ops);
			x = x.edges(false, "inner", 8, 1, ops);
			p = x.as_points_value(1, ops);
			out = replaceValues({NAN, exclude, target}, {target, NAN, NAN}, 1, false, NAN, false, ops);
		}
	} else if (!std::isnan(target)) {
		SpatRaster x = replaceValues({target}, {NAN}, 1, false, NAN, false, ops);
		x = x.edges(false, "inner", 8, 0, ops);
		p = x.as_points_value(1, ops);
		out = replaceValues({NAN, target}, {std::numeric_limits<double>::max(), NAN}, 1, false, NAN, false, ops);
		setNA = true;
	} else {
		out = edges(false, "inner", 8, 0, ops);
		p = out.as_points_value(1, ops);
	}
	if (p.empty()) {
		return out.init({0}, opt);
	}
	if (values) {
		std::vector<std::vector<double>> vv = extractXY(p[0], p[1], "", false, opt);
		out = out.distance_crds_vals(p[0], p[1], vv[0], method, true, setNA, unit, threshold, opt);				
	} else {
		out = out.distance_crds(p[0], p[1], method, true, setNA, unit, threshold, opt);
	}
	return out;
}



SpatRaster SpatRaster::direction(bool from, bool degrees, double target, double exclude, const std::string &method, SpatOptions &opt) {
	SpatRaster out = geometry(1);
	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return out;
	}

	SpatOptions ops(opt);
	size_t nl = nlyr();
	if (nl > 1) {
		out.source.resize(nl);
		std::vector<std::string> nms = getNames();
		if (ops.names.size() == nms.size()) {
			nms = opt.names;
		}
		for (size_t i=0; i<nl; i++) {
			std::vector<size_t> lyr = {i};
			SpatRaster r = subset(lyr, ops);
			ops.names = {nms[i]};
			r = r.direction(from, degrees, target, exclude, method, ops);
			out.source[i] = r.source[0];
		}
		if (!opt.get_filename().empty()) {
			out = out.writeRaster(opt);
		}
		return out;
	}

	if (!std::isnan(exclude)) {
		SpatOptions xopt(opt);
		SpatRaster x = replaceValues({exclude}, {NAN}, 1, false, NAN, false, xopt);
		out = x.edges(false, "inner", 8, target, ops);
	} else {
		out = edges(false, "inner", 8, target, ops);
	}
	SpatVector p = out.as_points(false, true, false, opt);
	if (p.empty()) {
		out.setError("no cells to compute direction from or to");
		return(out);
	}
	return direction_rasterize(p, from, degrees, target, exclude, method, opt);
}






inline double minCostDist(std::vector<double> &d) {
	d.erase(std::remove_if(d.begin(), d.end(),
		[](const double& v) { return std::isnan(v); }), d.end());
	std::sort(d.begin(), d.end());
	return d.empty() ? NAN : d[0];
}


inline void DxDxyCost(const double &lat, const size_t &row, double xres, double yres, const double &dir, double &dx,  double &dy, double &dxy, double distscale, const double mult=2) {
	double rlat = lat + (double)row * yres * dir;
	dx  = distance_lonlat(0, rlat, xres, rlat) / (mult * distscale);
	yres *= -dir;
	dy  = distance_lonlat(0, 0, 0, yres);
	dxy = distance_lonlat(0, rlat, xres, rlat+yres);
	dy = std::isnan(dy) ? NAN : dy / (mult * distscale);
	dxy = std::isnan(dxy) ? NAN : dxy / (mult * distscale);
}


void cost_dist(std::vector<double> &dist, std::vector<double> &dabove, std::vector<double> &v, std::vector<double> &vabove, std::vector<double> res, size_t nr, size_t nc, double lindist, bool geo, double lat, double latdir, bool global, bool npole, bool spole) {

	std::vector<double> cd;

	double dx, dy, dxy;
	if (geo) {
		DxDxyCost(lat, 0, res[0], res[1], latdir, dx, dy, dxy, lindist);
	} else {
		dx = res[0] * lindist / 2;
		dy = res[1] * lindist / 2;
		dxy = sqrt(dx*dx + dy*dy);
	}

	//top to bottom
    //left to right
	//first cell, no cell left of it
	if (!std::isnan(v[0])) {
		if (global) {
			cd = {dist[0], dabove[0] + (v[0]+vabove[0]) * dy,
				dist[nc-1] + (v[0] + v[nc-1]) * dx,
				dabove[nc-1] + dxy * (vabove[nc-1]+v[0])};
		} else {
			cd = {dist[0], dabove[0] + (v[0]+vabove[0]) * dy};
		}
		dist[0] = minCostDist(cd);
	}
	for (size_t i=1; i<nc; i++) { //first row, no row above it, use "above"
		if (!std::isnan(v[i])) {
			cd = {dist[i], dabove[i]+(vabove[i]+v[i])*dy, dabove[i-1]+(vabove[i-1]+v[i])*dxy, dist[i-1]+(v[i-1]+v[i])*dx};
			dist[i] = minCostDist(cd);
		}
	}
	if (npole) {
		double minp = *std::min_element(dist.begin(), dist.begin()+nc);
		minp += dy;
		for (size_t i=0; i<nc; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}


	for (size_t r=1; r<nr; r++) { //other rows
		if (geo) DxDxyCost(lat, r, res[0], res[1], latdir, dx, dy, dxy, lindist);
		size_t start=r*nc;
		if (!std::isnan(v[start])) {
			if (global) {
				cd = {dist[start-nc] + (v[start] + v[start-nc]) * dy, dist[start],
					dist[start+nc-1] + (v[start] + v[start+nc-1]) * dx,
					dist[start-1] + (v[start] + v[start-1]) * dxy};
			} else {
				cd = {dist[start-nc] + (v[start] + v[start-nc]) * dy, dist[start]};
			}
			dist[start] = minCostDist(cd);
		}
		size_t end = start+nc;
		for (size_t i=(start+1); i<end; i++) {
			if (!std::isnan(v[i])) {
				cd = {dist[i], dist[i-1]+(v[i]+v[i-1])*dx, dist[i-nc]+(v[i]+v[i-nc])*dy, dist[i-nc-1]+(v[i]+v[i-nc-1])*dxy};
				dist[i] = minCostDist(cd);
			}
		}
	}
	if (spole) {
		double minp = *std::min_element(dist.end()-nc, dist.end());
		minp += dy;
		size_t ds = dist.size();
		for (size_t i=ds-nc; i<ds; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	//right to left
	// first row, no need for first (last) cell (unless is global)
	if (geo) DxDxyCost(lat, 0, res[0], res[1], latdir, dx, dy, dxy, lindist);
	if (global) {
		size_t i=(nc-1);
		cd = {dist[i],
			dist[0] + (v[0] + v[i]) * dx,
			dabove[0] + dxy * (vabove[0]+v[i])};
		dist[i] = minCostDist(cd);
	}
	if (npole) {
		double minp = *std::min_element(dist.begin(), dist.begin()+nc);
		minp += dy;
		for (size_t i=0; i<nc; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	for (int i=(nc-2); i > -1; i--) { // other cells on first row
		if (!std::isnan(v[i])) {
			cd = {dabove[i]+(vabove[i]+v[i])*dy, dabove[i+1]+(vabove[i+1]+v[i])*dxy, dist[i+1]+(v[i+1]+v[i])*dx, dist[i]};
			dist[i] = minCostDist(cd);
		}
	}

	for (size_t r=1; r<nr; r++) { // other rows
		if (geo) DxDxyCost(lat, r, res[0], res[1], latdir, dx, dy, dxy, lindist);
		size_t start=(r+1)*nc-1;

		if (!std::isnan(v[start])) {
			if (global) {
				cd = { dist[start], dist[start-nc] + (v[start-nc]+v[start])* dy,
					dist[start-nc+1] + (v[start-nc+1] + v[start]) * dx,
					dist[start-(2*nc)+1] + (v[start-(2*nc)+1] + v[start]) * dxy
				};

			} else {
				cd = { dist[start], dist[start-nc] + (v[start-nc]+v[start])* dy };
			}
			dist[start] = minCostDist(cd);
		}

		size_t end=r*nc;
		start -= 1;
		for (size_t i=start; i>=end; i--) {
			if (!std::isnan(v[i])) {
				cd = { dist[i+1]+(v[i+1]+v[i])*dx, dist[i-nc+1]+(v[i]+v[i-nc+1])*dxy, dist[i-nc]+(v[i]+v[i-nc])*dy, dist[i]};
				dist[i] = minCostDist(cd);
			}
		}
	}
	if (spole) {
		double minp = *std::min_element(dist.end()-nc, dist.end());
		minp += dy;
		size_t ds = dist.size();
		for (size_t i=ds-nc; i<ds; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	size_t off = (nr-1) * nc;
	dabove = std::vector<double>(dist.begin()+off, dist.end());
	vabove = std::vector<double>(v.begin()+off, v.end());

}


void grid_dist(std::vector<double> &dist, std::vector<double> &dabove, std::vector<double> &v, std::vector<double> &vabove, std::vector<double> res, size_t nr, size_t nc, double lindist, bool geo, double lat, double latdir, bool global, bool npole, bool spole) {

	std::vector<double> cd;

	double dx, dy, dxy;
	if (geo) {
		DxDxyCost(lat, 0, res[0], res[1], latdir, dx, dy, dxy, lindist, 1);
	} else {
		dx = res[0] * lindist;
		dy = res[1] * lindist;
		dxy = sqrt(dx*dx + dy*dy);
	}

	//top to bottom
    //left to right
	//first cell, no cell left of it
	if (!std::isnan(v[0])) {
		if (global) {
			cd = {dist[0], dabove[0] + dy,
				dist[nc-1] + dx,
				dabove[nc-1] + dxy};
		} else {
			cd = {dist[0], dabove[0] + dy};
		}
		dist[0] = minCostDist(cd);
	}
	for (size_t i=1; i<nc; i++) { //first row, no row above it, use "above"
		if (!std::isnan(v[i])) {
			cd = {dist[i], dabove[i]+dy, dabove[i-1]+dxy, dist[i-1]+dx};
			dist[i] = minCostDist(cd);
		}
	}
	if (npole) {
		double minp = *std::min_element(dist.begin(), dist.begin()+nc);
		minp += dy;
		for (size_t i=0; i<nc; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}


	for (size_t r=1; r<nr; r++) { //other rows
		if (geo) DxDxyCost(lat, r, res[0], res[1], latdir, dx, dy, dxy, lindist, 1);
		size_t start=r*nc;
		if (!std::isnan(v[start])) {
			if (global) {
				cd = {dist[start-nc] + dy, dist[start],
					dist[start+nc-1] + dx,
					dist[start-1] + dxy};
			} else {
				cd = {dist[start-nc] + dy, dist[start]};
			}
			dist[start] = minCostDist(cd);
		}
		size_t end = start+nc;
		for (size_t i=(start+1); i<end; i++) {
			if (!std::isnan(v[i])) {
				cd = {dist[i], dist[i-1]+dx, dist[i-nc]+dy, dist[i-nc-1]+dxy};
				dist[i] = minCostDist(cd);
			}
		}
	}
	if (spole) {
		double minp = *std::min_element(dist.end()-nc, dist.end());
		minp += dy;
		size_t ds = dist.size();
		for (size_t i=ds-nc; i<ds; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	//right to left
	// first row, no need for first (last) cell (unless is global)
	if (geo) DxDxyCost(lat, 0, res[0], res[1], latdir, dx, dy, dxy, lindist, 1);
	if (global) {
		size_t i=(nc-1);
		cd = {dist[i],
			dist[0] + dx,
			dabove[0] + dxy};
		dist[i] = minCostDist(cd);
	}
	if (npole) {
		double minp = *std::min_element(dist.begin(), dist.begin()+nc);
		minp += dy;
		for (size_t i=0; i<nc; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	for (int i=(nc-2); i > -1; i--) { // other cells on first row
		if (!std::isnan(v[i])) {
			cd = {dabove[i]+dy, dabove[i+1]+dxy, dist[i+1]+dx, dist[i]};
			dist[i] = minCostDist(cd);
		}
	}

	for (size_t r=1; r<nr; r++) { // other rows
		if (geo) DxDxyCost(lat, r, res[0], res[1], latdir, dx, dy, dxy, lindist, 1);
		size_t start=(r+1)*nc-1;

		if (!std::isnan(v[start])) {
			if (global) {
				cd = { dist[start], dist[start-nc] + dy,
					dist[start-nc+1] + dx,
					dist[start-(2*nc)+1] + dxy
				};

			} else {
				cd = { dist[start], dist[start-nc] + dy };
			}
			dist[start] = minCostDist(cd);
		}

		size_t end=r*nc;
		start -= 1;
		for (size_t i=start; i>=end; i--) {
			if (!std::isnan(v[i])) {
				cd = { dist[i+1]+dx, dist[i-nc+1]+dxy, dist[i-nc]+dy, dist[i]};
				dist[i] = minCostDist(cd);
			}
		}
	}
	if (spole) {
		double minp = *std::min_element(dist.end()-nc, dist.end());
		minp += dy;
		size_t ds = dist.size();
		for (size_t i=ds-nc; i<ds; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	size_t off = (nr-1) * nc;
	dabove = std::vector<double>(dist.begin()+off, dist.end());
	vabove = std::vector<double>(v.begin()+off, v.end());
}

void block_is_same(bool& same, std::vector<double>& x,  std::vector<double>& y) {
	if (!same) return;
	for (size_t i=0; i<x.size(); i++) {
		if (!std::isnan(x[i]) && (x[i] != y[i])) {
			same = false;
			break;
		}
	}
}






SpatRaster SpatRaster::costDistanceRun(SpatRaster &old, bool &converged, double target, double m, bool lonlat, bool global, bool npole, bool spole, bool grid, SpatOptions &opt) {

	std::vector<double> res = resolution();

	SpatRaster first = geometry();
	SpatRaster second = first;
    std::vector<double> d, v, vv;
	if (!readStart()) {
		first.setError(getError());
		return(first);
	}
	opt.progressbar = false;
 	if (!first.writeStart(opt, filenames())) { return first; }

	size_t nc = ncol();
	std::vector<double> dabove(nc, NAN);
	std::vector<double> vabove(nc, 0);
	double lat = 0;
	if (old.hasValues()) {
		if (!old.readStart()) {
			first.setError(getError());
			return(first);
		}
		if (!first.writeStart(opt, filenames())) {
			readStop();
			old.readStop();
			return first;
		}

		for (size_t i = 0; i < first.bs.n; i++) {
			readBlock(v, first.bs, i);
			if (lonlat) {
				lat = yFromRow(first.bs.row[i]);
			}
			bool np = (i==0) && npole;
			bool sp = (i==first.bs.n-1) && spole;
			if (target != 0) {
				for (size_t j=0; j<v.size(); j++) {
					if (v[j] == target) {
						v[j] = 0;
					}
				}
			}
			old.readBlock(d, first.bs, i);
			if (grid) {
				grid_dist(d, dabove, v, vabove, res, first.bs.nrows[i], nc, m, lonlat, lat, -1, global, np, sp);
			} else {
				cost_dist(d, dabove, v, vabove, res, first.bs.nrows[i], nc, m, lonlat, lat, -1, global, np, sp);
			}
			if (!first.writeValues(d, first.bs.row[i], first.bs.nrows[i])) return first;
		}
	} else {
		converged = false;
		for (size_t i = 0; i < first.bs.n; i++) {
			if (lonlat) {
				lat = yFromRow(first.bs.row[i]);
			}
			bool np = (i==0) && npole;
			bool sp = (i==first.bs.n-1) && spole;
			readBlock(v, first.bs, i);
			d.clear();
			d.resize(v.size(), NAN);
			for (size_t j = 0; j < v.size(); j++) {
				if (v[j] == target) {
					v[j] = 0;
					d[j] = 0;
				} else if ((!grid) && (v[j] < 0)) {
					readStop();
					first.writeStop();
					first.setError("negative friction values not allowed");
					return first;
				}
			}
			if (grid) {
				grid_dist(d, dabove, v, vabove, res, first.bs.nrows[i], nc, m, lonlat, lat, -1, global, np, sp);
			} else {
				cost_dist(d, dabove, v, vabove, res, first.bs.nrows[i], nc, m, lonlat, lat, -1, global, np, sp);
			}
			if (!first.writeValuesRect(d, first.bs.row[i], first.bs.nrows[i], 0, nc)) return first;
		}
	}
	first.writeStop();

	if (!first.readStart()) {
		return(first);
	}

	dabove = std::vector<double>(nc, NAN);
	vabove = std::vector<double>(nc, 0);
  	if (!second.writeStart(opt, filenames())) {
		readStop();
		first.readStop();
		return second;
	}
	for (int i = second.bs.n; i>0; i--) {
		if (lonlat) {
			lat = yFromRow(second.bs.row[i-1] + second.bs.nrows[i-1] - 1);
		}
		bool sp = (i==1) && spole; //! reverse order
		bool np = (i==(int)second.bs.n) && npole;
        readBlock(v, second.bs, i-1);
		if (target != 0) {
			for (size_t j=0; j<v.size(); j++) {
				if (v[j] == target) {
					v[j] = 0;
				}
			}
		}
		first.readBlock(d, second.bs, i-1);
		std::reverse(v.begin(), v.end());
		std::reverse(d.begin(), d.end());
		if (grid) {
			grid_dist(d, dabove, v, vabove, res, second.bs.nrows[i-1], nc, m, lonlat, lat, 1, global, np, sp);
		} else {
			cost_dist(d, dabove, v, vabove, res, second.bs.nrows[i-1], nc, m, lonlat, lat, 1, global, np, sp);
		}
		std::reverse(d.begin(), d.end());
		if (converged) {
			old.readBlock(v, second.bs, i-1);
			block_is_same(converged, d, v);
		}
		if (!second.writeValuesRect(d, second.bs.row[i-1], second.bs.nrows[i-1], 0, nc)) return second;
	}
	second.writeStop();
	first.readStop();
	if (old.hasValues()) {
		old.readStop();
	}
	readStop();
	return(second);
}

SpatRaster SpatRaster::costDistance(double target, double m, size_t maxiter, bool grid, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (!hasValues()) {
		out.setError("cannot compute distance for a raster with no values");
		return out;
	}

	std::string filename = opt.get_filename();
	if (!filename.empty()) {
		if (file_exists(filename) && (!opt.get_overwrite())) {
			out.setError("output file exists. You can use 'overwrite=TRUE' to overwrite it");
			return(out);
		}
	}

	SpatOptions ops(opt);
	if (nlyr() > 1) {
		std::vector<size_t> lyr = {0};
		out = subset(lyr, ops);
		out = out.costDistance(target, m, maxiter, grid, opt);
		out.addWarning("distance computations are only done for the first input layer");
		return out;
	}

	bool lonlat = is_lonlat();
	bool global = is_global_lonlat();
	int polar = ns_polar();
	bool npole = (polar == 1) || (polar == 2);
	bool spole = (polar == -1) || (polar == 2);

	double scale;
	if (!lonlat) {
		scale = source[0].srs.to_meter();
		scale = std::isnan(scale) ? 1 : scale;
		scale /= m;
	} else {
		scale = m;
	}

	std::vector<double> res = resolution();

	size_t i = 0;
	bool converged=false;
	while (i < maxiter) {
		out = costDistanceRun(out, converged, target, scale, lonlat, global, npole, spole, grid, ops);
		if (out.hasError()) return out;
		if (converged) break;
		converged = true;
		i++;
	}
	if (!filename.empty()) {
		out = out.writeRaster(opt);
	}
	if (i == maxiter) {
		out.addWarning("distance algorithm did not converge");
	}
	return(out);
}






std::vector<double> broom_dist_planar(std::vector<double> &v, std::vector<double> &above, std::vector<double> res, size_t nr, size_t nc, double lindist) {

	double dx = res[0] * lindist;
	double dy = res[1] * lindist;
	double dxy = sqrt(dx * dx + dy *dy);

	std::vector<double> dist(v.size(), 0);

    //left to right
	//first cell, no cell left of it
	if ( std::isnan(v[0]) ) {
		dist[0] = above[0] + dy;
	}
	//first row, no row above it, use "above"
	for (size_t i=1; i<nc; i++) {
		if (std::isnan(v[i])) {
			dist[i] = std::min(std::min(above[i] + dy, above[i-1] + dxy), dist[i-1] + dx);
		}
	}

	//other rows
	for (size_t r=1; r<nr; r++) {
		size_t start=r*nc;
		if (std::isnan(v[start])) {
			dist[start] = dist[start-nc] + dy;
		}

		size_t end = start+nc;
		for (size_t i=(start+1); i<end; i++) {
			if (std::isnan(v[i])) {
				dist[i] = std::min(std::min(dist[i-1] + dx, dist[i-nc] + dy), dist[i-nc-1] + dxy);
			}
		}
	}

	//right to left
	//first cell
	if ( std::isnan(v[nc-1])) {
		dist[nc-1] = std::min(dist[nc-1], above[nc-1] + dy);
	}

	// other cells on first row
	for (int i=(nc-2); i > -1; i--) {
		if (std::isnan(v[i])) {
			dist[i] = std::min(std::min(std::min(dist[i+1] + dx, above[i+1] + dxy), above[i] + dy), dist[i]);
		}
	}

	// other rows
	for (size_t r=1; r<nr; r++) {
		size_t start=(r+1)*nc-1;
		if (std::isnan(v[start])) {
			dist[start] = std::min(dist[start], dist[start-nc] + dy);
		}
		for (size_t i=start-1; i>=(r*nc); i--) {
			if (std::isnan(v[i])) {
				dist[i] = std::min(std::min(std::min(dist[i], dist[i+1] + dx), dist[i-nc] + dy), dist[i-nc+1] + dxy);
			}
		}
	}
	size_t off = (nr-1) * nc;
	above = std::vector<double>(dist.begin()+off, dist.end());
	return dist;
}

/*
inline double minNArm(const double &a, const double &b) {
	if (std::isnan(a)) return b;
	if (std::isnan(b)) return a;
	return std::min(a, b);
}
*/

inline void DxDxy(const double &lat, const size_t &row, const double &xres, double yres, const double &dir, const double &scale, double &dx,  double &dy, double &dxy) {
	double rlat = lat + (double)row * yres * dir;
	dx  = distance_lonlat(0, rlat, xres, rlat) / scale;
	yres *= -dir;
	dy  = distance_lonlat(0, rlat, 0, rlat+yres);
	dxy = distance_lonlat(0, rlat, xres, rlat+yres);
	dy = std::isnan(dy) ? std::numeric_limits<double>::infinity() : dy / scale;
	dxy = std::isnan(dxy) ? std::numeric_limits<double>::infinity() : dxy / scale;
}


void broom_dist_geo(std::vector<double> &dist, std::vector<double> &v, std::vector<double> &above, std::vector<double> res, size_t nr, size_t nc, double lat, double latdir, double scale, bool npole, bool spole) {

	double dx, dy, dxy;
	//top to bottom
    //left to right
	DxDxy(lat, 0, res[0], res[1], latdir, scale, dx, dy, dxy);
	//first cell, no cell left of it
	if ( std::isnan(v[0]) ) {
		dist[0] = std::min(above[0] + dy, dist[0]);
	}
	//first row, no row above it, use "above"
	for (size_t i=1; i<nc; i++) {
		if (std::isnan(v[i])) {
			dist[i] = std::min(std::min(std::min(above[i] + dy, above[i-1] + dxy), dist[i-1] + dx), dist[i]);
		}
	}
	if (npole) {
		double minp = *std::min_element(dist.begin(), dist.begin()+nc);
		minp += dy;
		for (size_t i=0; i<nc; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	//other rows
	for (size_t r=1; r<nr; r++) {
		DxDxy(lat, r, res[0], res[1], latdir, scale, dx, dy, dxy);
		size_t start=r*nc;
		if (std::isnan(v[start])) {
			dist[start] = std::min(dist[start], dist[start-nc] + dy);
		}
		size_t end = start+nc;
		for (size_t i=(start+1); i<end; i++) {
			if (std::isnan(v[i])) {
				dist[i] = std::min(dist[i], std::min(std::min(dist[i-1] + dx, dist[i-nc] + dy), dist[i-nc-1] + dxy));
			}
		}
	}
	if (spole) {
		double minp = *std::min_element(dist.end()-nc, dist.end());
		minp += dy;
		size_t ds = dist.size();
		for (size_t i=ds-nc; i<ds; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	//right to left
	DxDxy(lat, 0, res[0], res[1], latdir, scale, dx, dy, dxy);
	 //first cell
	if ( std::isnan(v[nc-1])) {
		dist[nc-1] = std::min(dist[nc-1], above[nc-1] + dy);
	}

	// other cells on first row
	for (int i=(nc-2); i > -1; i--) {
		if (std::isnan(v[i])) {
			dist[i] = std::min(std::min(std::min(dist[i+1] + dx, above[i+1] + dxy), above[i] + dy), dist[i]);
		}
	}
	if (npole) {
		double minp = *std::min_element(dist.begin(), dist.begin()+nc);
		minp += dy;
		for (size_t i=0; i<nc; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}
	// other rows
	for (size_t r=1; r<nr; r++) {
		DxDxy(lat, r, res[0], res[1], latdir, scale, dx, dy, dxy);

		size_t start=(r+1)*nc-1;
		if (std::isnan(v[start])) {
			dist[start] = std::min(dist[start], dist[start-nc] + dy);
		}
		for (size_t i=start-1; i>=(r*nc); i--) {
			if (std::isnan(v[i])) {
				dist[i] = std::min(std::min(std::min(dist[i], dist[i+1] + dx), dist[i-nc] + dy), dist[i-nc+1] + dxy);
			}
		}
	}
	if (spole) {
		double minp = *std::min_element(dist.end()-nc, dist.end());
		minp += dy;
		size_t ds = dist.size();
		for (size_t i=ds-nc; i<ds; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	size_t off = (nr-1) * nc;
	above = std::vector<double>(dist.begin()+off, dist.end());
}


void broom_dist_geo_global(std::vector<double> &dist, std::vector<double> &v, std::vector<double> &above, std::vector<double> res, size_t nr, size_t nc, double lat, double latdir, double scale, bool npole, bool spole) {

//	double dy = distance_lonlat(0, 0, 0, res[1]);
	double dx, dy, dxy;
	size_t stopnc = nc - 1;
	//top to bottom
    //left to right
	DxDxy(lat, 0, res[0], res[1], latdir, scale, dx, dy, dxy);
	//first cell, no cell left of it
	if ( std::isnan(v[0]) ) {
		dist[0] = std::min(std::min(std::min(above[0] + dy, above[stopnc] + dxy), dist[stopnc] + dx), dist[0]);
	}
	//first row, no row above it, use "above"
	for (size_t i=1; i<nc; i++) {
		if (std::isnan(v[i])) {
			dist[i] = std::min(std::min(std::min(above[i] + dy, above[i-1] + dxy),
								dist[i-1] + dx), dist[i]);
		}
	}
	if (npole) {
		double minp = *std::min_element(dist.begin(), dist.begin()+nc);
		minp += dy;
		for (size_t i=0; i<nc; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	//other rows
	for (size_t r=1; r<nr; r++) {
		DxDxy(lat, r, res[0], res[1], latdir, scale, dx, dy, dxy);
		size_t start=r*nc;
		if (std::isnan(v[start])) {
			dist[start] = std::min(std::min(std::min(dist[start-nc] + dy, dist[start-1] + dxy),
									dist[start+stopnc] + dx), dist[start]);
		}

		size_t end = start+nc;
		for (size_t i=(start+1); i<end; i++) {
			if (std::isnan(v[i])) {
				dist[i] = std::min(std::min(std::min(dist[i-1] + dx, dist[i-nc] + dy),
						dist[i-nc-1] + dxy), dist[i]);
			}
		}
	}

	if (spole) {
		double minp = *std::min_element(dist.end()-nc, dist.end());
		minp += dy;
		size_t ds = dist.size();
		for (size_t i=ds-nc; i<ds; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	//right to left
	DxDxy(lat, 0, res[0], res[1], latdir, scale, dx, dy, dxy);
	 //first cell
	if ( std::isnan(v[stopnc])) {
		dist[stopnc] = std::min(std::min(std::min(dist[stopnc], above[stopnc] + dy), above[0] + dxy), dist[0] + dx);
	}

	// other cells on first row
	for (int i=(nc-2); i > -1; i--) {
		if (std::isnan(v[i])) {
			dist[i] = std::min(std::min(std::min(dist[i+1] + dx, above[i+1] + dxy),
						above[i] + dy), dist[i]);
		}
	}
	if (npole) {
		double minp = *std::min_element(dist.begin(), dist.begin()+nc);
		minp += dy;
		for (size_t i=0; i<nc; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	// other rows
	for (size_t r=1; r<nr; r++) {
		DxDxy(lat, r, res[0], res[1], latdir, scale, dx, dy, dxy);

		size_t start=(r+1)*nc-1;
		if (std::isnan(v[start])) {
			//dist[start] = std::min(dist[start], dist[start-nc] + dy);
			dist[start] = std::min(std::min(std::min(dist[start], dist[start-nc] + dy), dist[start-nc-stopnc] + dxy)
							, dist[start-stopnc] + dx);
		}
		size_t end = r*nc;
		for (size_t i=start-1; i>=end; i--) {
			if (std::isnan(v[i])) {
			//	dist[i] = std::min(std::min(std::min(dist[i], dist[i+1] + dx), dist[i-nc] + dy), dist[i-nc+1] + dxy);
				dist[i] = std::min(std::min(std::min(dist[i], dist[i+1] + dx), dist[i-nc] + dy),
						dist[i-nc+1] + dxy);
			}
		}
	}
	if (spole) {
		double minp = *std::min_element(dist.end()-nc, dist.end());
		minp += dy;
		size_t ds = dist.size();
		for (size_t i=ds-nc; i<ds; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	size_t off = (nr-1) * nc;
	above = std::vector<double>(dist.begin()+off, dist.end());

}


SpatRaster SpatRaster::gridDistance(double m, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (!hasValues()) {
		out.setError("cannot compute distance for a raster with no values");
		return out;
	}

	SpatOptions ops(opt);
	size_t nl = nlyr();
	if (nl > 1) {
		out.source.resize(nl);
		for (size_t i=0; i<nl; i++) {
			std::vector<size_t> lyr = {i};
			SpatRaster r = subset(lyr, ops);
			r = r.gridDistance(m, ops);
			out.source[i] = r.source[0];
		}
		if (!opt.get_filename().empty()) {
			out = out.writeRaster(opt);
		}
		return out;
	}

	SpatRaster first = out.geometry();

	std::vector<double> res = resolution();
	size_t nc = ncol();
	bool lonlat = is_lonlat();
    std::vector<double> d, v;
	std::vector<double> above(nc, std::numeric_limits<double>::infinity());

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	if (lonlat) {
		bool global = is_global_lonlat();
		int polar = ns_polar();
		bool npole = (polar == 1) || (polar == 2);
		bool spole = (polar == -1) || (polar == 2);
		SpatRaster second = first;
		if (!first.writeStart(ops, filenames())) { return first; }
		for (size_t i = 0; i < first.bs.n; i++) {
			readBlock(v, first.bs, i);
			d.resize(v.size(), std::numeric_limits<double>::infinity());
			for (size_t j=0; j<d.size(); j++) {
				if (!std::isnan(v[j])) d[j] = 0;
			}
			double lat = yFromRow(first.bs.row[i]);
			bool np = (i==0) && npole;
			bool sp = (i==first.bs.n-1) && spole;
			if (global) {
				broom_dist_geo_global(d, v, above, res, first.bs.nrows[i], nc, lat, -1, m, np, sp);
			} else {
				broom_dist_geo(d, v, above, res, first.bs.nrows[i], nc, lat, -1, m, np, sp);
			}
			if (!first.writeValues(d, first.bs.row[i], first.bs.nrows[i])) return first;
		}
		first.writeStop();

		if (!first.readStart()) {
			readStop();
			return(first);
		}
		if (!second.writeStart(ops, filenames())) {
			readStop();
			return second;
		}
		above = std::vector<double>(ncol(), std::numeric_limits<double>::infinity());
		for (int i = second.bs.n; i>0; i--) {
			readBlock(v, second.bs, i-1);
			first.readBlock(d, second.bs, i-1);
			std::reverse(v.begin(), v.end());
			std::reverse(d.begin(), d.end());
			double lat = yFromRow(second.bs.row[i-1] + second.bs.nrows[i-1] - 1);
			bool sp = (i==1) && spole; //! reverse order
			bool np = (i==(int)second.bs.n) && npole;
			if (global) {
				broom_dist_geo_global(d, v, above, res, second.bs.nrows[i-1], nc, lat, 1, m, np, sp);
			} else {
				broom_dist_geo(d, v, above, res, second.bs.nrows[i-1], nc, lat, 1, m, np, sp);
			}
			std::reverse(d.begin(), d.end());
			if (!second.writeValuesRect(d, second.bs.row[i-1], second.bs.nrows[i-1], 0, nc)) return second;
		}
		second.writeStop();

		if (!second.readStart()) {
			readStop();
			return(second);
		}
		if (!out.writeStart(opt, filenames())) {
			readStop();
			return out;
		}
		above = std::vector<double>(ncol(), std::numeric_limits<double>::infinity());
		for (size_t i = 0; i < out.bs.n; i++) {
			readBlock(v, out.bs, i);
			second.readBlock(d, out.bs, i);
			double lat = yFromRow(first.bs.row[i]);
			bool np = (i==0) && npole;
			bool sp = (i==out.bs.n-1) && spole;
			if (global) {
				broom_dist_geo_global(d, v, above, res, out.bs.nrows[i], nc, lat, -1, m, np, sp);
			} else {
				broom_dist_geo(d, v, above, res, out.bs.nrows[i], nc, lat, -1, m, np, sp);
			}
			if (!out.writeValues(d, out.bs.row[i], out.bs.nrows[i])) return out;
		}
		second.readStop();
		out.writeStop();
		readStop();

	} else {
		double scale = source[0].srs.to_meter() / m;
		scale = std::isnan(scale) ? 1 : scale;

		if (!first.writeStart(ops, filenames())) { return first; }
		std::vector<double> vv;
		for (size_t i = 0; i < first.bs.n; i++) {
			readBlock(v, first.bs, i);
			d = broom_dist_planar(v, above, res, first.bs.nrows[i], nc, scale);
			if (!first.writeValues(d, first.bs.row[i], first.bs.nrows[i])) return first;
		}
		first.writeStop();

		if (!first.readStart()) {
			out.setError(first.getError());
			return(out);
		}
		above = std::vector<double>(ncol(), std::numeric_limits<double>::infinity());
		if (!out.writeStart(opt, filenames())) {
			readStop();
			return out;
		}
		for (int i = out.bs.n; i>0; i--) {
			readBlock(v, out.bs, i-1);
			std::reverse(v.begin(), v.end());
			d = broom_dist_planar(v, above, res, out.bs.nrows[i-1], nc, scale);
			first.readBlock(vv, out.bs, i-1);
			std::transform(d.rbegin(), d.rend(), vv.begin(), vv.begin(), [](double a, double b) {return std::min(a,b);});
			if (!out.writeValuesRect(vv, out.bs.row[i-1], out.bs.nrows[i-1], 0, nc)) return out;
		}
		out.writeStop();
		readStop();
	}
	return(out);
}


std::vector<double> do_edge(const std::vector<double> &d, const size_t nrow, const size_t ncol, const bool classes, const bool inner, const unsigned dirs, double falseval) {

	size_t n = nrow * ncol;
	std::vector<double> val(n, falseval);

	int r[8] = { -1,0,0,1 , -1,-1,1,1};
	int c[8] = { 0,-1,1,0 , -1,1,-1,1};

	if (!classes) {
		if (inner) { // inner
			for (size_t i = 1; i < (nrow-1); i++) {
				for (size_t j = 1; j < (ncol-1); j++) {
					size_t cell = i*ncol+j;
					val[cell] = NAN;
					if ( !std::isnan(d[cell])) {
						val[cell] = falseval;
						for (size_t k=0; k< dirs; k++) {
							if ( std::isnan(d[cell + r[k] * ncol + c[k]])) {
								val[cell] = 1;
								break;
							}
						}
					}
				}
			}

		} else { //outer
			for (size_t i = 1; i < (nrow-1); i++) {
				for (size_t j = 1; j < (ncol-1); j++) {
					size_t cell = i*ncol+j;
					val[cell] = falseval;
					if (std::isnan(d[cell])) {
						val[cell] = NAN;
						for (size_t k=0; k < dirs; k++) {
							if ( !std::isnan(d[cell+ r[k] * ncol + c[k] ])) {
								val[cell] = 1;
								break;
							}
						}
					}
				}
			}
		}
	} else { // by class
		for (size_t i = 1; i < (nrow-1); i++) {
			for (size_t j = 1; j < (ncol-1); j++) {
				size_t cell = i*ncol+j;
				double test = d[cell+r[0]*ncol+c[0]];
				val[cell] = std::isnan(test) ? NAN : falseval;
				for (size_t k=1; k<dirs; k++) {
					double v = d[cell+r[k]*ncol +c[k]];
					if (std::isnan(test)) {
						if (!std::isnan(v)) {
							val[cell] = 1;
							break;
						}
					} else if (test != v) {
						val[cell] = 1;
						break;
					}
				}
			}
		}

	}
	return(val);
}



void addrowcol(std::vector<double> &v, size_t nr, size_t nc, bool rowbefore, bool rowafter, bool cols) {

	if (rowbefore) {
		v.insert(v.begin(), v.begin(), v.begin()+nc);
		nr++;
	}
	if (rowafter) {
		v.insert(v.end(), v.end()-nc, v.end());
		nr++;
	}
	if (cols) {
		for (size_t i=0; i<nr; i++) {
			size_t j = i*(nc+2);
			v.insert(v.begin()+j+nc, v[j+nc-1]);
			v.insert(v.begin()+j, v[j]);
		}
	}
}


void striprowcol(std::vector<double> &v, size_t nr, size_t nc, bool rows, bool cols) {
	if (rows) {
		v.erase(v.begin(), v.begin()+nc);
		v.erase(v.end()-nc, v.end());
		nr -= 2;
	}
	if (cols) {
		nc -= 2;
		for (size_t i=0; i<nr; i++) {
			size_t j = i*nc;
			v.erase(v.begin()+j);
			v.erase(v.begin()+j+nc);
		}
	}
}


SpatRaster SpatRaster::edges(bool classes, std::string type, unsigned directions, double falseval, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (nlyr() > 1) {
		std::vector<size_t> lyr = {0};
		SpatOptions ops(opt);
		out = subset(lyr, ops);
		out = out.edges(classes, type, directions, falseval, opt);
		out.addWarning("boundary detection is only done for the first layer");
		return out;
	}
	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return out;
	}


	if ((directions != 4) && (directions != 8)) {
		out.setError("directions should be 4 or 8");
		return(out);
	}
	if ((type != "inner") && (type != "outer")) {
		out.setError("directions should be 'inner' or 'outer'");
		return(out);
	}
	bool inner = type == "inner";

	size_t nc = ncol();
	size_t nr = nrow();


	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	opt.minrows = 2;
 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v;
		//bool before = false;
		//bool after = false;
		if (i == 0) {
			if (out.bs.n == 1) {
				readValues(v, out.bs.row[i], out.bs.nrows[i], 0, nc);
				addrowcol(v, nr, nc, true, true, true);
			} else {
				readValues(v, out.bs.row[i], out.bs.nrows[i]+1, 0, nc);
				addrowcol(v, nr, nc, true, false, true);
				//after = true;
			}
		} else {
			//before = true;
			if (i == out.bs.n) {
				readValues(v, out.bs.row[i]-1, out.bs.nrows[i]+1, 0, nc);
				addrowcol(v, nr, nc, false, true, true);
			} else {
				readValues(v, out.bs.row[i]-1, out.bs.nrows[i]+2, 0, nc);
				addrowcol(v, nr, nc, false, false, true);
				//after = true;
			}
		}
		//before, after,
		std::vector<double> vv = do_edge(v, out.bs.nrows[i]+2, nc+2, classes, inner, directions, falseval);
		striprowcol(vv, out.bs.nrows[i]+2, nc+2, true, true);
		if (!out.writeBlock(vv, i)) return out;
	}
	out.writeStop();
	readStop();

	return(out);
}



SpatRaster SpatRaster::buffer(double d, double background, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return out;
	}

	if (d <= 0) {
		out.setError("buffer should be > 0");
		return out;
	}

	if (background == 1) {
		out.setError("the background value cannot be 1");
		return out;
	}

	SpatOptions ops(opt);
	size_t nl = nlyr();
	if (nl > 1) {
		std::vector<std::string> nms = getNames();
		if (ops.names.size() == nms.size()) {
			nms = opt.names;
		}
		out.source.resize(nl);
		for (size_t i=0; i<nl; i++) {
			std::vector<size_t> lyr = {i};
			SpatRaster r = subset(lyr, ops);
			ops.names = {nms[i]};
			r = r.buffer(d, background, ops);
			out.source[i] = r.source[0];
		}
		if (!opt.get_filename().empty()) {
			out = out.writeRaster(opt);
		}
		return out;
	}


	if (!is_lonlat()) {
		if (!std::isnan(background)) {
			out = proximity(NAN, NAN, false, "", true, d, true, ops);
			if (background == 0) {
				out = out.isnotnan(false, opt);
			} else {
				out = out.replaceValues({NAN}, {background}, 1, false, NAN, false, opt);
			}
		} else {
			out = proximity(NAN, NAN, false, "", true, d, true, opt);
		}
	} else {
		SpatRaster e = edges(false, "inner", 8, NAN, ops);
		SpatVector p = e.as_points(false, true, false, ops);
		p = p.buffer({d}, 10, "", "", NAN, false);
		p = p.aggregate(true);
		out = out.rasterize(p, "", {1}, background, false, "", false, false, true, opt);
		if (background == 0) {
			out.setValueType(3);
		}
		//out = out.disdir_vector_rasterize(p, false, true, false, false, NAN, NAN, "m", ops);
		//out = out.arith(d, "<=", false, opt);
	}
	if (source[0].srs.is_empty()) {
		out.addWarning("unknown CRS. Results may be wrong");
	} 
	return out;
}



SpatRaster SpatRaster::rst_area(bool mask, std::string unit, bool transform, int rcmax, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (out.source[0].srs.wkt.empty()) {
		addWarning("unknown CRS. Results can be wrong");
		transform = false;
	}

	std::vector<std::string> f {"m", "km", "ha"};
	if (std::find(f.begin(), f.end(), unit) == f.end()) {
		out.setError("invalid unit");
		return out;
	}

	if (opt.names.empty()) {
		opt.names = {"area"};
	}
	bool lonlat = is_lonlat();

	SpatOptions mopt(opt);
	if (mask) {
		if (!hasValues()) {
			mask = false;
		} else {
			mopt.filenames = opt.filenames;
			opt.filenames = {""};
		}
	}

	SpatOptions xopt(mopt);

	if (lonlat) {
		bool disagg = false;
		SpatExtent extent = getExtent();
		if ((out.ncol() == 1) && ((extent.xmax - extent.xmin) > 180)) {
			disagg = true;
			std::vector<size_t> fact = {1,2};
			out = out.disaggregate(fact, xopt);
		}
		SpatExtent e = {extent.xmin, extent.xmin+out.xres(), extent.ymin, extent.ymax};
		SpatRaster onecol = out.crop(e, "near", false, xopt);
		SpatVector p = onecol.as_polygons(false, false, false, false, false, 0, xopt);
		if (p.hasError()) {
			out.setError(p.getError());
			return out;
		}
		std::vector<double> a = p.area(unit, true, {});
		size_t nc = out.ncol();
		if (disagg) {
			if (!out.writeStart(xopt, filenames())) { return out; }
		} else {
			if (!out.writeStart(opt, filenames())) { return out; }
		}
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v;
			v.reserve(out.bs.nrows[i] * nc);
			size_t r = out.bs.row[i];
			for (size_t j=0; j<out.bs.nrows[i]; j++) {
				v.insert(v.end(), nc, a[r]);
				r++;
			}
			if (!out.writeBlock(v, i)) return out;
		}
		out.writeStop();
		if (disagg) {
			SpatRaster tmp = out.to_memory_copy(xopt);
			std::vector<size_t> fact = {1,2};
			opt.overwrite=true;
			out = tmp.aggregate(fact, "sum", true, opt);
		}

	} else {
		if (transform) {
			transform = can_transform(source[0].srs.wkt, "EPSG:4326");
		}

		if (transform) {
			bool resample = false;
//			SpatRaster empty = out.geometry(1);
			size_t rcx = std::max(rcmax, 10);
			size_t frow = 1, fcol = 1;
			SpatRaster target = out.geometry(1);
			if ((nrow() > rcx) || (ncol() > rcx)) {
				resample = true;
				frow = (nrow() / rcx) + 1;
				fcol = (ncol() / rcx) + 1;
				out = out.aggregate({frow, fcol}, "mean", false, xopt);
				xopt.ncopies *= 5;
				if (!out.writeStart(xopt, filenames())) { return out; }
			} else {
				opt.ncopies *= 5;
				if (!out.writeStart(opt, filenames())) { return out; }
			}
			SpatRaster empty = out.geometry(1);
			SpatExtent extent = out.getExtent();
			double dy = out.yres() / 2;
			for (size_t i = 0; i < out.bs.n; i++) {
				double ymax = out.yFromRow(out.bs.row[i]) + dy;
				double ymin = out.yFromRow(out.bs.row[i] + out.bs.nrows[i]-1) - dy;
				SpatExtent e = {extent.xmin, extent.xmax, ymin, ymax};
				SpatRaster chunk = empty.crop(e, "near", false, xopt);
				SpatVector p = chunk.as_polygons(false, false, false, false, false, 0, xopt);
				std::vector<double> v = p.area(unit, true, {});
				if (!out.writeBlock(v, i)) return out;
				out.writeStop();
			}
			if (resample) {
				double divr = (double)(frow*fcol);
				out = out.arith(divr, "/", false, false, xopt);
				out = out.warper(target, "", "bilinear", false, false, true, opt);
			}
		} else {
			if (!out.writeStart(opt, filenames())) { return out; }
			double u = unit == "m" ? 1 : unit == "km" ? 1000000 : 10000;
			double m = out.source[0].srs.to_meter();
			double a = std::isnan(m) ? 1 : m;
			a *= xres() * yres() / u;
			for (size_t i = 0; i < out.bs.n; i++) {
				std::vector<double> v(out.bs.nrows[i]*ncol(), a);
				if (!out.writeBlock(v, i)) return out;
			}
			out.writeStop();
		}
	}

	if (mask) {
		out = out.mask(*this, false, NAN, NAN, mopt);
	}
	return(out);
}


std::vector<std::vector<double>> SpatRaster::sum_area(std::string unit, bool transform, bool by_value, SpatOptions &opt) {

	if (source[0].srs.wkt.empty()) {
		addWarning("unknown CRS. Results can be wrong");
		transform = false;
	}

	std::vector<std::string> f {"m", "km", "ha"};
	if (std::find(f.begin(), f.end(), unit) == f.end()) {
		setError("invalid unit");
		return {{NAN}};
	}

	if (transform) { //avoid very large polygon objects
		if (!is_lonlat()) {
			transform = can_transform(getSRS("wkt"), "+proj=longlat");
		}
	}
	BlockSize bs = getBlockSize(opt);
	if (!readStart()) {
		return {{NAN}};
	}
	size_t nc = ncol();
	size_t nl = nlyr();
	std::vector<double> out(nl, 0);
	std::vector<std::map<double, double>> m;
	if (by_value) {
		m.resize(nl);
	}

	if (is_lonlat()) {
		SpatRaster x = geometry(1);
		SpatExtent extent = x.getExtent();
		if ((nc == 1) && ((extent.xmax - extent.xmin) > 180)) {
			std::vector<size_t> fact= {1,2};
			x = x.disaggregate(fact, opt);
		}
		SpatExtent e = {extent.xmin, extent.xmin+x.xres(), extent.ymin, extent.ymax};
		SpatRaster onecol = x.crop(e, "near", false, opt);
		SpatVector p = onecol.as_polygons(false, false, false, false, false, 0, opt);
		std::vector<double> ar = p.area(unit, true, {});
		if (!hasValues()) {
			out.resize(1);
			for (size_t i=0; i<ar.size(); i++) {
				out[0] += ar[i] * (double)nc;
			}
		} else {
			for (size_t i=0; i<bs.n; i++) {
				std::vector<double> v;
				readValues(v, bs.row[i], bs.nrows[i], 0, ncol());
				size_t blockoff = bs.nrows[i] * nc;
				for (size_t lyr=0; lyr<nl; lyr++) {
					size_t lyroff = lyr * blockoff;
					if (by_value) {
						for (size_t j=0; j<bs.nrows[i]; j++) {
							size_t row = bs.row[i] + j;
							size_t offset = lyroff + j * nc;
							size_t n = offset + nc;
							for (size_t k=offset; k<n; k++) {
								if (!std::isnan(v[k])) {
									if (m[lyr].find(v[k]) == m[lyr].end()) {
										m[lyr][v[k]] = ar[row];
									} else {
										m[lyr][v[k]] += ar[row];
									}
								}
							}
						}
					} else {
						for (size_t j=0; j<bs.nrows[i]; j++) {
							size_t row = bs.row[i] + j;
							size_t offset = lyroff + j * nc;
							size_t n = offset + nc;
							for (size_t k=offset; k<n; k++) {
								if (!std::isnan(v[k])) out[lyr] += ar[row];
							}
						}
					}
				}
			}
		}
	} else if (transform) {
		opt.set_memfrac(std::max(0.1, opt.get_memfrac()/2));
		SpatRaster x = geometry(1);
		SpatExtent extent = x.getExtent();
		double dy = x.yres() / 2;
		SpatOptions popt(opt);
		if (!hasValues()) {
			out.resize(1);
			for (size_t i=0; i<bs.n; i++) {
				double ymax = x.yFromRow(bs.row[i]) + dy;
				double ymin = x.yFromRow(bs.row[i] + bs.nrows[i]-1) - dy;
				SpatExtent e = {extent.xmin, extent.xmax, ymin, ymax};
				SpatRaster onechunk = x.crop(e, "near", false, popt);
				SpatVector p = onechunk.as_polygons(false, false, false, false, false, 0, popt);
				p = p.project("EPSG:4326", false);
				std::vector<double> v = p.area(unit, true, 	{});
				out[0] += accumulate(v.begin(), v.end(), 0.0);
			}
		} else {
			for (size_t i=0; i<bs.n; i++) {
				double ymax = x.yFromRow(bs.row[i]) + dy;
				double ymin = x.yFromRow(bs.row[i] + bs.nrows[i]-1) - dy;
				SpatExtent e = {extent.xmin, extent.xmax, ymin, ymax};
				SpatRaster onechunk = x.crop(e, "near", false, popt);
				SpatVector p = onechunk.as_polygons(false, false, false, false, false, 0, popt);
				p = p.project("EPSG:4326", false);
				std::vector<double> ar = p.area(unit, true, {});
				std::vector<double> v;
				readValues(v, bs.row[i], bs.nrows[i], 0, ncol());
				size_t blockoff = bs.nrows[i] * nc;
				for (size_t lyr=0; lyr<nl; lyr++) {
					size_t lyroff = lyr * blockoff;
					if (by_value) {
						for (size_t j=0; j<bs.nrows[i]; j++) {
							size_t row = bs.row[i] + j;
							size_t offset = lyroff + j * nc;
							size_t n = offset + nc;
							for (size_t k=offset; k<n; k++) {
								if (!std::isnan(v[k])) {
									if (m[lyr].find(v[k]) == m[lyr].end()) {
										m[lyr][v[k]] = ar[row];
									} else {
										m[lyr][v[k]] += ar[row];
									}
								}
							}
						}
					} else {
						for (size_t j=0; j<bs.nrows[i]; j++) {
							size_t row = bs.row[i] + j;
							size_t offset = lyroff + j * nc;
							size_t n = offset + nc;
							for (size_t k=offset; k<n; k++) {
								if (!std::isnan(v[k])) out[lyr] += ar[row];
							}
						}
					}
				}
			}
		}
	} else {
		double adj = unit == "m" ? 1 : unit == "km" ? 1000000 : 10000;
		double unit = source[0].srs.to_meter();
		unit = std::isnan(unit) ? 1 : unit;
		double ar = xres() * yres() * unit * unit / adj;
		if (!hasValues()) {
			out.resize(1);
			out[0] = ncell() * ar;
		} else {
			for (size_t i=0; i<bs.n; i++) {
				std::vector<double> v;
				readValues(v, bs.row[i], bs.nrows[i], 0, ncol());
				size_t blockoff = bs.nrows[i] * nc;
				for (size_t lyr=0; lyr<nl; lyr++) {
					size_t offset = lyr * blockoff;
					size_t n = offset + blockoff;
					if (by_value) {
						for (size_t k=offset; k<n; k++) {
							if (!std::isnan(v[k])) {
								if (m[lyr].find(v[k]) == m[lyr].end()) {
									m[lyr][v[k]] = ar;
								} else {
									m[lyr][v[k]] += ar;
								}
							}
						}
					} else {
						for (size_t k=offset; k<n; k++) {
							if (!std::isnan(v[k])) {
								out[lyr] += ar;
							}
						}
					}
				}
			}
		}
	}
	readStop();
	if (by_value) {
		std::vector<std::vector<double>> dout(nl);
		for (size_t i=0; i<nl; i++) {
			for (auto& it:m[i]) {
				dout[i].push_back(it.first);
				dout[i].push_back(it.second);
			}
		}
		return dout;
	} 
	return {out};
}




std::vector<std::vector<double>> SpatRaster::sum_area_group(SpatRaster group, std::string unit, bool transform, bool by_value, SpatOptions &opt) {

	if (source[0].srs.wkt.empty()) {
		addWarning("unknown CRS. Results can be wrong");
		transform = false;
	}

	if (!(hasValues() && group.hasValues())) {
		setError("raster has no values");
		return {{NAN}};		
	}

	std::vector<std::string> f {"m", "km", "ha"};
	if (std::find(f.begin(), f.end(), unit) == f.end()) {
		setError("invalid unit");
		return {{NAN}};
	}

	BlockSize bs = getBlockSize(opt);
	if (!readStart()) {
		return {{NAN}};
	}
	if (!group.readStart()) {
		return {{NAN}};
	}
	
	size_t nc = ncol();
	size_t nl = nlyr();
	std::vector<std::map<double, std::map<double, double>>> m(nl);

	if (is_lonlat()) {
		SpatRaster x = geometry(1);
		SpatExtent extent = x.getExtent();
		if ((nc == 1) && ((extent.xmax - extent.xmin) > 180)) {
			std::vector<size_t> fact= {1,2};
			x = x.disaggregate(fact, opt);
		}
		SpatExtent e = {extent.xmin, extent.xmin+x.xres(), extent.ymin, extent.ymax};
		SpatRaster onecol = x.crop(e, "near", false, opt);
		SpatVector p = onecol.as_polygons(false, false, false, false, false, 0, opt);
		std::vector<double> ar = p.area(unit, true, {});
		for (size_t i=0; i<bs.n; i++) {
			std::vector<double> v, g;
			readValues(v, bs.row[i], bs.nrows[i], 0, ncol());
			group.readValues(g, bs.row[i], bs.nrows[i], 0, ncol());
			size_t blockoff = bs.nrows[i] * nc;
			for (size_t lyr=0; lyr<nl; lyr++) {
				size_t lyroff = lyr * blockoff;
				if (by_value) {
					for (size_t j=0; j<bs.nrows[i]; j++) {
						size_t row = bs.row[i] + j;
						size_t rowoff = j * nc;
						size_t offset = lyroff + rowoff;
						for (size_t k=0; k<nc; k++) {
							double vk = v[k+offset];
							double gk = g[k+rowoff];
							if (!(std::isnan(vk) || std::isnan(gk))) {
								if (m[lyr].find(vk) == m[lyr].end() || 
										m[lyr][vk].find(gk) == m[lyr][vk].end()) {	
									m[lyr][vk][gk] = ar[row];
								} else {
									m[lyr][vk][gk] += ar[row];
								}
							}
						}
					}
				} else {
					for (size_t j=0; j<bs.nrows[i]; j++) {
						size_t row = bs.row[i] + j;
						size_t rowoff = j * nc;
						size_t offset = lyroff + rowoff;
						for (size_t k=0; k<nc; k++) {
							double vk = v[k+offset];
							double gk = g[k+rowoff];
							if (!(std::isnan(vk) || std::isnan(gk))) {
								vk = 0;
								if (m[lyr].find(vk) == m[lyr].end() || 
										m[lyr][vk].find(gk) == m[lyr][vk].end()) {	
									m[lyr][vk][gk] = ar[row];
								} else {
									m[lyr][vk][gk] += ar[row];
								}
							}
						}
					}
				}
			}
		}
	} else {
		transform = transform && can_transform(getSRS("wkt"), "+proj=longlat");
		if (transform) {
			opt.set_memfrac(std::max(0.1, opt.get_memfrac()/2));
			SpatRaster x = geometry(1);
			SpatExtent extent = x.getExtent();
			double dy = x.yres() / 2;
			SpatOptions popt(opt);
			for (size_t i=0; i<bs.n; i++) {
				double ymax = x.yFromRow(bs.row[i]) + dy;
				double ymin = x.yFromRow(bs.row[i] + bs.nrows[i]-1) - dy;
				SpatExtent e = {extent.xmin, extent.xmax, ymin, ymax};
				SpatRaster onechunk = x.crop(e, "near", false, popt);
				SpatVector p = onechunk.as_polygons(false, false, false, false, false, 0, popt);
				p = p.project("EPSG:4326", false);
				std::vector<double> ar = p.area(unit, true, {});
				std::vector<double> v, g;
				readValues(v, bs.row[i], bs.nrows[i], 0, ncol());
				group.readValues(g, bs.row[i], bs.nrows[i], 0, ncol());
				size_t blockoff = bs.nrows[i] * nc;
				for (size_t lyr=0; lyr<nl; lyr++) {
					size_t lyroff = lyr * blockoff;
					if (by_value) {
						for (size_t j=0; j<bs.nrows[i]; j++) {
							size_t row = bs.row[i] + j;
							size_t rowoff = j * nc;
							size_t offset = lyroff + rowoff;
							for (size_t k=0; k<nc; k++) {
								double vk = v[k+offset];
								double gk = g[k+rowoff];							
								if (!(std::isnan(vk) || std::isnan(gk))) {
									if (m[lyr].find(vk) == m[lyr].end() ||
											m[lyr][vk].find(gk) == m[lyr][vk].end()) {	
										m[lyr][vk][gk] = ar[row];
									} else {
										m[lyr][vk][gk] += ar[row];
									}
								}
							}
						}
					} else {
						for (size_t j=0; j<bs.nrows[i]; j++) {
							size_t row = bs.row[i] + j;
							size_t rowoff = j * nc;
							size_t offset = lyroff + rowoff;
							for (size_t k=0; k<nc; k++) {
								double vk = v[k+offset];
								double gk = g[k+rowoff];							
								if (!(std::isnan(vk) || std::isnan(gk))) {
									vk = 0;
									if (m[lyr].find(vk) == m[lyr].end() ||
											m[lyr][vk].find(gk) == m[lyr][vk].end()) {	
										m[lyr][vk][gk] = ar[row];
									} else {
										m[lyr][vk][gk] += ar[row];
									}
								}
							}
						}
					}
				}
			}
		} else {
			double adj = unit == "m" ? 1 : unit == "km" ? 1000000 : 10000;
			double unit = source[0].srs.to_meter();
			unit = std::isnan(unit) ? 1 : unit;
			double ar = xres() * yres() * unit * unit / adj;
			for (size_t i=0; i<bs.n; i++) {
				std::vector<double> v, g;
				readValues(v, bs.row[i], bs.nrows[i], 0, ncol());
				group.readValues(g, bs.row[i], bs.nrows[i], 0, ncol());
				size_t blockoff = bs.nrows[i] * nc;
				for (size_t lyr=0; lyr<nl; lyr++) {
					size_t offset = lyr * blockoff;
					if (by_value) {
						for (size_t k=0; k<blockoff; k++) {
							double vk = v[offset+k];
							double gk = g[blockoff+k];
							if (!(std::isnan(vk) || std::isnan(gk))) {
								if (m[lyr].find(vk) == m[lyr].end() || 
										m[lyr][vk].find(gk) == m[lyr][vk].end()) {	
									m[lyr][vk][gk] = ar;
								} else {
									m[lyr][vk][gk] += ar;
								}
							}
						}
					} else {
						for (size_t k=0; k<blockoff; k++) {
							double vk = v[offset+k];
							double gk = g[blockoff+k];
							if (!(std::isnan(vk) || std::isnan(gk))) {
								vk = 0;
								if (m[lyr].find(vk) == m[lyr].end() ||
										m[lyr][vk].find(gk) == m[lyr][vk].end()) {	
									m[lyr][vk][gk] = ar;
								} else {
									m[lyr][vk][gk] += ar;
								}
							}
						}
					}
				}
			}
		}
	}
	
	readStop();
	group.readStop();
	std::vector<std::vector<double>> out(nl);
	for (size_t i=0; i<nl; i++) {
		for (auto& it1:m[i]) {
			std::map<double, double> gm = it1.second;
			for (auto& it2:gm) {
				out[i].push_back((double)i);
				out[i].push_back(it1.first);
				out[i].push_back(it2.first);
				out[i].push_back(it2.second);
			}
		}
	} 
	return out;
}

	


size_t get_k(const std::vector<double> &r, std::default_random_engine &generator, std::uniform_int_distribution<> &U) {
	double dmin = 0;
	size_t k = 0;
	for (size_t j=0; j<8; j++) {
		if (r[j] > dmin) {
			dmin = r[j];
			k = j + 1;
		} else if (r[j] == dmin) {
			if (U(generator)) {
				dmin = r[j];
				k = j + 1;
			}
		}
	}
	return k;
}


void do_flowdir(std::vector<double> &val, std::vector<double> &d, size_t nrow, size_t ncol, double dx, double dy, unsigned seed, bool before, bool after) {

	if (!before) {
	//	val.resize(val.size() + ncol, NAN);
		std::vector<double> rna(ncol, NAN);
		d.insert(d.begin(), rna.begin(), rna.end());
		nrow++;
	}
	if (!after) {
	//	val.resize(val.size() + ncol, NAN);
		d.resize(d.size()+ncol, NAN);
		nrow++;
	}

	std::vector<double> r(8);
	std::vector<double> p = {0, 1, 2, 4, 8, 16, 32, 64, 128}; // pow(2, j)
	//std::vector<double> p2 = {0, 1, 2, 3, 4, 5, 6, 7, 8}; 
	double dxy = sqrt(dx * dx + dy * dy);

	std::default_random_engine generator(seed);
	std::uniform_int_distribution<> U(0, 1);

	size_t nc1 = ncol - 1;
	for (size_t row=1; row<(nrow-1); row++) {
		//first col
		size_t i = row * ncol;
		if (std::isnan(d[i])) {
			val.push_back( NAN );
		} else {
			r[0] = (d[i] - d[i+1]) / dx;
			r[1] = (d[i] - d[i+1+ncol]) / dxy;
			r[2] = (d[i] - d[i+ncol]) / dy;
			r[3] = NAN;
			r[4] = NAN;
			r[5] = NAN;
			r[6] = (d[i] - d[i-ncol]) / dy;
			r[7] = (d[i] - d[i+1-ncol]) / dxy;
			size_t k = get_k(r, generator, U);
			val.push_back( p[k] );
		}
		for (size_t col=1; col<nc1; col++) {
			i = row * ncol + col;
			if (!std::isnan(d[i])) {
				r[0] = (d[i] - d[i+1]) / dx;
				r[1] = (d[i] - d[i+1+ncol]) / dxy;
				r[2] = (d[i] - d[i+ncol]) / dy;
				r[3] = (d[i] - d[i-1+ncol]) / dxy;
				r[4] = (d[i] - d[i-1]) / dx;
				r[5] = (d[i] - d[i-1-ncol]) / dxy;
				r[6] = (d[i] - d[i-ncol]) / dy;
				r[7] = (d[i] - d[i+1-ncol]) / dxy;
				size_t k = get_k(r, generator, U);
				val.push_back( p[k] );
			} else {
				val.push_back( NAN );
			}
		}	
		//last col
		i = row * ncol + nc1;
		if (!std::isnan(d[i])) {
			r[0] = NAN;
			r[1] = NAN;
			r[2] = (d[i] - d[i+ncol]) / dy;
			r[3] = (d[i] - d[i-1+ncol]) / dxy;
			r[4] = (d[i] - d[i-1]) / dx;
			r[5] = (d[i] - d[i-1-ncol]) / dxy;
			r[6] = (d[i] - d[i-ncol]) / dy;
			r[7] = NAN;
			size_t k = get_k(r, generator, U);
			val.push_back( p[k] );
		} else {
			val.push_back( NAN );
		}
	}
/*	
	if (!before) {
		if (!after) {
			val = std::vector<double>(val.begin()+ncol, val.end()-ncol);
		} else {
			val = std::vector<double>(val.begin()+ncol, val.end());
		}
	} else if (!after) {
		val = std::vector<double>(val.begin(), val.end()-ncol);
	}
*/
//	if (!after) {
//		val.resize(val.size() + ncol, NAN);
//	}

}


void do_TRI(std::vector<double> &val, std::vector<double> const &d, size_t nrow, size_t ncol, bool before, bool after) {
	if (!before) {
		val.resize(val.size() + ncol, NAN);
	}
	for (size_t row=1; row< (nrow-1); row++) {
		val.push_back(NAN);
		for (size_t col=1; col< (ncol-1); col++) {
			size_t i = row * ncol + col;
			val.push_back(
				(fabs(d[i-1-ncol]-d[i]) + fabs(d[i-1]-d[i]) + fabs(d[i-1+ncol]-d[i]) +  fabs(d[i-ncol]-d[i]) + fabs(d[i+ncol]-d[i]) +  fabs(d[i+1-ncol]-d[i]) + fabs(d[i+1]-d[i]) +  fabs(d[i+1+ncol]-d[i])) / 8
			);
		}
		val.push_back(NAN);
	}
	if (!after) {
		val.resize(val.size() + ncol, NAN);
	}
}

inline double pow2(double x) {
	return pow(x, 2);
}

void do_TRI_riley(std::vector<double> &val, std::vector<double> const &d, size_t nrow, size_t ncol, bool before, bool after) {
	if (!before) {
		val.resize(val.size() + ncol, NAN);
	}
	for (size_t row=1; row< (nrow-1); row++) {
		val.push_back(NAN);
		for (size_t col=1; col< (ncol-1); col++) {
			size_t i = row * ncol + col;
			val.push_back(
				sqrt(pow2(d[i-1-ncol]-d[i]) + pow2(d[i-1]-d[i]) + pow2(d[i-1+ncol]-d[i]) + 
				pow2(d[i-ncol]-d[i]) + pow2(d[i+ncol]-d[i]) + pow2(d[i+1-ncol]-d[i]) + 
				pow2(d[i+1]-d[i]) + pow2(d[i+1+ncol]-d[i]))
			);
		}
		val.push_back(NAN);
	}
	if (!after) {
		val.resize(val.size() + ncol, NAN);
	}
}

void do_TRI_rmsd(std::vector<double> &val, std::vector<double> const &d, size_t nrow, size_t ncol, bool before, bool after) {
	if (!before) {
		val.resize(val.size() + ncol, NAN);
	}
	for (size_t row=1; row< (nrow-1); row++) {
		val.push_back(NAN);
		for (size_t col=1; col< (ncol-1); col++) {
			size_t i = row * ncol + col;
			val.push_back(
				sqrt((pow2(d[i-1-ncol]-d[i]) + pow2(d[i-1]-d[i]) + pow2(d[i-1+ncol]-d[i]) + 
				pow2(d[i-ncol]-d[i]) + pow2(d[i+ncol]-d[i]) + pow2(d[i+1-ncol]-d[i]) + 
				pow2(d[i+1]-d[i]) + pow2(d[i+1+ncol]-d[i]))/8)
			);
		}
		val.push_back(NAN);
	}
	if (!after) {
		val.resize(val.size() + ncol, NAN);
	}
}


void do_TPI(std::vector<double> &val, const std::vector<double> &d, const size_t nrow, const size_t ncol, bool before, bool after) {
	if (!before) {
		val.resize(val.size() + ncol, NAN);
	}
	for (size_t row=1; row< (nrow-1); row++) {
		val.push_back(NAN);
		for (size_t col=1; col< (ncol-1); col++) {
			size_t i = row * ncol + col;
			val.push_back( d[i] - (d[i-1-ncol] + d[i-1] + d[i-1+ncol] + d[i-ncol]
			+ d[i+ncol] + d[i+1-ncol] + d[i+1] + d[i+1+ncol]) / 8 );
		}
		val.push_back(NAN);
	}
/*
	if (expand) {
		for (size_t i=1; i < (ncol-1); i++) {
			val[i+add] = d[i] - (d[i-1] + d[i-1+ncol] + d[i+ncol] + d[i+1] + d[i+1+ncol]) / 5;
			size_t j = i+(nrow-1) * ncol;
			val[j+add] = d[j] - (d[j-1-ncol] + d[j-1] + d[j-ncol] + d[j+1-ncol] + d[j+1]) / 5;
		}
		for (size_t row=1; row< (nrow-1); row++) {
			size_t i = row * ncol;
			val[i+add] = d[i] - (d[i-ncol] + d[i+ncol] + d[i+1-ncol] + d[i+1] + d[i+1+ncol]) / 5;
			i += ncol - 1;
			val[i+add] = d[i] - (d[i-ncol] + d[i] + d[i+ncol] + d[i-ncol]) / 5;
		}
		size_t i = 0;
		val[i+add] = d[i] - (d[i+ncol] + d[i+1] + d[i+1+ncol]) / 3;
		i = ncol-1;
		val[i+add] = d[i] - (d[i+ncol] + d[i-1] + d[i-1+ncol]) / 3;
		i = (nrow-1)*ncol;
		val[i+add] = d[i] - (d[i-ncol] + d[i+1] + d[i+1-ncol]) / 3;
		i = (nrow*ncol)-1;
		val[i+add] = d[i] - (d[i-ncol] + d[i-1] + d[i-1-ncol]) / 3;
	}
*/
	if (!after) {
		val.resize(val.size() + ncol, NAN);
	}


}



void do_roughness(std::vector<double> &val, const std::vector<double> &d, size_t nrow, size_t ncol, bool before, bool after) {

	if (!before) {
		val.resize(val.size() + ncol, NAN);
	}

	int incol = (int)ncol;
	int a[9] = { -1-incol, -1, -1+incol, -incol, 0, incol, 1-incol, 1, 1+incol };
	double min, max, v;
	for (size_t row=1; row< (nrow-1); row++) {
		val.push_back(NAN);
		for (size_t col=1; col< (ncol-1); col++) {
			size_t i = row * ncol + col;
			min = d[i + a[0]];
			max = d[i + a[0]];
			for (size_t j = 1; j < 9; j++) {
				v = d[i + a[j]];
				if (v > max) {
					max = v;
				} else if (v < min) {
					min = v;
				}
			}
			val.push_back(max - min);
		}
		val.push_back(NAN);
	}
	if (!after) {
		val.resize(val.size() + ncol, NAN);
	}
}


#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif




void to_degrees(std::vector<double>& x, size_t start) {
	double adj = 180 / M_PI;
	for (size_t i=start; i<x.size(); i++) {
		x[i] *= adj;
	}
}


void do_slope(std::vector<double> &val, const std::vector<double> &d, size_t ngb, size_t nrow, size_t ncol, double dx, double dy, bool geo, std::vector<double> &gy, bool degrees, bool before, bool after) {

	size_t start = val.size();
	if (!before) {
		val.resize(start + ncol, NAN);
	}

	std::vector<double> ddx;
	if (geo) {
		ddx.resize(nrow);
		for (size_t i=0; i<nrow; i++) {
			ddx[i] = distance_lonlat(-dx, gy[i], dx, gy[i]) / 2 ;
		}
	}

	if (ngb == 4) {
		if (geo) {
			double xwi[2] = {-1,1};
			double xw[2] = {0,0};
			double yw[2] = {-1,1};

			for (size_t i=0; i<2; i++) {
				yw[i] = yw[i] / (2 * dy);
			}

			for (size_t row=1; row< (nrow-1); row++) {
				for (size_t k=0; k<2; k++) {
					xw[k] = xwi[k] / (-2 * ddx[row]);
				}
				val.push_back(NAN); // first col
				for (size_t col=1; col< (ncol-1); col++) {
					size_t i = row * ncol + col;

					double zx = d[i-1] * xw[0] + d[i+1] * xw[1];
					double zy = d[i-ncol] * yw[0] + d[i+ncol] * yw[1];
					val.push_back( atan(sqrt( pow(zy, 2) + pow(zx, 2) )));
				}
				val.push_back(NAN); // last col
			}
		} else { // ngb == 8

			double xw[2] = {-1,1};
			double yw[2] = {-1,1};
			for (size_t i=0; i<2; i++) {
				xw[i] /= -2 * dx;
				yw[i] /=  2 * dy;
			}
			for (size_t row=1; row< (nrow-1); row++) {
				val.push_back(NAN);
				for (size_t col=1; col< (ncol-1); col++) {
					size_t i = row * ncol + col;
					double zx = d[i-1] * xw[0] + d[i+1] * xw[1];
					double zy = d[i-ncol] * yw[0] + d[i+ncol] * yw[1];
					val.push_back(  atan( sqrt( pow(zy, 2) + pow(zx, 2)  )));
				}
				val.push_back(NAN);
			}
		}
	} else {

		if (geo) {

			double xwi[6] = {-1,-2,-1,1,2,1};
			double xw[6] = {0,0,0,0,0,0};
			double yw[6] = {-1,1,-2,2,-1,1};

			for (size_t i=0; i<6; i++) {
				yw[i] = yw[i] / (8 * dy);
			}

			for (size_t row=1; row< (nrow-1); row++) {
				for (size_t k=0; k<6; k++) {
					xw[k] = xwi[k] / (8 * ddx[row]);
				}
				val.push_back(NAN);
				for (size_t col=1; col< (ncol-1); col++) {
					size_t i = row * ncol + col;
					double zx = d[i-1-ncol] * xw[0] + d[i-1] * xw[1] + d[i-1+ncol] * xw[2]
					   + d[i+1-ncol] * xw[3] + d[i+1] * xw[4] + d[i+1+ncol] * xw[5];

					double zy = d[i-1-ncol] * yw[0] + d[i-1+ncol] * yw[1] + d[i-ncol] * yw[2]
							+ d[i+ncol] * yw[3] + d[i+1-ncol] * yw[4] + d[i+1+ncol] * yw[5];
					val.push_back( atan(sqrt( pow(zy, 2) + pow(zx, 2))));
				}
				val.push_back(NAN);
			}
		} else {

			double xw[6] = {-1,-2,-1,1,2,1};
			double yw[6] = {-1,1,-2,2,-1,1};
			for (size_t i=0; i<6; i++) {
				xw[i] /= -8 * dx;
				yw[i] /= 8 * dy;
			}
			for (size_t row=1; row< (nrow-1); row++) {
				val.push_back(NAN);
				for (size_t col=1; col< (ncol-1); col++) {
					size_t i = row * ncol + col;
					double zx = d[i-1-ncol] * xw[0] + d[i-1] * xw[1] + d[i-1+ncol] * xw[2]
							+ d[i+1-ncol] * xw[3] + d[i+1] * xw[4] + d[i+1+ncol] * xw[5];
					double zy = d[i-1-ncol] * yw[0] + d[i-1+ncol] * yw[1] + d[i-ncol] * yw[2]
							+ d[i+ncol] * yw[3] + d[i+1-ncol] * yw[4] + d[i+1+ncol] * yw[5];
					val.push_back( atan(sqrt( pow(zy, 2) + pow(zx, 2) )));
				}
				val.push_back(NAN);
			}
		}
	}
	if (!after) {
		val.resize(val.size() + ncol, NAN);
	}

	if (degrees) {
		to_degrees(val, start);
	}
}


double dmod(double x, double n) {
	return(x - n * std::floor(x/n));
}



void do_aspect(std::vector<double> &val, const std::vector<double> &d, size_t ngb, size_t nrow, size_t ncol, double dx, double dy, bool geo, std::vector<double> &gy, bool degrees, bool before, bool after) {

	size_t start = val.size();
	if (!before) {
		val.resize(start + ncol, NAN);
	}
	std::vector<double> ddx;
	if (geo) {
		ddx.resize(nrow);
		for (size_t i=0; i<nrow; i++) {
			ddx[i] = distance_lonlat(-dx, gy[i], dx, gy[i]) / 2 ;
		}
	}
	double zy, zx;

	//double const pi2 = M_PI / 2;
	double const twoPI = 2 * M_PI;
	double const halfPI = M_PI / 2;

	if (ngb == 4) {
		if (geo) {
			double xwi[2] = {-1,1};
			double xw[2] = {0,0};
			double yw[2] = {-1,1};

			for (size_t i=0; i<2; i++) {
				yw[i] = yw[i] / (2 * dy);
			}

			for (size_t row=1; row < (nrow-1); row++) {
				val.push_back(NAN);
				for (size_t k=0; k<2; k++) {
					xw[k] = xwi[k] / (-2 * ddx[row]);
				}
				for (size_t col=1; col< (ncol-1); col++) {
					size_t i = row * ncol + col;
					zx = d[i-1] * xw[0] + d[i+1] * xw[1];
					zy = d[i-ncol] * yw[0] + d[i+ncol] * yw[1];
					zx = atan2(zy, zx);
					val.push_back( std::fmod( halfPI - zx, twoPI) );
				}
				val.push_back(NAN);
			}
		} else {

			double xw[2] = {-1,1};
			double yw[2] = {-1,1};
			for (size_t i=0; i<2; i++) {
				xw[i] = xw[i] / (-2 * dx);
				yw[i] = yw[i] / (2 * dy);
			}
			for (size_t row=1; row< (nrow-1); row++) {
				val.push_back(NAN);
				for (size_t col=1; col< (ncol-1); col++) {
					size_t i = row * ncol + col;
					zx = d[i-1] * xw[0] + d[i+1] * xw[1];
					zy = d[i-ncol] * yw[0] + d[i+ncol] * yw[1];
					zx = atan2(zy, zx);
					val.push_back( dmod( halfPI - zx, twoPI));
				}
				val.push_back(NAN);
			}
		}
	} else { //	(ngb == 8) {

		if (geo) {

			double xwi[6] = {-1,-2,-1,1,2,1};
			double xw[6] = {0,0,0,0,0,0};
			double yw[6] = {-1,1,-2,2,-1,1};

			for (size_t i=0; i<6; i++) {
				yw[i] /= (8 * dy);
			}

			for (size_t row=1; row < (nrow-1); row++) {
				val.push_back(NAN);
				for (size_t k=0; k<6; k++) {
					xw[k] = xwi[k] / (-8 * ddx[row]);
				}
				for (size_t col=1; col< (ncol-1); col++) {
					size_t i = row * ncol + col;

					zx = d[i-1-ncol] * xw[0] + d[i-1] * xw[1] + d[i-1+ncol] * xw[2]
							+ d[i+1-ncol] * xw[3] + d[i+1] * xw[4] + d[i+1+ncol] * xw[5];
					zy = d[i-1-ncol] * yw[0] + d[i-1+ncol] * yw[1] + d[i-ncol] * yw[2]
							+ d[i+ncol] * yw[3] + d[i+1-ncol] * yw[4] + d[i+1+ncol] * yw[5];

					zx = atan2(zy, zx);
					val.push_back( dmod( halfPI - zx, twoPI) );
				}
				val.push_back(NAN);
			}

		} else {

			double xw[6] = {-1,-2,-1,1,2,1};
			double yw[6] = {-1,1,-2,2,-1,1};
			for (size_t i=0; i<6; i++) {
				xw[i] /= -8 * dx;
				yw[i] /= 8 * dy;
			}
			for (size_t row=1; row< (nrow-1); row++) {
				val.push_back(NAN);
				for (size_t col=1; col< (ncol-1); col++) {
					size_t i = row * ncol + col;
					zx = d[i-1-ncol] * xw[0] + d[i-1] * xw[1] + d[i-1+ncol] * xw[2]
						+ d[i+1-ncol] * xw[3] + d[i+1] * xw[4] + d[i+1+ncol] * xw[5];
					zy = d[i-1-ncol] * yw[0] + d[i-1+ncol] * yw[1] + d[i-ncol] * yw[2]
							+ d[i+ncol] * yw[3] + d[i+1-ncol] * yw[4] + d[i+1+ncol] * yw[5];
					zx = atan2(zy, zx);
					val.push_back( dmod( halfPI - zx, twoPI) );
				}
				val.push_back(NAN);
			}
		}
	}
	if (!after) {
		val.resize(val.size() + ncol, NAN);
	}

	if (degrees) {
		to_degrees(val, start);
	}
}


SpatRaster SpatRaster::terrain(std::vector<std::string> v, unsigned neighbors, bool degrees, unsigned seed, SpatOptions &opt) {

//TPI, TRI, aspect, flowdir, slope, roughness
	//std::sort(v.begin(), v.end());
	//v.erase(std::unique(v.begin(), v.end()), v.end());

	SpatRaster out = geometry(v.size());
	out.setNames(v);

	if (nlyr() > 1) {
		out.setError("terrain needs a single layer object");
		return out;
	}

	bool aspslope = false;
	std::vector<std::string> f {"TPI", "TRI", "TRIriley", "TRIrmsd", "aspect", "flowdir", "slope", "roughness"};
	for (size_t i=0; i<v.size(); i++) {
		if (std::find(f.begin(), f.end(), v[i]) == f.end()) {
			out.setError("unknown terrain variable: " + v[i]);
			return(out);
		}
		if (v[i] == "slope" || v[i] == "aspect") {
			aspslope=true;
		}
	}
	if ((v.size() == 1) && (v[0] == "flowdir")) {
		out.setValueType(1);
	}

	if ((neighbors != 4) && (neighbors != 8)) {
		out.setError("neighbors should be 4 or 8");
		return out;
	}
	bool lonlat = is_lonlat();
	double yr = yres();

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	opt.minrows = 3;
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	size_t nc = ncol();

	if (nrow() < 3 || nc < 3) {
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> val(out.bs.nrows[i] * nc, NAN);
			if (!out.writeBlock(val, i)) return out;
		}
		return out;
	}


	std::vector<double> y;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> d;
		bool before= false;
		bool after = false;
		size_t rrow = out.bs.row[i];
		size_t rnrw = out.bs.nrows[i];
		if (i > 0) {
			rrow--;
			rnrw++;
			before=true;
		}
		if ((out.bs.row[i] + out.bs.nrows[i]) < nrow()) {
			rnrw++;
			after = true;
		}
		readValues(d, rrow, rnrw, 0, nc);
		if (lonlat && aspslope) {
			std::vector<int64_t> rows(rnrw);
			std::iota(rows.begin(), rows.end(), rrow);
			y = yFromRow(rows);
			yr = distance_lonlat(0, 0, 0, yres());
		}
		std::vector<double> val;
		val.reserve(out.bs.nrows[i] * ncol() * v.size());
		for (size_t j =0; j<v.size(); j++) {
			if (v[j] == "slope") {
				do_slope(val, d, neighbors, rnrw, nc, xres(), yr, lonlat, y, degrees, before, after);
			} else if (v[j] == "aspect") {
				do_aspect(val, d, neighbors, rnrw, nc, xres(), yr, lonlat, y, degrees, before, after);
			} else if (v[j] == "flowdir") {
				double dx = xres();
				double dy = yres();
				if (lonlat) {
					double yhalf = yFromRow((size_t) nrow()/2);
					dx = distance_lonlat(0, yhalf, dx, yhalf);
					dy = distance_lonlat(0, 0, 0, dy);
				}
				do_flowdir(val, d, rnrw, nc, dx, dy, seed, before, after);
			} else if (v[j] == "roughness") {
				do_roughness(val, d, rnrw, nc, before, after);
			} else if (v[j] == "TPI") {
				do_TPI(val, d, rnrw, nc, before, after);
			} else if (v[j] == "TRI") {
				do_TRI(val, d, rnrw, nc, before, after);
			} else if (v[j] == "TRIriley") {
				do_TRI_riley(val, d, rnrw, nc, before, after);
			} else if (v[j] == "TRIrmsd") {
				do_TRI_rmsd(val, d, rnrw, nc, before, after);			
			} else {
				out.setError("?"); return out;
			}
		}
		if (!out.writeBlock(val, i)) return out;
	}
	out.writeStop();
	readStop();
	return out;
}




SpatRaster SpatRaster::hillshade(SpatRaster aspect, std::vector<double> angle, std::vector<double> direction, bool normalize, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if ((nlyr() != 1) || (aspect.nlyr() != 1)) {
		out.setError("slope and aspect should have one layer");
		return out;
	}
	if (angle.empty() || direction.empty()) {
		out.setError("you must provide a value for aspect and direction");
		return out;		
	}
	
	std::vector<std::string> nms;

	if ((angle.size() > 1) || (direction.size() > 1)) {
		recycle(angle, direction);
		recycle(direction, angle);	
		//nms = opt.names;
		SpatOptions ops(opt);
		ops.ncopies *= angle.size();
		size_t nl = angle.size();
		out.source.resize(nl);
		if (ops.names.size() == nl) {
			nms = opt.names;
		} else {
			nms.reserve(nl);
			for (size_t i=0; i<nl; i++) {
				std::string nmi = "hs_" + double_to_string(angle[i]) + "_" + double_to_string(direction[i]);
				nms.push_back(nmi);
			}
		}
				
		for (size_t i=0; i<nl; i++) {
			ops.names = {nms[i]};
			SpatRaster r = hillshade(aspect, {angle[i]}, {direction[i]}, normalize, ops);
			out.source[i] = r.source[0];
		}
		if (!opt.get_filename().empty()) {
			out = out.writeRaster(opt);
		}
		return out;
	} else {
		if ((opt.names.size() == 1) && !opt.names[0].empty()) {
			out.setNames(  {opt.names[0] } );
		} else {
			out.setNames( {"hillshade"} );
		}
	}

	double dir = direction[0] * M_PI / 180.0;
	double zen = (90.0 - angle[0]) * M_PI/180.0;
	double coszen = cos(zen);
	double sinzen = sin(zen);

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!aspect.readStart()) {
		out.setError(getError());
		return(out);
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> slp;
		std::vector<double> asp;
		readBlock(slp, out.bs, i);
		aspect.readBlock(asp, out.bs, i);
		if (normalize) {
			for (size_t i=0; i<slp.size(); i++) {
				slp[i] = cos(slp[i]) * coszen + sin(slp[i]) * sinzen * cos(dir-asp[i]);
				if (slp[i] < 0) {
					slp[i] = 0;
				} else {
					slp[i] *= 255;
				}
			}
		} else {
			for (size_t i=0; i<slp.size(); i++) {
				slp[i] = cos(slp[i]) * coszen + sin(slp[i]) * sinzen * cos(dir-asp[i]);
			}
		}
		if (!out.writeBlock(slp, i)) return out;
	}
	out.writeStop();
	readStop();
	aspect.readStop();
	return out;
}



