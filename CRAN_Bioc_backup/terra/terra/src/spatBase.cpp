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
#include "string_utils.h"
#include "math_utils.h"


SpatOptions::SpatOptions() {}


SpatOptions::SpatOptions(const SpatOptions &opt) {
	tempdir = opt.tempdir;
	memfrac = opt.memfrac;
	memmax = opt.memmax;
	memmin = opt.memmin;
	parallel = opt.parallel;
	todisk = opt.todisk;
	tolerance = opt.tolerance;

	def_datatype = opt.def_datatype;
	def_filetype = opt.def_filetype;
	filenames = {""};
	overwrite = false;
	progress = opt.progress;
	ncopies = opt.ncopies;
	verbose = opt.verbose;
	def_verbose = opt.def_verbose;
	statistics = opt.statistics;
	steps = opt.steps;
	minrows = opt.minrows;
	names = opt.names;
	//ncdfcopy = opt.ncdfcopy;
	gdal_options = opt.gdal_options;
	overwrite = opt.overwrite;
	hasNAflag = false;
	NAflag = NAN;
	datatype_set = opt.datatype_set;
	datatype = opt.datatype;
	filetype = opt.filetype;
	tmpfile = opt.tmpfile + "_2";
}

SpatOptions SpatOptions::deepCopy() {
	return *this;
}


//SpatOptions SpatOptions::deepCopy(const SpatOptions &opt) {
//	return SpatOptions(opt);
//}

//void SpatOptions::set_def_bandorder(std::string d) { def_bandorder = d; }
//std::string SpatOptions::get_def_bandorder() { return def_bandorder; }
//void SpatOptions::set_bandorder(std::string d) { bandorder = d; }
//std::string SpatOptions::get_bandorder() {if (bandorder != "") {return bandorder;} else {return def_datatype;}}

void SpatOptions::set_def_datatype(std::string d) {
#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 7
	std::vector<std::string> ss = {"INT1U", "INT2U", "INT4U", "INT8U", "INT2S", "INT4S", "INT8S", "FLT4S", "FLT8S"} ;
#else 
	std::vector<std::string> ss = {"INT1U", "INT2U", "INT4U", "INT8U", "INT1S", "INT2S", "INT4S", "INT8S", "FLT4S", "FLT8S"};
#endif
	if (is_in_vector(d, ss)) def_datatype = d;
}
std::string SpatOptions::get_def_datatype() { return def_datatype; }

void SpatOptions::set_datatype(std::string d) {
#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 7
	std::vector<std::string> ss = {"INT1U", "INT2U", "INT4U", "INT8U", "INT2S", "INT4S", "INT8S", "FLT4S", "FLT8S"} ;
#else 
	std::vector<std::string> ss = {"INT1U", "INT2U", "INT4U", "INT8U", "INT1S", "INT2S", "INT4S", "INT8S", "FLT4S", "FLT8S"};
#endif
	if (is_in_vector(d, ss)) {
		datatype = d;
		datatype_set = TRUE;
	} else {
		msg.addWarning(d + " is not a valid datatype");
	}
}
std::string SpatOptions::get_datatype() {if (datatype.empty()) {return def_datatype;} else {return datatype;}}

void SpatOptions::set_def_filetype(std::string d) { def_filetype = d; }
std::string SpatOptions::get_def_filetype() { return def_filetype;}

void SpatOptions::set_filetype(std::string d) { filetype = d; }
std::string SpatOptions::get_filetype() { return filetype;}

bool SpatOptions::get_overwrite() { return overwrite; }
void SpatOptions::set_overwrite(bool b) { overwrite = b; }

//bool SpatOptions::get_append() { return append; }
//void SpatOptions::set_append(bool b) { append = b; }

int SpatOptions::get_statistics() { return statistics; }
void SpatOptions::set_statistics(int s) { if ((s> 0) && (s<7)) statistics = s; }

//bool SpatOptions::get_ncdfcopy() { return ncdfcopy;}
//void SpatOptions::set_ncdfcopy(bool x) { ncdfcopy = x; }

void SpatOptions::set_def_verbose(bool v) { def_verbose = v; }
bool SpatOptions::get_def_verbose() { return def_verbose; }
bool SpatOptions::get_verbose() { return verbose; }
void SpatOptions::set_verbose(bool v) { verbose = v; }

bool SpatOptions::has_NAflag(double &flag) {
	flag = NAflag;
	return hasNAflag;
}

double SpatOptions::get_NAflag() {
	return NAflag;
}

void SpatOptions::set_NAflag(double flag) {
	NAflag = flag;
	hasNAflag = true;
}

size_t SpatOptions::get_progress() { return progress; }

void SpatOptions::set_progress(size_t p) {
	progress = p;
}

bool SpatOptions::show_progress(size_t n) {
	return ((progress > 0) & (progress <= n));
}


//void SpatOptions::set_filename(std::string f) {
//	f = lrtrim_copy(f);
//	filenames = {f};
//}

void SpatOptions::set_filenames(std::vector<std::string> f) {
	for (size_t i=0; i<f.size(); i++) {
		f[i] = lrtrim_copy(f[i]);
	}
	filenames = f;
}

std::string SpatOptions::get_filename() {
	if (!filenames.empty() ) {
		return filenames[0];
	} else {
		return "";
	}
}


std::vector<std::string> SpatOptions::get_filenames() {
	if (!filenames.empty() ) {
		return filenames;
	} else {
		return {""};
	}
}

std::string SpatOptions::get_tempdir() { return tempdir; }

void SpatOptions::set_tempdir(std::string d) {
	// check if exists?
	tempdir = d;
}

double SpatOptions::get_memfrac() { return memfrac; }

void SpatOptions::set_memfrac(double d) {
	// allowing very high values for testing purposes
	if ((d >= 0) && (d <= 100)) {
		memfrac = d;
	}
}

double SpatOptions::get_memmax() { return memmax; }

void SpatOptions::set_memmax(double d) {
	if (std::isnan(d) || (d <= 0)) {
		memmax = -1;
	} else {
		memmax = d * 1024 * 1024 * 1024 / 8;
	}
}

double SpatOptions::get_memmin() { return memmin; }

void SpatOptions::set_memmin(double d) {
	if (std::isnan(d) || (d <= 0)) {
		memmin = 1024 * 1024 * 1024 / 8;
	} else {
		memmin = d * 1024 * 1024 * 1024 / 8;
	}
}

double SpatOptions::get_tolerance() { return tolerance; }

void SpatOptions::set_tolerance(double d) {
	if (d > 0) {
		tolerance = d;
	}
}


bool SpatOptions::get_todisk() { return todisk; }
void SpatOptions::set_todisk(bool b) { todisk = b; }


void SpatOptions::set_steps(size_t n) { steps = std::max((size_t)1, n); }
size_t SpatOptions::get_steps(){ return steps; }

void SpatOptions::set_ncopies(size_t n) { ncopies = std::max((size_t)1, n); }
size_t SpatOptions::get_ncopies(){ return ncopies; }



void SpatOptions::set_offset(std::vector<double> d) { offset = d ; }
std::vector<double> SpatOptions::get_offset() {return offset;}

void SpatOptions::set_scale(std::vector<double> d) {scale=d;}
std::vector<double> SpatOptions::get_scale(){return scale;}


bool extent_operator(std::string oper) {
	std::vector<std::string> f {"==", "!=", ">", "<", ">=", "<="};
	return (std::find(f.begin(), f.end(), oper) != f.end());
}

bool SpatExtent::compare(SpatExtent e, std::string oper, double tolerance) {

	if (!extent_operator(oper)) {
		return false;  // not very useful
	}

	//double xr = (xmax - xmin) / tolerance;
	//double yr = (ymax - ymin) / tolerance;

	bool e1 = fabs(xmax - e.xmax) <= tolerance;
	bool e2 = fabs(xmin - e.xmin) <= tolerance;
	bool e3 = fabs(ymax - e.ymax) <= tolerance;
	bool e4 = fabs(ymin - e.ymin) <= tolerance;
	bool equal = (e1 && e2 && e3 && e4);
	if (oper == "==") {
		return equal;
	} else if (oper == "!=") {
		return (!equal);
	}
	if (oper == "<" || oper == "<=") {
		bool c1 = xmax < e.xmax;
		bool c2 = xmin > e.xmin;
		bool c3 = ymax < e.ymax;
		bool c4 = ymin > e.ymin;
		bool smaller = (c1 && c2 && c3 && c4);
		if (oper == "<") {
			return smaller;
		} else {
			return (equal || smaller);
		}
	}
	if (oper == ">" || oper == ">=") {
		bool c1 = xmax > e.xmax;
		bool c2 = xmin < e.xmin;
		bool c3 = ymax > e.ymax;
		bool c4 = ymin < e.ymin;
		bool larger = (c1 && c2 && c3 && c4);
		if (oper == ">") {
			return larger;
		} else {
			return (equal || larger);
		}
	}
	return false;
}

SpatExtent SpatExtent::round(int n) {
	double xn = roundn(xmin, n);
	double xx = roundn(xmax, n);
	double yn = roundn(ymin, n);
	double yx = roundn(ymax, n);
	SpatExtent e(xn, xx, yn, yx);
	return e;
}


SpatExtent SpatExtent::floor() {
	double xn = std::floor(xmin);
	double xx = std::ceil(xmax);
	double yn = std::floor(ymin);
	double yx = std::ceil(ymax);
	SpatExtent e(xn, xx, yn, yx);
	return e;
}

SpatExtent SpatExtent::ceil() {
	double xn = std::ceil(xmin);
	double xx = std::floor(xmax);
	double yn = std::ceil(ymin);
	double yx = std::floor(ymax);
	SpatExtent e(xn, xx, yn, yx);
	return e;
}

SpatExtent SpatRaster::getExtent() {
	if (source.empty()) {
		SpatExtent e;
		return e;
	} else {
		return source[0].extent;
	}
}



void SpatRaster::setExtent(SpatExtent e) {
	for (size_t i=0; i<nsrc(); i++) {
		source[i].extent = e;
		source[i].extset = true;
	}
}


void SpatRaster::setExtent(SpatExtent ext, bool keepRes, bool expand, std::string snap) {

	if (!snap.empty()) {
		ext = align(ext, snap);
	}
	if (!expand) {
		ext = ext.intersect(getExtent());
	}

	if (keepRes) {
		std::vector<double> res = resolution();
		double xrs = res[0];
		double yrs = res[1];
		unsigned nc = std::max(1.0, round( (ext.xmax - ext.xmin) / xrs ));
		unsigned nr = std::max(1.0, round( (ext.ymax - ext.ymin) / yrs ));
		ext.xmax = ext.xmin + nc * xrs;
		ext.ymax = ext.ymin + nr * yrs;
		for (size_t i=0; i<nsrc(); i++) {
			source[i].extent = ext;
			source[i].extset = true;
			source[i].nrow = nr;
			source[i].ncol = nc;
		}
	} else {
		for (size_t i=0; i<nsrc(); i++) {
			source[i].extent = ext;
			source[i].extset = true;
		}
	}
}


SpatExtent SpatExtent::align(double d, std::string snap) {
    std::vector<double> e = asVector();
	if (d == 0) {
		SpatExtent out = *this;
		return(out);
	}
	d = d < 0 ? -d : d;


	for (size_t i=0; i<4; i++) {
		double x = d * trunc(e[i] / d);
		if ((i == 0) | (i == 2)) {
			if (x > e[i]) {
				x -= d;
			}
		} else {
			if (x < e[i]) {
				x += d;
			}
		}
		e[i] = x;
	}
	SpatExtent out(e[0], e[1], e[2], e[3]);
	return(out)	;
}


SpatExtent SpatRaster::align(SpatExtent e, std::string snap) {

	snap = is_in_set_default(snap, std::vector<std::string> {"near", "in", "out"}, "near", true);
	std::vector<double> res = resolution();
	std::vector<double> orig = origin();

	// snap points to cell boundaries
	double xmn, xmx, ymn, ymx;
	if (snap == "near") {
		xmn = round((e.xmin-orig[0]) / res[0]) * res[0] + orig[0];
		xmx = round((e.xmax-orig[0]) / res[0]) * res[0] + orig[0];
		ymn = round((e.ymin-orig[1]) / res[1]) * res[1] + orig[1];
		ymx = round((e.ymax-orig[1]) / res[1]) * res[1] + orig[1];
	} else if (snap == "out") {
		xmn = std::floor((e.xmin-orig[0]) / res[0]) * res[0] + orig[0];
		xmx = std::ceil((e.xmax-orig[0]) / res[0]) * res[0] + orig[0];
		ymn = std::floor((e.ymin-orig[1]) / res[1]) * res[1] + orig[1];
		ymx = std::ceil((e.ymax-orig[1]) / res[1]) * res[1] + orig[1];
	} else { //if (snap == "in") {
		xmn = std::ceil((e.xmin-orig[0]) / res[0]) * res[0] + orig[0];
		xmx = std::floor((e.xmax-orig[0]) / res[0]) * res[0] + orig[0];
		ymn = std::ceil((e.ymin-orig[1]) / res[1]) * res[1] + orig[1];
		ymx = std::floor((e.ymax-orig[1]) / res[1]) * res[1] + orig[1];
		if (xmn > xmx) std::swap(xmn, xmx);
		if (ymn > ymx) std::swap(ymn, ymx);
	}


	if (xmn == xmx) {
		if (xmn < e.xmin) {
			xmx = xmx + res[0];
		} else {
			xmn = xmn - res[0];
		}
	}
	if (ymn == ymx) {
		if (ymn < e.ymin) {
			ymx = ymx + res[1];
		} else {
			ymn = ymn - res[1];
		}
	}
	return SpatExtent(xmn, xmx, ymn, ymx);
}


std::vector<double> SpatRaster::origin() {
	std::vector<double> r = resolution();
	SpatExtent extent = getExtent();
	double x = extent.xmin - r[0] * (round(extent.xmin / r[0]));
	double y = extent.ymax - r[1] * (round(extent.ymax / r[1]));
	if (is_equal((r[0] + x), abs(x))) {
		x = fabs(x);
	}
	if (is_equal((r[1] + y), abs(y))) {
		y = fabs(y);
	}
	std::vector<double> out {x, y};
	return out;
}


bool SpatRaster::compare_geom(SpatRaster &x, bool lyrs, bool crs, double tol, bool warncrs, bool ext, bool rowcol, bool res) {

	tol = tol < 0 ? 0 : (tol > 0.5 ? 0.5 : tol);

	if (ext) {
		SpatExtent extent = getExtent();
		double res = std::max(xres(), yres());
		if (extent.compare(x.getExtent(), "!=", tol * res)) {
			setError("extents do not match");
			return false;
		}
	}
	if (rowcol) {
		if (! ((nrow() == x.nrow()) && (ncol() == x.ncol())) ) {
			setError("number of rows and/or columns do not match");
			return false;
		}
	}
	if (res) {
		if (! ((is_equal_relative(x.xres(), xres(), 0.0001)) && (is_equal_relative(x.yres(), yres(), 0.0001)))) {
			setError("resolution does not match");
			return false;
		}
	}

	if (lyrs) {
		if (!(nlyr() == x.nlyr())) {
			setError("number of layers does not match" + std::to_string(nlyr()) + " != " + std::to_string(x.nlyr()));
			return false;
		}
	}

	if (crs) {
		if (!source[0].srs.is_equal(x.source[0].srs)) {
			if (warncrs) {
				addWarning("CRS do not match");
			} else {
				setError("CRS do not match");
				return false;
			}
		}
	}
	return true;
}


bool SpatCategories::combine(SpatCategories &x) {
	bool ok = d.rbind(x.d);
	if (!ok) {
		return(false);
	}
	d = d.unique();
	std::vector<long> ids = d.getI(0);
	size_t n = ids.size();
	std::sort(ids.begin(), ids.end());
	ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
	if (ids.size() < n) {
		return false;
	}
	return true;
}


bool SpatCategories::concatenate(SpatCategories &x) {

	std::vector<long> ids = d.getI(0);
	std::vector<long> xids = x.d.getI(0);
	std::vector<std::string> labs = d.as_string(index);
	std::vector<std::string> xlabs = x.d.as_string(x.index);
	size_t n = ids.size() * xids.size();
	std::vector<long> id1, id2;
	std::vector<std::string> news;
	id1.reserve(n);
	id2.reserve(n);
	news.reserve(n);
	std::string nm = d.names[index] + "_" + x.d.names[index];

	for (size_t i=0; i<ids.size(); i++) {
		for (size_t j=0; j<xids.size(); j++) {
			id1.push_back(ids[i]);
			id2.push_back(xids[j]);
			news.push_back(labs[i] + "_" + xlabs[j]);
		}
	}
	std::vector<long> id(n);
	std::iota(id.begin(), id.end(), 0);
	SpatDataFrame dd;
	dd.add_column(id, "ID");
	dd.add_column(news, nm);
	dd.add_column(id1, "idx");
	dd.add_column(id2, "idy");
	d = dd;
	return true;
}


#ifdef useRcpp

void SpatProgress::init(size_t n, int nmin) {

	if ((nmin <= 0) || ((int)n < nmin)) {
		show = false;
		return;
	} 

	show = true;

	std::string bar = "|---------|---------|---------|---------|";
	Rcpp::Rcout << "\r" << bar << "\r";
	R_FlushConsole();

	nstep = n;
	step = 0;
	size_t width = bar.size();

	double increment = (double) width / double(nstep);

	steps.resize(0);
	steps.reserve(nstep+1);
	for (size_t i=0; i<nstep; i++) {
		int val = round(i * increment);
		steps.push_back(val);
	}
	steps.push_back(width);
}

void SpatProgress::stepit() {
	if (show) {
		if (step < nstep) {
			int n = steps[step+1] - steps[step];
			if (n > 0) {
				for (int i=0; i<n; i++) Rcpp::Rcout << "=";
			}
		}
		step++;
		R_FlushConsole();
	}
}

void SpatProgress::finish() {
	if (show) {
		Rcpp::Rcout << "\r                                          \r";
		step++;
		R_FlushConsole();
	}
}

void SpatProgress::interrupt() {
	if (show) {
		Rcpp::Rcout << "\r                                          \r";
		R_FlushConsole();
	}
}

#else 

void SpatProgress::init(size_t n, int nmin) {}
void SpatProgress::stepit() {}
void SpatProgress::interrupt() {}

#endif
