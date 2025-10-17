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
#include "file_utils.h"
#include "string_utils.h"
#include "math_utils.h"
#include "recycle.h"


bool SpatRaster::writeValuesMem(std::vector<double> &vals, size_t startrow, size_t nrows) {

	//if (source[0].has_scale_offset[0]) {
	//	for (double &d : vals) d = d * source[0].scale[0] + source[0].offset[0];
	//}

	if (vals.size() == size()) {
		source[0].values = std::move(vals);
		return true;
	}

	if (nlyr() == 1) {
		source[0].values.insert(source[0].values.end(), vals.begin(), vals.end());
		return true;
	}

	if (source[0].values.empty()) { // && startrow != 0 && startcol != 0) {
		source[0].values = std::vector<double>(size(), NAN);
	}

	size_t nc = ncell();
	size_t ncols = ncol();
	size_t chunk = nrows * ncols;
	for (size_t i=0; i<nlyr(); i++) {
		size_t off1 = i * chunk;
		size_t off2 = startrow * ncols + i * nc;
		std::copy( vals.begin()+off1, vals.begin()+off1+chunk, source[0].values.begin()+off2 );
	}
	return true;
}

bool SpatRaster::writeValuesMemRect(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols) {

	if (source[0].values.empty()) { // && startrow != 0 && startcol != 0) {
		source[0].values = std::vector<double>(size(), NAN);
	}

	size_t nc = ncell();
	size_t chunk = nrows * ncols;

	for (size_t i=0; i<nlyr(); i++) {
		unsigned off = i*chunk;
		for (size_t r=0; r<nrows; r++) {
			size_t start1 = r * ncols + off;
			size_t start2 = (startrow+r)*ncol() + i*nc + startcol;
			std::copy(vals.begin()+start1, vals.begin()+start1+ncols, source[0].values.begin()+start2);
		}
	}
	return true;
}



void SpatRaster::fill(double x) {
	if (source[0].driver == "gdal") {
		#ifdef useGDAL
		fillValuesGDAL(x);
		#endif
	} else {
		source[0].values.resize(size(), x);
	}

}



bool SpatRaster::isSource(std::string filename) {
	std::vector<std::string> ff = filenames();
	for (size_t i=0; i<ff.size(); i++) {
		if (ff[i] == filename) {
			return true;
		}
	}
	return false;
}


SpatRaster SpatRaster::writeRaster(SpatOptions &opt) {

	SpatRaster out = geometry_opt(nlyr(), true, true, true, true, true, opt);
	if (!hasValues()) {
		out.setError("there are no cell values");
		return out;
	}

	// recursive writing of layers
	std::vector<std::string> fnames = opt.get_filenames();
	std::string msg;

	size_t nl = nlyr();
	if (fnames.size() > 1) {
		if (fnames.size() != nl) {
			out.setError("the number of filenames should either be one, or equal to the number of layers");
			return out;
		} else {
			bool overwrite = opt.get_overwrite();
			std::string errmsg;
			if (!can_write(fnames, filenames(), overwrite, errmsg)) {
				out.setError(errmsg);
				return(out);
			}
			for (unsigned i=0; i<nl; i++) {
				opt.set_filenames({fnames[i]});
				SpatRaster out = subset({i}, opt);
				if (out.hasError()) {
					return out;
				}
				fnames[i] = out.source[0].filename;
			}
			SpatRaster out(fnames, {-1}, {""}, false, {}, {}, {}, false, false, {});
			return out;
		}
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	opt.ncopies = 2;
	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v;
		readBlock(v, out.bs, i);
		if (!out.writeBlock(v, i)) {
			readStop();
			out.writeStop();
			return out;
		}
	}
	out.writeStop();
	readStop();
	return out;
}


SpatRaster SpatRaster::writeTempRaster(SpatOptions &opt) {
	SpatOptions xopt(opt);
	std::string fname = tempFile(xopt.get_tempdir(), xopt.tmpfile, "_temp_raster.tif");
	xopt.set_filenames({fname});
	return writeRaster(xopt);
}


bool SpatRaster::writeStart(SpatOptions &opt, const std::vector<std::string> srcnames) {

	if (opt.names.size() == nlyr()) {
		setNames(opt.names);
	}

	std::vector<std::string> fnames = opt.get_filenames();
	if (fnames.size() > 1) {
		addWarning("only the first filename supplied is used");
	}
	std::string filename = fnames[0];
	if (filename.empty()) {
		if (!canProcessInMemory(opt)) {
			//std::string extension = ".tif";
			//filename = tempFile(opt.get_tempdir(), opt.pid, extension);
			std::string driver;
			if (!getTempFile(filename, driver, opt)) {
				return false;
			}
			opt.set_filenames({filename});
			//opt.gdal_options = {"COMPRESS=NONE"};
		}
	}

	if (!opt.tags.empty()) {
		user_tags = std::vector<std::vector<std::string>>(3);
		for (size_t i=0; i<opt.tags.size(); i++) {
			std::vector<std::string> s = strsplit(opt.tags[i], "_#_");
			if (s.size() == 3) {
				addTag(s[0], s[1], s[2]);
			}
		}
	} 
	
	size_t nl = nlyr();
	bs = getBlockSize(opt);
	if (!filename.empty()) {
		// open GDAL filestream
		#ifdef useGDAL
		if (! writeStartGDAL(opt, srcnames) ) {
			return false;
		}
		#else
		setError("GDAL is not available");
		return false;
		#endif
	} else if ((nl == 1) && (bs.n > 1)) {
		source[0].values.reserve(ncell());
	}

	if (source[0].open_write) {
		addWarning("file was already open");
	}
	source[0].open_write = true;
	source[0].filename = filename;
	//bs = getBlockSize(opt);
    #ifdef useRcpp
	if (opt.verbose) {
		std::vector<double> mems = mem_needs(opt);
		double gb = 1073741824 / 8;
		//{memneed, memavail, frac, csize, inmem} ;
		// << "max vect size : " << roundn(mems.max_size() / gb, 2) << " GB" << std::endl;
		Rcpp::Rcout<< "memory avail. : " << roundn(mems[1] / gb, 2) << " GB" << std::endl;
		Rcpp::Rcout<< "memory allow. : " << roundn(mems[2] * mems[1] / gb, 2) << " GB" << std::endl;
		Rcpp::Rcout<< "memory needed : " << roundn(mems[0] / gb, 3) << " GB" << "  (" << opt.ncopies << " copies)" << std::endl;
		std::string inmem = mems[4] < 0.5 ? "false" : "true";
		Rcpp::Rcout<< "in memory     : " << inmem << std::endl;
		Rcpp::Rcout<< "block size    : " << mems[3] << " rows" << std::endl;
		Rcpp::Rcout<< "n blocks      : " << bs.n << std::endl;
		Rcpp::Rcout<< "pb            : " << opt.get_progress() << std::endl << std::endl;
	}

	if (opt.progressbar) {
		pbar.init(bs.n, opt.get_progress());
		progressbar = true;
	} else {
		progressbar = false;
	}
	#endif
	return true;
}

#ifdef useRcpp

static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}
bool checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

#endif


bool SpatRaster::writeValues(std::vector<double> &vals, size_t startrow, size_t nrows) {
	bool success = true;

	if (!source[0].open_write) {
		setError("cannot write (no open file)");
		return false;
	}

	if ((startrow + nrows) > nrow()) {
		setError("incorrect start and/or nrows value");
		return false;
	}

	size_t nv = nrows * ncol() * nlyr();
	if (vals.size() != nv) {
		if (vals.size() > nv) {
			setError("too many values for writing: " + std::to_string(vals.size()) + " > " + std::to_string(nv));
		} else {
			setError("too few values for writing: " + std::to_string(vals.size()) + " < " + std::to_string(nv));
		}
		return false;
	}

	if (source[0].driver == "gdal") {
		#ifdef useGDAL
		success = writeValuesGDAL(vals, startrow, nrows, 0, ncol());
		#else
		setError("GDAL is not available");
		return false;
		#endif
	} else {
		success = writeValuesMem(vals, startrow, nrows);
	}

//return success;
#ifdef useRcpp
	if (checkInterrupt()) {
		pbar.interrupt();
		setError("interrupted");
		return(false);
	}
	if (progressbar) {
		pbar.stepit();
	}
#endif
	return success;
}


bool SpatRaster::writeValuesRect(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols) {
	bool success = true;

	if (!source[0].open_write) {
		setError("cannot write (no open file)");
		return false;
	}

	if ((startrow + nrows) > nrow()) {
		setError("incorrect start and/or nrows value");
		return false;
	}

	if (source[0].driver == "gdal") {
		#ifdef useGDAL

		success = writeValuesGDAL(vals, startrow, nrows, startcol, ncols);
		#else
		setError("GDAL is not available");
		return false;
		#endif
	} else {
		success = writeValuesMemRect(vals, startrow, nrows, startcol, ncols);
	}

#ifdef useRcpp
	if (checkInterrupt()) {
		pbar.interrupt();
		setError("aborted");
		return(false);
	}
	if (progressbar) {
		pbar.stepit();
	}
#endif
	return success;
}


bool SpatRaster::writeValuesRectRast(SpatRaster &r, SpatOptions& opt) {
	bool success = true;

	if (!compare_geom(r, false, false, opt.get_tolerance(), false, false, false, true)) {
		return(false);
	}
	double hxr = xres() / 2;
	double hyr = yres() / 2;

	SpatExtent e = r.getExtent();
	int64_t row1  = rowFromY(e.ymax - hyr);
	int64_t row2  = rowFromY(e.ymin + hyr);
	int64_t col1  = colFromX(e.xmin + hxr);
	int64_t col2  = colFromX(e.xmax - hxr);
	if ((row1 < 0) || (row2 < 0) || (col1 < 0) || (col2 < 0)) {
		setError("block outside raster");
		return(false);		
	}
	size_t ncols = col2-col1+1;
	size_t nrows = row2-row1+1;
	size_t startrow = row1;
	size_t startcol = col1;
	
	if ((startrow + nrows) > nrow()) {
		setError("incorrect start row and/or nrows value");
		return false;
	}
	if ((startcol + ncols) > ncol()) {
		setError("incorrect start col and/or ncols value");
		return false;
	}
	if (!source[0].open_write) {
		setError("cannot write (no open file)");
		return false;
	}
	std::vector<double> vals = r.getValues(-1, opt);
	recycle(vals, ncols * nrows * nlyr());

	if ((nrows * ncols * nlyr()) != vals.size()) {
		setError("incorrect row/col size");
		return false;
	}


	if (source[0].driver == "gdal") {
		#ifdef useGDAL
		success = writeValuesGDAL(vals, startrow, nrows, startcol, ncols);
		#else
		setError("GDAL is not available");
		return false;
		#endif
	} else {
		success = writeValuesMemRect(vals, startrow, nrows, startcol, ncols);
	}

#ifdef useRcpp
	if (checkInterrupt()) {
		pbar.interrupt();
		setError("aborted");
		return(false);
	}
	if (progressbar) {
		pbar.stepit();
	}
#endif
	return success;
}



/*
bool SpatRaster::writeValues2(std::vector<std::vector<double>> &vals, size_t startrow, size_t nrows) {
    std::vector<double> vv = flatten(vals);
    return writeValues(vv, startrow, nrows, 0, ncol());
}
*/

bool SpatRaster::writeStop(){
	if (!source[0].open_write) {
		setError("cannot close a file that is not open");
		return false;
	}
	source[0].open_write = false;
	bool success = true;
	source[0].memory = false;
	if (source[0].driver=="gdal") {
		#ifdef useGDAL
		success = writeStopGDAL();
		//source[0].hasValues = true;
		#else
		return false;
		#endif
	} else {
   		source[0].setRange();
		//source[0].driver = "memory";
		source[0].memory = true;
		if (!source[0].values.empty()) {
			source[0].hasValues = true;
		}
	}

#ifdef useRcpp
	if (progressbar) {
		pbar.finish();
	}
/*
	if (progressbar) {
		pbar->increment();
		pbar->cleanup();
		delete pbar;
	}
*/
#endif

	return success;
}


#ifdef useRcpp
bool SpatRaster::setValuesRcpp(Rcpp::NumericVector &v, SpatOptions &opt) {
	SpatRaster g = geometry(nlyr(), true, true, true);
	source = g.source;
	source[0].hasValues = true;
	source[0].memory = true;
	//source[0].names = getNames();
	source[0].driver = "memory";

	if (v.size() < g.size()) {
		std::vector<double> vv = Rcpp::as<std::vector<double> >(v);
		*this = g.init(vv, opt);
		return (!hasError());
	} else if (v.size() == g.size()) {
		source[0].values = Rcpp::as<std::vector<double> >(v);
		source[0].setRange();
	} else {
		setError("incorrect number of values");
		return false;
	}
	return true;

}
#endif


bool SpatRaster::setValues(std::vector<double> &v, SpatOptions &opt) {

	SpatRaster g = geometry(nlyr(), true, true, true);
	source = g.source;
	source[0].hasValues = true;
	source[0].memory = true;
	//source[0].names = getNames();
	source[0].driver = "memory";

	if (v.size() < g.size()) {
		*this = g.init(v, opt);
		return (!hasError());
	} else if (v.size() == g.size()) {
		source[0].values = v;
		source[0].setRange();
	} else {
		setError("incorrect number of values");
		return false;
	}
	return true;
}

void SpatRaster::setRange(SpatOptions &opt, bool force) {

	for (size_t i=0; i<nsrc(); i++) {
		if (source[i].hasRange[0] && (!force)) continue;
		if (source[i].memory) {
			source[i].setRange();
		} else {
			SpatRaster r(source[i]);
			SpatDataFrame x = r.global("range", true, opt);
			source[i].range_min = x.getD(0);
			source[i].range_max = x.getD(1);
			source[i].hasRange = std::vector<bool>(source[i].hasRange.size(), true);
		}
	}
}


void SpatRasterSource::setRange() {
	range_min.resize(nlyr);
	range_max.resize(nlyr);
	hasRange.resize(nlyr);
	if (nlyr==1) {
		minmax(values.begin(), values.end(), range_min[0], range_max[0], NAN);
		hasRange[0] = true;
		return;
	}
	size_t nc = ncol * nrow;
	if (values.size() == (nc * nlyr)) {
		for (size_t i=0; i<nlyr; i++) {
			size_t start = nc * i;
			minmax(values.begin()+start, values.begin()+start+nc, range_min[i], range_max[i], NAN);
			hasRange[i] = true;
		}
	}
}



/*
bool SpatRaster::writeStartFs(std::string filename, std::string format, std::string datatype, bool overwrite,  fstream& f) {

	lrtrim(filename);
	if (filename == "") {
		if (!canProcessInMemory(4)) {
			filename = "random_file_name.grd";
		}
	}

	if (filename == "") {
		source.driver = {"memory"};

	} else {
		// if (!overwrite) check if file exists
		string ext = getFileExt(filename);
		lowercase(ext);
		if (ext == ".grd") {
			source.driver = {"raster"};
			string fname = setFileExt(filename, ".gri");
			f.open(fname, ios::out | ios::binary);
			(*fs).open(fname, ios::out | ios::binary);
			fs = &f;
		} else {
			source.driver = {"gdal"} ;
			// open GDAL filestream
		}
	}

	source.filename = {filename};
	bs = getBlockSize(4);
	return true;
}
*/


bool SpatRaster::writeDelim(std::string filename, std::string delim, bool cell, bool xy, SpatOptions &opt) {

	if (!hasValues()) {
		setError("there are no cell values");
		return false;
	}
	if (!readStart()) {
		setError(getError());
		return(false);
	}

	std::ofstream f;
	f.open(filename);
	if (!f.is_open()) {
		setError("could not open the csv file for writing");
		return false;
	}
	std::vector<std::string> nms = getNames();
	if (xy | cell) {
		std::vector<std::string> add;
		if (xy) {
			add.push_back("x");
			add.push_back("y");
		} 
		if (cell) {
			add.push_back("cell");
		}
		nms.insert(nms.begin(), add.begin(), add.end());
	}

	std::string s = concatenate(nms, delim);
	f << s << std::endl;

	BlockSize bs = getBlockSize(opt);
	for (size_t i=0; i<bs.n; i++) {
		std::vector<double> v;
		readBlock(v, bs, i);
		//s = get_delim_string(v, delim);
		//f << s << std::endl;
	}
	f.close();
	readStop();
	return true;
}

