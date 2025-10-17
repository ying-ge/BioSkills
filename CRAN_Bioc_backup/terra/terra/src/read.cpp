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

#include "spatRasterMultiple.h"

bool SpatRaster::readStart() {

	//if (!valid_sources(true, true)) {
	//	return false;
	//}

	for (size_t i=0; i<nsrc(); i++) {
		if (source[i].open_read) {
			addWarning("source already open for reading");
			continue;
		}
		if (source[i].memory) {
			source[i].open_read = true;
		} else if (source[i].is_multidim) {
			if (!readStartMulti(i)) {
				readStop();
				return false;
			}
		} else {
			if (!readStartGDAL(i)) {
				readStop();
				return false;
			}
		}
	}
	return true;
}

bool SpatRaster::readStop() {
	for (size_t i=0; i<nsrc(); i++) {
		if (source[i].open_read) {
			if (source[i].memory) {
				source[i].open_read = false;
			} else if (source[i].is_multidim) {
				readStopMulti(i);
			} else {
				readStopGDAL(i);
			}
		}
	}
	return true;
}


// 2D BSQ
void SpatRaster::readBlock2(std::vector<std::vector<double>> &v, BlockSize bs, size_t i) {
	std::vector<double> x;
	readValues(x, bs.row[i], bs.nrows[i], 0, ncol());
	v.resize(nlyr());
	size_t off = bs.nrows[i] * ncol();
	for (size_t i=0; i<nlyr(); i++) {
		v[i] = std::vector<double>(x.begin()+(i*off), x.begin()+((i+1)*off));
	}
}

// BIP
void SpatRaster::readBlockIP(std::vector<double> &x, BlockSize bs, size_t i) {
	readValues(x, bs.row[i], bs.nrows[i], 0, ncol());
	std::vector<double> v(x.size());
	size_t off = bs.nrows[i] * ncol();
	size_t nl = nlyr();
	for (size_t i=0; i<nl; i++) {
		std::vector<double> lyr = std::vector<double>(x.begin()+(i*off), x.begin()+((i+1)*off));
		for (size_t j=0; j<off; j++){
			size_t jj = j * nl + i;
			v[jj] = lyr[j];
		}
	}
	x = std::move(v);
}


void SpatRaster::readChunkMEM(std::vector<double> &out, size_t src, size_t row, size_t nrows, size_t col, size_t ncols){

	size_t nl = source[src].nlyr;

	if (source[src].hasWindow) {
		row += source[src].window.off_row;
		col += source[src].window.off_col;
		size_t endrow = row + nrows;
		size_t endcol = col + ncols;
		size_t nc = source[src].window.full_ncol;
		double ncells = source[src].window.full_nrow * nc;

		for (size_t lyr=0; lyr < nl; lyr++) {
			size_t add = ncells * lyr;
			for (size_t r = row; r < endrow; r++) {
				size_t off = add + r * nc;
				out.insert(out.end(), source[src].values.begin()+off+col, source[src].values.begin()+off+endcol);
			}
		}
			/*
			else if (source[0].window.expanded) {
				unsigned add = ncells * lyr;
				std::vector<double> v1(source[0].window.expand[0] * ncols, NAN);
				out.insert(out.end(), v1.begin(), v1.end());
				v1.resize(source[0].window.expand[1], NAN);
				std::vector<double> v2(source[0].window.expand[2], NAN);
				for (size_t r = wrow; r < endrow; r++) {
					unsigned a = add + r * source[0].window.full_ncol;
					out.insert(out.end(), v1.begin(), v1.end());
					out.insert(out.end(), source[src].values.begin()+a+wcol, source[src].values.begin()+a+endcol);
					out.insert(out.end(), v2.begin(), v2.end());
				}
				v1.resize(source[0].window.expand[3] * ncols, NAN);
				out.insert(out.end(), v1.begin(), v1.end());
			}
			*/

	} else { //	no window
		size_t nc = ncol();
		if (row==0 && nrows==nrow() && col==0 && ncols==nc) {
			out.insert(out.end(), source[src].values.begin(), source[src].values.end());
		} else {
			double ncells = ncell();
			if (col==0 && ncols==nc) {
				for (size_t lyr=0; lyr < nl; lyr++) {
					size_t add = ncells * lyr;
					size_t a = add + row * nc;
					size_t b = a + nrows * nc;
					out.insert(out.end(), source[src].values.begin()+a, source[src].values.begin()+b);
				}
			} else {
				size_t endrow = row + nrows;
				size_t endcol = col + ncols;
				for (size_t lyr=0; lyr < nl; lyr++) {
					size_t add = ncells * lyr;
					for (size_t r = row; r < endrow; r++) {
						size_t a = add + r * nc;
						out.insert(out.end(), source[src].values.begin()+a+col, source[src].values.begin()+a+endcol);
					}
				}
			}
		}
	}
}



std::vector<double> SpatRaster::readValuesR(size_t row, size_t nrows, size_t col, size_t ncols){

	std::vector<double> out;

	if (((row + nrows) > nrow()) || ((col + ncols) > ncol())) {
		setError("invalid rows/columns");
		return out;
	}


	//row = std::min(std::max(size_t(0), row), nrow()-1);
	//col = std::min(std::max(size_t(0), col), ncol()-1);
	//nrows = std::max(size_t(1), std::min(nrows, nrow()-row));
	//ncols = std::max(size_t(1), std::min(ncols, ncol()-col));
	if ((nrows==0) | (ncols==0)) {
		return out;
	}

	if (!hasValues()) {
		out.resize(nrows * ncols * nlyr(), NAN);
		addWarning("raster has no values");
		return out; // or NAs?
	}


	unsigned n = nsrc();

	out.reserve(nrows * ncols * nlyr());
	for (size_t src=0; src<n; src++) {
		if (source[src].memory) {
			readChunkMEM(out, src, row, nrows, col, ncols);
		} else {
			// read from file
			#ifdef useGDAL

/*
				if (source[0].window.expanded) {
					std::vector<double> gout;
					readChunkGDAL(gout, src, source[0].window.off_row, nrows, source[0].window.off_col, ncols);

					size_t rrow = row + source[0].window.off_row;
					size_t rcol = col + source[0].window.off_col;
					unsigned endrow = rrow + nrows;
					unsigned endcol = rcol + ncols;
					unsigned ncells = source[0].window.full_nrow * source[0].window.full_ncol;
					unsigned nl = source[src].nlyr;

					for (size_t lyr=0; lyr < nl; lyr++) {
						unsigned add = ncells * lyr;
						std::vector<double> v1(source[0].window.expand[0] * ncols, NAN);
						out.insert(out.end(), v1.begin(), v1.end());
						v1.resize(source[0].window.expand[1], NAN);
						std::vector<double> v2(source[0].window.expand[2], NAN);
						for (size_t r = rrow; r < endrow; r++) {
							unsigned a = add + r * source[0].window.full_ncol;
							out.insert(out.end(), v1.begin(), v1.end());
							out.insert(out.end(), gout.begin()+a+rcol, gout.begin()+a+endcol);
							out.insert(out.end(), v2.begin(), v2.end());
						}
						v1.resize(source[0].window.expand[3] * ncols, NAN);
						out.insert(out.end(), v1.begin(), v1.end());
					}
				}
				*/
			readChunkGDAL(out, src, row, nrows, col, ncols);
			#endif // useGDAL
		}
	}
	return out;
}

void SpatRaster::readValues(std::vector<double> &out, size_t row, size_t nrows, size_t col, size_t ncols){

	if (((row + nrows) > nrow()) || ((col + ncols) > ncol())) {
		setError("invalid rows/columns");
		return;
	}

	if ((nrows==0) | (ncols==0)) {
		return;
	}

	out.resize(0);
	if (!hasValues()) {
		out.resize(nrows * ncols * nlyr(), NAN);
		addWarning("raster has no values");
		return; // or NAs?
	}

	unsigned n = nsrc();
	out.reserve(nrows * ncols * nlyr());
	for (size_t src=0; src<n; src++) {
		if (source[src].memory) {
			readChunkMEM(out, src, row, nrows, col, ncols);
		} else {
			// read from file
			#ifdef useGDAL
			readChunkGDAL(out, src, row, nrows, col, ncols);
			#endif // useGDAL
		}
	}
	return;
}


void SpatRaster::readValuesWhileWriting(std::vector<double> &out, size_t row, size_t nrows, size_t col, size_t ncols){

	if (((row + nrows) > nrow()) || ((col + ncols) > ncol())) {
		setError("invalid rows/columns");
		return;
	}

	if ((nrows==0) | (ncols==0)) {
		return;
	}

	unsigned n = nsrc();
	out.resize(0);
	out.reserve(nrows * ncols * nlyr());

	for (size_t src=0; src<n; src++) {
		if (source[src].memory) {
			readChunkMEM(out, src, row, nrows, col, ncols);
		} else {
			// read from file
			#ifdef useGDAL
			readChunkGDAL(out, src, row, nrows, col, ncols);
			#endif // useGDAL
		}
	}
	return;
}



bool SpatRaster::readAll() {
	if (!hasValues()) {
		return true;
	}

	size_t row =0, col=0, nrows=nrow(), ncols=ncol();
	readStart();
	size_t n = nsrc();
	for (size_t src=0; src<n; src++) {
		if (!source[src].memory) {
			readChunkGDAL(source[src].values, src, row, nrows, col, ncols);
			source[src].memory = true;
			source[src].extset = false;
			source[src].flipped = false;
			source[src].filename = "";
			std::iota(source[src].layers.begin(), source[src].layers.end(), 0);			
		}
		if (src > 0) {
			if (!source[0].combine_sources(source[src])) {
				setError("could not combine sources");
				return false;
			}
			source[src].values.resize(0);
		}
	}
	readStop();
	if (n > 1) source.resize(1);

	source[0].hasWindow = false;
	return true;
}



std::vector<double> SpatRaster::getValues(long lyr, SpatOptions &opt) {

	std::vector<double> out;

	bool hw = false;
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].hasWindow) {
			hw = true;
			break;
		}
	}
	if (hw) {
		if (!readStart()) return out;
		if (lyr < 0) { // default; read all
			readValues(out, 0, nrow(), 0, ncol());
		} else {
			unsigned lyrnr = lyr;
			SpatRaster sub = subset({lyrnr}, opt);
			sub.readValues(out, 0, nrow(), 0, ncol());
		}
		readStop();
		return out;
	}

	if (lyr < 0) { // default; read all
		unsigned n = nsrc();
		for (size_t src=0; src<n; src++) {
			if (source[src].memory) {
				out.insert(out.end(), source[src].values.begin(), source[src].values.end());
			} else {
				#ifdef useGDAL
				std::vector<double> fvals = readValuesGDAL(src, 0, nrow(), 0, ncol());
				out.insert(out.end(), fvals.begin(), fvals.end());
				#endif // useGDAL
			}
		}
	} else { // read one lyr
		std::vector<size_t> sl = findLyr(lyr);
		unsigned src=sl[0];
		if (source[src].memory) {
			size_t start = sl[1] * ncell();
			out = std::vector<double>(source[src].values.begin()+start, source[src].values.begin()+start+ncell());
		} else {
			#ifdef useGDAL
			out = readValuesGDAL(src, 0, nrow(), 0, ncol(), sl[1]);
			#endif // useGDAL
		}
	}
	return out;
}


bool SpatRaster::getValuesSource(size_t src, std::vector<double> &out) {

	unsigned n = nsrc();
	if (src > n) {
		return false;
	}

	bool hw = false;
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].hasWindow) {
			hw = true;
			break;
		}
	}
	if (hw) {
		SpatRaster sub = SpatRaster(source[src]);
		if (!readStart()) return false;
		sub.readValues(out, 0, nrow(), 0, ncol());
		readStop();
		return true;
	}


	if (source[src].memory) {
		out = std::vector<double>(source[src].values.begin(), source[src].values.end());
	} else {
		#ifdef useGDAL
		out = readValuesGDAL(src, 0, nrow(), 0, ncol());
		#endif // useGDAL
	}
	return true;
}


SpatRasterStack SpatRasterCollection::read_into(SpatRaster &tmp, size_t row, size_t nrows) {

	size_t n = size();
	double nan = NAN;
	std::vector<double> v(nan, nrows * tmp.ncol() * n);
	SpatRasterStack out;

	SpatExtent e = tmp.getExtent();
	e.ymax = tmp.source[0].extent.ymax - row * tmp.yres(); 	
	e.ymin = tmp.source[0].extent.ymax - (row + nrows) * tmp.yres(); 
	SpatOptions ops;
	for (size_t i=0; i<n; i++) {
		SpatExtent ee = ds[i].getExtent();
		if ((ee.ymax > e.ymin) & (ee.ymin < e.ymax)) {
			if (!tmp.compare_geom(ds[i], false, false, ops.get_tolerance(), false, false, false, true)) {
				out.setError(tmp.msg.error);
				return(out);
			}
			SpatRaster x = ds[i].crop(e, "near", true, ops);
			out.ds.push_back(x);
		}
	}
	return out;
}


