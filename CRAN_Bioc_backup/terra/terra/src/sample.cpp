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

#include <vector>
#include "spatRaster.h"
#include "distance.h"
#include "recycle.h"
#include <random>
#include <unordered_set>
#include "string_utils.h"
#include "geodesic.h"
#include "table_utils.h"
#include "sort.h"

std::vector<std::vector<double>> SpatRaster::sampleStratifiedCells(double size, bool each, bool replace, unsigned seed, SpatOptions &opt) {

	std::vector<std::vector<double>> out;

	if (!hasValues()) {
		setError("raster has no values");
		return out;
	}

	if (nlyr() > 1) {
		SpatRaster r = subset({0}, opt);
		addWarning("only the first layer of the raster is used");		
		return r.sampleStratifiedCells(size, each, replace, seed, opt);
	}

	std::default_random_engine gen1(seed);

	size_t szz=size;
	std::set<double> add_set;
	if (!each) {
		std::vector<std::vector<double>> u = unique(false, NAN, true, opt);
		size_t nuv = u[0].size();
		szz = floor(size / nuv);
		size_t d = size - szz * nuv;
		if (d > 0) {
			std::shuffle(u[0].begin(), u[0].end(), gen1);
			u[0].erase(u[0].begin()+d, u[0].end());
			add_set = std::set<double>(u[0].begin(), u[0].end());
		} 
	} 

	if (!readStart()) {
		return(out);
	}

	opt.ncopies *= 2;
	BlockSize bs = getBlockSize(opt);
	std::vector<double> vals, vcell, vwght, outvals, outcell;
	
	for (size_t i=0; i<bs.n; i++) {
		std::vector<double> v;
		readBlock(v, bs, i);
		std::vector<double> cells(v.size());
		std::iota(cells.begin(), cells.end(), bs.row[i] * ncol());
		std::vector<std::size_t> pm = sort_order_nan_a(v);
		permute(v, pm);
		permute(cells, pm);
		std::map<double, size_t> tab = table(v);
		std::vector<std::vector<double>> tv = table2vector2(tab);
		size_t start = 0;
		for (size_t j=0; j<tv[0].size(); j++) {
			size_t size_j = add_set.count(tv[0][j]) ? szz + 1 : szz;
			if (size_j == 0) continue;
			std::vector<size_t> z;
			z.resize(tv[1][j]);
			std::iota(z.begin(), z.end(), 0);
			if (tv[1][j] > size_j) {
				std::shuffle(z.begin(), z.end(), gen1);
				z.erase(z.begin()+size_j, z.end());
			}
			double weight_j = tv[1][j] / z.size();
			for (size_t k=0; k<z.size(); k++) {
				vals.push_back(tv[0][j]);
				vwght.push_back(weight_j);
				vcell.push_back(cells[start+z[k]]);
			}
			start += tv[1][j];		
		}
	}
	readStop();

	if (bs.n == 1) {
		out.push_back(vcell);
		out.push_back(vals);
	} else {
		std::vector<std::size_t> pm = sort_order_a(vals);
		permute(vals, pm);
		permute(vwght, pm);
		permute(vcell, pm);
		std::map<double, size_t> tab = table(vals);
		std::vector<std::vector<double>> tv = table2vector2(tab);
		size_t start = 0;
				
		std::mt19937 gen2(seed);
		for (size_t j=0; j<tv[0].size(); j++) {
			size_t size_j = add_set.count(tv[0][j]) ? szz + 1 : szz;
			if (size_j == 0) continue;
			size_t end = start + tv[1][j];
			std::vector<size_t> z;
			z.resize(tv[1][j]);
			std::iota(z.begin(), z.end(), 0);
			if (tv[1][j] > size_j) {
				std::discrete_distribution<int> dist(vwght.begin()+start, vwght.begin()+end);
				std::vector<size_t> Z;
				if (replace) {
					Z.resize(size_j);
					for(auto & i: Z) i = dist(gen2);
				} else {
					std::unordered_set<size_t> z;
					while (z.size() < size_j) {
						z.insert(dist(gen2));
					}
					Z = std::vector<size_t>(z.begin(), z.end());
				}
				for (size_t k=0; k<z.size(); k++) {
					outvals.push_back(tv[0][j]);
					outcell.push_back(vcell[Z[k]+start]);
				}
			} else {
				outvals.insert(outvals.end(), vals.begin()+start, vals.begin()+end);
				outcell.insert(outcell.end(), vcell.begin()+start, vcell.begin()+end);				
			}
			start = end;		
		}
		out.push_back(outcell);
		out.push_back(outvals);
	}
	return(out);
}

	


void get_nx_ny(double size, size_t &nx, size_t &ny) {
	double nxy = nx * ny;
	if (std::isfinite(size) && (!std::isnan(size)) && (size < nxy)) {
		double f = sqrt(size / nxy );
		double fnx = nx * f;
		double fny = ny * f;
		double s = fnx * fny;
		f = size / s;
		nx = std::max((size_t)1, (size_t) std::ceil(fnx * f));
		ny = std::max((size_t)1, (size_t) std::ceil(fny * f));
	}
}


void getSampleRowCol(std::vector<size_t> &oldrow, std::vector<size_t> &oldcol, size_t nrows, size_t ncols, size_t snrow, size_t sncol) {

	double rf = nrows / (double)(snrow);
	double cf = ncols / (double)(sncol);
	//double rstart = std::floor(0.5 * rf);
	//double cstart = std::floor(0.5 * cf);
	double rstart = 0.5 * rf;
	double cstart = 0.5 * cf;
	oldcol.reserve(sncol);
	for (size_t i =0; i<sncol; i++) {
        oldcol.push_back(i * cf + cstart);
	}
	oldrow.reserve(snrow);
	for (size_t i =0; i<snrow; i++) {
        oldrow.push_back(i * rf + rstart);
	}
}


std::vector<double> SpatRaster::readSample(size_t src, size_t srows, size_t scols) {

	size_t nl = source[src].nlyr;
	std::vector<size_t> oldcol, oldrow;
	std::vector<double>	out;
	getSampleRowCol(oldrow, oldcol, nrow(), ncol(), srows, scols);

	out.reserve(srows*scols);
	if (source[src].hasWindow) {
		size_t offrow = source[src].window.off_row;
		size_t offcol = source[src].window.off_col;
		size_t fncol = source[src].window.full_ncol;
		size_t oldnc = fncol * source[src].window.full_nrow;
		for (size_t lyr=0; lyr<nl; lyr++) {
			size_t off1 = lyr * oldnc;
			for (size_t r=0; r<srows; r++) {
				size_t off2 = off1 + (oldrow[r]+offrow) * fncol;
				for (size_t c=0; c<scols; c++) {
					size_t oldcell = off2 + oldcol[c] + offcol;
					out.push_back(source[src].values[oldcell]);
				}
			}
		}
	} else {
		size_t oldnc = ncell();
		for (size_t lyr=0; lyr<nl; lyr++) {
			size_t off = lyr * oldnc;
			for (size_t r=0; r<srows; r++) {
				size_t oldc = off + oldrow[r] * ncol();
				for (size_t c=0; c<scols; c++) {
					size_t oldcell = oldc + oldcol[c];
					out.push_back(source[src].values[oldcell]);
				}
			}
		}
	}
	return out;
}


SpatRaster SpatRaster::sampleRegularRaster(double size, bool overview) {

	if (size >= ncell()) {
		return( *this );
	}
	if (size < 0.5) {
		SpatRaster out;
		out.setError("sample size must be > 0");
		return out;
	}
	
	double f = std::min(1.0, sqrt(size / ncell()));
	size_t nr = std::min((size_t)ceil(nrow() * f), nrow());
	size_t nc = std::min((size_t)ceil(ncol() * f), ncol());
	if ((nc == ncol()) && (nr == nrow())) {
		return( *this );
	}
	SpatRaster out = geometry(nlyr(), true);
	out.source[0].nrow = nr;
	out.source[0].ncol = nc;

	std::vector<int> vt = getValueType(true);
	if (vt.size() == 1) {
		out.setValueType(vt[0]);
	}

	if (!source[0].hasValues) return (out);

	std::vector<double> v;
	for (size_t src=0; src<nsrc(); src++) {
		if (source[src].memory) {
			v = readSample(src, nr, nc);
		//} else if (source[src].driver == "raster") {
		//	v = readSampleBinary(src, nr, nc);
		} else {
		    #ifdef useGDAL
			v = readGDALsample(src, nr, nc, overview);
			#endif
		}
		if (hasError()) return out;
		out.source[0].values.insert(out.source[0].values.end(), v.begin(), v.end());
	}
	out.source[0].memory = true;
	out.source[0].hasValues = true;
	out.source[0].setRange();

	return out;
}


SpatRaster SpatRaster::sampleRowColRaster(size_t nr, size_t nc, bool warn) {

	SpatRaster out = geometry(nlyr(), true);
	if ((nr == 0) || (nc ==0)) {
		out.setError("number of rows and columns must be > 0");
	}

	if (nr > nrow()) {
		if (warn) out.addWarning("number of rows cannot be larger than nrow(x)");
		nr = nrow();
	}
	if (nc > ncol()) {
		if (warn) out.addWarning("number of rows cannot be larger than nrow(x)");
		nc = ncol();
	}

	if ((nc == ncol()) && (nr == nrow())) {
		return( *this );
	}

	out.source[0].nrow = nr;
	out.source[0].ncol = nc;

	std::vector<int> vt = getValueType(true);
	if (vt.size() == 1) {
		out.setValueType(vt[0]);
	}

	if (!source[0].hasValues) return (out);

	std::vector<double> v;
	for (size_t src=0; src<nsrc(); src++) {
		if (source[src].memory) {
			v = readSample(src, nr, nc);
		//} else if (source[src].driver == "raster") {
		//	v = readSampleBinary(src, nr, nc);
		} else {
		    #ifdef useGDAL
			v = readGDALsample(src, nr, nc, false);
			#endif
		}
		if (hasError()) return out;
		out.source[0].values.insert(out.source[0].values.end(), v.begin(), v.end());
	}
	out.source[0].memory = true;
	out.source[0].hasValues = true;
	out.source[0].setRange();

	return out;
}



std::vector<std::vector<double>> SpatRaster::sampleRegularValues(double size, SpatOptions &opt) {

	std::vector<std::vector<double>> out;
	if (!source[0].hasValues) return (out);


	size_t nr = nrow();
	size_t nc = ncol();
	get_nx_ny(size, nc, nr);
/*	
	if (size < ncell()) {
		double f = sqrt(size / ncell());
		double nr1 = nrow() * f;
		double nc1 = ncol() * f;
		double s = nr1 * nc1;
		f = size / s;
		nr = std::max((size_t)1, (size_t) std::ceil(nr1 * f));
		nc = std::max((size_t)1, (size_t) std::ceil(nc1 * f));
	}
*/	

	size_t nsize = nc * nr;
	std::vector<double> v;
	if ((size >= ncell()) || ((nc == ncol()) && (nr == nrow()))) {
		v = getValues(-1, opt) ;
		if (hasError()) return out;
		for (size_t i=0; i<nlyr(); i++) {
			size_t offset = i * nsize;
			std::vector<double> vv(v.begin()+offset, v.begin()+offset+nsize);
			out.push_back(vv);
		}
		return out;
	}

	for (size_t src=0; src<nsrc(); src++) {
		if (source[src].memory) {
			v = readSample(src, nr, nc);
		//} else if (source[src].driver == "raster") {
		//	v = readSampleBinary(src, nr, nc);
		} else {
		    #ifdef useGDAL
			v = readGDALsample(src, nr, nc, false);
			#endif
		}
		if (hasError()) return out;
		for (size_t i=0; i<source[src].nlyr; i++) {
			size_t offset = i * nsize;
			std::vector<double> vv(v.begin()+offset, v.begin()+offset+nsize);
			out.push_back(vv);
		}
	}
	return out;
}



std::vector<double> SpatRaster::sampleRowCol(size_t nr, size_t nc) {

	std::vector<double> out;
	if (!source[0].hasValues) return (out);
	if ((nr == 0) || (nc ==0)) {
		return(out);
	}
	nr = std::min(nr, nrow());
	nc = std::min(nc, ncol());
	
	if ((nr == nrow()) && (nc == ncol())) {
		out.resize(ncell());
		std::iota(out.begin(), out.end(), 0);
		return out;
	}

	std::vector<int64_t> rows, cols;

	double d = nrow()/nr;
	for (size_t i=0; i<nr; i++) {
		rows.push_back((i + .5) * d);
	}
	d = ncol()/nc;
	for (size_t i=0; i<nc; i++) {
		cols.push_back((i + .5) * d);
	}
	return cellFromRowColCombine(rows, cols);
}


std::vector<std::vector<double>> SpatRaster::sampleRowColValues(size_t nr, size_t nc, SpatOptions &opt) {
	std::vector<std::vector<double>> out;

	if ((nr >= nrow()) && (nc >= ncol())) {
		std::vector<double> v = getValues(-1, opt);
		size_t nsize = ncell();
		if (hasError()) return out;
		for (size_t i=0; i<nlyr(); i++) {
			size_t offset = i * nsize;
			std::vector<double> vv(v.begin()+offset, v.begin()+offset+nsize);
			out.push_back(vv);
		}
		return out;
	}
	
	std::vector<double> cells = sampleRowCol(nr, nc);
	if (cells.empty()) return out;
	return extractCell(cells, opt);
}


std::vector<size_t> sample_replace(size_t size, size_t N, unsigned seed){
	std::default_random_engine gen(seed);
	std::uniform_int_distribution<> U(0, N-1);
	std::vector<size_t> sample;
	sample.reserve(size);
	for (size_t i=0; i<size; i++) {
		sample.push_back( U(gen) );
	}
	return sample;
}

std::vector<size_t> sample_replace_weights(size_t size, size_t N, std::vector<double> prob, unsigned seed){
	std::discrete_distribution<int> dist(std::begin(prob), std::end(prob));
	std::mt19937 gen;
	gen.seed(seed);
	std::vector<size_t> sample(size);
	for(auto & i: sample) i = dist(gen);
	return sample;
}


std::vector<size_t> sample_no_replace(size_t size, size_t N, unsigned seed){
	size_t one = 1;
	size = std::max(one, std::min(size, N));
	std::vector<size_t> sample;
	if (size == N) {
		sample.resize(size);
		std::iota(sample.begin(), sample.end(), 0);
		return sample;
	}
	std::default_random_engine gen(seed);

	if (size >= .66 * N) {
		sample.resize(N);
		std::iota(std::begin(sample), std::end(sample), 0);
		std::shuffle(sample.begin(), sample.end(), gen);
		if (size < N) {
			sample.erase(sample.begin()+size, sample.end());
		}
		return sample;
	}

	std::uniform_real_distribution<> U( 0, std::nextafter(1.0, std::numeric_limits<double>::max() ) );

	sample.reserve(size);
	for (size_t i=0; i<N; i++) {
		if ( ((N-i) * U(gen)) <  (size - sample.size()) ) {
			sample.push_back(i);
    		if (sample.size() == size ) {
				break;
			}
		}
	}
	return sample;
}


std::vector<size_t> sample_no_replace_weights(size_t size, size_t N, std::vector<double> prob, unsigned seed){
	size_t one = 1;
	size = std::max(one, std::min(size, N));
	std::vector<size_t> sample;
	std::default_random_engine gen(seed);
	if (size == N) {
		sample.resize(size);
		std::iota(sample.begin(), sample.end(), 0);
		std::shuffle(sample.begin(), sample.end(), gen);
		return sample;
	}

	std::uniform_int_distribution<> U(0, std::numeric_limits<int>::max());
	std::unordered_set<size_t> sampleset;

	size_t isize = size;
	if (size > (0.8 * N)) {
		isize = N - size;
		for (double &d : prob) d = 1-d;
		size_t ssize = isize * (1.1 + isize / N);
		size_t cnt=0;
		while (sampleset.size() < isize) {
			seed = U(gen);
			std::vector<size_t> s = sample_replace_weights(ssize, N, prob, seed);
			for (size_t i=0; i<s.size(); i++) {
				sampleset.insert(s[i]);
			}
			cnt++;
			if (cnt > 10) break;
		}
		std::vector<size_t> invsamp;
		invsamp.insert(invsamp.begin(), sampleset.begin(), sampleset.end());
		std::sort(invsamp.begin(), invsamp.end());
		invsamp.push_back(N+1);
		size_t j=0;
		sample.reserve(size);
		for (size_t i=0; i<N; i++) {
			if (i != invsamp[j]) {
				sample.push_back(i);
			} else {
				j++;
			}
		}
		std::shuffle(sample.begin(), sample.end(), gen);

	} else {
		size_t ssize = size * (1.1 + size / N);
		size_t cnt=0;
		while (sampleset.size() < size) {
			seed = U(gen);
			std::vector<size_t> s = sample_replace_weights(ssize, N, prob, seed);
			for (size_t i=0; i<s.size(); i++) {
				sampleset.insert(s[i]);
			}
			cnt++;
			if (cnt > 10) break;
		}
		sample.insert(sample.begin(), sampleset.begin(), sampleset.end());
		if (sample.size() > size) {
			sample.resize(size);
		};
	}

	return(sample);
}


std::vector<size_t> sample(size_t size, size_t N, bool replace, std::vector<double> prob, unsigned seed){
	if ((size == 0) || (N == 0)) {
		std::vector<size_t> s;
		return s;
	}
	bool w = prob.size() == N;
	if (replace) {
		if (N == 1) {
			std::vector<size_t> s(size,0);
			return s;
		}
		if (w) {
			return sample_replace_weights(size, N, prob, seed);
		} else {
			return sample_replace(size, N, seed);
		}
	} else {
		if (N == 1) {
			std::vector<size_t> s(1,0);
			return s;
		}
		if (w) {
			return sample_no_replace_weights(size, N, prob, seed);
		} else {
			return sample_no_replace(size, N, seed);
		}
	}
}



std::vector<std::vector<double>> SpatRaster::sampleRandomValues(double size, bool replace, unsigned seed) {

	double nc = ncell();
	std::vector<size_t> cells;
	std::vector<double> w;
	if (replace) {
		cells = sample((size_t)size, nc, false, w, seed);
	} else {
		cells = sample((size_t)size, nc, true, w, seed);
	}

	std::vector<double> dcells(cells.begin(), cells.end());
	SpatOptions opt;
	std::vector<std::vector<double>> d = extractCell(dcells, opt);
	return d;
}


SpatRaster SpatRaster::sampleRandomRaster(double size, bool replace, unsigned seed) {

	unsigned nsize;
	unsigned nr = nrow();
	unsigned nc = ncol();
	if (size < ncell()) {
		double f = sqrt(size / ncell());
		nr = std::ceil(nrow() * f);
		nc = std::ceil(ncol() * f);
	}
	SpatRaster out = geometry(nlyr(), true);
	out.source[0].nrow = nr;
	out.source[0].ncol = nc;
	if (!source[0].hasValues) return (out);

	nsize = nr * nc;
	std::vector<std::vector<double>> vv = sampleRandomValues(nsize, replace, seed);

	for (size_t i=0; i<vv.size(); i++) {
		out.source[0].values.insert(out.source[0].values.end(), vv[i].begin(), vv[i].end());
	}
	out.source[0].memory = true;
	out.source[0].hasValues = true;
	out.source[0].setRange();

	return out;
}

std::vector<size_t> SpatExtent::test_sample(size_t size, size_t N, bool replace, std::vector<double> w, unsigned seed) {
	return sample(size, N, replace, w, seed);
}


std::vector<std::vector<double>> SpatExtent::sampleRandom(size_t size, bool lonlat, unsigned seed){
	std::vector<std::vector<double>> out(2);
	if (size == 0) return out;
	std::default_random_engine gen(seed);


	if (lonlat) {
		double d = (ymax - ymin) / 1000.0;
		std::vector<double> r = seq(ymin, ymax, d);
		std::vector<double> w;
		w.reserve(r.size());
		for (size_t i=0; i<r.size(); i++) {
			double ww = std::abs(cos(M_PI * r[i]/180.0));
			w.push_back(ww );
		}

		std::vector	<size_t> x = sample(size, r.size(), true, w, seed);
		std::vector <double> lat, lon;
		lat.reserve(size);
		lon.reserve(size);
		std::uniform_real_distribution<> U1(-0.5, 0.5);

		double dx = 0.5 * d;
		for (size_t i=0; i<x.size(); i++) {
			double v = r[x[i]] + dx * U1(gen);
			lat.push_back(v);
		}
		std::uniform_real_distribution<> U2(xmin, xmax);
		for (size_t i=0; i<size; i++) {
			lon.push_back(U2(gen));
		}
		out[0] = lon;
		out[1] = lat;

	} else {
		std::vector <double> x, y;
		x.reserve(size);
		y.reserve(size);
		std::uniform_real_distribution<> runifx(xmin, xmax);
		std::uniform_real_distribution<> runify(ymin, ymax);
		for (size_t i=0; i<size; i++) {
			x.push_back(runifx(gen));
			y.push_back(runify(gen));
		}
		out[0] = x;
		out[1] = y;
	}
	return out;
}




std::vector<std::vector<double>> SpatExtent::sampleRegular(size_t size, bool lonlat) {
	std::vector<std::vector<double>> out(2);
	if (size == 0) return out;

	double r1 = xmax - xmin;
	double r2 = ymax - ymin;

	if (lonlat) {
		double halfy = ymin + r2/2;
		// beware that -180 is the same as 180; and that latitude can only go from -90:90 therefore:
		double dx = distance_lonlat(xmin, halfy, xmin + 1, halfy) * std::min(180.0, r1);
		double dy = distance_lonlat(0, ymin, 0, ymax);

		double ratio = dy/dx;
		double n = sqrt(size);
		double ny = n * ratio;
		double nx = n / ratio;
		double s = nx * ny;
		ratio = size / s;
		ny = std::max((size_t)1, (size_t) std::ceil(ny * ratio));
		nx = std::max((size_t)1, (size_t) std::ceil(nx * ratio));

		double x_i = r1 / nx;
		double y_i = r2 / ny;

		std::vector<double> lat, lon, w, xi;
		lat.reserve(ny);
		lat.push_back(ymin+0.5*y_i);
		for (size_t i=1; i<ny; i++) {
			lat.push_back(lat[i-1] + y_i);
		}

		w.reserve(lat.size());
		for (size_t i=0; i<lat.size(); i++) {
			w.push_back(cos(M_PI * lat[i] / 180.0));
		}

		double nwsumw = w.size() / accumulate(w.begin(), w.end(), 0.0);
		xi.reserve(w.size());
		for (size_t i=0; i<w.size(); i++) {
			xi.push_back(x_i / (w[i] * nwsumw));
		}
		bool global = (xmax - xmin) > 355; // needs refinement
		if (global) {
			xmax -= 0.000001;
			for (size_t i=0; i<lat.size(); i++) {
				size_t n = std::max(1, (int)(360.0/xi[i]));
				double step = 360.0 / n;
				std::vector<double> x = seq(xmin+0.5*step, xmax, step);
				std::vector<double> y(x.size(), lat[i]);
				out[0].insert(out[0].end(), x.begin(), x.end());
				out[1].insert(out[1].end(), y.begin(), y.end());
			}

		} else {
			double halfx = xmin + (xmax - xmin)/2;
			for (size_t i=0; i<lat.size(); i++) {
				std::vector<double> x = seq(halfx, xmax, xi[i]);
				double start = halfx-xi[i];
				if (start > xmin) {
					std::vector <double> x2 = seq(start, xmin, -xi[i]);
					x.insert(x.end(), x2.begin(), x2.end());
				}
				std::vector<double> y(x.size(), lat[i]);
				out[0].insert(out[0].end(), x.begin(), x.end());
				out[1].insert(out[1].end(), y.begin(), y.end());
			}
		}
	} else {
		double ratio = r1/r2;
		double ny = sqrt(size / ratio);
		double nx = size / ny;
		double s = nx * ny;
		ratio = size / s;
		ny = std::max((size_t)1, (size_t) std::ceil(ny * ratio));
		nx = std::max((size_t)1, (size_t) std::ceil(nx * ratio));

		double x_i = r1 / nx;
		double y_i = r2 / ny;

		std::vector<double> x, y;
		x.reserve(nx);
		y.reserve(ny);
		x.push_back(xmin+0.5*x_i);
		for (size_t i=1; i<nx; i++) {
			x.push_back(x[i-1] + x_i);
		}
		y.reserve(ny);
		y.push_back(ymin+0.5*y_i);
		for (size_t i=1; i<ny; i++) {
			y.push_back(y[i-1] + y_i);
		}
		rep(x, ny);
		rep_each(y, nx);
		out[0] = x;
		out[1] = y;
	}

	return out;
}


/*
std::vector<size_t> SpatRaster::sampleCells(double size, std::string method, bool replace, unsigned seed) {

	std::default_random_engine gen(seed);
	std::vector<size_t> out;
	if ((size >= ncell()) & (!replace)) {
		out.resize(ncell());
		std::iota(out.begin(), out.end(), 0);
		if (method == "random") {
			std::shuffle(out.begin(), out.end(), gen);
		}
		return out;
	}

	if (method == "random") {

	} else if (method == "regular") {

	} else { //method == "stratified"

	} // else "Cluster"
	return out;
}
*/

SpatVector SpatVector::sample(size_t n, std::string method, unsigned seed) {

	std::string gt = type();
	SpatVector out;
	if (gt == "points") {
		out.setError("only implemented for lines and polygons");
		return out;
	}
	if (n == 0) {
		out.srs = srs;
		return out;
	}
	bool lonlat = is_lonlat();
	bool random = (method == "random");

	if (gt == "lines") {
		
		std::vector<double> x, y;
		std::vector<double> steps;
		steps.reserve(n);		
		std::vector<double> pp = length();
		double p = std::accumulate(pp.begin(), pp.end(), 0.0); 

		if (random) {
			std::default_random_engine gen(seed);
			std::uniform_real_distribution<> U2(0, p);
			for (size_t i=0; i<n; i++) {
				steps.push_back(U2(gen));
			}
			std::sort(steps.begin(), steps.end());
		} else {
			double d = p/n;
			for (size_t i=0; i<n; i++) {
				steps.push_back((i + .5) * d);
			}
		}
		size_t k = 0;
		double length = 0;
		double oldlength = 0;
		if (lonlat) {
			x.resize(n);
			y.resize(n);
			struct geod_geodesic gd;
			double a = 6378137;
			double f = 1 / 298.257223563;
			geod_init(&gd, a, f);
			double azi1, azi2, dist;
			for (size_t g=0; g<geoms.size(); g++) {
				for (size_t i=0; i<geoms[g].parts.size(); i++) {
					for (size_t j=1; j<geoms[g].parts[i].x.size(); j++) {
						geod_inverse(&gd, geoms[g].parts[i].y[j-1], geoms[g].parts[i].x[j-1], 
											 geoms[g].parts[i].y[j], geoms[g].parts[i].x[j], &dist, &azi1, &azi2);
						length += dist;
						while (length > steps[k]) {
							geod_direct(&gd, geoms[g].parts[i].y[j-1], geoms[g].parts[i].x[j-1], azi1, 
										steps[k]-oldlength, &y[k], &x[k], &azi2);
							k++;
							if (k == n) { goto endloop; }
						}
						oldlength = length;
					}
				}
			}
		} else {
			x.reserve(n);
			y.reserve(n);

			for (size_t g=0; g < geoms.size(); g++) {
				for (size_t i=0; i < geoms[g].parts.size(); i++) {
					for (size_t j=1; j<geoms[g].parts[i].x.size(); j++) {
						length += sqrt(pow(geoms[g].parts[i].x[j-1] - geoms[g].parts[i].x[j], 2) +
									   pow(geoms[g].parts[i].y[j-1] - geoms[g].parts[i].y[j], 2));
						while (length > steps[k]) {
							double bearing = direction_plane(geoms[g].parts[i].x[j-1], geoms[g].parts[i].y[j-1], 
															geoms[g].parts[i].x[j], geoms[g].parts[i].y[j], false);
							double distance = steps[k]-oldlength;
							x.push_back(geoms[g].parts[i].x[j-1] + distance * sin(bearing));
							y.push_back(geoms[g].parts[i].y[j-1] + distance * cos(bearing));
							k++;
							if (k == n) { goto endloop; }
						}
						oldlength = length;
					}
				}
			}
		}
		endloop:
		out = SpatVector(x, y, points, "");
		out.srs = srs;
		return out;
	}

/*
	if (strata != "") {
		// should use
		// SpatVector a = aggregate(strata, false);
		// but get nasty self-intersection precision probs.

		int i = where_in_vector(strata, get_names());
		if (i < 0) {
			out.setError("cannot find field");
			return out;
		}
		SpatDataFrame uv;
		std::vector<int> idx = df.getIndex(i, uv);
		for (size_t i=0; i<uv.nrow(); i++) {
			std::vector<int> g;
			g.resize(0);
			for (size_t j=0; j<idx.size(); j++) {
				if (i == (size_t)idx[j]) {
					g.push_back(j);
				}
			}
			SpatVector s = subset_rows(g);
			s = s.sample(n, "random", false, "", seed);
			for (long &v : s.df.iv[0]) v = v+i;
			out = out.append(s, true);
		}
		return out;
	}
*/

	std::vector<double> a = area("m", true, {});
	if (hasError()) {
		out.setError(getError());
		return out;
	}
	double suma = accumulate(a.begin(), a.end(), 0.0);

/*
	if (by_geom) {
		std::vector<double> pa;
		pa.reserve(a.size());
		for (size_t i=0; i<a.size(); i++) {
			pa.push_back(a[i] / suma);
		}
		std::vector<std::vector<double>> pxy(2);

		std::vector<size_t> nsamp(size());
		for (size_t i=0; i<size(); i++) {
			if (pa[i] > 0) {
				SpatGeom g = getGeom(i);
				SpatVector ve(g.extent, "");
				ve.srs = srs;
				double vea = ve.area()[0];
				if (random) {
					double m = vea / a[i];
					m = std::max(2.0, std::min(m*m, 100.0));
					size_t ssize = pa[i] * n * m;
					pxy = g.extent.sampleRandom(ssize, lonlat, seed);
				} else {
					size_t ssize = std::round(pa[i] * n * vea / a[i]);
					pxy = g.extent.sampleRegular(ssize, lonlat);
				}
				SpatVector vpnt(pxy[0], pxy[1], points, "");
				SpatVector vpol(g);
				vpnt = vpnt.intersect(vpol);
				if (random) {
					size_t psize = pa[i] * n;
					if (vpnt.size() > psize) {
						std::vector<int> rows(psize);
						std::iota(rows.begin(), rows.end(), 0);
						vpnt = vpnt.subset_rows(rows);
					}
				}
				nsamp[i] = vpnt.size();
				if (out.size() == 0) {
					out = vpnt;
				} else {
					out = out.append(vpnt, true);
				}
			}
		}
		std::vector<long> id(size());
		std::iota(id.begin(), id.end(), 1);
		rep_each_vect(id, nsamp);
		SpatDataFrame df;
		df.add_column(id, "pol.id");
		out.df = df;
	} else {
*/
	std::vector<std::vector<double>> pxy(2);

	SpatVector ve(extent, "");
	ve.srs = srs;
	double vea = ve.area("m", true, {})[0];
	if (random) {
		double m = vea / suma;
	// the larger the sample size, the fewer extra samples needed
		double smx = sqrt(std::max(9.0, 100.0 - n));
		m = std::max(smx, std::min(m*m, 100.0));
		size_t ssize = n * m;
		pxy = extent.sampleRandom(ssize, lonlat, seed);
	} else {
		size_t ssize = std::round(n * vea / suma);
		pxy = extent.sampleRegular(ssize, lonlat);
	}
	out = SpatVector(pxy[0], pxy[1], points, "");
	out = intersect(out, true);
	if (random) {
		if (out.size() > n) {
			std::vector<long> rows(out.size());
			std::iota(rows.begin(), rows.end(), 0);
			std::default_random_engine gen(seed);
			std::shuffle(rows.begin(), rows.end(), gen);
			rows.resize(n);
			out = out.subset_rows(rows);
		}
	}
	//std::vector<long> id(out.size(), 1);
	//SpatDataFrame df;
	//df.add_column(id, "pol.id");
	//out.df = df;
//	}
	out.srs = srs;

	return out;
}


SpatVector SpatVector::sample_geom(std::vector<size_t> n, std::string method, unsigned seed) {

	SpatVector out;
	if (n.size() != size()) {
		out.setError("length of samples does not match number of geoms");
		return out;
	}
	if (n.empty()) {
		out.srs = srs;
		return out;
	}

	for (size_t i=0; i<size(); i++) {
		if (n[i] == 0) {
			continue;
		}
		int j = i;
		SpatVector v = subset_rows(j);
		SpatVector p = v.sample(n[i], method, seed+i);
		out = out.append(p, true);
	}
	out.srs = srs;
	return out;
}

/*
std::vector<double> sample(size_t size, size_t N, bool replace, std::vector<double> prob, unsigned seed){
    // Sample "size" elements from [1, N]
	std::vector<double> result;
	std::default_random_engine gen(seed);

	bool weights = false;
	if (prob.size() == N) {
		weights = true;
		// should check for neg numbers
		double minw = *min_element(prob.begin(),prob.end());
		double maxw = *max_element(prob.begin(),prob.end()) - minw;
		for (double& d : prob)  d = (d - minw) / maxw;
	}

	if (replace) {
		//std::vector<double> samples;
		result.reserve(size);
		std::uniform_int_distribution<> distribution(0, N-1);
		if (weights) {
			std::uniform_real_distribution<> wdist(0, 1);
			size_t cnt = 0;
			while (cnt < size) {
				double w = wdist(gen);
				double v = distribution(gen);
				if (prob[v] >= w) {
					result.push_back(v);
					cnt++;
				}
			}
		} else {
			for (size_t i=0; i<size; i++) {
				result.push_back( distribution(gen) );
			}
		}
	} else {//without replacement

		size = std::max(size_t(1), std::min(size, N)); // k <= N

		std::uniform_int_distribution<> distribution(1, N);
		std::unordered_set<unsigned> samples;

		if (weights) {
			std::uniform_int_distribution<> wdist(0, N-1);
			size_t cnt = 0;
			size_t r = 0;
			while (cnt < size) {
				double w = wdist(gen)/N;
				double v = distribution(gen);
				if (prob[v] >= w) {
					if (!samples.insert(v).second) {
						samples.insert(r);
						cnt++;
					}
				}
				r++;
				r = r%(N-1);
			}
		} else {
			for (size_t r = N - size; r < N; ++r) {
				unsigned v = distribution(gen) - 1;
				if (!samples.insert(v).second) samples.insert(r);
			}
		}
		result = std::vector<double>(samples.begin(), samples.end());
		std::shuffle(result.begin(), result.end(), gen);
	}
    return result;
}
*/

