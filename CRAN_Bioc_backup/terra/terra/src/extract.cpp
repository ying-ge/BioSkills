// Copyright (c) 2018-2025  Robert J. Hijmans
//
// This file is part of the "spat" library
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

#include <functional>

#include "spatRasterMultiple.h"
#include "vecmath.h"
#include "vecmathse.h"
#include "geosphere.h"
#include "sort.h"
#include "math_utils.h"


void SpatRaster::readRowColBlock(size_t src, std::vector<std::vector<double>> &out, size_t outstart, std::vector<int64_t> &rows, std::vector<int64_t> &cols, SpatOptions &opt) {

	std::vector<std::vector<double>> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return;
	}
	
	SpatRaster rs(source[src]);
	size_t nl = rs.nlyr();
	size_t nc = rs.ncol();
	size_t n = rows.size();

	BlockSize bs = getBlockSize(opt);

/*
	std::vector<std::pair<int64_t, int64_t>> rowcol;
	rowcol.reserve(rows.size());
	for (size_t i=0; i < rows.size(); i++) {
		rowcol.push_back(std::make_pair(rows[i], cols[i]));
	}
	std::vector<std::size_t> pm = sort_order_a(rowcol);
	permute(rows, pm);
	permute(cols, pm);
*/

//for (size_t i=0; i<bs.n; i++) Rcpp::Rcout << bs.row[i] << " " << bs.nrows[i] << "; ";
//Rcpp::Rcout << std::endl;

	std::vector<int64_t> urows = vunique(rows);
	std::vector<bool> useblock(bs.n, false);
	size_t jj = 0;
	for (size_t i=0; i<bs.n; i++) {
		int64_t st = bs.row[i];
		int64_t ed = bs.row[i] + bs.nrows[i] - 1;
		for (size_t j=jj; j<urows.size(); j++) {
			if ((urows[j] >= st) && (urows[j] <= ed)) {
				useblock[i] = true;
				bs.row[i] = urows[j];
				for (size_t k=j; k<urows.size(); k++) {
					if (urows[k] > ed) {
						jj = k;
						break;
					} else {
						bs.nrows[i] = urows[k]-urows[j]+1;
					}
				}
				break;
			}
		}
	}
	bs.nrows[bs.n-1] = urows[urows.size()-1] - bs.row[bs.n-1] + 1;

//for (size_t i=0; i<bs.n; i++) Rcpp::Rcout << bs.row[i] << " " << bs.nrows[i] << "; ";
//Rcpp::Rcout << std::endl;


	if (!rs.readStart()) {
		setError(getError());
		return;
	}

	size_t outend  = outstart + nl;

//	std::vector<std::vector<double>> out(nl, std::vector<double>(n, NAN));
	for (size_t k=outstart; k<outend; k++) {
		out[k] = std::vector<double>(n, NAN);
	}

	for (size_t i=0; i<bs.n; i++) {
		if (!useblock[i]) continue;
		std::vector<double> v;
		rs.readBlock(v, bs, i);
		int64_t rstart = bs.row[i];
		int64_t rend = bs.row[i] + bs.nrows[i];
		size_t off1 = bs.nrows[i] * nc;
		for (size_t j=0; j<n; j++) {
//			if (rows[j] >= rend) break; // if rows are sorted
			if ((rows[j] >= rstart) && (rows[j] < rend)) {
				size_t cell = (rows[j]-rstart) * nc + cols[j];
				for (size_t lyr=0; lyr<nl; lyr++) {
					out[outstart + lyr][j] = v[cell+lyr * off1]; 
				}
			}
		}
	}
	rs.readStop();
/*
	pm = sort_order_a(pm);
	for (size_t i=0; i<out.size(); i++) {
		permute(out[i], pm);
	}
*/	
//	return out;
}

/*
std::vector<double> SpatRaster::readRowColBlockFlat(size_t src, std::vector<int64_t> &rows, std::vector<int64_t> &cols) {

	std::vector<double> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return errout;
	}
	std::vector<std::vector<double>> v = readRowColBlock(src, rows, cols);
	
	size_t nr = v[0].size();
	size_t nl = v.size();
	std::vector<double> out;
	out.reserve(nl * nr);

	for (size_t i=0; i<nr; i++) {
		for (size_t j=0; j<nl; j++) {
			out.push_back(v[j][i]);
		}
	}
	return out;
}
*/


std::vector<double> circ_dist(double xres, double yres, double d, size_t nrows, size_t ncols, std::vector<size_t> &dim, bool lonlat, double ymean) {
	
	size_t nx, ny;
	std::string crs;
	if (lonlat) {
		deg2rad(ymean);
		double xr = xres;
		double yr = yres;
		deg2rad(xr);
		deg2rad(yr);
		double dx = distance_cos(0, ymean, xr, ymean);
		double dy = distance_cos(0, ymean-0.5*yr, 0, ymean+0.5*yr);
		nx = 1 + 2 * floor(d/dx);
		ny = 1 + 2 * floor(d/dy);
		crs = "+proj=longlat";
	} else {
		nx = 1 + 2 * floor(d/xres);
		ny = 1 + 2 * floor(d/yres);
		crs = "+proj=utm +zone=1";
	}
	nx = std::min(nx, ncols);
	ny = std::min(ny, nrows);
	if ((nx == 1) || (ny == 1)) {
		dim = {1, 1};
		std::vector<double> out{1};
		return out;
	} 
	dim = {ny, nx};

	SpatRaster x({ny, nx, 1}, {0., nx * xres, 0., ny * yres}, crs);
	std::vector<double> v(nx*ny, NAN);
	v[v.size()/2] = 1;
	SpatOptions opt;
	x.setValues(v, opt);
	x = x.distance(NAN, NAN, false, "m", false, "cosine", false, -1, opt);	

	std::vector<double> out;
	x.getValuesSource(0, out);				
	out[out.size()/2] = 1;

/*
	for (size_t i=0; i<ny; i++) {
		for (size_t j=0; j<nx; j++) {
			Rcpp::Rcout << out[i*nx+j] << " ";
		}
		Rcpp::Rcout << std::endl;
	}
*/
	return out;
}


std::vector<std::vector<double>> SpatRaster::extractBuffer(const std::vector<double> &x, const std::vector<double> &y, double b, SpatOptions &opt) {

	std::vector<std::vector<double>> out;
	
	if (!hasValues()) {
		setError("the raster has no values");
		return out;
	}
	if (nlyr() > 1) {
		setError("can only use a search_radius for one layer at a time");
		return out;
	}

	std::vector<size_t> dim;
	std::vector<double> cd;
	double ymean = 0;
	bool lonlat = is_lonlat();
	if (lonlat) {
		ymean = vmean(y, true);
	}
	
	std::vector<double> cb = circ_dist(xres(), yres(), b, nrow(), ncol(), dim, lonlat, ymean);
	bool docb = false;
	std::vector<bool> adj(cb.size(), false);
	if (cb.size() > 1) {
		cd.reserve(cb.size() * .67);
		for (size_t i=0; i<cb.size(); i++) {
			if (cb[i] <= b) {
				adj[i] = true;
				cd.push_back(cb[i]);
			}
		}
		docb = true;
	} else {
		
	}
	
    std::vector<double> cells = cellFromXY(x, y);
	std::vector<std::vector<double>> v = extractCell(cells, opt);

	std::vector<double> bestcell = cells;
	std::vector<double> bestdist(cells.size());

	std::vector<std::size_t> pm = sort_order_a(cd);
	permute(cd, pm);
	
	if (docb) {
		size_t n = x.size();
		for (size_t i=0; i<n; i++) {
			if (std::isnan(v[0][i]) && (!std::isnan(cells[i]))) {
				std::vector<double> acells = adjacentMat({cells[i]}, adj, dim, false);
				permute(acells, pm);
				std::vector<std::vector<double>> vv = extractCell(acells, opt);
	// take the first nearest. Instead could average over the cells with same distance
				for (size_t j=0; j<vv[0].size(); j++) {
					if (!std::isnan(vv[0][j])) {
						v[0][i] = vv[0][j];
						bestdist[i] = cd[j];
						bestcell[i] = acells[j];
						break;
					}
				}
			}
		}
	}
	
	out.push_back(v[0]);
	out.push_back(bestdist);
	out.push_back(bestcell);
	return out;
	
}


void SpatRaster::fourCellsFromXY(std::vector<double> &out, const std::vector<double> &x, const std::vector<double> &y) {

 	size_t n = x.size();
	SpatExtent e = getExtent();

	out.resize(0);
	out.reserve(4*n);

	double xmin = e.xmin;
	double xmax = e.xmax;
	double xr = xres();
	double ymin = e.ymin;
	double ymax = e.ymax;
	double yr = yres();
	double nc = ncol();
	int64_t mxr = nrow()-1;
	int64_t mxc = ncol()-1;
	int64_t r1, r2, c1, c2;
	std::vector<double> bad = {NAN, NAN, NAN, NAN};
	for (size_t i = 0; i < n; i++) {
		if (y[i] < ymin || y[i] > ymax || x[i] < xmin || x[i] > xmax) {
			out.insert(out.end(), bad.begin(), bad.end());
			continue;
		}
		if (y[i] == ymin) {
			r1 = mxr;
			r2 = mxr;
		} else {
			double p = (ymax - y[i]) / yr;
			r1 = trunc(p);
			if ((p - r1) > 0.5) {
				r2 = r1 == mxr ? mxr : r1 + 1;
			} else {
				r2 = r1;
				r1 = r1 == 0 ? 0 : r1 - 1;
			}
		}
		if (x[i] == xmax) {
			c1 = mxc;
			c2 = mxc;
		} else {
			double p = (x[i] - xmin) / xr;
			c1 = trunc(p);
			if ((p - c1) > 0.5) {
				c2 = c1 == mxc ? mxc : c1 + 1;
			} else {
				c2 = c1;
				c1 = c2 == 0 ? 0 : c2 - 1;
			}
		}
		out.push_back(r1 * nc + c1);
		out.push_back(r1 * nc + c2);
		out.push_back(r2 * nc + c1);
		out.push_back(r2 * nc + c2);
	}
}


std::vector<double> return_NAN(bool weights) {
	if (weights) {
		return std::vector<double>(4, NAN);
	}
	return std::vector<double>(1, NAN);
}

std::vector<double> bilinearInt(const double& x, const double& y,
                   const double& x1, const double& x2, const double& y1, const double& y2,
                   double& v11, double& v12, double& v21, double& v22, bool weights) {

	bool n1 = std::isnan(v11);
	bool n2 = std::isnan(v12);
	bool n3 = std::isnan(v21);
	bool n4 = std::isnan(v22);

	double dx = (x2 - x1);
	bool intx = dx > 0;
	double dy = (y1 - y2);
	bool inty = dy > 0;
	double w11, w12, w21, w22;

	if (std::isnan(x) || std::isnan(y) || (n1 && n2 && n3 && n4)) {
		return return_NAN(weights);
	}


	if (weights) {
		v11 = 1;
		v12 = 1;
		v21 = 1;
		v22 = 1;
	}


	if (intx && inty) {
		double d = dx * dy;
		if (!(n1 || n2 || n3 || n4)){
			w11 = v11 * ((x2 - x) * (y - y2)) / d;
			w12 = v12 * ((x - x1) * (y - y2)) / d;
			w21 = v21 * ((x2 - x) * (y1 - y)) / d;
			w22 = v22 * ((x - x1) * (y1 - y)) / d;
		} else if (!(n1 || n2 || n3)){
			w11 = v11 * ((x2 - x) * (y - y2)) / d;
			w12 = v12 * ((x - x1) * (y - y2)) / d;
			w21 = v21 * ((y1 - y)) / dy;
			w22 = 0;
		} else if (!(n1 || n2 || n4)){
			w11 = v11 * ((x2 - x) * (y - y2)) / d;
			w12 = v12 * ((x - x1) * (y - y2)) / d;
			w21 = 0;
			w22 = v22 * ((y1 - y)) / dy;
		} else if (!(n1 || n2 || n3)){
			w11 = v11 * ((x2 - x) * (y - y2)) / d;
			w12 = v12 * ((x - x1) * (y - y2)) / d;
			w21 = v21 * ((y1 - y)) / dy;
			w22 = 0;
		} else 	if (!(n1 || n3 || n4)){
			w11 = v11 * ((y - y2)) / dy;
			w12 = 0;
			w21 = v21 * ((x2 - x) * (y1 - y)) / d;
			w22 = v22 * ((x - x1) * (y1 - y)) / d;
		} else if (!(n2 || n3 || n4)){
			w11 = 0;
			w12 = v12 * ((y - y2)) / dy;
			w21 = v21 * ((x2 - x) * (y1 - y)) / d;
			w22 = v22 * ((x - x1) * (y1 - y)) / d;
		} else if (!(n1 || n2 )){
			w11 = v11 * ((x2 - x)) / dx;
			w12 = v12 * ((x - x1)) / dx;
			w21 = 0;
			w22 = 0;
		} else if (!(n1 || n3)){
			w11 = v11 * ((y - y2)) / dy;
			w12 = 0;
			w21 = v21 * ((y1 - y)) / dy;
			w22 = 0;
		} else if (!(n1 || n4)){
			w11 = v11 * ((y - y2)) / dy;
			w12 = 0;
			w21 = 0;
			w22 = v22 * ((y1 - y)) / dy;
		} else if (!(n2 || n3)){
			w11 = 0;
			w12 = v12 * ((y - y2)) / dy;
			w21 = v21 * ((y1 - y)) / dy;
			w22 = 0;
		} else if (!(n2 || n4)){
			w11 = 0;
			w12 = v12 * ((y - y2)) / dy;
			w21 = 0;
			w22 = v22 * ((y1 - y)) / dy;
		} else if (!(n3 || n4)){
			w11 = 0;
			w12 = 0;
			w21 = v21 * ((x2 - x)) / dx;
			w22 = v22 * ((x - x1)) / dx;
		} else if (!n1){
			w11 = v11;
			w12 = 0;
			w21 = 0;
			w22 = 0;
		} else if (!n2){
			w11 = 0;
			w12 = v12;
			w21 = 0;
			w22 = 0;
		} else if (!n3){
			w11 = 0;
			w12 = 0;
			w21 = v21;
			w22 = 0;
		} else if (!n4){
			w11 = 0;
			w12 = 0;
			w21 = 0;
			w22 = v22;
		} else {
			return return_NAN(weights);
		}
	} else if (intx) {
		w21 = 0.0;
		w22 = 0.0;
		if (!(n1 || n2)) {
			w11 = v11 * (x2 - x) / dx;
			w12 = v12 * (x - x1) / dx;
		} else if (!n1) {
			w11 = v11;
			w12 = 0.0;
		} else if (!n2){
			w11 = 0.0;
			w12 = v12;
		} else {
			return return_NAN(weights);
		}
	} else if (inty) {
		w12 = 0.0;
		w22 = 0.0;
		if (!(n1 || n3)) {
			w11 = v11 * (y - y2) / dy;
			w21 = v21 * (y1 - y) / dy;
		} else if (!n1) {
			w11 = v11;
			w21 = 0;
		} else if (!n3) {
			w11 = 0;
			w21 = v21;
		} else{
			return return_NAN(weights);
		}
	} else {
		w11 = v11;
		w21 = 0.0;
		w12 = 0.0;
		w22 = 0.0;
	}

	if (weights) {
		return std::vector<double>{ w11, w12, w21, w22 };
	}
	return std::vector<double>{ w11 + w12 + w21 + w22 };

}


void SpatRaster::bilinearValues(std::vector<std::vector<double>> &out, const std::vector<double> &x, const std::vector<double> &y, SpatOptions &opt) {


	std::vector<double> four;
	fourCellsFromXY(four, x, y);
	std::vector<std::vector<double>> xy = xyFromCell(four);
	std::vector<std::vector<double>> v = extractCell(four, opt);
	size_t n = x.size();
	out.resize(nlyr(), std::vector<double>(n));

	for (size_t i=0; i<n; i++) {
		size_t ii = i * 4;
		for (size_t j=0; j<nlyr(); j++) {
//			Rcpp::Rcout << xy[0][ii] << " " << xy[0][ii+1] << " " << xy[0][ii+2] << " " << xy[0][ii+3] << std::endl;
//			Rcpp::Rcout << xy[1][ii] << " " << xy[1][ii+1] << " " << xy[1][ii+2] << " " << xy[1][ii+3] << std::endl;
			
			std::vector<double> value = bilinearInt(x[i], y[i], xy[0][ii], xy[0][ii+1], xy[1][ii], xy[1][ii+3], 
												v[j][ii], v[j][ii+1], v[j][ii+2], v[j][ii+3], false);
			out[j][i] = value[0];
		}
	}
}



std::vector<double> SpatRaster::bilinearCells(const std::vector<double> &x, const std::vector<double> &y) {
	std::vector<double> four;
	fourCellsFromXY(four, x, y);
	std::vector<std::vector<double>> xy = xyFromCell(four);
//	std::vector<std::vector<double>> v = extractCell(four);
	size_t n = x.size();
	std::vector<double> res;
	res.reserve(n * 8);
	double v1=1, v2=1, v3=1, v4=1;
	
    for (size_t i=0; i<n; i++) {
        size_t ii = i * 4;
		//size_t j=0;
        std::vector<double> w = bilinearInt(x[i], y[i], xy[0][ii], xy[0][ii+1], xy[1][ii], xy[1][ii+3], v1, v2, v3, v4, true);
		res.insert(res.end(), four.begin()+ii, four.begin()+ii+4);
		res.insert(res.end(), w.begin(), w.end());
    }
	return res;
}


double bilinear(const std::vector<double> &v, const  std::vector<double> &e, const double &dxdy, const double &x, const double &y) {
	// values
	// v[0] v[1]
	// v[2] v[3]

	// coordinates
	//           e[3] (ymax)
	// (xmin)e[0]  e[1] (xmax)
	//           e[2] (ymin)

    double dx1 = x - e[0];
    double dx2 = e[1] - x;
    double dy1 = y - e[2];
    double dy2 = e[3] - y;
    return (v[2] * dx2 * dy2 + v[3] * dx1 * dy2 + v[0] * dx2 * dy1 + v[1] * dx1 * dy1) / dxdy;
}



std::vector<double> SpatRaster::line_cells(SpatGeom& g) {

	unsigned nrows = nrow();
	unsigned ncols = ncol();
	SpatExtent extent = getExtent();

	double xmin = extent.xmin;
	double ymax = extent.ymax;
	double rx = xres();
	double ry = yres();
	std::vector<double> out;

	unsigned np = g.size();
	for (size_t prt=0; prt<np; prt++) {
        SpatPart p = g.getPart(prt);
        double miny = vmin(p.y, true);
        double maxy = vmax(p.y, true);

        double minrow = rowFromY(miny);
        double maxrow = rowFromY(maxy);
        if (minrow > nrows || maxrow < 0) {
            return(out);
        }
        size_t startrow = minrow < 0 ? 0 : minrow;
        size_t endrow = maxrow >= nrows ? (nrows-1) : maxrow;
        unsigned n = p.x.size();
        out.reserve(2*(startrow-endrow+1));

        for (size_t row=startrow; row<endrow; row++) {
            double y = ymax - (row+0.5) * ry;
            unsigned rowcell = ncols * row;
            for (size_t i=1; i<n; i++) {
                size_t j = i-1;
                if (((p.y[i] < y) && (p.y[j] >= y)) || ((p.y[j] < y) && (p.y[i] >= y))) {
                    double col = ((p.x[i] - xmin + (y-p.y[i])/(p.y[j]-p.y[i]) * (p.x[j]-p.x[i])) + 0.5 * rx ) / rx;
                    if ((col >= 0) & (col < ncols)) {
                        out.push_back(rowcell + col);
                    }
                }
            }
        }
	}
	return(out);
}




std::vector<double> SpatRaster::polygon_cells(SpatGeom& g) {

// does not deal with holes yet.

	unsigned nrows = nrow();
	unsigned ncols = ncol();
	SpatExtent extent = getExtent();

	double xmin = extent.xmin;
	double ymax = extent.ymax;
	double rx = xres();
	double ry = yres();
	std::vector<double> out;

	unsigned np = g.size();
	for (size_t prt=0; prt<np; prt++) {

        SpatPart p = g.getPart(prt);
        double miny = vmin(p.y, true);
        double maxy = vmax(p.y, true);
        double minrow = rowFromY(miny);
        double maxrow = rowFromY(maxy);
        if (minrow > nrows || maxrow < 0) {
            return(out);
        }
        size_t startrow = minrow < 0 ? 0 : minrow;
        size_t endrow = maxrow >= nrows ? (nrows-1) : maxrow;
        unsigned n = p.x.size();
        out.reserve(5*(startrow-endrow+1));

        std::vector<unsigned> nCol(n);
        for (size_t row=0; row<nrows; row++) {
            double y = ymax - (row+0.5) * ry;
            // find nodes.
            unsigned nodes = 0;
            size_t j = n-1;
            for (size_t i=0; i<n; i++) {
                if (((p.y[i] < y) && (p.y[j] >= y)) || ((p.y[j] < y) && (p.y[i] >= y))) {
                //	nCol[nodes++]=(int)  (((pX[i] - xmin + (y-pY[i])/(pY[j]-pY[i]) * (pX[j]-pX[i])) + 0.5 * rx ) / rx);
                    double nds = ((p.x[i] - xmin + (y-p.y[i])/(p.y[j]-p.y[i]) * (p.x[j]-p.x[i])) + 0.5 * rx ) / rx;
                    nds = nds < 0 ? 0 : nds;
                    nds = nds > ncols ? ncols : nds;
                    nCol[nodes] = (unsigned) nds;
                    nodes++;
                }
                j = i;
            }

            // now remove the holes?

            std::sort(nCol.begin(), nCol.begin()+nodes);
            unsigned rowcell = ncols * row;

            // fill  cells between node pairs.
            for (size_t i=0; i < nodes; i+=2) {
                if (nCol[i+1] > 0 && nCol[i] < ncols) { // surely should be >= 0?
                    for (size_t col = nCol[i]; col < nCol[i+1]; col++) {
                        out.push_back(col + rowcell);
                    }
                }
            }
        }
	}
	return(out);
}


/*
idw
        bool lonlat = could_be_lonlat();
        //bool globalLonLat = is_global_lonlat();
        //size_t n = x.size();

        if (method == "idw") {
            std::function<std::vector<double>(std::vector<double>&,std::vector<double>&,double,double)> distFun;
 //           std::vector<double> distance_plane(std::vector<double> &x1, std::vector<double> &y1, std::vector<double> &x2, std::vector<double> &y2);
            if (lonlat) {
                distFun = distance_lonlat_vd;
            } else {
                distFun = distance_plane_vd;
            }
*/
/*
                cxy = xyFromCell(cells);
                d = distFun(cxy[0], cxy[1], x[i], y[i]);
                v = extractCell(cells);

                double a=0, b=0;
                for (size_t j=0; j<4; j++) {
                    a += v[j] * d[j];
                    b += d[j];
                }
                out[i] = a / b;
*/


// <layer<values>>
std::vector<std::vector<double>> SpatRaster::extractXY(const std::vector<double> &x, const std::vector<double> &y, const std::string &method, const bool &cells, SpatOptions &opt) {

    unsigned nl = nlyr();
    unsigned np = x.size();
	if (!hasValues()) {
		std::vector<std::vector<double>> out(nl+cells, std::vector<double>(np, NAN));
		return out;
	}
	std::vector<std::vector<double>> out;
	
    if (method == "bilinear") {
		bilinearValues(out, x, y, opt);
		if (cells) {
			std::vector<double> cell = cellFromXY(x, y);
			out.push_back(cell);
		}
	} else {
        std::vector<double> cell = cellFromXY(x, y);
        out = extractCell(cell, opt);
		if (cells) {
			out.push_back(cell);
		}
	}

    return out;
}


std::vector<double> SpatRaster::extractXYFlat(const std::vector<double> &x, const std::vector<double> &y, const std::string & method, const bool &cells, SpatOptions &opt) {


// <layer<values>>
	std::vector<std::vector<double>> e = extractXY(x, y, method, cells, opt);
	std::vector<double> out = e[0];
	for (size_t i=1; i<e.size(); i++) {
		out.insert(out.end(), e[i].begin(), e[i].end());
	}
	return out;
}

/*


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
*/

/*

std::vector<double> SpatRaster::extractXYFlat(const std::vector<double> &x, const std::vector<double> &y, const std::string & method, const bool &cells) {

    unsigned nl = nlyr();
    unsigned np = x.size();
	if (!hasValues()) {
		std::vector<double> out(nl * np, NAN);
		return out;
	}
	std::vector<double> out;
    if (method == "bilinear") {
		std::vector<std::vector<double>> bil = bilinearValues(x, y);
		if (cells) {
			std::vector<double> cell = cellFromXY(x, y);
			bil.push_back(cell);
		}
		out = flatten(bil);
	} else {
        std::vector<double> cell = cellFromXY(x, y);
		if (cells) {
			std::vector<std::vector<double>> xout;
			xout = extractCell(cell);
			xout.push_back(cell);
			out = flatten(xout);
		} else {
			out = extractCellFlat(cell);
		}
	}
	return out;
}
*/


// <geom<layer<values>>>
std::vector<std::vector<std::vector<double>>> SpatRaster::extractVector(SpatVector v, bool touches, bool small, std::string method, bool cells, bool xy, bool weights, bool exact, SpatOptions &opt) {

	if (!source[0].srs.is_same(v.srs, true)) {
		v = v.project(getSRS("wkt"), false);
		addWarning("transforming vector data to the CRS of the raster");
	}

	std::string gtype = v.type();
	if (gtype == "points") weights = false;

	if (exact) weights = false;
    unsigned nl = nlyr();
    unsigned ng = v.size();
    std::vector<std::vector<std::vector<double>>> out(ng, std::vector<std::vector<double>>(nl + cells + 2*xy + (weights || exact)));

	if (!hasValues()) {
		setError("raster has no values");
		return out;
	}
/*
	#if GDAL_VERSION_MAJOR < 3
	if (weights) {
		setError("extract with weights not supported for your GDAL version");
		return out;
	}
	#endif
*/
	std::vector<std::vector<double>> srcout;
	if (gtype == "points") {
		if (method != "bilinear") method = "simple";
		SpatDataFrame vd = v.getGeometryDF();
		if (vd.nrow() == ng) {  // single point geometry
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			srcout = extractXY(x, y, method, cells, opt);
			for (size_t i=0; i<ng; i++) {
				for (size_t j=0; j<nl; j++) {
					out[i][j].push_back( srcout[j][i] );
				}
				if (cells) {
					out[i][nl].push_back( srcout[nl][i] );
				}
				if (xy) {
					out[i][nl+cells].push_back(x[i]);
					out[i][nl+cells+1].push_back(y[i]);
				}
			}
		} else { //multipoint
//			Rcpp::Rcout << "multipoint" << std::endl;
			for (size_t i=0; i<ng; i++) {
				SpatVector vv = v.subset_rows(i);
				SpatDataFrame vd = vv.getGeometryDF();
				std::vector<double> x = vd.getD(0);
				std::vector<double> y = vd.getD(1);
				//srcout = extractXY(x, y, method, cells);

				/*
				for (size_t j=0; j<nl; j++) {
					out[i][j] = srcout[j];
				}
				if (cells) {
					out[i][nl] = srcout[nl];
				}
				if (xy) {
					out[i][nl+cells]   = x;
					out[i][nl+cells+1] = y;
				}
				*/
			}
		}
	} else {
	    SpatRaster r = geometry(1);
		//std::vector<double> feats(1, 1) ;
        for (size_t i=0; i<ng; i++) {
            SpatGeom g = v.getGeom(i);
            SpatVector p(g);
			p.srs = v.srs;
			std::vector<double> cell, wgt;
			if (weights) {
				if (gtype == "lines") {
					rasterizeLinesLength(cell, wgt, p, opt);
				} else {
					rasterizeCellsWeights(cell, wgt, p, opt);
				}
			} else if (exact) {
				if (gtype == "lines") {
					rasterizeLinesLength(cell, wgt, p, opt);
				} else {
					rasterizeCellsExact(cell, wgt, p, opt);
				}
			} else {
				cell = rasterizeCells(p, touches, small, opt);
            }
			srcout = extractCell(cell, opt);
            for (size_t j=0; j<nl; j++) {
                out[i][j] = srcout[j];
            }
			if (cells) {
				out[i][nl] = cell;
			}
			if (xy) {
				std::vector<std::vector<double>> crds = xyFromCell(cell);
				out[i][nl+cells]   = crds[0];
				out[i][nl+cells+1] = crds[1];
			}
			if (weights || exact) {
				out[i][nl + cells + 2*xy] = wgt;
			}
        }
	}
	return out;
}

/*
std::vector<double> SpatRaster::extractVectorFlat(SpatVector v, std::string fun, bool narm, bool touches, std::string method, bool cells, bool xy, bool weights, bool exact, SpatOptions &opt) {

	std::vector<double> flat;
	std::string gtype = v.type();
	if (gtype == "points") {
		weights = false;
		exact = false;
	}
	if (exact) weights = false;

    unsigned nl = nlyr();
    unsigned ng = v.size();

	if (!hasValues()) {
		setError("raster has no values");
		return flat;
	}

    std::vector<std::vector<std::vector<double>>> out;
	if (gtype != "points") {
		out.resize(ng, std::vector<std::vector<double>>(nl + cells + 2*xy + (weights||exact)));
	}
	std::vector<std::vector<double>> srcout;
	if (gtype == "points") {
		if (method != "bilinear") method = "simple";
			SpatDataFrame vd = v.getGeometryDF();
			//if (vd.nrow() == ng) {  // single point geometry
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			std::vector<std::vector<double>> xycells;
			if (xy) {
				std::vector<double> cellxy = cellFromXY(x, y);
				xycells = xyFromCell(cellxy);
			}
			if (!cells & !xy) {
				return( extractXYFlat(x, y, method, cells));
			} else {
				srcout = extractXY(x, y, method, cells);
				nl += cells;
				flat.reserve(ng * nl);
				for (size_t i=0; i<ng; i++) {
					//flat.push_back( i+1 );//no id for points
					for (size_t j=0; j<nl; j++) {
						flat.push_back( srcout[j][i] );
					}
					if (xy) {
						flat.push_back(xycells[0][i]);
						flat.push_back(xycells[1][i]);
					}
				}
			}
			return flat;
		*/
		/*
		} else { // multipoint
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			if (!cells & !xy & !weights) {
				return( extractXYFlat(x, y, method, cells));
			}
			for (size_t i=0; i<ng; i++) {
				SpatVector vv = v.subset_rows(i);
				SpatDataFrame vd = vv.getGeometryDF();
				std::vector<double> x = vd.getD(0);
				std::vector<double> y = vd.getD(1);
				srcout = extractXY(x, y, method, cells);

				out.push_back(srcout);

				if (cells) {
					out[i][nl] = srcout[nl];
				}
				if (xy) {
					out[i][nl+cells]   = x;
					out[i][nl+cells+1] = y;
				}
			}
		}
		*/
/*	} else {
	    SpatRaster r = geometry(1);
		//std::vector<double> feats(1, 1) ;
        for (size_t i=0; i<ng; i++) {
            SpatGeom g = v.getGeom(i);
            SpatVector p(g);
			p.srs = v.srs;
			std::vector<double> cell, wgt;
			if (weights) {
				if (gtype == "lines") {
					rasterizeLinesLength(cell, wgt, p, opt);
				} else {
					rasterizeCellsWeights(cell, wgt, p, opt);
				}
			} else if (exact) {
				if (gtype == "lines") {
					rasterizeLinesLength(cell, wgt, p, opt);
				} else {
					rasterizeCellsExact(cell, wgt, p, opt);
				}
			} else {
				cell = rasterizeCells(p, touches, opt);
            }
			srcout = extractCell(cell);
            for (size_t j=0; j<nl; j++) {
                out[i][j] = srcout[j];
            }
			if (cells) {
				out[i][nl] = cell;
			}
			if (xy) {
				std::vector<std::vector<double>> crds = xyFromCell(cell);
				out[i][nl+cells]   = crds[0];
				out[i][nl+cells+1] = crds[1];
			}
			if (weights || exact) {
				out[i][nl + cells + 2*xy] = wgt;
			}
        }
	}

	size_t fsize = 0;
	for (size_t i=0; i<out.size(); i++) { // geoms
		fsize += (out[i].size()+1) * nl;
	}
	flat.reserve(fsize);

	for (size_t i=0; i<out.size(); i++) { // geoms
		for (size_t j=0; j<out[i][0].size(); j++) { // cells
			flat.push_back(i+1);
			for (size_t k=0; k<out[i].size(); k++) { // layers
				flat.push_back(out[i][k][j]);
			}
		}
	}
	return flat;
}
*/


std::vector<double> SpatRaster::extractVectorFlat(SpatVector v, std::vector<std::string> funs, bool narm, bool touches, bool small, std::string method, bool cells, bool xy, bool weights, bool exact, SpatOptions &opt) {

	if (!source[0].srs.is_same(v.srs, true)) {
		v = v.project(getSRS("wkt"), false);
		addWarning("transforming vector data to the CRS of the raster");
//		addWarning("CRS of raster and vector data do not match");
	}

	std::vector<double> flat;
	std::string gtype = v.type();
	if (gtype == "points") {
		weights = false;
		exact = false;
	}
	if (exact) weights = false;

    unsigned nl = nlyr();
    unsigned ng = v.size();

	if (!hasValues()) {
		setError("raster has no values");
		return flat;
	}

 	if (gtype == "points") {
		if (method != "bilinear") method = "simple";
		SpatDataFrame vd = v.getGeometryDF();
			//if (vd.nrow() == ng) {  // single point geometry
		std::vector<double> x = vd.getD(0);
		std::vector<double> y = vd.getD(1);
		std::vector<std::vector<double>> xycells;
		if (xy) {
			std::vector<double> cellxy = cellFromXY(x, y);
			xycells = xyFromCell(cellxy);
		}
		if (!cells & !xy) {
			return( extractXYFlat(x, y, method, cells, opt));
		} else {
			std::vector<std::vector<double>> srcout = extractXY(x, y, method, cells, opt);
			nl += cells;
			flat.reserve(ng * nl);
			for (size_t i=0; i<ng; i++) {
				//flat.push_back( i+1 );//no id for points
				for (size_t j=0; j<nl; j++) {
					flat.push_back( srcout[j][i] );
				}
				if (xy) {
					flat.push_back(xycells[0][i]);
					flat.push_back(xycells[1][i]);
				}
			}
		}
		return flat;
		/*
		} else { // multipoint
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			if (!cells & !xy & !weights) {
				return( extractXYFlat(x, y, method, cells));
			}
			for (size_t i=0; i<ng; i++) {
				SpatVector vv = v.subset_rows(i);
				SpatDataFrame vd = vv.getGeometryDF();
				std::vector<double> x = vd.getD(0);
				std::vector<double> y = vd.getD(1);
				srcout = extractXY(x, y, method, cells);

				out.push_back(srcout);

				if (cells) {
					out[i][nl] = srcout[nl];
				}
				if (xy) {
					out[i][nl+cells]   = x;
					out[i][nl+cells+1] = y;
				}
			}
		}
		*/
	} 

	std::vector<std::vector<std::vector<double>>> out;

	SpatRaster r = geometry(1);
	//std::vector<double> feats(1, 1) ;
	std::vector<std::function<double(std::vector<double>&, size_t, size_t)>> efuns;
	std::vector<std::function<double(std::vector<double>&, std::vector<double>&, size_t, size_t)>> wfuns;
	bool havefun = false;
	if (!funs[0].empty()) {
		if (weights | exact) {
			wfuns.resize(funs.size());
			for (size_t i=0; i<funs.size(); i++) {
				if (!getseWfun(wfuns[i], funs[i], narm)) {
					setError(funs[i] + " is not a valid function");
					return flat;
				}
			}
		} else {
			efuns.resize(funs.size());
			for (size_t i=0; i<funs.size(); i++) {
				if (!getseFun(efuns[i], funs[i], narm)) {
					setError(funs[i] + " is not a valid function");
					return flat;
				}
			}
		}
		havefun = true;
		flat.reserve(nl*ng);
	} else {
		out.resize(ng);
	}	
	for (size_t i=0; i<ng; i++) {
		SpatGeom g = v.getGeom(i);
		SpatVector p(g);
		p.srs = v.srs;
		std::vector<double> cell, wgt;
		if (weights) {
			if (gtype == "lines") {
				rasterizeLinesLength(cell, wgt, p, opt);
			} else {
				rasterizeCellsWeights(cell, wgt, p, opt);
			}
		} else if (exact) {
			if (gtype == "lines") {
				rasterizeLinesLength(cell, wgt, p, opt);
			} else {
				rasterizeCellsExact(cell, wgt, p, opt);
			}
		} else {
			cell = rasterizeCells(p, touches, small, opt);
		}
		
		if (havefun) {
			std::vector<std::vector<double>> cvals = extractCell(cell, opt);
			if (weights | exact) {
				for (size_t j=0; j<funs.size(); j++) {
					for (size_t k=0; k<nl; k++) {
						flat.push_back( wfuns[j](cvals[k], wgt, 0, cvals[k].size()) );
					}
				}
			} else {
				for (size_t j=0; j<funs.size(); j++) {
					for (size_t k=0; k<nl; k++) {
						flat.push_back( efuns[j](cvals[k], 0, cvals[k].size()) );
					}						
				}
			}			
		} else {
			out[i] = extractCell(cell, opt);
			if (cells) {
				out[i].push_back(cell);
			}
			if (xy) {
				std::vector<std::vector<double>> crds = xyFromCell(cell);
				out[i].push_back(crds[0]);
				out[i].push_back(crds[1]);
			}
			if (weights || exact) {
				out[i].push_back(wgt);
			}
		}
	}

	if (havefun) return flat;
	
	size_t fsize = 0;
	for (size_t i=0; i<out.size(); i++) { // geoms
		fsize += (out[i].size()+1) * nl;
	}
	flat.reserve(fsize);

	for (size_t i=0; i<out.size(); i++) { // geoms
		for (size_t j=0; j<out[i][0].size(); j++) { // cells
			flat.push_back(i+1);
			for (size_t k=0; k<out[i].size(); k++) { // layers
				flat.push_back(out[i][k][j]);
			}
		}
	}
	return flat;
}




std::vector<std::vector<double>> SpatRaster::extractCell(std::vector<double> &cell, SpatOptions opt) {


	std::vector<double> wcell;
	std::vector<std::vector<int64_t>> rc, wrc;
	rc = rowColFromCell(cell);

	size_t n = cell.size();
	if (!hasValues()) {
		std::vector<std::vector<double>> out(nlyr(), std::vector<double>(n, NAN));
		return out;
	}
	
	unsigned ns = nsrc();
	unsigned lyr = 0;
	size_t nc;

	std::vector<std::vector<double>> out(nlyr());
	for (size_t src=0; src<ns; src++) {
		unsigned slyrs = source[src].layers.size();
		bool win = source[src].hasWindow;
		if (win) {
			nc = source[src].window.full_ncol * source[src].window.full_nrow;
			wrc = rc;
			wcell.reserve(cell.size());
			for (size_t i=0; i<cell.size(); i++) {
				if ((wrc[0][i] < 0) || (wrc[1][i] <0)) {
					wcell.push_back( NAN );
				} else {
					wrc[0][i] += source[src].window.off_row;
					wrc[1][i] += source[src].window.off_col;
					wcell.push_back( wrc[0][i] * source[src].window.full_ncol + wrc[1][i] );
				}
			}
		} else {
			nc = ncell();
		}
		if (source[src].memory) {
			for (size_t i=0; i<slyrs; i++) {
				out[lyr] = std::vector<double>(n, NAN);
				size_t j = i * nc;
				if (win) {
					for (size_t k=0; k<n; k++) {
						if (!is_NA(wcell[k]) && wcell[k] >= 0 && wcell[k] < nc) {
							out[lyr][k] = source[src].values[j + wcell[k]];
						}
					}
				} else {
					for (size_t k=0; k<n; k++) {
						if (!is_NA(cell[k]) && cell[k] >= 0 && cell[k] < nc) {
							out[lyr][k] = source[src].values[j + cell[k]];
						}
					}
				}
				lyr++;
			}
		} else {
			#ifdef useGDAL
			size_t pos = source[src].filename.find("https://");
			if ((pos != std::string::npos) && (rc[0].size() > 200)) {
				if (win) {
					readRowColBlock(src, out, lyr, wrc[0], wrc[1], opt);
				} else {
					readRowColBlock(src, out, lyr, rc[0], rc[1], opt);
				}
			} else {
				if (win) {
					readRowColGDAL(src, out, lyr, wrc[0], wrc[1]);
				} else {
					readRowColGDAL(src, out, lyr, rc[0], rc[1]);
				}
				if (hasError()) return out;
			}
			lyr += slyrs;
			#else 
			out.resize(slyrs);
			for (size_t i=0; i<slyrs; i++) {
				out[lyr] = std::vector<double>(n, NAN);
			}
			#endif
			if (hasError()) return out;
		}
	}
	return out;
}





//std::vector<std::vector<double>> SpatRaster::extractRowCol(std::vector<int64_t> &row, std::vector<int64_t> &col) {


/*

std::vector<double> SpatRaster::extractCellFlat(std::vector<double> &cell) {

	std::vector<double> wcell;
	std::vector<std::vector<int64_t>> rc, wrc;
	rc = rowColFromCell(cell);

	size_t n  = cell.size();
	std::vector<double> out(nlyr() * n, NAN);

	unsigned ns = nsrc();
//	unsigned lyr = 0;
	size_t nc;
	size_t off = 0;
	for (size_t src=0; src<ns; src++) {
		unsigned slyrs = source[src].layers.size();
		bool win = source[src].hasWindow;
		if (win) {
			nc = source[src].window.full_ncol * source[src].window.full_nrow;
			wrc = rc;
			wcell.reserve(cell.size());
			for (size_t i=0; i<cell.size(); i++) {
				if ((wrc[0][i] < 0) || (wrc[1][i] <0)) {
					wcell.push_back( NAN );
				} else {
					wrc[0][i] += source[src].window.off_row;
					wrc[1][i] += source[src].window.off_col;
					wcell.push_back( wrc[0][i] * source[src].window.full_ncol + wrc[1][i] );
				}
			}
		} else {
			nc = ncell();
		}

		if (source[src].memory) {
			for (size_t i=0; i<slyrs; i++) {
				size_t j = i * nc;
				size_t off2 = off + i*n;
				if (win) {
					for (size_t k=0; k<n; k++) {
						if (!is_NA(wcell[k]) && wcell[k] >= 0 && wcell[k] < nc) {
							out[off2+k] = source[src].values[j + wcell[k]] ;
						}
					}
				} else {
					for (size_t k=0; k<n; k++) {
						if (!is_NA(cell[k]) && cell[k] >= 0 && cell[k] < nc) {
							out[off2+k] = source[src].values[j + cell[k]];
						}
					}
				}
				//lyr++;
			}
		} else {
			//if (source[0].driver == "raster") {
			//	srcout = readCellsBinary(src, cell);
			//} else {
			#ifdef useGDAL
			std::vector<double> g;
			size_t pos = source[0].filename.find("https://");
			if ((pos != std::string::npos) && (rc[0].size() > 200)) {
				if (win) {
					g = readRowColBlockFlat(src, wrc[0], wrc[1]);
				} else {
					g = readRowColBlockFlat(src, rc[0], rc[1]);
				}
			} else {	
				if (win) {
					g = readRowColGDALFlat(src, wrc[0], wrc[1]);
				} else {
					g = readRowColGDALFlat(src, rc[0], rc[1]);
				}
			}
			for (size_t i=0; i<slyrs; i++) {
				size_t j = i * n;
				size_t off2 = off + j;
				for (size_t k=0; k<n; k++) {
					out[off2+k] = g[j + k];
				}
			}
			#endif
			if (hasError()) {
				return out;
			}
		}
		off += source[src].nlyr;
	}
	return out;
}
*/

std::vector<double> SpatRaster::vectCells(SpatVector v, bool touches, bool small, std::string method, bool weights, bool exact, SpatOptions &opt) {

	std::string gtype = v.type();
	if (gtype != "polygons") weights = false;
	std::vector<double> out, cells, wghts;
	if (gtype == "points") {
		SpatDataFrame vd = v.getGeometryDF();
		//std::vector<long> id = vd.getI(0);
		if (method == "bilinear") {
			return bilinearCells(vd.getD(0), vd.getD(1));
		} else {
			return cellFromXY(vd.getD(0), vd.getD(1));
			//cells = cellFromXY(vd.getD(0), vd.getD(1));
			//out.insert(out.end(), id.begin(), id.end());
			//out.insert(out.end(), cells.begin(), cells.end());
		}
	} else {
		unsigned ng = v.size();
		SpatRaster r = geometry(1);
		std::vector<double> feats(1, 1) ;
        for (size_t i=0; i<ng; i++) {
            SpatGeom g = v.getGeom(i);
            SpatVector p(g);
			p.srs = v.srs;
			if (weights) {
				std::vector<double> cnr, wght;
				rasterizeCellsWeights(cnr, wght, p, opt);
				std::vector<double> id(cnr.size(), i);
				out.insert(out.end(), id.begin(), id.end());
				cells.insert(cells.end(), cnr.begin(), cnr.end());
				wghts.insert(wghts.end(), wght.begin(), wght.end());
			} else if (exact) {
				std::vector<double> cnr, wght;
				rasterizeCellsExact(cnr, wght, p, opt);
				std::vector<double> id(cnr.size(), i);
				out.insert(out.end(), id.begin(), id.end());
				cells.insert(cells.end(), cnr.begin(), cnr.end());
				wghts.insert(wghts.end(), wght.begin(), wght.end());
			} else {
				std::vector<double> geomc = rasterizeCells(p, touches, small, opt);
				std::vector<double> id(geomc.size(), i);
				out.insert(out.end(), id.begin(), id.end());
				cells.insert(cells.end(), geomc.begin(), geomc.end());
			}
        }
		if (weights || exact) {
			out.insert(out.end(), cells.begin(), cells.end());
			out.insert(out.end(), wghts.begin(), wghts.end());
		} else {
			out.insert(out.end(), cells.begin(), cells.end());
		}
	}
	return out;
}


std::vector<double> SpatRaster::extCells(SpatExtent ext) {

	std::vector<double> out;
	ext = align(ext, "near");
	ext = ext.intersect(getExtent());
	if (!ext.valid()) {
		return(out);
	}
	double resx = xres() / 2;
	double resy = yres() / 2;
	std::vector<double> e = ext.asVector();
	e[0] += resx;
	e[1] -= resx;
	e[2] += resy;
	e[3] -= resy;
	std::vector<double> ex = {e[0], e[1]};
	std::vector<double> ey = {e[3], e[2]};
	std::vector<int64_t> r = rowFromY(ey);
	std::vector<int64_t> c = colFromX(ex);
	int64_t nc = ncol();
	out.reserve((r[1]-r[0]) * (c[1]-c[0]));
	for (int64_t i=r[0]; i <= r[1]; i++) {
		for (int64_t j=c[0]; j <= c[1]; j++) {
			out.push_back(i*nc+j);
		}
	}
	return out;
}



std::vector<std::vector<std::vector<double>>> SpatRasterStack::extractXY(std::vector<double> &x, std::vector<double> &y, std::string method, SpatOptions &opt) {
	unsigned ns = nsds();
	std::vector<std::vector<std::vector<double>>> out(ns);
	bool cells = false;
	for (size_t i=0; i<ns; i++) {
		SpatRaster r = getsds(i);
		out[i] = r.extractXY(x, y, method, cells, opt);
	}
	return out;
}

std::vector<std::vector<std::vector<double>>> SpatRasterStack::extractCell(std::vector<double> &cell, SpatOptions &opt) {
	unsigned ns = nsds();
	std::vector<std::vector<std::vector<double>>> out(ns);
	for (size_t i=0; i<ns; i++) {
		SpatRaster r = getsds(i);
		out[i] = r.extractCell(cell, opt);
	}
	return out;
}


// this is rather inefficient (repeated rasterization)
std::vector<std::vector<std::vector<std::vector<double>>>> SpatRasterStack::extractVector(SpatVector v, bool touches, bool small, std::string method, SpatOptions &opt) {
	unsigned ns = nsds();
	std::vector<std::vector<std::vector<std::vector<double>>>> out(ns);
	for (size_t i=0; i<ns; i++) {
		SpatRaster r = getsds(i);
		out[i] = r.extractVector(v, touches, small, method, false, false, false, false, opt);
	}
	return out;
}



/*
        } else if (method == "oldbilinear") {

// this is much too slow
			SpatExtent extent = getExtent();

			double ymax = extent.ymax;
            double xmin = extent.xmin;
            double yrs = yres();
            double xrs = xres();

            //SpatOptions opt;
            SpatRaster g = geometry();
			std::vector<unsigned> f = {2,2};
            SpatRaster gd = g.disaggregate(f, opt);

            double dyrs = gd.yres();
            double dxrs =  gd.xres();
            std::vector<double> d, cells(4);

            std::vector<std::vector<double> > cxy;

            std::vector<double> rc(4);
			unsigned nr = nrow();
			unsigned nc = ncol();
            unsigned mnr = nr-1;
            unsigned mnc = nc-1;

            // needs row-wise adjustment for lonlat
            double dxdy = xres() * yres();

            for (size_t i=0; i<n; i++) {
                long row_d = ((ymax - y[i]) / dyrs);
                long col_d = ((x[i] - xmin) / dxrs);
                unsigned rq = row_d % 2;
                unsigned cq = col_d % 2;

                double row1 = std::floor((ymax - y[i]) / yrs);
                double col1 = std::floor((x[i] - xmin) / xrs);
                if ((row1 < 0) || (row1 > mnr)) { continue; }

                double row2 = (rq == 0) ? row1-1 : row1+1;
                row2 = row2 < 0 ? row1+1 : row2==nr ? row1-1 : row2;
                double col2;
                if (globalLonLat) {
                    if ((col1 < -1) | (col1 > nc)) { continue; }
                    col1 = col1 < 0 ? mnc : col1 > mnc ? 0 : col1;
                    col2 = (cq == 0) ? col1-1 : col1 + 1;
                    col2 = col2 < 0 ? mnc : col2 > mnc ? 0 : col2;
                } else {
                    if ((col1 < 0) | (col1 > mnc)) { continue; }
                    col2 = (cq == 0) ? col1-1 : col1 + 1;
                    col2 = col2 < 0 ? col1+1 : col2 == nc ? col1-1 : col2;
                }
                cells[0] = nc * row1 + col1;
                cells[1] = nc * row1 + col2;
                cells[2] = nc * row2 + col1;
                cells[3] = nc * row2 + col2;
                std::sort(cells.begin(), cells.end());
                std::vector<std::vector<double>> xy = xyFromCell(cells);
                std::vector<std::vector<double>> v = extractCell(cells);
                std::vector<double> e = {xy[0][0], xy[0][1], xy[1][2], xy[1][0]};
                for (size_t j=0; j<nl; j++) {
                    out[j][i] = bilinear(v[j], e, dxdy, x[i], y[i]);
                }
            }
        }
		*/


/*
std::vector<double> SpatRaster::extractCell(std::vector<double> &cell) {

	unsigned n  = cell.size();
	unsigned nc = ncell();

	std::vector<double> out;
	if (!hasValues()) {
		out = std::vector<double>(n * nlyr(), NAN)
		return out;
	}

	unsigned ns = nsrc();
	for (size_t src=0; src<ns; src++) {
		unsigned slyrs = source[src].layers.size();
		std::vector<double> srcout;
		if (source[src].memory) {
			srcout = std::vector<double>(n * slyrs, NAN)
			std::vector<size_t> off1(slyrs);
			std::vector<size_t> off2(slyrs);
			for (size_t i=0; i<slyrs; i++) {
				off1[i] = i * n;
				off2[i] = i * nc;
			}
			for (size_t i=0; i<n; i++) {
				if (!is_NA(cell[i]) && cell[i] >= 0 && cell[i] < nc) {
					for (size_t j=0; j<slyrs; j++) {
						out[off1[j]+i] = source[src].values[off2[j] + cell[i]];
					}
				}
			}
		} else {
			#ifdef useGDAL
			std::vector<std::vector<int64_t>> rc = rowColFromCell(cell);
			srcout = readRowColGDAL(src, rc[0], rc[1]);
			#endif
			if (hasError()) return out;
			//}
		}
		out.insert(out.end(), srcout.begin(), srcout.end());
	}
	return out;
}
*/



/*
double distInt(double d, double pd1, double pd2, double v1, double v2) {
  double result = (v2 * pd1 + v1 * pd2) / d;
  return result;
}

inline double rowColToCell(unsigned ncols, unsigned row, unsigned col) {
  return row * ncols + col;
}

double linearInt(const double& d, const double& x, const double& x1, const double& x2, const double& v1, const double& v2) {
	double result = (v2 * (x - x1) + v1 * (x2 - x)) / d;
	return result;
}

// ok but cannot handle NA
double bilinearInt(const double& x, const double& y, const double& x1, const double& x2, const double& y1, const double& y2, const double& v11, const double& v21, const double& v12, const double& v22) {
  double d = x2-x1;
  double h1 =  linearInt(d, x, x1, x2, v11, v21);
  double h2 =  linearInt(d, x, x1, x2, v12, v22);
  d = y2-y1;
  double v =  linearInt(d, y, y1, y2, h1, h2);
  return v;
}

double bilinearIntold(const double& x, const double& y,
                   const double& x1, const double& x2, const double& y1, const double& y2,
                   const double& v11, const double& v21, const double& v12, const double& v22) {
	double d = x2-x1;
	double h1=NAN;
	double h2=NAN;
	if (!std::isnan(v11) && !std::isnan(v21)) {
		h1 =  linearInt(d, x, x1, x2, v11, v21);
	} else if (!std::isnan(v11)) {
		h1 = v11;
	} else if (!std::isnan(v21)) {
		h1 = v21;
	}

	if (!std::isnan(v12) && !std::isnan(v22)) {
		h2 =  linearInt(d, x, x1, x2, v12, v22);
	} else if (!std::isnan(v12)) {
		h2 = v12;
	} else if (!std::isnan(v22)) {
		h2 = v22;
	}
	if (!std::isnan(h1) && !std::isnan(h2)) {
		d = y2-y1;
		double v = linearInt(d, y, y1, y2, h1, h2);
		return v;
	} else if (!std::isnan(h1)) {
		return h1;
	} else if (!std::isnan(h2)) {
		return h2;
	}
	return NAN;
}

double bilinear_geo(double x, double y, double x1, double x2, double y1, double y2, double halfyres, std::vector<double> vv) {

    double a = 6378137.0;
    double f = 1/298.257223563;
	double hy = y1 - halfyres;
	double d = distance_lonlat(x1, hy, x2, hy, a, f);

    std::vector<double> dist(4);
	double pd1 = distance_lonlat(x, hy, x1, hy, a, f);
	double pd2 = distance_lonlat(x, hy, x2, hy, a, f);
	double h1 = distInt(d, pd1, pd2, vv[0], vv[1]);
	double h2 = distInt(d, pd1, pd2, vv[2], vv[3]);
	d = y2 - y1;
	double v =  linearInt(d, y, y1, y2, h1, h2);

	return v;
}
*/
