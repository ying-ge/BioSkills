
#include "ogr_spatialref.h"

#include "spatRaster.h"
#include "file_utils.h"

#include "gdal_alg.h"
#include "ogrsf_frmts.h"

//#include "spatFactor.h"
#include "recycle.h"
#include "sort.h"
#include "gdalio.h"


SpatRaster SpatRaster::rasterizePoints(std::vector<double>&x, std::vector<double> &y, std::string fun, std::vector<double> &values, bool narm, double background, SpatOptions &opt) {

	SpatRaster out = geometry(1, false, false, false);
	
	if (!out.writeStart(opt, filenames())) {
		return out;
	}
	
	if (y.size() != x.size()) {
		out.setError("number of x and y coordinates do not match");
		return out;
	}
	if ((fun == "count") && (values.size() != x.size()) && (!values.empty())) {
		out.setError("number of values does not match the number of points");
		return out;
	} else if (values.size() != x.size()) {
		out.setError("number of values does not match the number of points");
		return out;
	}

	size_t nc = ncol();
	std::vector<double> cells = cellFromXY(x, y, -9);
	// order for multiple chunks, but also to remove NAs (-9)
	std::vector<std::size_t> so = sort_order_a(cells);
	permute(cells, so);
	permute(values, so);

	size_t cellcnt = 0;
	for (size_t i=0; i<cells.size(); i++) {
		if (cells[i] < 0) {
			cellcnt++;
		} else {
			break;
		}
	}

	if (fun == "count") {
		bool dotest = (!values.empty()) && narm;
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), 0);

			for (size_t j=cellcnt; j<cells.size(); j++) {
				if (dotest && std::isnan(values[j])) continue;
				if (cells[j] <= cmax) {
					size_t k = cells[j] - cmin;
					v[k]++;
				} else {
					cellcnt = j;
					break;
				}
			}
			if (background != 0) {
				for (size_t j=0; j<v.size(); j++) {
					if (v[j] == 0) {
						v[j] = background;
					}
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	} else if (fun == "mean") {
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), 0);
			std::vector<double> cnt = v;
			for (size_t j=cellcnt; j<cells.size(); j++) {
				if (narm && std::isnan(values[j])) continue;
				if (cells[j] <= cmax) {
					size_t k = cells[j] - cmin;
					v[k] += values[j];
					cnt[k]++;
				} else {
					cellcnt = j;
					break;
				}
			}
			if (background != 0) {
				for (size_t j=0; j<v.size(); j++) {
					if (cnt[j] == 0) {
						v[j] = background;
					} else {
						v[j] /= cnt[j];						
					}
				}
			} else {
				for (size_t j=0; j<cnt.size(); j++) {
					if (cnt[j] > 0) {
						v[j] /= cnt[j];
					}
				}				
			}			
			
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	} else if (fun == "sum") {
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), background);
			std::vector<bool> newcell(out.bs.nrows[i] * out.ncol(), true);
			for (size_t j=cellcnt; j<cells.size(); j++) {
				if (narm && std::isnan(values[j])) continue;
				if (cells[j] <= cmax) {
					size_t k = cells[j] - cmin;
					if (newcell[k]) {
						v[k] = values[j];
						newcell[k] = false;
					} else {
						v[k] += values[j];
					}
				} else {
					cellcnt = j;
					break;
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	} else if (fun == "min") {
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), background);
			std::vector<bool> newcell(out.bs.nrows[i] * out.ncol(), true);
			for (size_t j=cellcnt; j<cells.size(); j++) {
				if (narm && std::isnan(values[j])) continue;
				if (cells[j] <= cmax) {
					size_t k = cells[j] - cmin;
					if (newcell[k]) {
						v[k] = values[j];
						newcell[k] = false;
					} else {
						v[k] = std::min(v[k], values[j]);
					}
				} else {
					cellcnt = j;
					break;
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	} else if (fun == "max") {
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), background);
			std::vector<bool> newcell(out.bs.nrows[i] * out.ncol(), true);
			for (size_t j=cellcnt; j<cells.size(); j++) {
				if (narm && std::isnan(values[j])) continue;
				if (cells[j] <= cmax) {
					size_t k = cells[j] - cmin;
					if (newcell[k]) {
						v[k] = values[j];
						newcell[k] = false;
					} else {
						v[k] = std::max(v[k], values[j]);
					}
				} else {
					cellcnt = j;
					break;
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	} else if (fun == "prod") {
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), background);
			std::vector<bool> newcell(out.bs.nrows[i] * out.ncol(), true);
			for (size_t j=cellcnt; j<cells.size(); j++) {
				if (narm && std::isnan(values[j])) continue;
				if (cells[j] <= cmax) {
					size_t k = cells[j] - cmin;
					if (newcell[k]) {
						v[k] = values[j];
						newcell[k] = false;
					} else {
						v[k] *= values[j];
					}
				} else {
					cellcnt = j;
					break;
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	} else if (fun == "pa") {
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), background);
			for (size_t j=cellcnt; j<cells.size(); j++) {
				if (narm && std::isnan(values[j])) continue;
				if (cells[j] <= cmax) {
					size_t k = cells[j] - cmin;
					v[k] = 1;
				} else {
					cellcnt = j;
					break;
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	} else if (fun == "first") {
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), background);
			std::vector<bool> newcell(out.bs.nrows[i] * out.ncol(), true);;
			for (size_t j=cellcnt; j<cells.size(); j++) {
				if (narm && std::isnan(values[j])) continue;
				if (cells[j] <= cmax) {
					size_t k = cells[j] - cmin;
					if (newcell[k]) {
						v[k] = values[j];
						newcell[k] = false;
					} 
				} else {
					cellcnt = j;
					break;
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	} else { // "last"
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), background);
			for (size_t j=cellcnt; j<cells.size(); j++) {
				if (narm && std::isnan(values[j])) continue;
				if (cells[j] <= cmax) {
					size_t k = cells[j] - cmin;
					v[k] = values[j];
				} else {
					cellcnt = j;
					break;
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	}
	out.writeStop();
	return out;
}


SpatRaster SpatRaster::rasterizePoints(SpatVector &x, std::string fun, std::vector<double> &values, bool narm, double background, SpatOptions &opt) {
	if (values.empty()) {
		values = std::vector<double>(x.nrow(), 1);
	}
	std::vector<std::vector<double>> pxy = x.coordinates();
	return rasterizePoints(pxy[0], pxy[1], fun, values, narm, background, opt);
}


SpatRaster SpatRaster::rasterizeGeom(SpatVector x, std::string unit, std::string fun, SpatOptions &opt) {

	if (x.type() == "points") {
		std::vector<double> v;
		return rasterizePoints(x, "count", v, false, 0.0, opt);

	} else {	

		SpatRaster out = geometry(1, false, false, false);
		SpatOptions ops(opt);

		std::vector<std::string> ss {"m", "km"};
		if (std::find(ss.begin(), ss.end(), unit) == ss.end()) {
			out.setError("invalid unit (not 'm' or 'km')");
			return out;
		}
		if ((x.type() == "lines")) {
			ss = {"count", "length", "crosses"};
			if (std::find(ss.begin(), ss.end(), fun) == ss.end()) {
				out.setError("invalid value for 'fun' (not 'count', 'crosses', or 'length')");
				return out;
			}
		} else {
			ss = {"area", "count"};
			if (std::find(ss.begin(), ss.end(), fun) == ss.end()) {
				out.setError("invalid value for 'fun' (not 'area' or 'count')");
				return out;
			}
		}

		SpatRaster empty = out.geometry();
		SpatExtent e = out.getExtent();
		double rsy = out.yres() / 2;

		double m = unit == "m" ? 1 : 1000;
		if (!x.is_lonlat()) {
			double tom = x.srs.to_meter();
			tom = std::isnan(tom) ? 1 : tom;
			m *= tom;
		}
		out.setNames({fun});
		opt.ncopies = std::max(opt.ncopies, (size_t)4) * 8;
		if (!out.writeStart(opt, filenames())) {
			return out;
		}
		for (size_t i=0; i < out.bs.n; i++) {
			e.ymax = yFromRow(out.bs.row[i]) + rsy;
			e.ymin = yFromRow(out.bs.row[i] + out.bs.nrows[i] - 1) - rsy;
			SpatRaster tmp = empty.crop(e, "near", false, ops);

			SpatVector p = tmp.as_polygons(true, false, false, false, false, 0, ops);
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), 0);

			if (fun == "crosses") {
				std::vector<int> r = p.relate(x, "crosses", true, true);
				size_t nx = x.size();
				for (size_t j=0; j< r.size(); j++) {
					size_t k= j / nx;
					v[k] += r[j];
				}
			} else {
				std::vector<long> cell(p.size());
				std::iota(cell.begin(), cell.end(), 0);
				p.df.add_column(cell, "cell");
				p = p.intersect(x, true);
				std::vector<double> stat;
				if (x.type() == "lines") {
					stat = p.length();
				} else {
					stat = p.area("m", false, {});
				}
				if (fun == "count") {
					for (size_t j=0; j<stat.size(); j++) {
						size_t k = p.df.iv[0][j];
						v[k]++;
					}
				} else {
					for (size_t j=0; j<stat.size(); j++) {
						size_t k = p.df.iv[0][j];
						v[k] += (stat[j] / m);
					}
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
		out.writeStop();
		return(out);

	}
}


SpatRaster SpatRaster::hardCopy(SpatOptions &opt) {
	SpatRaster out = geometry(nlyr(), true, true);
	if (!hasValues()) {
		out.addWarning("raster has no values");
		return out;
	}
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
		readBlock(v, out.bs, i);
		if (!out.writeBlock(v, i)) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}

/* gdalcopy
			GDALDatasetH hSrcDS = GDALOpenEx(out.source[0].filename.c_str(), GDAL_OF_RASTER | GDAL_OF_UPDATE, NULL, NULL, NULL);
			if(hSrcDS == NULL) {
				out.setError("cannot open dataset";)
				return false;
			}
			GDALDriverH hDriver = GDALGetDatasetDriver(hSrcDS);
			GDALDatasetH hDstDS = GDALCreateCopy( hDriver, filename.c_str(), hSrcDS, FALSE, NULL, NULL, NULL );
			GDALClose(hSrcDS);
			if(hDstDS == NULL) {
				out.setError("cannot create dataset";
				return false;
			}
			GDALClose(hDstDS);
*/

bool SpatRaster::getDSh(GDALDatasetH &rstDS, SpatRaster &out, std::string &filename, std::string &driver, double &naval, bool update, double background, SpatOptions &opt) {

	filename = opt.get_filename();
	SpatOptions ops(opt);
	ops.ncopies += 4;
	if (filename.empty()) {
		if (canProcessInMemory(ops)) {
			driver = "MEM";
		} else {
			filename = tempFile(opt.get_tempdir(), opt.tmpfile, ".tif");
			opt.set_filenames({filename});
			driver = "GTiff";
		}
	} else {
		driver = opt.get_filetype();
		getGDALdriver(filename, driver);
		if (driver.empty()) {
			out.setError("cannot guess file type from filename");
			return false;
		}
		std::string msg;
		if (!can_write({filename}, filenames(), opt.get_overwrite(), msg)) {
			out.setError(msg);
			return false;
		}
	}

	if (opt.names.size() == nlyr()) {
		out.setNames(opt.names);
	}


	if (update) {
		out = hardCopy(opt);
		//size_t ns = source.size();
		if (!out.open_gdal(rstDS, 0, true, opt)) {
			return false;
		}
	} else if (!out.create_gdalDS(rstDS, filename, driver, true, background, source[0].has_scale_offset, source[0].scale, source[0].offset, opt)) {
		out.setError("cannot create dataset");
		return false;
	}

	GDALRasterBandH hBand = GDALGetRasterBand(rstDS, 1);
	GDALDataType gdt = GDALGetRasterDataType(hBand);
	getNAvalue(gdt, naval);
	int hasNA;
	double naflag = GDALGetRasterNoDataValue(hBand, &hasNA);
	naval = hasNA ? naflag : naval;
	return true;
}


bool SpatRaster::getDShMEM(GDALDatasetH &rstDS, SpatRaster &out, double &naval, double background, SpatOptions &opt) {

	SpatOptions ops(opt);
	if (opt.names.size() == nlyr()) {
		out.setNames(opt.names);
	}
	if (!out.create_gdalDS(rstDS, "", "MEM", true, background, source[0].has_scale_offset, source[0].scale, source[0].offset, ops)) {
		out.setError("cannot create dataset");
		return false;
	}
	GDALRasterBandH hBand = GDALGetRasterBand(rstDS, 1);
	GDALDataType gdt = GDALGetRasterDataType(hBand);
	getNAvalue(gdt, naval);
	int hasNA;
	double naflag = GDALGetRasterNoDataValue(hBand, &hasNA);
	naval = hasNA ? naflag : naval;

	return true;
}




SpatRaster SpatRaster::rasterizeLyr(SpatVector x, double value, double background, bool touches, bool update, SpatOptions &opt) {

// not working well in some cases. See #552
	std::string gtype = x.type();
	SpatRaster out;
	out.setNames({"ID"});

	if ( !hasValues() ) update = false;
	if (update) { // all lyrs
		out = geometry();
	} else {
		out = geometry(1);
	}

	GDALDataset *vecDS = x.write_ogr("", "lyr", "Memory", false, true, std::vector<std::string>());
	if (x.hasError()) {
		out.setError(x.getError());
		return out;
	}

	OGRLayer *poLayer = vecDS->GetLayer(0);
#if GDAL_VERSION_MAJOR <= 2 && GDAL_VERSION_MINOR <= 2
	OGRLayerH hLyr = poLayer;
#else
	OGRLayerH hLyr = poLayer->ToHandle(poLayer);
#endif
    std::vector<OGRLayerH> ahLayers;
	ahLayers.push_back( hLyr );

	std::string driver, filename;
	GDALDatasetH rstDS;
	double naval;
	
	if (!opt.datatype_set) {
		if ((value < -16777216) || (value > 16777216)) {
			opt.datatype = "FLT8S";
		}
	}

	
	if (!getDSh(rstDS, out, filename, driver, naval, update, background, opt)) {
		return out;
	}
	if (std::isnan(value)) {
		// passing NULL instead may also work.
		value = naval;
	}

	std::vector<int> bands(out.nlyr());
	std::iota(bands.begin(), bands.end(), 1);
	std::vector<double> values(out.nlyr(), value);

	char** papszOptions = NULL;
	CPLErr err;
	if (touches) {
		papszOptions = CSLSetNameValue(papszOptions, "ALL_TOUCHED", "TRUE");
	}
	err = GDALRasterizeLayers(rstDS, static_cast<int>(bands.size()), &(bands[0]),
			1, &(ahLayers[0]), NULL, NULL, &(values[0]), papszOptions, NULL, NULL);

	CSLDestroy(papszOptions);

//	for (size_t i=0; i<ahGeometries.size(); i++) {
//		OGR_G_DestroyGeometry(ahGeometries[i]);
//	}
	GDALClose(vecDS);

	if ( err != CE_None ) {
		out.setError("rasterization failed");
		GDALClose(rstDS);
		return out;
	}

	if (driver == "MEM") {
		if (!out.from_gdalMEM(rstDS, false, true)) {
			out.setError("rasterization failed (mem)");
		}
	}

	GDALRasterBandH band = GDALGetRasterBand(rstDS, 1);
	double adfMinMax[2];
	GDALComputeRasterMinMax(band, false, adfMinMax);
	GDALSetRasterStatistics(band, adfMinMax[0], adfMinMax[1], -9999, -9999);

	GDALClose(rstDS);
	if (driver != "MEM") {
		out = SpatRaster(filename, {-1}, {""}, {}, {}, false, false, {});
	}
	return out;
}

#include "vecmath.h"
SpatRaster SpatRaster::rasterize(SpatVector x, std::string field, std::vector<double> values,
	double background, bool touches, std::string fun, bool weights, bool update, bool minmax, SpatOptions &opt) {
	
	std::string gtype = x.type();
	bool ispol = gtype == "polygons";
	if (weights) update = false;

	if (weights && ispol) {
		SpatOptions sopts(opt);
		SpatRaster wout = geometry(1);
		field = "";
		unsigned agx = 1000 / ncol();
		agx = std::max((unsigned)10, agx);
		unsigned agy = 1000 / nrow();
		agy = std::max((unsigned)10, agy);
		//unsigned agx = 100;
		//unsigned agy = 100;
		wout = wout.disaggregate({agx, agy}, sopts);
		double f = agx * agy;
		wout = wout.rasterize(x, field, {1/f}, background, touches, fun, false, false, false, sopts);
		wout = wout.aggregate({agx, agy}, "sum", true, opt);
		return wout;
	}

//	Rcpp::Rcout << "x" << std::endl;

	bool add = fun == "sum";

	SpatRaster out;
	if ( !hasValues() ) update = false;
	if (update) {
		out = geometry();
	} else {
		out = geometry(1);
	}
	if (field.empty()) {
		out.setNames({"layer"});
	} else {
		out.setNames({field});
	}

	size_t nGeoms = x.size();
	if (nGeoms == 0) {
		if (update) {
			out = *this;
		} else {
			out = out.init({background}, opt);
		} 
		return out;
	}

	if (ispol && touches && add) {
		add = false;
		out.addWarning("you cannot use 'sum' and 'touches' at the same time");
	}

	if (!field.empty()) {
		int i = x.df.get_fieldindex(field);
		if (i < 0) {
			out.setError("field " + field + " not found");
			return out;
		}
		std::string dt = x.df.get_datatype(field);
		if (dt == "double") {
			values = x.df.getD(i);
		} else if (dt == "long") {
			values = x.df.as_double(i);
			out.setValueType(1);
		} else if (dt == "bool") {
			values = x.df.as_double(i);
			out.setValueType(3);
		} else if (dt == "time") {
			// tbd
			values = x.df.as_double(i);
		} else {
			std::vector<std::string> sv = x.df.as_string(i);
			SpatFactor f(sv);
			values.resize(f.v.size());
			for (size_t j=0; j<values.size(); j++) {
				values[j] = f.v[j];
			}
			if (!add && !update) {
				std::vector<long> u(f.labels.size());
				std::iota(u.begin(), u.end(), 0);
				std::vector<std::string> nms = getNames();
				out.setLabels(0, u, f.labels, field);
			}
			if (add) {
				add = false;
				addWarning("cannot add factors");
			}
			out.setValueType(1);
		}
	}

	if (values.size() != nGeoms) {
		recycle(values, nGeoms);
	}

	if (!opt.datatype_set) {
		double vmn = vmin(values, true);
		double vmx = vmax(values, true);
		if ((vmn < -16777216) || (vmx > 16777216)) {
			opt.set_datatype("FLT8S");
		}
	}


	GDALDataset *vecDS = x.write_ogr("", "lyr", "Memory", false, true, std::vector<std::string>());
	if (x.hasError()) {
		out.setError(x.getError());
		return out;
	}
    std::vector<OGRGeometry *> ogrGeoms;
	ogrGeoms.reserve(nGeoms);

	OGRLayer *poLayer = vecDS->GetLayer(0);
	poLayer->ResetReading();

	OGRFeature *poFeature;
	while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
		OGRGeometry *poGeometry = poFeature->GetGeometryRef();
		if (poGeometry != NULL) {
			OGRGeometry *copyGeom = poGeometry->clone();
			ogrGeoms.push_back( copyGeom );
		}
		OGRFeature::DestroyFeature( poFeature );
	}
	GDALClose(vecDS);

	std::string errmsg, driver, filename;
	GDALDatasetH rstDS;
	double naval;
	if (add) {	background = 0;	}
	std::vector<int> bands(out.nlyr());
	std::iota(bands.begin(), bands.end(), 1);
	rep_each(values, out.nlyr());
	
	SpatRaster temp = out;
  	if (!out.writeStart(opt, filenames())) {
		return out;
	}
	
	bool hasError = false;
	SpatExtent e = temp.getExtent();
	SpatRaster tmp;
	SpatOptions topt(opt);

	char** papszOptions = NULL;

	for (size_t i = 0; i < out.bs.n; i++) {

		if (out.bs.n > 1) {
			double halfres = temp.yres() / 2;
			e.ymax = temp.yFromRow(out.bs.row[i]) + halfres;
			e.ymin = temp.yFromRow(out.bs.row[i] + out.bs.nrows[i] - 1) - halfres;
			if (update) {
				tmp = crop(e, "near", false, topt);	
			} else {
				tmp = temp.crop(e, "near", false, topt);
			}
		} else {
			if (update) {
				tmp = hardCopy(topt);
			} else {
				tmp = temp;
			}
		}
		if (update) {
			std::string filename;
			if (!tmp.getDSh(rstDS, tmp, filename, driver, naval, update, background, topt)) {
				return tmp;
			}
		} else if (!tmp.getDShMEM(rstDS, tmp, naval, background, opt)) {
			return tmp;
		}
		
	
		if (i==1) for (double &d : values) d = std::isnan(d) ? naval : d;

		CPLErr err;
		if (ispol && touches && (nGeoms > 1)) {
			// first to get the touches
			if (i == 0) papszOptions = CSLSetNameValue(papszOptions, "ALL_TOUCHED", "TRUE");
			err = GDALRasterizeGeometries(rstDS,
					static_cast<int>(bands.size()), &(bands[0]),
					static_cast<int>(ogrGeoms.size()),
					(OGRGeometryH *) &(ogrGeoms[0]),
					NULL, NULL, &(values[0]), papszOptions, NULL, NULL);

			if ( err != CE_None ) {
				tmp.setError("rasterization failed");
				GDALClose(rstDS);
				hasError = true; break;
			}
			//GDALFlushCache(rstDS);
			// second time to fix the internal area
			err = GDALRasterizeGeometries(rstDS,
					static_cast<int>(bands.size()), &(bands[0]),
					static_cast<int>(ogrGeoms.size()),
					(OGRGeometryH *) &(ogrGeoms[0]),
					NULL, NULL, &(values[0]), NULL, NULL, NULL);

		} else {
			if (i == 0)  {
				if (touches) {
					papszOptions = CSLSetNameValue(papszOptions, "ALL_TOUCHED", "TRUE");
				} else if (add) {
					papszOptions = CSLSetNameValue(papszOptions, "MERGE_ALG", "ADD");
				}
			}
			err = GDALRasterizeGeometries(rstDS,
					static_cast<int>(bands.size()), &(bands[0]),
					static_cast<int>(ogrGeoms.size()),
					(OGRGeometryH *) &(ogrGeoms[0]),
					NULL, NULL, &(values[0]), papszOptions, NULL, NULL);

		}

		if ( err != CE_None ) {
			tmp.setError("rasterization failed");
			GDALClose(rstDS);
			hasError = true; break;
		}

		if (!tmp.from_gdalMEM(rstDS, false, true)) {
			tmp.setError("rasterization failed");
			GDALClose(rstDS);
			hasError = true; break;
		}
		GDALClose(rstDS);

		std::vector<double> v = tmp.getValues(-1, topt);
		if (!out.writeBlock(v, i)) return out;
	}
	
	CSLDestroy(papszOptions);	
	for (size_t i=0; i<ogrGeoms.size(); i++) {
		OGR_G_DestroyGeometry(ogrGeoms[i]);
	}
	if (update) readStop();
	out.writeStop();
	
	if (hasError) return tmp; 
		
	return out;
}


std::vector<double> SpatRaster::rasterizeCells(SpatVector &v, bool touches, bool small, SpatOptions &opt) {
// note that this is only for lines and polygons
    SpatOptions ropt(opt);
	SpatRaster r = geometry(1);
	SpatExtent e = getExtent();
	SpatExtent ev = v.getExtent();
	if (ev.xmin >= ev.xmax) {
		double xr = 0.1 * xres();
		ev.xmin -= xr;
		ev.xmax += xr;
	}
	if (ev.ymin >= ev.ymax) {
		double yr = 0.1 * yres();
		ev.ymin -= yr;
		ev.ymax += yr;
	}

	e = e.intersect(ev);
	if ( !e.valid() ) {
		std::vector<double> out(1, NAN);
		return out;
	}

	SpatRaster rc = r.crop(e, "out", false, ropt);
	std::vector<double> feats(1, 1) ;
    SpatRaster rcr = rc.rasterize(v, "", feats, NAN, touches, "", false, false, false, ropt);
	SpatVector pts = rcr.as_points(false, true, false, ropt);
	std::vector<double> cells;
	if (pts.empty()) {
		if (small) {
			pts = v.as_points(false, true);
			SpatDataFrame vd = pts.getGeometryDF();
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			cells = r.cellFromXY(x, y);	
			cells.erase(std::remove_if(cells.begin(), cells.end(),
						[](const double& value) { return std::isnan(value); }), cells.end());
			std::sort( cells.begin(), cells.end() );
			cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
			if (cells.empty()) {
				cells.resize(1, NAN);
			}
		} else {
			cells.resize(1, NAN);			
		}
	} else {
		SpatDataFrame vd = pts.getGeometryDF();
		std::vector<double> x = vd.getD(0);
		std::vector<double> y = vd.getD(1);
		cells = r.cellFromXY(x, y);
//		cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
		if (cells.empty()) {
			cells.resize(1, NAN);
		}
	}
	return cells;
}

void SpatRaster::rasterizeCellsWeights(std::vector<double> &cells, std::vector<double> &weights, SpatVector &v, SpatOptions &opt) {
// note that this is only for polygons
    SpatOptions ropt(opt);
	//opt.progress = nrow()+1;
	SpatRaster r = geometry(1);
	//std::vector<unsigned> fact = {10, 10};
	SpatExtent e = getExtent();
	SpatExtent ve = v.getExtent();
	e = e.intersect(ve);
	if ( !e.valid() ) {
		return;
	}
	bool cropped = false;
	SpatRaster rc = r.crop(v.extent, "out", false, ropt);
	if ( ((ncol() > 1000) && ((ncol() / rc.ncol()) > 1.5))
			|| ((nrow() > 1000) && ((nrow() / rc.nrow()) > 1.5) )) {
		cropped = true;
		r = rc;
	}
	std::vector<double> feats;
	r = r.rasterize(v, "", feats, NAN, false, "", true, false, false, ropt);
	std::vector<std::vector<double>> cv = r.cells_notna(ropt);

	if (cv[0].empty()) {
		weights.resize(1);
		weights[0] = NAN;
		cells.resize(1);
		cells[0] = NAN;
	} else {
		weights = cv[1];
		if (cropped) {
			cv = r.xyFromCell(cv[0]);
			cells = cellFromXY(cv[0], cv[1]);
		} else {
			cells = cv[0];
		}
	}
	return;
}

void SpatRaster::rasterizeCellsExact(std::vector<double> &cells, std::vector<double> &weights, SpatVector &v, SpatOptions &opt) {

	SpatOptions ropt(opt);
	opt.progress = nrow()+1;
	SpatRaster r = geometry(1);
	r = r.crop(v.extent, "out", false, ropt);

//	if (r.ncell() < 1000) {
		std::vector<double> feats(1, 1) ;
		r = r.rasterize(v, "", feats, NAN, true, "", false, false, false, ropt);

		SpatVector pts = r.as_points(true, true, false, ropt);
		if (pts.empty()) {
			weights.resize(1);
			weights[0] = NAN;
			cells.resize(1);
			cells[0] = NAN;
		} else {
			SpatDataFrame vd = pts.getGeometryDF();
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			cells = cellFromXY(x, y);

			SpatVector rv = r.as_polygons(false, false, false, true, false, 0, ropt);
			std::vector<double> csize = rv.area("m", true, {});
			rv.df.add_column(csize, "area");
			rv.df.add_column(cells, "cells");
			rv = rv.crop(v);
			weights = rv.area("m", true, {});
			for (size_t i=0; i<weights.size(); i++) {
				weights[i] /= rv.df.dv[0][i];
			}
			cells = rv.df.dv[1];
		}
//	}

/*
// the below would need to remove all cells already included in the above
// because touches=false includes partly overlapped cells [#346]

	else {
		std::vector<double> feats(1, 1) ;
		SpatVector vv = v.as_lines();
		SpatRaster b = r.rasterize(vv, "", feats, NAN, true, false, false, false, false, opt);
		SpatVector pts = b.as_points(true, true, opt);
		if (pts.nrow() > 0) {
			SpatDataFrame vd = pts.getGeometryDF();
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			cells = cellFromXY(x, y);

			SpatVector bv = b.as_polygons(false, false, false, true, opt);
			std::vector<double> csize = bv.area("m", true, {});
			bv.df.add_column(csize, "cellsize");
			bv.df.add_column(cells, "cellnr");
			bv = bv.crop(v);
			weights = bv.area("m", true, {});
			for (size_t i=0; i<weights.size(); i++) {
				weights[i] /= bv.df.dv[0][i];
			}
			cells = bv.df.dv[1];
		}
		// touches = false
		r = r.rasterize(v, "", feats, NAN, false, false, false, false, false, opt);
		pts = r.as_points(true, true, opt);
		if (pts.nrow() > 0) {
			SpatDataFrame vd = pts.getGeometryDF();
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			std::vector<double> cells2 = cellFromXY(x, y);
			cells.insert(cells.end(), cells2.begin(), cells2.end());
			weights.resize(weights.size() + cells2.size(), 1);
		}

		if (cells.size() == 0) {
			weights.resize(1);
			weights[0] = NAN;
			cells.resize(1);
			cells[0] = NAN;
		}
	}
*/

}


void SpatRaster::rasterizeLinesLength(std::vector<double> &cells, std::vector<double> &weights, SpatVector &v, SpatOptions &opt) {

	if (v.type() != "lines") {
		setError("expected lines");
		return;
	}

/*
	double m = 1;
	if (!v.is_lonlat()) {
		double tom = v.srs.to_meter();
		tom = std::isnan(tom) ? 1 : tom;
		m *= tom;
	}
*/

	SpatOptions xopt(opt);
	xopt.ncopies = std::max(xopt.ncopies, (size_t)4) * 8;
	SpatRaster x = geometry(1);

	SpatExtent ev = v.getExtent();
	x = x.crop(ev, "out", false, xopt);
	BlockSize bs = x.getBlockSize(xopt);

	SpatExtent e = x.getExtent();
	double rsy = x.yres() / 2;
	for (size_t i=0; i < bs.n; i++) {
		e.ymax = yFromRow(bs.row[i]) + rsy;
		e.ymin = yFromRow(bs.row[i] + bs.nrows[i] - 1) - rsy;
		SpatRaster tmp = x.crop(e, "near", false, xopt);
		std::vector<double> cell(tmp.ncell());
		std::iota(cell.begin(), cell.end(), 0);
		std::vector<std::vector<double>> xy = tmp.xyFromCell(cell);
		cell = cellFromXY(xy[0], xy[1]);
		SpatVector p = tmp.as_polygons(true, false, false, false, false, 0, xopt);
		p.df.add_column(cell, "cell");
		p = p.intersect(v, true);
		if (p.nrow() > 1) {
			cells.insert(cells.end(), p.df.dv[0].begin(), p.df.dv[0].end());
			std::vector<double> w = p.length();
			double sm = std::accumulate(w.begin(), w.end(), 0.0);
			for (double &d : w) d /= sm;
			weights.insert(weights.end(), w.begin(), w.end());
		}
	}
}


