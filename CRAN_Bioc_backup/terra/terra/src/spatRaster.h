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

#include <fstream>
#include <numeric>
#include "spatVector.h"

#ifdef useGDAL
#include "gdal_priv.h"
#endif

#ifdef useRcpp
#include <Rcpp.h>
#endif


class SpatCategories {
	public:
		virtual ~SpatCategories(){}
		SpatDataFrame d;
		int index = 0;	
		bool combine(SpatCategories &x);
		bool concatenate(SpatCategories &x);
};


class SpatWindow {
	public:
		virtual ~SpatWindow(){}

		SpatExtent full_extent;
		size_t full_ncol, full_nrow, off_row, off_col;
		bool expanded = false;
		std::vector<size_t> expand;
};




class SpatRasterSource {
    private:
//		std::ofstream ofs;
	public:
#ifdef useGDAL
		GDALDataset* gdalconnection;
#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 4
		std::shared_ptr<GDALMDArray> m_array;
#endif
#endif
		bool open_read=false;
		bool open_write=false;

		SpatRasterSource();
		virtual ~SpatRasterSource(){}

//		void fsopen(std::string filename);
//		bool fswrite(std::vector<double> &v);
//		void fsclose();

		size_t ncol, nrow;
		size_t nlyr;
		size_t nlyrfile = 0;
		SpatExtent extent;
		bool extset=false;
		bool rotated=false;
		bool flipped=false;
		bool hasWindow=false;
		SpatWindow window;
	
		bool is_multidim = false;
		std::string m_arrayname;
		size_t m_ndims;
		std::vector<size_t> m_dims;
		std::vector<std::string> m_names;
//		std::vector<double> m_dimstart;
//		std::vector<double> m_dimend;
		std::vector<size_t> m_size;
		std::vector<size_t> m_order;
		std::vector<size_t> m_subset;
		bool m_hasNA = false;
		double m_missing_value;
		
		
		std::vector<std::vector<std::string>> bmdata;
		std::vector<std::string> smdata;
		
		//std::vector<std::string> crs = std::vector<std::string>(2, "");
		SpatSRS srs;
		std::vector<size_t> layers;
		// layer names 
		std::vector<std::string> names;
		// data source (sds) has one "variable name" / long_name
		std::string source_name;
		std::string source_name_long;
		
		std::vector<int64_t> time;
		std::string timestep = "seconds";
		std::string timezone = "";
		bool hasTime = false;
		std::vector<double> depth;
		std::string depthname = "depth";
		std::string depthunit = "";
		bool hasDepth = false;

		std::vector<std::string> unit;
		bool hasUnit = false;

		//std::vector< std::vector<double> values;
        std::vector<double> values;
        //std::vector<int64_t> ivalues;
        //std::vector<bool> bvalues;

//		unsigned char datatype;

		std::vector<int> blockrows;
		std::vector<int> blockcols;
		std::vector<bool> hasRange;
		std::vector<double> range_min;
		std::vector<double> range_max;
//		std::vector<bool> hasAttributes;
//		std::vector<SpatDataFrame> atts;
//		std::vector<int> attsIndex;
		std::vector<bool> hasCategories;
		std::vector<SpatCategories> cats;
		std::vector<unsigned char> valueType;
		// 0:double; 1:int; 3:bool

		//std::vector<std::string> dataType;

		std::vector<bool> hasColors;
		std::vector<SpatDataFrame> cols;
		SpatDataFrame legend;

		bool memory=true;
		bool hasValues=false;
		std::string filename;
		std::string driver;
		std::string dtype; 
		std::vector<std::string> open_ops;
		std::vector<std::string> open_drivers;
		
		// user set for reading:
		bool hasNAflag = false;
		double NAflag = NAN;

		std::vector<bool> has_scale_offset;
		std::vector<double> scale;
		std::vector<double> offset;

//		std::vector<SpatRasterSource> subset(std::vector<unsigned> lyrs);
		SpatRasterSource subset(std::vector<size_t> lyrs);
//		void getValues(std::vector<double> &v, unsigned lyr, SpatOptions &opt);
		void appendValues(std::vector<double> &v, size_t lyr);
		
		void setRange();
		void resize(size_t n);
		void reserve(size_t n);
		bool in_order(bool all);
		bool combine_sources(const SpatRasterSource &x);
		bool combine(SpatRasterSource &x);
		
		bool parameters_changed = false;		
		
		void set_names_time_ncdf(std::vector<std::string> metadata, std::vector<std::vector<std::string>> bandmeta, std::string &msg);
		void set_names_time_grib(std::vector<std::vector<std::string>> bandmeta, std::string &msg);
		void set_names_time_tif(std::vector<std::vector<std::string>> bandmeta, std::string &msg);
		
		std::vector<std::map<std::string, std::string>> lyrTags;
		void addLyrTag(size_t slyr, std::string name, std::string value);
		
};


class BlockSize {
	public:
		virtual ~BlockSize(){}
		std::vector<size_t> row;
		std::vector<size_t> nrows;
		size_t n;
};


class SpatRaster {

    private:
		std::string copy_driver = "";
		std::string copy_filename = "";
		std::vector<std::string> gdal_options;
		bool compute_stats = true;
		bool gdal_stats = false;
		bool gdal_approx = true;
		bool gdal_minmax = true;
	protected:
		SpatExtent window;

	public:

#ifdef useRcpp
		SpatProgress pbar;
		bool progressbar = false;
#endif

////////////////////////////////////////////////////
// properties and property-like methods for entire object
////////////////////////////////////////////////////
		
		std::vector<SpatRasterSource> source;

		BlockSize bs;
		//BlockSize getBlockSize(unsigned n, double frac, unsigned steps=0);
		BlockSize getBlockSize(SpatOptions &opt);
		std::vector<double> mem_needs(SpatOptions &opt);

		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		void addWarning(std::string s) { msg.addWarning(s); }
		void setMessage(std::string s) { msg.setMessage(s); }
		bool hasError() { return msg.has_error; }
		bool hasWarning() { return msg.has_warning; }
		std::vector<std::string> getWarnings() { return msg.getWarnings();}
		std::string getError() { return msg.getError();}
		std::string getMessage() { return msg.getMessage();}

		std::vector<std::vector<std::string>> user_tags;
		bool addTag(std::string name, std::string value, std::string domain);
		bool removeTag(std::string name, std::string domain);
		std::string getTag(std::string name, std::string domain);
		std::vector<std::vector<std::string>> getTags();

		void addLyrTags(std::vector<size_t> lyrs, std::vector<std::string> names, std::vector<std::string> values);

		bool removeLyrTags();
		bool removeLyrTag(size_t lyr, std::string name);
		std::string getLyrTag(size_t lyr, std::string name);
		std::vector<std::string> getLyrTags(std::vector<size_t> lyrs);
		std::vector<std::map<std::string, std::string>> getAllLyrTags();
		//double NA = std::numeric_limits<double>::quiet_NaN();

		size_t ncol();
		size_t nrow();
		SpatExtent getExtent();
		void setExtent(SpatExtent e);
//		void setExtent(SpatExtent ext, bool keepRes=false, std::string snap="");  // also set it for sources?
		void setExtent(SpatExtent ext, bool keepRes, bool no_expand, std::string snap);

		SpatVector dense_extent(bool inside, bool geobounds);

		//std::vector<std::string> getCRS();
		//void setCRS(std::vector<std::string> _crs);

		std::string getSRS(std::string x);
		bool setSRS(std::string crs);

		bool rgb=false;
		std::string rgbtype;
		std::vector<int> rgblyrs;
		bool setRGB(int r, int g, int b, int alpha, std::string type);
		std::vector<int> getRGB();
		void removeRGB();

		std::vector<bool> is_multidim();
		std::vector<std::vector<std::string>> dim_names();
		std::vector<std::vector<size_t>> dim_order();
		std::vector<std::vector<size_t>> dim_size();

/*
#ifdef useGDAL	
		bool setSRS(OGRSpatialReference *poSRS, std::string &msg) {
#endif 
*/

		bool is_lonlat();
		bool could_be_lonlat();
		bool is_global_lonlat();
		int ns_polar();

		std::vector<bool> is_flipped();

		std::vector<double> resolution();
		SpatRaster setResolution(double xres, double yres);
		double ncell() { return (double)nrow() * (double)ncol(); }
		double size() { return (double)ncol() * (double)nrow() * (double)nlyr() ; }

		std::vector<bool> is_rotated();

		double xres();
		double yres();
		std::vector<double> origin();
		size_t nlyr();

		// only no values allowed with a single SpatRasterSource
		bool hasValues();
		std::vector<double> getValues(long lyr, SpatOptions &opt);
		
		bool getValuesSource(size_t src, std::vector<double> &out);				
		bool setValues(std::vector<double> &v, SpatOptions &opt);

#ifdef useRcpp
		bool setValuesRcpp(Rcpp::NumericVector &v, SpatOptions &opt);
#endif

		bool replaceCellValues(std::vector<double> &cells, std::vector<double> &v, bool bylyr, SpatOptions &opt);
		bool replaceCellValuesLayer(std::vector<size_t> layers, std::vector<double> &cells, std::vector<double> &v, bool bylyr, SpatOptions &opt);

		void setRange(SpatOptions &opt, bool force);
		
////////////////////////////////////////////////////
// property like methods for RasterSources
////////////////////////////////////////////////////
		std::vector<std::string> filenames();
		bool isSource(std::string filename);
		std::vector<bool> inMemory();

////////////////////////////////////////////////////
// property like methods for layers
////////////////////////////////////////////////////

		std::vector<bool> hasRange();
		std::vector<double> range_min();
		std::vector<double> range_max();

		std::vector<int> getValueType(bool unique);
		bool setValueType(unsigned char d);


		std::vector<std::string> getNames();
		bool setNames(std::vector<std::string> names, bool make_valid=false);

		std::vector<std::string> getSourceNames();
		bool setSourceNames(std::vector<std::string>);
		std::vector<std::string> getLongSourceNames();
		bool setLongSourceNames(std::vector<std::string>);

		bool hasTime();
		std::vector<int64_t> getTime();
		std::string getTimeStep();
		std::string getTimeZone();
		std::vector<std::string> getTimeStr(bool addstep, std::string timesep);
		bool setTime(std::vector<int64_t> time, std::string step, std::string zone);
		
		bool hasDepth();
		std::vector<double> getDepth();
		bool setDepth(std::vector<double> depths);
		std::string getDepthName();
		bool setDepthName(std::string name);
		std::string getDepthUnit();
		bool setDepthUnit(std::string unit);
		
		bool hasUnit();
		std::vector<std::string> getUnit();
		bool setUnit(std::vector<std::string> units);

		bool setNAflag(std::vector<double> flag);
		std::vector<double> getNAflag();

		std::vector<std::vector<std::string>> getMetadata(bool layers);

		std::vector<bool> isMD();

////////////////////////////////////////////////////
// constructors
////////////////////////////////////////////////////

		SpatRaster();
		SpatRaster(size_t nr, size_t nc, size_t nl, SpatExtent ext, std::string crs);
		SpatRaster(std::vector<size_t> rcl, std::vector<double> ext, std::string crs);
		SpatRaster(std::vector<std::string> fname, std::vector<int> subds, std::vector<std::string> subdsname, bool multi, std::vector<std::string> drivers, std::vector<std::string> options, std::vector<int> dims, bool noflip, bool guessCRS, std::vector<std::string> domains);
		SpatRaster(std::string fname, std::vector<int> subds, std::vector<std::string> subdsname, std::vector<std::string> drivers, std::vector<std::string> options, bool noflip, bool guessCRS, std::vector<std::string> domains);
		SpatRaster(SpatRasterSource &s);
		virtual ~SpatRaster(){}

		void setSource(SpatRasterSource &s);
		void setSources(std::vector<SpatRasterSource> &s);
		//SpatRaster(const SpatRaster& x);


        SpatRaster deepCopy();
		SpatRaster hardCopy(SpatOptions &opt);
        SpatRaster geometry(size_t nlyrs=0, bool properties=false, bool time=true, bool units=false, bool tags=false);
		SpatRaster geometry_opt(long nlyrs, bool properties, bool time, bool units, bool tags, bool datatype, SpatOptions &opt);

		bool constructFromFile(std::string fname, std::vector<int> subds, std::vector<std::string> subdsname, std::vector<std::string> drivers, std::vector<std::string> options, bool noflip, bool guessCRS, std::vector<std::string> domains);
		bool constructFromFileMulti(std::string fname, std::vector<int> subds, std::vector<std::string> subname, std::vector<std::string> drivers, std::vector<std::string> options, std::vector<int> dims, bool noflip, bool guessCRS, std::vector<std::string> domains);
		bool constructFromSDS(std::string filename, std::vector<std::string> meta, std::vector<int> subds, std::vector<std::string> subdsname, std::vector<std::string> options, std::string driver, bool noflip, bool guessCRS, std::vector<std::string> domains);

		
		//SpatRaster fromFiles(std::vector<std::string> fname, std::vector<int> subds, std::vector<std::string> subdsname, std::string drivers, std::vector<std::string> options);
		
//		bool constructFromNCDFsds(std::string filename, std::vector<std::string> meta, std::vector<int> subds, std::vector<std::string> subdsname);

		void addSource(SpatRaster &x, bool warn, SpatOptions &opt);	
		void checkDepth(SpatRaster &x);	
		void checkTime(SpatRaster &x);	
		SpatRaster combineSources(SpatRaster &x, bool warn);
		void combine(SpatRaster &x);
		
		SpatRaster subsetSource(size_t snr);
		SpatRaster subset(std::vector<size_t> lyrs, SpatOptions &opt);
		SpatRaster replace(SpatRaster x, size_t layer, SpatOptions &opt);
////////////////////////////////////////////////////
// helper methods
////////////////////////////////////////////////////


		std::vector<std::string> getAllFiles();

		void gdalogrproj_init(std::string path);

		bool compare_geom(SpatRaster &x, bool lyrs, bool crs, double tol, bool warncrs=false, bool ext=true, bool rowcol=true, bool res=false);
		bool compare_origin(std::vector<double> x, double tol);
		bool shared_basegeom(SpatRaster &x, double tol, bool test_overlap);

		std::vector<double> cellFromXY (std::vector<double> x, std::vector<double> y, double missing=NAN);
		double cellFromXY(double x, double y, double missing=NAN);
		std::vector<double> cellFromRowCol(std::vector<int64_t> row, std::vector<int64_t> col);
		double cellFromRowCol(int64_t row, int64_t col);
		std::vector<double> cellFromRowColCombine(std::vector<int64_t> row, std::vector<int64_t> col);
		double cellFromRowColCombine(int64_t row, int64_t col);
		std::vector<double> yFromRow(const std::vector<int64_t> &row);
		double yFromRow(int64_t row);
		void yFromRow(std::vector<double> &y);
		std::vector<double> xFromCol(const std::vector<int64_t> &col);
		double xFromCol(int64_t col);
		void xFromCol(std::vector<double> &x);

		std::vector<int64_t> colFromX(const std::vector<double> &x);
		int64_t colFromX(double x);
		std::vector<int64_t> rowFromY(const std::vector<double> &y);
		int64_t rowFromY(double y);
		void xyFromCell( std::vector<std::vector<double>> &xy );
		std::vector<std::vector<double>> xyFromCell( std::vector<double> &cell);
		std::vector<std::vector<double>> xyFromCell( double cell);
		std::vector<std::vector<int64_t>> rowColFromCell(std::vector<double> &cell);
		std::vector<int64_t> rowColFromY(std::vector<double> &y);
		std::vector<std::vector<int64_t>> rowColFromExtent(SpatExtent e);
		std::vector<std::vector<double>> coordinates(bool narm, bool nall, SpatOptions &opt);
		
        std::vector<size_t> sourcesFromLyrs(std::vector<size_t> lyrs);
		int sourceFromLyr(size_t lyr);
		std::vector<size_t> findLyr(size_t lyr);
		std::vector<size_t> getBands();

        std::vector<size_t> nlyrBySource();
        std::vector<size_t> lyrsBySource();
        size_t nsrc();

		SpatRaster makeCategorical(long layer, SpatOptions &opt);
		bool createCategories(size_t layer, SpatOptions &opt);
		std::vector<bool> hasCategories();
		bool setCategories(size_t layer, SpatDataFrame d, size_t index);
		bool removeCategories(long layer);
		std::vector<SpatCategories> getCategories();
		SpatCategories getLayerCategories(size_t layer);
		std::vector<std::string> getLabels(size_t layer);
		bool setLabels(size_t layer, std::vector<long> value, std::vector<std::string> labels, std::string name);
		int getCatIndex(size_t layer);
		bool setCatIndex(size_t layer, int index);

		bool hasLegend();
		bool setLegend(SpatDataFrame x);
		SpatDataFrame getLegend();

		
		bool hasScaleOffset();
		bool setScaleOffset(std::vector<double> sc, std::vector<double> of);
		std::vector<std::vector<double>> getScaleOffset();

		//bool setAttrIndex(size_t layer, int i);
		//std::vector<int> getAttrIndex();
		//void createAttributes(size_t layer);
		//std::vector<bool> hasAttributes();
		//void setAttributes(size_t layer, SpatDataFrame df);
		//std::vector<SpatDataFrame> getAttributes();
		//SpatDataFrame getLayerAttributes(size_t layer);
		
		std::vector<bool> hasColors();
		std::vector<SpatDataFrame> getColors();
		bool setColors(size_t layer, SpatDataFrame cols);
		bool removeColors(size_t layer);

		double valuesCell(double);
		double valuesCell(int, int);
		std::vector<double> valuesCell(std::vector<double>);
		std::vector<double> valuesRow(int);

////////////////////////////////////////////////////
// read and write
////////////////////////////////////////////////////

		bool valid_sources(bool files=true, bool rotated=true);
		bool readStart();
		std::vector<double> readValuesR(size_t row, size_t nrows, size_t col, size_t ncols);
		void readValues(std::vector<double> &out, size_t row, size_t nrows, size_t col, size_t ncols);
		void readValuesWhileWriting(std::vector<double> &out, size_t row, size_t nrows, size_t col, size_t ncols);
		void readChunkMEM(std::vector<double> &out, size_t src, size_t row, size_t nrows, size_t col, size_t ncols);

		void readBlock(std::vector<double> &v, BlockSize bs, size_t i){ // inline
			readValues(v, bs.row[i], bs.nrows[i], 0, ncol());
		}

		void readBlock2(std::vector<std::vector<double>> &v, BlockSize bs, size_t i);
		void readBlockIP(std::vector<double> &x, BlockSize bs, size_t i);		
		std::vector<double> readExtent(SpatExtent e);
		bool readStop();

		bool readAll();

		bool writeStart(SpatOptions &opt, std::vector<std::string> srcnames);

		bool writeBlock(std::vector<double> &v, size_t i){ // inline
			// for debugging?
			// if (bs.row.size() <= i) {
			//    setError("invalid block number"); return false;	
			// }
			return writeValues(v, bs.row[i], bs.nrows[i]);
		}

		bool writeValues(std::vector<double> &vals, size_t startrow, size_t nrows);
		bool writeValuesRect(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols);
		bool writeValuesRectRast(SpatRaster &r, SpatOptions& opt);
		
		//bool writeValues2(std::vector<std::vector<double>> &vals, size_t startrow, size_t nrows);
		bool writeStop();
		bool writeHDR(std::string filename);
		std::string make_vrt(std::vector<std::string> filenames, std::vector<std::string> options, SpatOptions &opt);
		bool write_aux_json(std::string filename);

		//bool writeStartGDAL(std::string filename, std::string driver, std::string datatype, bool overwrite, SpatOptions &opt);
		bool writeStartGDAL(SpatOptions &opt, const std::vector<std::string> &srcnames);		
		bool fillValuesGDAL(double fillvalue);
		bool writeValuesGDAL(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols);
		bool writeStopGDAL();
		
		bool writeStartMulti(SpatOptions &opt, const std::vector<std::string> &srcnames);
		bool writeValuesMulti(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols);
		bool writeStopMulti();
		
		bool getTempFile(std::string &filename, std::string &driver, SpatOptions& opt);

		bool readStartMulti(size_t src);
		bool readStopMulti(size_t src);
		bool readChunkMulti(std::vector<double> &data, size_t src, size_t row, size_t nrows, size_t col, size_t ncols);
		std::vector<double> readValuesMulti(size_t src, size_t row, size_t nrows, size_t col, size_t ncols, int lyr);
		std::vector<double> readSampleMulti(size_t src, size_t srows, size_t scols, bool overview);
		bool readRowColMulti(size_t src, std::vector<std::vector<double>> &out, size_t outstart, std::vector<int64_t> &rows, const std::vector<int64_t> &cols);

		//bool writeStartBinary(std::string filename, std::string datatype, std::string bandorder, bool overwrite);
		//bool writeValuesBinary(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols);

		bool writeValuesMem(std::vector<double> &vals, size_t startrow, size_t nrows);
		bool writeValuesMemRect(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols);

		// binary (flat) source
		//std::vector<double> readValuesBinary(size_t src, size_t row, size_t nrows, size_t col, size_t ncols);
		//std::vector<double> readSampleBinary(size_t src, size_t srows, size_t scols);
		//std::vector<std::vector<double>> readCellsBinary(size_t src, std::vector<double> cells);

		// gdal source
		std::vector<double> readValuesGDAL(size_t src, size_t row, size_t nrows, size_t col, size_t ncols, int lyr = -1);
		std::vector<double> readGDALsample(size_t src, size_t srows, size_t scols, bool overview);

		void readRowColGDAL(size_t src, std::vector<std::vector<double>> &out, size_t outstart, std::vector<int64_t> &rows, const std::vector<int64_t> &cols);

//		std::vector<double> readRowColGDALFlat(size_t src, std::vector<int64_t> &rows, const std::vector<int64_t> &cols);

//		std::vector<double> readRowColBlockFlat(size_t src, std::vector<int64_t> &rows, std::vector<int64_t> &cols);
//		std::vector<std::vector<double>> readRowColBlock(size_t src, std::vector<int64_t> &rows, std::vector<int64_t> &cols);

		void readRowColBlock(size_t src, std::vector<std::vector<double>> &out, size_t outstart, std::vector<int64_t> &rows, std::vector<int64_t> &cols, SpatOptions &opt);


		bool readStartGDAL(size_t src);
		bool readStopGDAL(size_t src);
		void readChunkGDAL(std::vector<double> &data, size_t src, size_t row, size_t nrows, size_t col, size_t ncols);

		bool setWindow(SpatExtent x);
		bool removeWindow();
		std::vector<bool> hasWindow();

		void openFS(std::string const &filename);

		SpatRaster writeRaster(SpatOptions &opt);
		SpatRaster writeRasterM(SpatOptions &opt);
		SpatRaster writeTempRaster(SpatOptions &opt);
		bool writeDelim(std::string filename, std::string delim, bool cell, bool xy, SpatOptions &opt);
		bool update_meta(bool names, bool crs, bool ext, SpatOptions &opt);

		//SpatRaster writeRasterGDAL(std::string filename, std::string format, std::string datatype, bool overwrite, SpatOptions &opt);
		//SpatRaster writeRasterBinary(std::string filename, std::string datatype, std::string bandorder, bool overwrite);
		//bool checkFormatRequirements(const std::string &driver, std::string &filename);

		bool canProcessInMemory(SpatOptions &opt);
		size_t chunkSize(SpatOptions &opt);

		void fill(double x);

		SpatRaster sources_to_disk(std::vector<std::string> &tmpfs, bool unique, SpatOptions &opt);
		bool sources_from_file();

		std::vector<int> getFileBlocksize();

////////////////////////////////////////////////////
// main methods
////////////////////////////////////////////////////

		SpatRaster collapse_sources();
		void collapse();
		
		SpatRaster rectify(std::string method, SpatRaster aoi, unsigned useaoi, bool snap, SpatOptions &opt);
		
        std::vector<double> adjacent(std::vector<double> cells, std::string directions, bool include);
        std::vector<double> adjacentMat(std::vector<double> cells, std::vector<bool> mat, std::vector<size_t> dim, bool include);
 		SpatRaster aggregate(std::vector<size_t> fact, std::string fun, bool narm, SpatOptions &opt);
		SpatExtent align(SpatExtent e, std::string snap);
		SpatRaster rst_area(bool mask, std::string unit, bool transform, int rcmax, SpatOptions &opt);
		std::vector<std::vector<double>> sum_area(std::string unit, bool transform, bool by_value, SpatOptions &opt);
		std::vector<std::vector<double>> sum_area_group(SpatRaster group, std::string unit, bool transform, bool by_value, SpatOptions &opt);
		SpatRaster surfaceArea(SpatOptions &opt);

		SpatRaster roll(size_t n, std::string fun, std::string type, bool circular, bool narm, SpatOptions &opt);

		SpatRaster arith(SpatRaster x, std::string oper, bool falseNA, SpatOptions &opt);
		SpatRaster arith(double x, std::string oper, bool reverse, bool falseNA, SpatOptions &opt);
		SpatRaster arith(std::vector<double> x, std::string oper, bool reverse, bool falseNA, SpatOptions &opt);
		SpatRaster arith_m(std::vector<double> x, std::string oper, std::vector<size_t> dim, bool reverse, SpatOptions &opt);

		SpatRaster apply(std::vector<size_t> ind, std::string fun, bool narm, std::vector<std::string> nms, std::vector<int64_t> time, std::string timestep, std::string timezone, SpatOptions &opt);
	
		SpatRaster rapply(SpatRaster x, double first, double last, std::string fun, bool clamp, bool narm, bool circular, SpatOptions &opt);
		std::vector<std::vector<double>> rappvals(SpatRaster x, double first, double last, bool clamp, bool all, double fill, size_t startrow, size_t nrows, bool circular);
		SpatRaster fill_range(long limit, bool circular, SpatOptions &opt);

		SpatVector as_polygons(bool round, bool dissolve, bool values, bool narm, bool nall, int digits, SpatOptions &opt);
		SpatVector polygonize(bool round, bool values, bool narm, bool aggregate, int digits, SpatOptions &opt);
		
		SpatVector as_lines(SpatOptions &opt);
		SpatVector as_points(bool values, bool narm, bool nall, SpatOptions &opt);
		std::vector<std::vector<double>> as_points_value(const double& target, SpatOptions &opt);
		std::vector<std::vector<double>> cells_notna(SpatOptions &opt);
		std::vector<double> cells_notna_novalues(SpatOptions &opt);


		SpatVector as_multipoints(bool narm, bool nall, SpatOptions &opt);
		SpatRaster atan_2(SpatRaster x, SpatOptions &opt);


		void bilinearValues(std::vector<std::vector<double>> &out, const std::vector<double> &x, const std::vector<double> &y, SpatOptions &opt);
		std::vector<double> bilinearCells(const std::vector<double> &x, const std::vector<double> &y);
		void fourCellsFromXY(std::vector<double> &out, const std::vector<double> &x, const std::vector<double> &y);

		SpatRaster buffer(double d, double background, SpatOptions &opt);
		SpatRaster clamp(std::vector<double> low, std::vector<double> high, bool usevalue, SpatOptions &opt);
		SpatRaster clamp_raster(SpatRaster &x, SpatRaster &y, std::vector<double> low, std::vector<double> high, bool usevalue, SpatOptions &opt);
		SpatRaster clamp_ts(bool min, bool max, SpatOptions &opt);

		SpatRaster combineCats(SpatRaster x, SpatOptions &opt);
		SpatRaster dropLevels();

		SpatRaster cover(SpatRaster x, std::vector<double> value, SpatOptions &opt);
		SpatRaster cover(std::vector<double> value, SpatOptions &opt);

		SpatRaster crop(SpatExtent e, std::string snap, bool expand, SpatOptions &opt);
		SpatRaster cropmask(SpatVector &v, std::string snap, bool touches, bool expand, SpatOptions &opt);
		SpatRaster cum(std::string fun, bool narm, SpatOptions &opt);
        SpatRaster disaggregate(std::vector<size_t> fact, SpatOptions &opt);
		SpatRaster proximity(double target, double exclude, bool keepNA, std::string unit, bool buffer, double maxdist, bool remove_zero, SpatOptions &opt);
		SpatRaster fillNA(double missing, double maxdist, int niter, SpatOptions &opt);
		
		SpatRaster distance(double target, double exclude, bool keepNA, std::string unit, bool remove_zero, std::string method, bool values, double threshold, SpatOptions &opt);
		SpatRaster nearest(double target, double exclude, bool keepNA, std::string unit, bool remove_zero, std::string method, SpatOptions &opt);

//		SpatRaster distance_spatvector(SpatVector p, std::string unit, const std::string& method, SpatOptions &opt);
//		SpatRaster distance_rasterize(SpatVector p, double target, double exclude, std::string unit, const std::string& method, SpatOptions &opt);
		SpatRaster distance_vector(SpatVector p, bool rasterize, std::string unit, const std::string& method, SpatOptions &opt);

		SpatRaster direction_rasterize(SpatVector p, bool from, bool degrees, double target, double exclude, const std::string& method, SpatOptions &opt);
		
		
		SpatRaster distance_crds(std::vector<double>& x, std::vector<double>& y, const std::string& method, bool skip, bool setNA, std::string unit,double threshold, SpatOptions &opt);
		SpatRaster distance_crds_vals(std::vector<double>& x, std::vector<double>& y, std::vector<double>& v, const std::string& method, bool skip, bool setNA, std::string unit, double threshold, SpatOptions &opt);

		SpatRaster dn_crds(std::vector<double>& x, std::vector<double>& y, const std::string& method, bool skip, bool setNA, std::string unit, SpatOptions &opt);

		SpatRaster direction(bool from, bool degrees, double target, double exclude, const std::string& method, SpatOptions &opt);
		SpatRaster direction_vector(SpatVector p, bool from, bool degrees, const std::string& method, SpatOptions &opt);
		
		SpatRaster clumps(int directions, bool zeroAsNA, SpatOptions &opt);
		SpatRaster patches(size_t directions, SpatOptions &opt);


		SpatRaster edges(bool classes, std::string type, unsigned directions, double falseval, SpatOptions &opt);
		SpatRaster extend(SpatExtent e, std::string snap, double fill, SpatOptions &opt);
		std::vector<std::vector<std::vector<double>>> extractVector(SpatVector v, bool touches, bool small, std::string method, bool cells, bool xy, bool weights, bool exact, SpatOptions &opt);
		std::vector<double> extractVectorFlat(SpatVector v, std::vector<std::string> funs, bool narm, bool touches, bool small, std::string method, bool cells, bool xy, bool weights, bool exact, SpatOptions &opt);
		std::vector<std::vector<double>> extractBuffer(const std::vector<double> &x, const std::vector<double> &y, double b, SpatOptions &opt);
		
		
//		std::vector<double> extract_interpolate(std::vector<double> x, std::vector<double> y, std::string algo);

		
		std::vector<double> vectCells(SpatVector v, bool touches, bool small, std::string method, bool weights, bool exact, SpatOptions &opt);
		std::vector<double> extCells(SpatExtent ext);

		std::vector<std::vector<double>> extractCell(std::vector<double> &cell, SpatOptions opt);
//		std::vector<double> extractCellFlat(std::vector<double> &cell);
	
		std::vector<std::vector<double>> extractXY(const std::vector<double> &x, const std::vector<double> &y, const std::string & method, const bool &cells, SpatOptions &opt);
		std::vector<double> extractXYFlat(const std::vector<double> &x, const std::vector<double> &y, const std::string & method, const bool &cells, SpatOptions &opt);
		
		SpatRaster flip(bool vertical, SpatOptions &opt);
		SpatRaster filler(SpatRaster x, SpatOptions &opt);
		
		SpatRaster focal(std::vector<unsigned> w, std::vector<double> m, double fillvalue, bool narm, bool naonly, bool naomit, std::string fun, bool expand, SpatOptions &opt);

		std::vector<double> focal_values(std::vector<unsigned> w, double fillvalue, int64_t row, int64_t nrows, SpatOptions &opt);
		std::vector<std::vector<double>> freq(bool bylayer, bool round, int digits, SpatOptions &opt);
		std::vector<size_t> count(double value, bool bylayer, bool round, int digits, SpatOptions &opt);
		
		bool get_aggregate_dims(std::vector<size_t> &fact, std::string &message);
		std::vector<size_t> get_aggregate_dims2(std::vector<size_t> fact);
		std::vector<std::vector<double> > get_aggregates(std::vector<double> &in, size_t nr, std::vector<size_t> dim);
//		std::vector<double> compute_aggregates(std::vector<double> &in, unsigned nr, std::vector<unsigned> dim, std::function<double(std::vector<double>&, bool)> fun, bool narm);

		SpatDataFrame mglobal(std::vector<std::string> funs, bool narm, SpatOptions &opt);

		SpatDataFrame global(std::string fun, bool narm, SpatOptions &opt);
		SpatDataFrame globalTF(std::string fun, SpatOptions &opt);
		SpatDataFrame global_weighted_mean(SpatRaster &weights, std::string fun, bool narm, SpatOptions &opt);

		SpatRaster gridDistance(double m, SpatOptions &opt);
		SpatRaster costDistanceRun(SpatRaster &old, bool &converged, double target, double m, bool lonlat, bool global, bool npole, bool spole, bool grid, SpatOptions &opt);
		SpatRaster costDistance(double target, double m, size_t maxiter, bool grid, SpatOptions &opt);

		SpatRaster init(std::string value, bool plusone, SpatOptions &opt);
		SpatRaster init(std::vector<double> values, SpatOptions &opt);
		
		SpatRaster is_in(std::vector<double> m, SpatOptions &opt);
		std::vector<std::vector<double>> is_in_cells(std::vector<double> m, bool keepvalue, SpatOptions &opt);

		std::vector<std::string> getDataType(bool unique, bool memtype);
		std::vector<std::string> dataType();

		SpatRaster isnot(bool falseNA, SpatOptions &opt);
		SpatRaster isnan(bool falseNA, SpatOptions &opt);
		SpatRaster isnotnan(bool falseNA, SpatOptions &opt);
		SpatRaster countnan(long n, SpatOptions &opt);

		SpatRaster isfinite(bool falseNA, SpatOptions &opt);
		SpatRaster isinfinite(bool falseNA, SpatOptions &opt);
		SpatRaster is_true(bool falseNA, SpatOptions &opt);
		SpatRaster is_false(bool falseNA, SpatOptions &opt);
		SpatRaster not_na(bool falseNA, SpatOptions &opt);

		SpatRaster allnan(bool falseNA, SpatOptions &opt);
		SpatRaster anynan(bool falseNA, SpatOptions &opt);
		SpatRaster nonan(bool falseNA, SpatOptions &opt);
		SpatRaster which(SpatOptions &opt);

		std::vector<std::vector<double>> layerCor(std::string fun, std::string use, bool asSample, SpatOptions &opt);

		std::vector<double> line_cells(SpatGeom& g);
		SpatRaster logic(SpatRaster x, std::string oper, SpatOptions &opt);
		SpatRaster logic(double x, std::string oper, SpatOptions &opt);
		SpatRaster logic(std::vector<double> x, std::string oper, SpatOptions &opt);

		SpatExtent ext_from_rc(int64_t r1, int64_t r2, int64_t c1, int64_t c2);
		SpatExtent ext_from_cell(double cell);

		std::vector<double> get_tiles_extent(SpatRaster x, bool expand, std::vector<int> buffer);
		std::vector<std::string> make_tiles(SpatRaster x, bool expand, std::vector<int> buffer, bool narm, std::string filename, SpatOptions &opt);
		
		std::vector<double> get_tiles_extent_vect(SpatVector x, bool expand, std::vector<int> buffer);
		std::vector<std::string> make_tiles_vect(SpatVector x, bool expand, std::vector<int> buffer, bool narm, std::string filename, SpatOptions &opt);

		SpatRaster mask(SpatRaster &x, bool inverse, double maskvalue, double updatevalue, SpatOptions &opt);
		SpatRaster mask(SpatRaster &x, bool inverse, std::vector<double> maskvalues, double updatevalue, SpatOptions &opt);
		SpatRaster mask(SpatOptions &opt);


		SpatRaster mask(SpatVector &x, bool inverse, double updatevalue, bool touches, SpatOptions &opt);
		SpatRaster math(std::string fun, SpatOptions &opt);
		SpatRaster math2(std::string fun, unsigned digits, SpatOptions &opt);


		SpatRaster separate(std::vector<double> classes, double keepvalue, double othervalue, bool round, int digits, SpatOptions &opt);

		SpatRaster modal(std::vector<double> add, std::string ties, bool narm, SpatOptions &opt);

        std::vector<double> polygon_cells(SpatGeom& g);
		SpatRaster quantile(std::vector<double> probs, bool narm, SpatOptions &opt);
		SpatRaster stretch(std::vector<double> minv, std::vector<double> maxv, std::vector<double> minq, std::vector<double> maxq, std::vector<double> smin, std::vector<double> smax, SpatOptions &opt);
		SpatRaster reverse(SpatOptions &opt);

		SpatRaster range(std::vector<double> add, bool narm, SpatOptions &opt);

		SpatRaster rasterizeLyr(SpatVector x, double value, double background, bool touches, bool update, SpatOptions &opt);

		SpatRaster rasterize(SpatVector x, std::string field, std::vector<double> values, double background, bool touches, std::string fun, bool weights, bool update, bool minmax, SpatOptions &opt);
		
		SpatRaster rasterizeWindow(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::string algo, std::vector<double> algops, SpatOptions &opt);

		std::vector<std::vector<double>> win_circle(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> win, SpatOptions &opt);
		std::vector<std::vector<double>> win_rect(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> win, SpatOptions &opt);
		
		std::vector<double> rasterizeCells(SpatVector &v, bool touches, bool small, SpatOptions &opt);
		//std::vector<std::vector<double>> rasterizeCellsWeights(SpatVector &v, bool touches);
		SpatRaster rasterizeGeom(SpatVector x, std::string unit, std::string count, SpatOptions &opt);
		SpatRaster rasterizePoints(std::vector<double>&x, std::vector<double> &y, std::string fun, std::vector<double> &values, bool narm, double background, SpatOptions &opt);
		SpatRaster rasterizePoints(SpatVector &x, std::string fun, std::vector<double> &values, bool narm, double background, SpatOptions &opt);
		void rasterizeCellsWeights(std::vector<double> &cells, std::vector<double> &weights, SpatVector &v, SpatOptions &opt); 
		void rasterizeCellsExact(std::vector<double> &cells, std::vector<double> &weights, SpatVector &v, SpatOptions &opt); 
		void rasterizeLinesLength(std::vector<double> &cells, std::vector<double> &weights, SpatVector &v, SpatOptions &opt);


		SpatRaster replaceValues(std::vector<double> from, std::vector<double> to, long nl, bool setothers, double others, bool keepcats, SpatOptions &opt);
		SpatRaster reclassify(std::vector<std::vector<double>> rcl, unsigned openclosed, bool lowest, bool others, double othersValue, bool bylayer, bool brackets, bool keepcats, SpatOptions &opt);
		SpatRaster reclassify(std::vector<double> rcl, size_t nc, unsigned openclosed, bool lowest, bool others, double othersValue, bool bylayer, bool brackets, bool keepcats, SpatOptions &opt);
		//SpatRaster classify_layers(std::vector<std::vector<double>> groups, std::vector<double> id, SpatOptions &opt);
		//SpatRaster classify_layers(std::vector<double> groups, size_t nc, std::vector<double> id, SpatOptions &opt);

		SpatRaster intersect(SpatRaster &x, SpatOptions &opt);
		
		std::vector<double> readSample(size_t src, size_t srows, size_t scols);
		SpatRaster rotate(bool left, SpatOptions &opt);

		std::vector<size_t> sampleCells(double size, std::string method, bool replace, unsigned seed);
		SpatRaster sampleRegularRaster(double size, bool overview);
		SpatRaster sampleRowColRaster(size_t nr, size_t nc, bool warn);
		SpatRaster sampleRandomRaster(double size, bool replace, unsigned seed);
		std::vector<std::vector<double>> sampleRegularValues(double size, SpatOptions &opt);
		std::vector<double> sampleRowCol(size_t nr, size_t nc);
		std::vector<std::vector<double>> sampleRowColValues(size_t nr, size_t nc, SpatOptions &opt);
		
		std::vector<std::vector<double>> sampleRandomValues(double size, bool replace, unsigned seed);
		std::vector<std::vector<double>> sampleStratifiedCells(double size, bool each, bool replace, unsigned seed, SpatOptions &opt);

		SpatRaster sort(bool decreasing, bool order, SpatOptions &opt);

		SpatRaster scale(std::vector<double> center, bool docenter, std::vector<double> scale, bool doscale, SpatOptions &opt);
		SpatRaster scale_linear(double smin, double smax, SpatOptions &opt);

		SpatRaster similarity(std::vector<double> x, SpatOptions &opt);

		SpatRaster terrain(std::vector<std::string> v, unsigned neighbors, bool degrees, unsigned seed, SpatOptions &opt);

    // watershed2 ecor 20210317; EC 20210702 
		SpatRaster watershed2(int pp_offset,SpatOptions &opt); 
		SpatRaster pitfinder2(SpatOptions &opt); 
		SpatRaster NIDP2(SpatOptions &opt); 
		SpatRaster flowAccu2(SpatOptions &opt); 
		SpatRaster flowAccu2_weight(SpatRaster weight,SpatOptions &opt);
	// END watershed2 
		
		SpatRaster hillshade(SpatRaster aspect, std::vector<double> angle, std::vector<double> direction, bool normalize, SpatOptions &opt);

		SpatRaster selRange(SpatRaster x, int z, int recycleby, SpatOptions &opt);
		SpatRaster selectHighest(size_t n, bool low, SpatOptions &opt);

		SpatRaster shift(double x, double y, SpatOptions &opt);
		SpatRaster summary(std::string fun, bool narm, SpatOptions &opt);
		SpatRaster summary_numb(std::string fun, std::vector<double> add, bool narm, SpatOptions &opt);
		std::vector<std::vector<double>> where(std::string what, bool values, SpatOptions &opt);

		SpatRaster transpose(SpatOptions &opt);
		SpatRaster trig(std::string fun, SpatOptions &opt);
		SpatRaster trim1(double value, size_t padding, SpatOptions &opt);
		SpatRaster trim2(double value, size_t padding, SpatOptions &opt);
		std::vector<std::vector<double>> unique(bool bylayer, double digits, bool narm, SpatOptions &opt);
		SpatRaster project1(std::string newcrs, std::string method, SpatOptions &opt);
		SpatRaster project2(SpatRaster &x, std::string method, SpatOptions &opt);
		void project3(SpatRaster &out, std::string method, SpatOptions &opt);


#ifdef useGDAL
		bool getDSh(GDALDatasetH &rstDS, SpatRaster &out, std::string &filename, std::string &driver, double &naval, bool update, double background, SpatOptions &opt);
		bool getDShMEM(GDALDatasetH &rstDS, SpatRaster &out, double &naval, double background, SpatOptions &opt);
		
		bool open_gdal(GDALDatasetH &hDS, int src, bool update, SpatOptions &opt);
		bool create_gdalDS(GDALDatasetH &hDS, std::string filename, std::string driver, bool fill, double fillvalue, std::vector<bool> has_so, std::vector<double> scale, std::vector<double> offset, SpatOptions& opt);
		bool from_gdalMEM(GDALDatasetH hDS, bool set_geometry, bool get_values);

		bool as_gdalvrt(GDALDatasetH &hVRT, SpatOptions &opt);
		//bool as_gdalmem(GDALDatasetH &hVRT);
		
#endif

		SpatRaster to_memory_copy(SpatOptions &opt);
		bool to_memory(SpatOptions &opt);

		SpatRaster weighted_mean(SpatRaster w, bool narm, SpatOptions &opt);
		SpatRaster weighted_mean(std::vector<double> w, bool narm, SpatOptions &opt);

		SpatRaster warper(SpatRaster x, std::string crs, std::string method, bool mask, bool align, bool resample, SpatOptions &opt);
		SpatRaster warper_by_util(SpatRaster x, std::string crs, std::string method, bool mask, bool align, bool resample, SpatOptions &opt);
		
		SpatRaster resample(SpatRaster x, std::string method, bool mask, bool agg, SpatOptions &opt);
		
		SpatRaster applyGCP(std::vector<double> fx, std::vector<double> fy, std::vector<double> tx, std::vector<double> ty, SpatOptions &opt);

		SpatDataFrame zonal(SpatRaster z, SpatRaster g, std::string fun, bool narm, SpatOptions &opt);
		SpatDataFrame zonal_weighted(SpatRaster x, SpatRaster w,  bool narm, SpatOptions &opt);

		SpatDataFrame zonal_poly(SpatVector x, std::string fun, bool weights, bool exact, bool touches, bool small, bool narm, SpatOptions &opt);
		SpatDataFrame zonal_poly_weighted(SpatVector x, SpatRaster w, bool weights, bool exact, bool touches, bool small, bool narm, SpatOptions &opt);
		std::vector<std::vector<double>> zonal_poly_table(SpatVector x, bool weights, bool exact, bool touches, bool small, bool narm, SpatOptions &opt);

//		SpatDataFrame zonal_old(SpatRaster x, std::string fun, bool narm, SpatOptions &opt);
		SpatRaster rgb2col(size_t r,  size_t g, size_t b, SpatOptions &opt);
		SpatRaster rgb2hsx(std::string type, SpatOptions &opt);	
		SpatRaster hsx2rgb(SpatOptions &opt);	

		SpatRaster viewshed(std::vector<double> obs, std::vector<double> vals, double curvcoef, int mode, double maxdist, int heightmode, SpatOptions &opt);
		SpatRaster sieveFilter(int threshold, int connections, SpatOptions &opt);	
		
//		SpatRaster panSharpen(SpatRaster pan, SpatOptions &opt);	
};

