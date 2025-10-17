#include <Rcpp.h>
#include "spatRasterMultiple.h"
#include "string_utils.h"
#include "math_utils.h"

#include "sort.h"

#include "gdal_priv.h"
#include "gdalio.h"
#include "ogr_spatialref.h"

//#define GEOS_USE_ONLY_R_API
#include <geos_c.h>

#if defined(HAVE_TBB) && !defined(__APPLE__)
#define USE_TBB
#endif


#if GDAL_VERSION_MAJOR >= 3
#include "proj.h"
#define projh
#if PROJ_VERSION_MAJOR >=6
# define PROJ_6
#endif
#if PROJ_VERSION_MAJOR > 7
# define PROJ_71
#else
# if PROJ_VERSION_MAJOR == 7
#  if PROJ_VERSION_MINOR >= 1
#   define PROJ_71
#  endif
# endif
#endif
#else
#if PROJ_VERSION_MAJOR >=8
#include "proj.h"
#else
#define ACCEPT_USE_OF_DEPRECATED_PROJ_API_H
#include <proj_api.h>
#endif
#endif



// [[Rcpp::export(name = ".have_TBB")]]
bool have_TBB() {
	#ifdef HAVE_TBB
		return true;
	#else 
		return false;
	#endif 
}


//from sf
#ifdef projh
// [[Rcpp::export]]
std::string proj_version() {
	std::stringstream buffer;
	buffer << PROJ_VERSION_MAJOR << "." << PROJ_VERSION_MINOR << "." << PROJ_VERSION_PATCH;
	return buffer.str();
}

#else

std::string proj_version() {
	int v = PJ_VERSION;
	std::stringstream buffer;
	buffer << v / 100 << "." << (v / 10) % 10 << "." << v % 10;
	return buffer.str();
}

#endif



// [[Rcpp::export]]
std::vector<unsigned char> hex2rgb(std::string s) {
	unsigned char r, g, b;
	s = s.erase(0,1); // remove the "#"
	sscanf(s.c_str(), "%02hhx%02hhx%02hhx", &r, &g, &b);
	std::vector<unsigned char> x = {r, g, b};
	return x;
}

/*
// [[Rcpp::export]]
std::vector<size_t> terra_order(std::vector<double> v) {
	return sort_order_a(v);
}

// [[Rcpp::export]]
std::vector<double> terra_permute(std::vector<double> v, std::vector<size_t> p) {
	permute(v, p);
	return v;
}
*/


// [[Rcpp::export]]
std::string rgb2hex(std::vector<unsigned char> x) {
	std::stringstream ss;
	ss << "#" << std::hex << std::setw(6) << (x[0] << 16 | x[1] << 8 | x[2] );
	std::string s = ss.str();
	//std::transform(s.begin(), s.end(), s.begin(), ::toupper);
	str_replace_all(s, " ", "0");
	return s;
}


// [[Rcpp::export(name = ".sameSRS")]]
bool sameSRS(std::string x, std::string y) {
	std::string msg;
	SpatSRS srs;
	if (!srs.set(x, msg)) return false;
	return srs.is_same(y, false);
}


// [[Rcpp::export(name = ".SRSinfo")]]
std::vector<std::string> getCRSname(std::string s) {
	OGRSpatialReference x;
	OGRErr erro = x.SetFromUserInput(s.c_str());
	if (erro != OGRERR_NONE) {
		return {"unknown", "", "", "", ""};
	}
	std::string node;
	if (x.IsGeographic()) {
		node = "geogcs";
	} else {
		node = "projcs";
	}

	const char *value;
	std::string name = "";
	value = x.GetAttrValue(node.c_str());
	if (value != NULL) {
		name = value;
	}
	std::string aname = "";
	value = x.GetAuthorityName(node.c_str());
	if (value != NULL) {
		aname = value;
	}

	std::string acode = "";
	value = x.GetAuthorityCode(node.c_str());
	if (value != NULL) {
		acode = value;
	}

	double west, south, east, north;
	west = -10000;
	east = -10000;
	south = -10000;
	north = -10000;

	std::string aoi="", box="";
	#if GDAL_VERSION_MAJOR >= 3
	if (x.GetAreaOfUse(&west, &south, &east, &north, &value)) {
		if (value != NULL) {
			if (west > -1000) {
				aoi	= value;
				box = std::to_string(west) + ", " + std::to_string(east) + ", " + std::to_string(south) + ", " + std::to_string(north);
			}
		}
	}
	#endif
	return {name, aname, acode, aoi, box};
}

// [[Rcpp::export(name = ".getLinearUnits")]]
double getLinearUnits(std::string s) {
	std::string msg;
	SpatSRS srs;
	if (!srs.set(s, msg)) return NAN;
	return srs.to_meter();
}




// [[Rcpp::export(name = ".geotransform")]]
std::vector<double> geotransform(std::string fname) {
	std::vector<double> out;
    GDALDataset *poDataset = static_cast<GDALDataset*>(GDALOpenEx( fname.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY, NULL, NULL, NULL ));

    if( poDataset == NULL )  {
		Rcpp::Rcout << "cannot read from: " + fname << std::endl;
		return out;
	}

	double gt[6];
	if( poDataset->GetGeoTransform( gt ) != CE_None ) {
		Rcpp::Rcout << "bad geotransform" << std::endl;
	}
	out = std::vector<double>(std::begin(gt), std::end(gt));
	GDALClose( (GDALDatasetH) poDataset );

	return out;
}

// [[Rcpp::export(name = ".gdal_setconfig")]]
void gdal_setconfig(std::string option, std::string value) {
	if (value.empty()) {
		CPLSetConfigOption(option.c_str(), NULL);
	} else {
		CPLSetConfigOption(option.c_str(), value.c_str());
	}
}

// [[Rcpp::export(name = ".gdal_getconfig")]]
std::string gdal_getconfig(std::string option) {
	const char * value = CPLGetConfigOption(option.c_str(), NULL);
	std::string out = "";
	if (value != NULL) {
		out = value;
	}
	return out;
}


// [[Rcpp::export(name = ".gdalinfo")]]
std::string ginfo(std::string filename, std::vector<std::string> options, std::vector<std::string> oo) {
	return gdalinfo(filename, options, oo);
}

// [[Rcpp::export(name = ".gdalmdinfo")]]
std::string gmdinfo(std::string filename, std::vector<std::string> options) {
	return gdalMDinfo(filename, options);
}


// [[Rcpp::export(name = ".sdinfo")]]
std::vector<std::vector<std::string>> sd_info(std::string filename) {
	std::vector<std::vector<std::string>> sd = sdinfo(filename);
	return sd;
}

// [[Rcpp::export(name = ".gdal_version")]]
std::string gdal_version() {
	const char* what = "RELEASE_NAME";
	const char* x = GDALVersionInfo(what);
	std::string s = (std::string) x;
	return s;
}

// [[Rcpp::export(name = ".geos_version")]]
std::string geos_version(bool runtime = false, bool capi = false) {
	if (runtime)
		return GEOSversion();
	else {
		if (capi)
			return GEOS_CAPI_VERSION;
		else
			return GEOS_VERSION;
	}
}

// [[Rcpp::export(name = ".metadata")]]
std::vector<std::string> metatdata(std::string filename) {
	std::vector<std::string> m = get_metadata(filename);
	return m;
}

// [[Rcpp::export(name = ".sdsmetadata")]]
std::vector<std::string> sdsmetatdata(std::string filename) {
	std::vector<std::string> m = get_metadata_sds(filename);
	return m;
}




// [[Rcpp::export(name = ".parsedsdsmetadata")]]
std::vector<std::vector<std::string>> sdsmetatdataparsed(std::string filename) {
	std::vector<std::string> m = sdsmetatdata(filename);
	std::vector<std::vector<std::string>> s = parse_metadata_sds(m);
	return s;
}


// [[Rcpp::export(name = ".gdaldrivers")]]
std::vector<std::vector<std::string>> gdal_drivers() {
	size_t n = GetGDALDriverManager()->GetDriverCount();
	std::vector<std::vector<std::string>> s(6, std::vector<std::string>(n));
    GDALDriver *poDriver;
    char **papszMetadata;
	for (size_t i=0; i<n; i++) {
	    poDriver = GetGDALDriverManager()->GetDriver((int)i);
		const char* ss = poDriver->GetDescription();
		if (ss != NULL ) s[0][i] = ss;
		ss = poDriver->GetMetadataItem( GDAL_DMD_LONGNAME );
		if (ss != NULL ) s[5][i] = ss;

		papszMetadata = poDriver->GetMetadata();
		bool rst = CSLFetchBoolean( papszMetadata, GDAL_DCAP_RASTER, FALSE);
		bool vct = CSLFetchBoolean( papszMetadata, GDAL_DCAP_VECTOR, FALSE);
		s[1][i] = std::to_string(rst);
		s[2][i] = std::to_string(vct);
		bool create = CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE);
		bool copy = CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATECOPY, FALSE);
		s[3][i] = std::to_string(create + copy);
		bool vsi = CSLFetchBoolean( papszMetadata, GDAL_DCAP_VIRTUALIO, FALSE);
		s[4][i] = std::to_string(vsi);
	}
	return s;
}


template <typename... Args>
inline void warningNoCall(const char* fmt, Args&&... args ) {
    Rf_warningcall(R_NilValue, "%s", tfm::format(fmt, std::forward<Args>(args)... ).c_str());
}

template <typename... Args>
inline void NORET stopNoCall(const char* fmt, Args&&... args) {
    throw Rcpp::exception(tfm::format(fmt, std::forward<Args>(args)... ).c_str(), false);
}

static void __err_warning(CPLErr eErrClass, int err_no, const char *msg) {
	switch ( eErrClass ) {
        case 0:
            break;
        case 1:
        case 2:
            warningNoCall("%s (GDAL %d)", msg, err_no);
            break;
        case 3:
            warningNoCall("%s (GDAL error %d)", msg, err_no);
            break;
        case 4:
            stopNoCall("%s (GDAL unrecoverable error %d)", msg, err_no);
            break;
        default:
            warningNoCall("%s (GDAL error class %d, #%d)", msg, eErrClass, err_no);
            break;
    }
    return;
}

static void __err_error(CPLErr eErrClass, int err_no, const char *msg) {
	switch ( eErrClass ) {
        case 0:
        case 1:
        case 2:
            break;
        case 3:
            warningNoCall("%s (GDAL error %d)", msg, err_no);
            break;
        case 4:
            stopNoCall("%s (GDAL unrecoverable error %d)", msg, err_no);
            break;
        default:
            stopNoCall("%s (GDAL unrecoverable error %d)", msg, err_no);
            break;
    }
    return;
}


static void __err_fatal(CPLErr eErrClass, int err_no, const char *msg) {
	switch ( eErrClass ) {
        case 0:
        case 1:
        case 2:
        case 3:
            break;
        case 4:
            stopNoCall("%s (GDAL unrecoverable error %d)", msg, err_no);
            break;
        default:
            break;
    }
    return;
}

static void __err_none(CPLErr eErrClass, int err_no, const char *msg) {
    return;
}



// [[Rcpp::export(name = ".set_gdal_warnings")]]
void set_gdal_warnings(int level) {
	if (level==4) {
		CPLSetErrorHandler((CPLErrorHandler)__err_none);
	} else if (level==1) {
		CPLSetErrorHandler((CPLErrorHandler)__err_warning);
	} else if (level==2) {
		CPLSetErrorHandler((CPLErrorHandler)__err_error);
	} else {
		CPLSetErrorHandler((CPLErrorHandler)__err_fatal);
	}
}

#include "common.h"

std::mt19937 my_rgen;

// [[Rcpp::export(name = ".seedinit")]]
void seed_init(uint32_t seed_val) {
  my_rgen.seed(seed_val);
}

/*
[[Rcpp::export(name = ".rnumb")]]
int rnumb() {
    std::uniform_int_distribution<std::mt19937::result_type> rand_nr(0, 9); 
	return  rand_nr(my_rgen);
}
*/

// [[Rcpp::export(name = ".gdalinit")]]
void gdal_init(std::string projpath, std::string datapath) {
	set_gdal_warnings(2);
    GDALAllRegister();
    OGRRegisterAll();
	CPLSetConfigOption("GDAL_MAX_BAND_COUNT", "9999999");
	CPLSetConfigOption("OGR_CT_FORCE_TRADITIONAL_GIS_ORDER", "YES");
	CPLSetConfigOption("GDAL_DATA", datapath.c_str());
	CPLSetConfigOption("CPL_VSIL_USE_TEMP_FILE_FOR_RANDOM_WRITE", "YES");
	//GDAL_NETCDF_IGNORE_XY_AXIS_NAME_CHECKS

	//GDALregistred = true;
#if GDAL_VERSION_MAJOR >= 3
 #ifdef PROJ_6
	if (!projpath.empty()) {
		const char *cp = projpath.c_str();
		proj_context_set_search_paths(PJ_DEFAULT_CTX, 1, &cp);
	}
 #endif
#endif
#ifdef PROJ_71
	#ifndef __EMSCRIPTEN__
		proj_context_set_enable_network(PJ_DEFAULT_CTX, 1);
	#endif
#endif
}

// [[Rcpp::export(name = ".precRank")]]
std::vector<double> percRank(std::vector<double> x, std::vector<double> y, double minc, double maxc, int tail) {

	std::vector<double> out;
	out.reserve(y.size());
	size_t nx = x.size();
	for (size_t i=0; i<y.size(); i++) {
		if (std::isnan(y[i]) ) {
			out.push_back( NAN );
		} else if ((y[i] < minc) || (y[i] > maxc )) {
			out.push_back( 0 );
		} else {
			size_t b = 0;
			size_t t = 0;
			for (size_t j=0; j<x.size(); j++) {
				if (y[i] > x[j]) {
					b++;
				} else if (y[i] == x[j]) {
					t++;
				} else {
				// y is sorted, so we need not continue
					break;
				}
			}
			double z = ((double)b + 0.5 * (double)t) / (double) nx;
			if (tail == 1) { // both
				if (z > 0.5) {
					z = 2 * (1 - z);
				} else {
					z = 2 * z;
				}
			} else if (tail == 2) { // high
				if (z < 0.5) {
					z = 1;
				} else {
					z = 2 * (1 - z);
				}
			} else { // if (tail == 3) { // low
				if (z > 0.5) {
					z = 1;
				} else {
					z = 2 * z;
				}
			}
			out.push_back(z);
		}
	}
	return(out);
}


// [[Rcpp::export(name = ".clearVSIcache")]]
void clearVSIcache(bool vsi) {
	//if (vsi)
	VSICurlClearCache();
}

// [[Rcpp::export(name = ".setGDALCacheSizeMB")]]
void setGDALCacheSizeMB(double x, bool vsi) {
	if (vsi) {
		int64_t v = int64_t(x * 1024 * 1024);
		CPLSetConfigOption("CPL_VSIL_CURL_CACHE_SIZE", std::to_string(v).c_str());
	} else {
		GDALSetCacheMax64(static_cast<int64_t>(x) * 1024 * 1024);
	}

}

// [[Rcpp::export(name = ".getGDALCacheSizeMB")]]
double getGDALCacheSizeMB(bool vsi) {
	if (vsi) {
		std::string out = gdal_getconfig("CPLGetConfigOption");
		Rcpp::Rcout << out << std::endl;
		if (out == "") return NAN;
		double v = -1; 
		try {
			v = stod(out) / (1024 * 1024);
		} catch(...){
			return(NAN);
		}
		return(v);
		
	} else {
		return static_cast<double>(GDALGetCacheMax64() / 1024 / 1024);
	}
}


// convert NULL-terminated array of strings to std::vector<std::string>
std::vector<std::string> charpp2vect(char **cp) {
	std::vector<std::string> out;
	if (cp == NULL) return out;
	size_t i=0;
	while (cp[i] != NULL) {
		out.push_back(cp[i]);
		i++;
	}
	return out;
}


// [[Rcpp::export(name = ".get_proj_search_paths")]]
std::vector<std::string> get_proj_search_paths() {
	std::vector<std::string> out;
#if GDAL_VERSION_NUM >= 3000300
	char **cp = OSRGetPROJSearchPaths();
	out = charpp2vect(cp);
	CSLDestroy(cp);
#else
	out.push_back("error: GDAL >= 3.0.3 required");
#endif
	return out;
}




// [[Rcpp::export(name = ".set_proj_search_paths")]]
bool set_proj_search_paths(std::vector<std::string> paths) {
	if (paths.empty()) {
		return false;
	}
#if GDAL_VERSION_NUM >= 3000000
	std::vector <char *> cpaths(paths.size()+1);
	for (size_t i = 0; i < paths.size(); i++) {
		cpaths[i] = (char *) (paths[i].c_str());
	}
	cpaths[cpaths.size()-1] = NULL;
	OSRSetPROJSearchPaths(cpaths.data());
	return true;
#else
	return false;
#endif
}


// [[Rcpp::export(name = ".PROJ_network")]]
std::string PROJ_network(bool enable, std::string url) {
	std::string s = "";
#ifdef PROJ_71
	if (enable) {
		proj_context_set_enable_network(PJ_DEFAULT_CTX, 1);
		if (url.size() > 5) {
			proj_context_set_url_endpoint(PJ_DEFAULT_CTX, url.c_str());
		}
		s = proj_context_get_url_endpoint(PJ_DEFAULT_CTX);
	} else { // disable:
		proj_context_set_enable_network(PJ_DEFAULT_CTX, 0);
	}
#endif
	return s;
}



// [[Rcpp::export(name = ".removeDriver")]]
void removeDriver(std::vector<std::string> d) {
	if ((d.size() == 0) || ((d.size() == 1) && (d[0] == ""))) {
		GDALAllRegister();	
	} else {
		for (size_t i=0; i<d.size(); i++) {
			GDALDriverH hDrv = GDALGetDriverByName(d[i].c_str());
			if (hDrv == NULL) {
				Rcpp::warning(d[i] + " is not a known driver\n");
			} else {
				GDALDeregisterDriver(hDrv);
			}
		}
	}
}


// [[Rcpp::export(name = ".pearson")]]
double pearson_cor(std::vector<double> x, std::vector<double> y, bool narm) {
 
	if (narm) {
		size_t n = x.size()-1;
		for (long i=n; i >= 0; i--) {
			if (std::isnan(x[i]) || std::isnan(y[i])) {
				x.erase(x.begin()+i);
				y.erase(y.begin()+i);
			}
		}	
		if (x.size() < 2) {
			return(NAN);
		}
	}
	size_t n = x.size();
	double xbar = accumulate(x.begin(), x.end(), 0.0) / (double)n;
	double ybar = accumulate(y.begin(), y.end(), 0.0) / (double)n;
	double numer = 0;
	for (size_t i=0; i<n; i++) {
		numer += (x[i]-xbar) * (y[i]-ybar);
	}
	double s1 = 0;
	double s2 = 0;
	for (size_t i=0; i<n; i++) {
		s1 += std::pow(x[i]-xbar, 2);
		s2 += std::pow(y[i]-ybar, 2);
	}
	return numer / std::sqrt(s1 * s2);
}



// [[Rcpp::export(name = ".weighted_pearson")]]
double weighted_pearson_cor(std::vector<double> x, std::vector<double> y, std::vector<double> weights, bool narm=true) {
  
	if (narm) {
		size_t n = x.size()-1;
		for (long i=n; i >= 0; i--) {
			if (std::isnan(x[i]) || std::isnan(y[i])) {
				x.erase(x.begin()+i);
				y.erase(y.begin()+i);
				weights.erase(weights.begin()+i);
			}
		}	
		if (x.size() < 2) {
			return(NAN);
		}
	}
	size_t n = x.size();
	double sw = accumulate(weights.begin(), weights.end(), 0.0);
	for (double &d : weights) d /= sw;
	double sxw = 0;
	double syw = 0;
	for (size_t i=0; i<n; i++) {
		sxw += x[i] * weights[i];
		syw += y[i] * weights[i];
	}
	for (size_t i=0; i<n; i++) {
		x[i] -= sxw;
		y[i] -= syw;
	}
	double vx = 0;
	double vy = 0;
	double vxy = 0;
	for (size_t i=0; i<n; i++) {
		vx += weights[i] * x[i] * x[i];
		vy += weights[i] * y[i] * y[i];
		vxy += weights[i] * x[i] * y[i];
	}
	return  vxy / std::sqrt(vx * vy);
}



// [[Rcpp::export(name = ".unique_symmetric_rows")]]
Rcpp::IntegerMatrix uniqueSymmetricRows(std::vector<size_t> x, std::vector<size_t> y) {
	size_t n = x.size();
	for (size_t i=0; i<n; i++) {
		if (x[i] > y[i]) {
			size_t tmp = x[i];
			x[i] = y[i];
			y[i] = tmp;
		}
	}
	sort_unique_2d(x, y);
	Rcpp::IntegerMatrix mat(x.size(), 2); 
    std::copy(x.begin(), x.end(), mat.begin());  	
    std::copy(y.begin(), y.end(), mat.begin()+x.size());  	
	return mat;
//	x.insert(x.end(), y.begin(), y.end());
///	return(x);
}





/*

# include "vecmathse.h"
// [[Rcpp::export(name = ".stattest1")]]
double stattest1(std::vector<double> x, std::string fun, bool narm) {
	if (!haveseFun(fun)) {
		Rcpp::Rcout << fun + " is not available" << std::endl;
		return NAN;
	}
	std::function<double(std::vector<double>&, size_t, size_t)> f;
	if (!getseFun(f, fun, narm)) {
		Rcpp::Rcout << "Unknown function" << std::endl;
		return NAN;
	}
	return f(x, 0, x.size());
}


# include "vecmath.h"
// [[Rcpp::export(name = ".stattest2")]]
double stattest2(std::vector<double> x, std::string fun, bool narm) {
	if (!haveFun(fun)) {
		Rcpp::Rcout << fun + " is not available" << std::endl;
		return NAN;
	}
	std::function<double(std::vector<double>&, bool)> f = getFun(fun);
	return f(x, narm);
}

*/


Rcpp::List get_output(std::vector<std::string> names, std::vector<long> sizes) {
	Rcpp::List L = Rcpp::List::create(Rcpp::Named("name") = names, Rcpp::Named("size") = sizes);
	return(L);
}


#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 4

// [[Rcpp::export(name = ".arnames")]]
std::vector<std::string> arnames(std::string filename, bool filter) {
// FROM GDAL 3.11 (while not widely available).
// * Author:   Even Rouault <even.rouault at spatialys.com>
// * Copyright (c) 2019, Even Rouault <even.rouault at spatialys.com>

    std::vector<std::string> ret;
    auto poDataset = std::unique_ptr<GDALDataset>(GDALDataset::Open(filename.c_str(), GDAL_OF_MULTIDIM_RASTER ));

    if (!poDataset) {
		ret.push_back("not a good dataset");
        return ret;
    }

	std::shared_ptr<GDALGroup> poRootGroup = poDataset->GetRootGroup();
    if (!poRootGroup) {
		ret.push_back("no root group");
        return ret;
    }

    std::list<std::shared_ptr<GDALGroup>> stackGroups;
    stackGroups.push_back(nullptr);  // nullptr means this
    while (!stackGroups.empty()) {
        std::shared_ptr<GDALGroup> groupPtr = std::move(stackGroups.front());
        stackGroups.erase(stackGroups.begin());
        const GDALGroup *poCurGroup = groupPtr ? groupPtr.get() : poRootGroup.get();
        for (const std::string &arrayName :  poCurGroup->GetMDArrayNames(nullptr)) {
            std::string osFullName = poCurGroup->GetFullName();
            if (!osFullName.empty() && osFullName.back() != '/') 
                osFullName += '/';
            osFullName += arrayName;
            ret.push_back(std::move(osFullName));
        }
        auto insertionPoint = stackGroups.begin();
        for (const auto &osSubGroup : poCurGroup->GetGroupNames(nullptr)) {
            auto poSubGroup = poCurGroup->OpenGroup(osSubGroup);
            if (poSubGroup)
                stackGroups.insert(insertionPoint, std::move(poSubGroup));
        }
    }
	if (filter) {
		ret = ncdf_filternames(ret);
	}
    return ret;
}


// [[Rcpp::export(name = ".dimfo")]]
Rcpp::List dimfo(std::string filename, std::string array_name) {

	std::vector<std::string> names;
	std::vector<long> sizes;

    auto poDataset = std::unique_ptr<GDALDataset>(GDALDataset::Open(filename.c_str(), GDAL_OF_MULTIDIM_RASTER ));
    if (!poDataset) {
		names.push_back("cannot open as md: " + filename);
		sizes.push_back(-99);
		return get_output(names, sizes);
    }

	std::shared_ptr<GDALGroup> poRootGroup = poDataset->GetRootGroup();
    if (!poRootGroup) {
		names.push_back("no root group: " + filename);
		sizes.push_back(-99);
		return get_output(names, sizes);
    }

	auto poVar = poRootGroup->OpenMDArray(array_name.c_str());
	if (!poVar)   {
		names.push_back("cannot open: " + array_name);
		sizes.push_back(-99);
		return get_output(names, sizes);
	}

	if (names.size() == 0) {
		for ( const auto &poDim: poVar->GetDimensions() ) {
			names.push_back(static_cast<std::string>(poDim->GetName()));
			sizes.push_back(static_cast<long>(poDim->GetSize()));
		}
	}
	return get_output(names, sizes);
}

#else

Rcpp::List dimfo(std::string filename, std::string array_name) {
	return get_output({"not supported with GDAL < 3.4"}, {-99});
}

std::vector<std::string> arnames(std::string filename, bool filter) {
	std::vector<std::string> out(1, "not supported with GDAL < 3.4");
	return out;
}

#endif

