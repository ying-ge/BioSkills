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

#include "spatVector.h"
#include "geodesic.h"
#include "recycle.h"
#include "geosphere.h"

#ifndef M_PI
#define M_PI  (3.1415926535897932384626433)
#endif

#ifndef M_2PI
#define M_2PI (M_PI * 2.0)
#endif

#ifndef M_hPI
#define M_hPI (M_PI / 2.0)
#endif

#ifndef WGS84_a
#define WGS84_a 6378137.0
#endif

#ifndef WGS84_f
#define WGS84_f 1/298.257223563
#endif


inline void normLon(double &lon) {
	lon = fmod(lon + 180, 360.) - 180;
}

inline void normLonRad(double &lon) {
	lon = fmod(lon + M_PI, M_2PI) - M_PI;
}




inline double get_sign(const double &x) {
	return (x > 0.0) ? 1.0 : (x < 0.0) ? -1.0 : 0;
}


double distance_geo(double lon1, double lat1, double lon2, double lat2) {
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, WGS84_a, WGS84_f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
  	return s12;
}

inline double distance_cos_r(double lon1, double lat1, double lon2, double lat2, double r = 6378137.) {
	return r * acos((sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1-lon2)));
}

inline double distance_hav_r(double lon1, double lat1, double lon2, double lat2, const double r = 6378137.) {
	double dLat = lat2-lat1;
	double dLon = lon2-lon1;
	double a = sin(dLat/2.) * sin(dLat/2.) + cos(lat1) * cos(lat2) * sin(dLon/2.) * sin(dLon/2.);
	return 2. * atan2(sqrt(a), sqrt(1. - a)) * 6378137.0;
}


/*
double distance_cosdeg(double lon1, double lat1, double lon2, double lat2, double r = 6378137.) {
	deg2rad(lon1);
	deg2rad(lon2);
	deg2rad(lat1);
	deg2rad(lat2);
	return r * acos((sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1-lon2)));
}
*/


void dest_geo(double slon, double slat, double sazi, double dist, double &dlon, double &dlat, double &dazi) {
	struct geod_geodesic g;
	geod_init(&g, WGS84_a, WGS84_f);
	geod_direct(&g, slat, slon, sazi, dist, &dlat, &dlon, &dazi);
}


double direction_geo(double lon1, double lat1, double lon2, double lat2) {
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, WGS84_a, WGS84_f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
	return( azi1) ;
}


double direction_cos(double& lon1, double& lat1, double& lon2, double& lat2) {
	if ((lon1 == lon2) && (lat1 == lat2)) return 0; // NAN?
	double dLon = lon2 - lon1;
    double y = sin(dLon)  * cos(lat2); 
    double x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon); 
    double azm = atan2(y, x);
	azm = fmod(azm+M_PI, M_PI);
	return azm > M_PI ? -(M_PI - azm) : azm;
}



double dist2track_geo(double lon1, double lat1, double lon2, double lat2, double plon, double plat, bool sign, double r=6378137) {
	double a = 1;
	double f = 0;
	struct geod_geodesic geod;
	geod_init(&geod, a, f);
	double d, b2, b3, azi;
	geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &b2, &azi);
	geod_inverse(&geod, lat1, lon1, plat, plon, &d, &b3, &azi);
	double toRad = M_PI / 180.;
	b2 *= toRad;
	b3 *= toRad;
	double xtr = asin(sin(b3-b2) * sin(d)) * r;
	return sign ? xtr : fabs(xtr);
}


inline double dist2track_cos(double lon1, double lat1, double lon2, double lat2, double plon, double plat, bool sign, double r=6378137) {
	double b2 = direction_cos(lon1, lat1, lon2, lat2);
	double b3 = direction_cos(lon1, lat1, plon, plat);
	double d = distance_cos_r(lon1, lat1, plon, plat, 1);
	double xtr = asin(sin(b3-b2) * sin(d)) * r;
	return sign ? xtr : fabs(xtr);
}


inline double dist2track_hav(double lon1, double lat1, double lon2, double lat2, double plon, double plat, bool sign, double r=6378137) {
	double b2 = direction_cos(lon1, lat1, lon2, lat2);
	double b3 = direction_cos(lon1, lat1, plon, plat);
	double d = distance_hav_r(lon1, lat1, plon, plat, 1);
	double xtr = asin(sin(b3-b2) * sin(d)) * r;
	return sign ? xtr : fabs(xtr);
}



double alongTrackDistance_geo(double lon1, double lat1, double lon2, double lat2, double plon, double plat, double r=6378137) {
	double a = 1;
	double f = 0;
	struct geod_geodesic geod;
	geod_init(&geod, a, f);
	double d, b2, b3, azi;
	geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &b2, &azi);
	geod_inverse(&geod, lat1, lon1, plat, plon, &d, &b3, &azi);
	deg2rad(b2);
	deg2rad(b3);
	double xtr = asin(sin(b3-b2) * sin(d));
	double bsign = get_sign(cos(b2-b3));
	return fabs(bsign * acos(cos(d) / cos(xtr)) * r);
}


double alongTrackDistance_cos(double lon1, double lat1, double lon2, double lat2, double plon, double plat, double r=6378137) {

	double tc = direction_cos(lon1, lat1, lon2, lat2); // * toRad
	double tcp = direction_cos(lon1, lat1, plon, plat); // * toRad
    double dp = distance_cos_r(lon1, lat1, plon, plat, 1);
	double xtr = asin(sin(tcp-tc) * sin(dp));

// +1/-1 for ahead/behind [lat1,lon1]
	double bearing = get_sign(cos(tc - tcp));
	double angle = cos(dp) / cos(xtr);

// Fixing limits for the angle between [-1, 1] to avoid NaNs from acos
	angle = angle > 1 ? 1 : angle < -1 ? -1 : angle;
	double dist = bearing * acos(angle) * r;
	return fabs(dist);
}



double alongTrackDistance_hav(double lon1, double lat1, double lon2, double lat2, double plon, double plat, double r=6378137) {

	double tc = direction_cos(lon1, lat1, lon2, lat2); // * toRad
	double tcp = direction_cos(lon1, lat1, plon, plat); // * toRad
    double dp = distance_hav_r(lon1, lat1, plon, plat, 1);
	double xtr = asin(sin(tcp-tc) * sin(dp));

// +1/-1 for ahead/behind [lat1,lon1]
	double bearing = get_sign(cos(tc - tcp));
	double angle = cos(dp) / cos(xtr);

// Fixing limits for the angle between [-1, 1] to avoid NaNs from acos
	angle = angle > 1 ? 1 : angle < -1 ? -1 : angle;
	double dist = bearing * acos(angle) * r;
	return fabs(dist);
}


// the alongTrackDistance is the length of the path along the great circle to the point of intersection
// there are two, depending on which node you start
// we want to use the min, but the max needs to be < segment length

double dist2segment_geo(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double notused) {

	double seglength = distance_geo(lon1, lat1, lon2, lat2);
	double trackdist1 = alongTrackDistance_geo(lon1, lat1, lon2, lat2, plon, plat);
	double trackdist2 = alongTrackDistance_geo(lon2, lat2, lon1, lat1, plon, plat);
	if ((trackdist1 >= seglength) || (trackdist2 >= seglength)) {
		double d1 = distance_geo(lon1, lat1, plon, plat);
		double d2 = distance_geo(lon2, lat2, plon, plat);
		return d1 < d2 ? d1 : d2;
	}
	return dist2track_geo(lon1, lat1, lon2, lat2, plon, plat, false);
}


double dist2segment_cos(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double r) {
	double seglength = distance_cos_r(lon1, lat1, lon2, lat2, r);
	double trackdist1 = alongTrackDistance_cos(lon1, lat1, lon2, lat2, plon, plat, r);
	double trackdist2 = alongTrackDistance_cos(lon2, lat2, lon1, lat1, plon, plat, r);
	if ((trackdist1 >= seglength) || (trackdist2 >= seglength)) {
		double d1 = distance_cos_r(lon1, lat1, plon, plat, r);
		double d2 = distance_cos_r(lon2, lat2, plon, plat, r);
		return d1 < d2 ? d1 : d2;
	}
	return dist2track_cos(lon1, lat1, lon2, lat2, plon, plat, false, r);
}

double dist2segment_hav(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double r) {
	double seglength = distance_hav_r(lon1, lat1, lon2, lat2, r);
	double trackdist1 = alongTrackDistance_hav(lon1, lat1, lon2, lat2, plon, plat, r);
	double trackdist2 = alongTrackDistance_hav(lon2, lat2, lon1, lat1, plon, plat, r);
	if ((trackdist1 >= seglength) || (trackdist2 >= seglength)) {
		double d1 = distance_hav_r(lon1, lat1, plon, plat, r);
		double d2 = distance_hav_r(lon2, lat2, plon, plat, r);
		return d1 < d2 ? d1 : d2;
	}
	return dist2track_hav(lon1, lat1, lon2, lat2, plon, plat, false, r);
}


// [[Rcpp::export]]
double dist2segmentPoint_geo(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double &ilon, double &ilat) {

	double seglength = distance_geo(lon1, lat1, lon2, lat2);
	double trackdist1 = alongTrackDistance_geo(lon1, lat1, lon2, lat2, plon, plat);
	double trackdist2 = alongTrackDistance_geo(lon2, lat2, lon1, lat1, plon, plat);
	if ((trackdist1 >= seglength) || (trackdist2 >= seglength)) {
		double d1 = distance_geo(lon1, lat1, plon, plat);
		double d2 = distance_geo(lat2, lat2, plon, plat);
		if (d1 < d2) {
			ilon = lon1;
			ilat = lat1;
			return d1;
		} else {
			ilon = lon2;
			ilat = lat2;
			return d2;
		}
	}
	double azi;
	double crossd = dist2track_geo(lon1, lat1, lon2, lat2, plon, plat, false);
	if (trackdist1 < trackdist2) {
		double bear = direction_geo(lon1, lat1, lon2, lat2);
		dest_geo(lon1, lat1, bear, trackdist1, ilon, ilat, azi);
	} else {
		double bear = direction_geo(lon2, lat2, lon1, lat1);
		dest_geo(lon2, lat2, bear, trackdist2, ilon, ilat, azi);
	}
	return crossd;
}



// [[Rcpp::export(name = "intermediate")]]
std::vector<std::vector<double>> intermediate(double lon1, double lat1, double lon2, double lat2, int n, double distance) {
	double a = 6378137.0;
	double f = 1/298.257223563;
	struct geod_geodesic geod;
	geod_init(&geod, a, f);
	double d, azi1, azi2;

	std::vector<std::vector<double>> out(2);
	if (n <= 0) {
		if (distance <= 0) {
			out[0] = {lon1, lon2};
			out[1] = {lon1, lon2};
			return out;
		} else {
			geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &azi1, &azi2);
			n = std::round(d / distance);
			if (n < 2) {
				out[0] = {lon1, lon2};
				out[1] = {lon1, lon2};
				return out;
			}
			distance = d / n;
		}
	} else if (n == 1) {
		out[0] = {lon1, lon2};
		out[1] = {lon1, lon2};
		return out;
	} else {
		geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &azi1, &azi2);
		//distance = d / n;
	}
	out[0].resize(n+1);
	out[1].resize(n+1);
	out[0][0] = lon1;
	out[1][0] = lat1;
	for (int i=1; i<n; i++) {
		geod_direct(&geod, lat1, lon1, azi1, distance*i, &out[1][i], &out[0][i], &azi2);
	}
	out[0][n] = lon2;
	out[1][n] = lat2;
	return out;
}




std::vector<bool> antipodal(std::vector<double> lon1, std::vector<double> lat1, std::vector<double> lon2, std::vector<double> lat2, double tol=1e-9) {
	recycle(lon1, lon2);
	recycle(lat1, lat2);
	std::vector<bool> out;
	out.reserve(lon1.size());
	double Pi180 = M_PI / 180.;
	for (size_t i=0; i<lon1.size(); i++){
		normLon(lon1[i]);
		normLon(lon2[i]);
		double diflon = fabs(lon1[i] - lon2[i]);
		double diflat = fabs(lat1[i] + lat2[i]);
		out.push_back(
			(diflat < tol) && ((cos(lat2[i] * Pi180) * fabs(fmod(diflon, 360.) - 180)) < tol)
		);
	}
	return out;
}


void antipodes(std::vector<double> &lon, std::vector<double> &lat) {
	size_t n=lon.size();
	for (size_t i=0; i<n; i++) {
		lon[i] = lon[i] + 180;
		normLon(lon[i]);
		lat[i] = -lat[i];
	}
}

