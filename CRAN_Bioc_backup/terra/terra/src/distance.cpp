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

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#include <vector>
//#include <math.h>
#include <cmath>
#include "geodesic.h"
#include "geosphere.h"
#include "recycle.h"
#include <functional>
#include <string>

// Convert degrees to radians
double toRad(double &deg) {
	return( deg * 0.0174532925199433 );
}

double toDeg(double &rad) {
	return( rad * 57.2957795130823 );
}




double distance_lonlat(const double &lon1, const double &lat1, const double &lon2, const double &lat2) {
	double a = 6378137.0;
	double f = 1/298.257223563;
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, a, f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
  	return s12;
}


// to circumvent overloading problem when assigning to std::function
double distLonlat(const double &lon1, const double &lat1, const double &lon2, const double &lat2) {
	double a = 6378137.0;
	double f = 1/298.257223563;
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, a, f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
  	return s12;
}


double distCosine(double lon1, double lat1, double lon2, double lat2) {
	const double r = 6378137;
	lon1 = toRad(lon1);
	lon2 = toRad(lon2);
	lat1 = toRad(lat1);
	lat2 = toRad(lat2);
	return r * acos((sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1-lon2)));
}

double distCosineRad(const double &lon1, const double &lat1, const double &lon2, const double &lat2) {
	const double r = 6378137;
	return r * acos((sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1-lon2)));
}


double distHaversine(double lon1, double lat1, double lon2, double lat2) {
	const double r = 6378137;
	lon1 = toRad(lon1);
	lon2 = toRad(lon2);
	lat1 = toRad(lat1);
	lat2 = toRad(lat2);

	double dLat = lat2-lat1;
	double dLon = lon2-lon1;
	double a = sin(dLat/2.) * sin(dLat/2.) + cos(lat1) * cos(lat2) * sin(dLon/2.) * sin(dLon/2.);
	return 2. * atan2(sqrt(a), sqrt(1.-a)) * r;
}

/*
double distHaversineRad(const double &lon1, const double &lat1, const double &lon2, const double &lat2) {
	const double r = 6378137;
	double dLat = lat2-lat1;
	double dLon = lon2-lon1;
	double a = sin(dLat/2.) * sin(dLat/2.) + cos(lat1) * cos(lat2) * sin(dLon/2.) * sin(dLon/2.);
	return 2. * atan2(sqrt(a), sqrt(1.-a)) * r;
}
*/



std::vector<double> distance_lonlat(std::vector<double> &lon1, std::vector<double> &lat1, std::vector<double> &lon2, std::vector<double> &lat2) {
	double a = 6378137.0;
	double f = 1/298.257223563;
	
	size_t n = std::max(std::max(std::max(lon1.size(), lat1.size()), lon2.size()), lat2.size());
    recycle(lon1, n);
    recycle(lon2, n);
    recycle(lat1, n);
    recycle(lat2, n);
	
	std::vector<double> r(n);
	double azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, a, f);
  	for (size_t i=0; i < n; i++) {
		geod_inverse(&g, lat1[i], lon1[i], lat2[i], lon2[i], &r[i], &azi1, &azi2);
	}
  	return r;
}


std::vector<double> distance_lonlat_vd(std::vector<double> &lon1, std::vector<double> &lat1, double lon2, double lat2) {
	std::vector<double> vlon2(lon1.size(), lon2);
	std::vector<double> vlat2(lat1.size(), lat2);
    return distance_lonlat(lon1, lat1, vlon2, vlat2);
}


std::vector<double> distance_lon(double &lon, std::vector<double> &lat) {
	double a = 6378137.0;
	double f = 1/298.257223563;
	
	size_t n = lat.size();
	std::vector<double> r(n);
	double azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, a, f);
  	for (size_t i=0; i < n; i++) {
		geod_inverse(&g, lat[i], 0, lat[i], lon, &r[i], &azi1, &azi2);
	}
  	return r;
}



double distance_plane(const double &x1, const double &y1, const double &x2, const double &y2) {
	return( sqrt(pow((x2-x1),2) + pow((y2-y1), 2)) );
}


std::vector<double> distance_plane(std::vector<double> &x1, std::vector<double> &y1, std::vector<double> &x2, std::vector<double> &y2) {
    recycle(x1, x2);
    recycle(y1, y2);
	std::vector<double> r (x1.size());
	size_t n = x1.size();
  	for (size_t i=0; i < n; i++) {
		r[i] = distance_plane(x1[i], y1[i], x2[i], y2[i]);
	}
  	return r;
}

std::vector<double> distance_plane_vd(std::vector<double> &x1, std::vector<double> &y1, double x2, double y2) {
	std::vector<double> vx2(x1.size(), x2);
	std::vector<double> vy2(y1.size(), y2);
  	return distance_plane(x1, y1, vx2, vy2);
}


/*
double distPlane(double x1, double y1, double x2, double y2) {
	return( sqrt(pow((x2-x1),2) + pow((y2-y1), 2)) );
}
*/



double direction_lonlat(double lon1, double lat1, double lon2, double lat2, bool degrees) {
	double a = 6378137.0;
	double f = 1/298.257223563;

	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, a, f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
	if (!degrees) {
		return(toRad(azi1));
	}
	return( azi1) ;
}

std::vector<double> direction_lonlat(std::vector<double> lon1, std::vector<double> lat1, std::vector<double> lon2, std::vector<double> lat2, bool degrees) {
	double a = 6378137.0;
	double f = 1/298.257223563;

// lonlat1 and lonlat2 should have the same length

	std::vector<double> azi1(lon1.size());
	double s12, azi2;
	struct geod_geodesic g;
	geod_init(&g, a, f);
	size_t n = lat1.size();
	if (degrees) {
		for (size_t i=0; i < n; i++) {
			geod_inverse(&g, lat1[i], lon1[i], lat2[i], lon2[i], &s12, &azi1[i], &azi2);
		}
	} else {
		for (size_t i=0; i < n; i++) {
			geod_inverse(&g, lat1[i], lon1[i], lat2[i], lon2[i], &s12, &azi1[i], &azi2);
			azi1[i] = toRad(azi1[i]);
		}
	}
  	return azi1;
}


void directionToNearest_lonlat(std::vector<double> &azi, std::vector<double> &lon1, std::vector<double> &lat1, std::vector<double> &lon2, std::vector<double> &lat2, bool& degrees, bool& from, const std::string &method) {

	if (method == "geo") {
		double a = 6378137.0;
		double f = 1/298.257223563;

		double azi1, azi2, s12, dist;
		size_t n = lon1.size();
		size_t m = lon2.size();

		azi.resize(n, NAN);

		struct geod_geodesic g;
		geod_init(&g, a, f);

		for (size_t i=0; i < n; i++) {
			if (std::isnan(lat1[i])) {
				azi[i] = NAN;
				continue;
			}
			geod_inverse(&g, lat1[i], lon1[i], lat2[0], lon2[0], &dist, &azi1, &azi2);
			size_t minj=0;
			azi[i] = azi1;
			for (size_t j=1; j<m; j++) {
				geod_inverse(&g, lat1[i], lon1[i], lat2[j], lon2[j], &s12, &azi1, &azi2);
				if (s12 < dist) {
					minj = j;
					dist = s12;
					azi[i] = azi1;
				}
			}
			if (from) {
				geod_inverse(&g, lat2[minj], lon2[minj], lat1[i], lon1[i], &s12, &azi1, &azi2);
				azi[i] = azi1;
			}
			if (!degrees) {
				azi[i] = toRad(azi[i]);
			}
		}
	} else {


//		double azi1, azi2, s12, dist;
		size_t n = lon1.size();
		size_t m = lon2.size();
		deg2rad(lat1);
		deg2rad(lon1);
		deg2rad(lat2);
		deg2rad(lon2);
		
		azi.resize(n, NAN);

		for (size_t i=0; i < n; i++) {
			if (std::isnan(lat1[i])) {
				azi[i] = NAN;
				continue;
			}
			azi[i] = direction_cos(lon1[i], lat1[i], lon2[0], lat2[0]);
			double dist = distance_cos(lon1[i], lat1[i], lon2[0], lat2[0]);
			size_t minj=0;
			for (size_t j=1; j<m; j++) {
				double d = distance_cos(lon1[i], lat1[i], lon2[j], lat2[j]);
				if (d < dist) {
					minj = j;
					dist = d;
					azi[i] = direction_cos(lon1[i], lat1[i], lon2[j], lat2[j]);
				}
			}
			if (from) {
				azi[i] = direction_cos(lon2[minj], lat2[minj], lon1[i], lat1[i]);
			}
			if (degrees) {
				azi[i] = toDeg(azi[i]);
			}
		}
	}	
}		
		



double direction_plane(double x1, double y1, double x2, double y2, bool degrees) {
	double a;
	double pi2 = M_PI * 2;
	a = fmod(atan2( x2 - x1, y2 - y1), pi2);
	a = a < 0 ? a + pi2 : a;
	return (degrees ? toDeg(a) : a);
}



std::vector<double> direction_plane(std::vector<double> x1, std::vector<double> y1, std::vector<double> x2, std::vector<double> y2, bool degrees) {
// xy1 and xy2 should have the same length
	std::vector<double> r (x1.size());
	//double a;
	size_t n = x1.size();
  	for (size_t i=0; i < n; i++) {
		r[i] = direction_plane(x1[i], y1[i], x2[i], y2[i], degrees);
	}
  	return r;
}



void directionToNearest_plane(std::vector<double> &r, const std::vector<double> &x1, const std::vector<double> &y1, const std::vector<double> &x2, const std::vector<double> &y2, bool& degrees, bool &from) {
	size_t n = x1.size();
	size_t m = x2.size();
	r.resize(n, NAN);
	double d, mind;
	size_t minj;
	for (size_t i = 0; i < n; i++) {
		r[i] = NAN;
		if (std::isnan(x1[i])) continue; // x2 must not be NAN
		mind = distance_plane(x1[i], y1[i], x2[0], y2[0]);
		minj = 0;
		for (size_t j = 1; j < m; j++) {
			d = distance_plane(x1[i], y1[i], x2[j], y2[j]);
			if (d < mind) {
				mind = d;
				minj = j;
			}
		}
		if (from) {
			r[i] = direction_plane(x2[minj], y2[minj], x1[i], y1[i], degrees);
		} else {
			r[i] = direction_plane(x1[i], y1[i], x2[minj], y2[minj], degrees);
		}
	}
}





std::vector<double> destpoint_lonlat(double longitude, double latitude, double  bearing, double distance) {
	double a = 6378137.0;
	double f = 1/298.257223563;

	struct geod_geodesic g;
	geod_init(&g, a, f);
	double lat2, lon2, azi2;
	geod_direct(&g, latitude, longitude, bearing, distance, &lat2, &lon2, &azi2);
	std::vector<double> out = { lon2, lat2, azi2 };
	return out;
}


std::vector<std::vector<double> > destpoint_lonlat(const std::vector<double> &longitude, const std::vector<double> &latitude, const std::vector<double> &bearing, const std::vector<double> &distance) {
	double a = 6378137.0;
	double f = 1/298.257223563;

	struct geod_geodesic g;
	geod_init(&g, a, f);
	size_t n = longitude.size();
	std::vector<std::vector<double> > out(3, std::vector<double>(n));
	double lat2, lon2, azi2;
	for (size_t i=0; i < n; i++) {
		geod_direct(&g, latitude[i], longitude[i], bearing[i], distance[i], &lat2, &lon2, &azi2);
		out[0][i] = lon2;
		out[1][i] = lat2;
		out[2][i] = azi2;
	}
	return out;
}


std::vector<std::vector<double> > destpoint_lonlat(const double &longitude, const double &latitude, const std::vector<double> &bearing, const double& distance, bool wrap) {
	double a = 6378137.0;
	double f = 1/298.257223563;

	struct geod_geodesic g;
	geod_init(&g, a, f);
	size_t n = bearing.size();
	std::vector<std::vector<double> > out(3, std::vector<double>(n));
	double lat2, lon2, azi2;
	if (wrap) {
		for (size_t i=0; i < n; i++) {
			geod_direct(&g, latitude, longitude, bearing[i], distance, &lat2, &lon2, &azi2);
			out[0][i] = lon2;
			out[1][i] = lat2;
			out[2][i] = azi2;
		}
	} else {
		for (size_t i=0; i < n; i++) {
			geod_direct(&g, latitude, 0, bearing[i], distance, &lat2, &lon2, &azi2);
			out[0][i] = lon2 + longitude;
			out[1][i] = lat2;
			out[2][i] = azi2;
		}
	}
	return out;
}



std::vector<double> destpoint_plane(double x, double y, double bearing, double distance) {
	bearing *= M_PI / 180;
	x += distance * sin(bearing);
	y += distance * cos(bearing);
	std::vector<double> out = {x, y};
	return(out);
}


std::vector<std::vector<double> > destpoint_plane(std::vector<double>  x, std::vector<double>  y, std::vector<double>  bearing, std::vector<double>  distance) {
	size_t n = x.size();
	std::vector<std::vector<double> > out;
	out.reserve(n);
	double xd, yd, b;
	for (size_t i=0; i < n; i++) {
		b = bearing[i] * M_PI / 180;
		xd = x[i] + distance[i] * sin(b);
		yd = y[i] + distance[i] * cos(b);
		out.push_back( {xd, yd });
	}
	return(out);
}



void distanceToNearest_lonlat(std::vector<double> &d, const std::vector<double> &lon1, const std::vector<double> &lat1, const std::vector<double> &lon2, const std::vector<double> &lat2, const double& adj_unit, std::string method) {

	int n = lon1.size();
	int m = lon2.size();

	if (method == "geo") {
		double a = 6378137.0;
		double f = 1/298.257223563;
		double azi1, azi2, s12;
		struct geod_geodesic g;
		geod_init(&g, a, f);

		for (int i=0; i < n; i++) {
			if (std::isnan(lat1[i])) {
				continue;
			}
			geod_inverse(&g, lat1[i], lon1[i], lat2[0], lon2[0], &d[i], &azi1, &azi2);
			for (int j=1; j<m; j++) {
				geod_inverse(&g, lat1[i], lon1[i], lat2[j], lon2[j], &s12, &azi1, &azi2);
				if (s12 < d[i]) {
					d[i] = s12;
				}
			}
			d[i] *= adj_unit;
		}
	} else {
		std::function<double(double, double, double, double)> dfun;
		if (method == "haversine") {
			dfun = distHaversine;
		} else {
			dfun = distCosine;			
		}
		
		for (int i=0; i < n; i++) {
			if (!std::isnan(lat1[i])) {
				d[i] = dfun(lon1[i], lat1[i], lon2[0], lat2[0]);
				for (int j=1; j<m; j++) {
					double s12 = dfun(lon1[i], lat1[i], lon2[j], lat2[j]);
					if (s12 < d[i]) {
						d[i] = s12;
					}
				}
				d[i] *= adj_unit;
			}
		}
	}

}





// input lon lat in radians
void distanceCosineToNearest_lonlat(std::vector<double> &d, const std::vector<double> &lon1, const std::vector<double> &lat1, const std::vector<double> &lon2, const std::vector<double> &lat2) {
	int n = lon1.size();
	int m = lon2.size();
	double s12;
 	for (int i=0; i < n; i++) {
		if (std::isnan(lat1[i])) {
			continue;
		}
		d[i] = distCosine(lat1[i], lon1[i], lat2[0], lon2[0]);
		for (int j=1; j<m; j++) {
			s12 = 6378137 * acos((sin(lat1[i]) * sin(lat2[j]) + cos(lat1[i]) * cos(lat2[j]) * cos(lon1[i]-lon2[j])));
			if (s12 < d[i]) {
				d[i] = s12;
			}
		}
	}
}




void distanceToNearest_plane(std::vector<double> &d, const std::vector<double> &x1, const  std::vector<double> &y1, const std::vector<double> &x2, const std::vector<double> &y2, const double& lindist) {
	int n = x1.size();
	int m = x2.size();
  	for (int i=0; i < n; i++) {
		if (std::isnan(x1[i])) continue;
		d[i] = sqrt(pow((x2[0]-x1[i]),2) + pow((y2[0]-y1[i]), 2));
		for (int j=1; j < m; j++) {
			double r = sqrt(pow((x2[j]-x1[i]),2) + pow((y2[j]-y1[i]), 2));
			if (r < d[i]) {
				d[i] = r;
			}
		}
		d[i] *= lindist;
	}
}



void nearest_lonlat(std::vector<long> &id, std::vector<double> &d, std::vector<double> &nlon, std::vector<double> &nlat, const std::vector<double> &lon1, const std::vector<double> &lat1, const std::vector<double> &lon2, const std::vector<double> &lat2, const std::string method) {
	
	size_t n = lon1.size();
	size_t m = lon2.size();

	nlon.resize(n);
	nlat.resize(n);
	id.resize(n);
	d.resize(n);
	
	if (method == "geo") {
		
		double a = 6378137.0;
		double f = 1/298.257223563;
		double azi1, azi2, s12;
		struct geod_geodesic g;
		geod_init(&g, a, f);

		for (size_t i=0; i < n; i++) {
			if (std::isnan(lat1[i])) {
				nlon[i] = NAN;
				nlat[i] = NAN;
				id[i] = -1;
				d[i] = NAN;
				continue;
			}
			geod_inverse(&g, lat1[i], lon1[i], lat2[0], lon2[0], &d[i], &azi1, &azi2);
			nlon[i] = lon2[0];
			nlat[i] = lat2[0];
			id[i] = 0;
			for (size_t j=1; j<m; j++) {
				geod_inverse(&g, lat1[i], lon1[i], lat2[j], lon2[j], &s12, &azi1, &azi2);
				if (s12 < d[i]) {
					d[i] = s12;
					id[i] = j;
					nlon[i] = lon2[j];
					nlat[i] = lat2[j];
				}
			}
		}
	} else {
		
		std::function<double(double, double, double, double)> dfun;
		if (method == "haversine") {
			dfun = distHaversine;
		} else {
			dfun = distCosine;			
		}	
		
		for (size_t i=0; i < n; i++) {
			if (std::isnan(lat1[i])) {
				nlon[i] = NAN;
				nlat[i] = NAN;
				id[i] = -1;
				d[i] = NAN;
				continue;
			}
			d[i] = dfun(lon1[i], lat1[i], lon2[0], lat2[0]);
			nlon[i] = lon2[0];
			nlat[i] = lat2[0];
			id[i] = 0;
			for (size_t j=1; j<m; j++) {
				double s12 = dfun(lon1[i], lat1[i], lon2[j], lat2[j]);
				if (s12 < d[i]) {
					d[i] = s12;
					id[i] = j;
					nlon[i] = lon2[j];
					nlat[i] = lat2[j];
				}
			}
		}	
	} 
}

void nearest_lonlat_self(std::vector<long> &id, std::vector<double> &d, std::vector<double> &nlon, std::vector<double> &nlat, const std::vector<double> &lon, const std::vector<double> &lat, const std::string method) {

	size_t n = lon.size();
	if (n <= 1) {
		nlon = lon;
		nlat = lat;
		if (nlon.size() == 1) {
			id.resize(1);
			id[0] = 0;
		}
		return;
	}

	nlon.resize(n);
	nlat.resize(n);
	id.resize(n);
	d.resize(n);

	if (method == "geo") {
		
		double a = 6378137.0;
		double f = 1/298.257223563;
		double azi1, azi2, s12;
		struct geod_geodesic g;
		geod_init(&g, a, f);
		for (size_t i=0; i < n; i++) {
			if (std::isnan(lat[i])) {
				id[i] = -1;
				d[i] = NAN;
				nlon[i] = NAN;
				nlat[i] = NAN;
				continue;
			}
			if (i>0) {
				geod_inverse(&g, lat[i], lon[i], lat[0], lon[0], &d[i], &azi1, &azi2);
				nlon[i] = lon[0];
				nlat[i] = lat[0];
				id[i] = 0;
			} else {
				geod_inverse(&g, lat[1], lon[1], lat[0], lon[0], &d[i], &azi1, &azi2);
				nlon[i] = lon[1];
				nlat[i] = lat[1];
				id[i] = 1;
			}

			for (size_t j=1; j < n; j++) {
				if (j == i) continue;
				geod_inverse(&g, lat[i], lon[i], lat[j], lon[j], &s12, &azi1, &azi2);
				if (s12 < d[i]) {
					d[i] = s12;
					id[i] = j;
					nlon[i] = lon[j];
					nlat[i] = lat[j];
				}
			}
		}

	} else {

		std::function<double(double, double, double, double)> dfun;
		if (method == "haversine") {
			dfun = distHaversine;
		} else if (method == "cosine") {
			dfun = distCosine;			
		} else {
			dfun = distLonlat;
		}

		for (size_t i=0; i < n; i++) {
			if (std::isnan(lat[i])) {
				id[i] = -1;
				d[i] = NAN;
				nlon[i] = NAN;
				nlat[i] = NAN;
				continue;
			}
			if (i>0) {
				d[i] = dfun(lon[i], lat[i], lon[0], lat[0]);
				nlon[i] = lon[0];
				nlat[i] = lat[0];
				id[i] = 0;
			} else {
				d[i] = dfun(lon[1], lat[1], lon[0], lat[0]);
				nlon[i] = lon[1];
				nlat[i] = lat[1];
				id[i] = 1;
			}

			for (size_t j=1; j < n; j++) {
				if (j == i) continue;
				double s12 = dfun(lon[i], lat[i], lon[j], lat[j]);

				if (s12 < d[i]) {
					d[i] = s12;
					id[i] = j;
					nlon[i] = lon[j];
					nlat[i] = lat[j];
				}
			}
		}
	} 
}

