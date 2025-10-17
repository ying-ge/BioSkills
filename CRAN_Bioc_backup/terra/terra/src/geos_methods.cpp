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


#include <numeric>
#include "geos_spat.h"
#include "distance.h"
#include "recycle.h"
#include "string_utils.h"

void callbck(void *item, void *userdata) { // callback function for tree selection
	std::vector<size_t> *ret = (std::vector<size_t> *) userdata;
	ret->push_back(*((size_t *) item));
}

SpatVector SpatVector::allerretour() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	SpatVector out = vect_from_geos(g, hGEOSCtxt, type());
	geos_finish(hGEOSCtxt);
	return out;
}

SpatVectorCollection SpatVector::bienvenue() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	SpatVectorCollection out = coll_from_geos(g, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	return out;
}

/*
std::vector<std::string> SpatVector::wkt() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<std::string> out::crop
	out.reserve(g.size());
	char * wkt;
	for (size_t i = 0; i < g.size(); i++) {
		wkt = GEOSGeomToWKT_r(hGEOSCtxt, g[i].get());
		out.push_back(wkt);
	}
	geos_finish(hGEOSCtxt);
	return out;
}
*/

std::vector<std::string> SpatVector::wkt() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<std::string> out;
	out.reserve(g.size());
	GEOSWKTWriter* writer = GEOSWKTWriter_create_r(hGEOSCtxt);
	for (size_t i = 0; i < g.size(); i++) {
		char *wkt = GEOSGeomToWKT_r(hGEOSCtxt, g[i].get());
		out.push_back(wkt);
	}
	GEOSWKTWriter_destroy_r(hGEOSCtxt, writer);
	geos_finish(hGEOSCtxt);
	return out;
}



std::vector<std::string> SpatVector::wkb() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<std::string> out;
	out.reserve(g.size());
	GEOSWKBWriter* writer = GEOSWKBWriter_create_r(hGEOSCtxt);
	size_t len=0;
	for (size_t i=0; i<g.size(); i++) {
		unsigned char *wkb = GEOSWKBWriter_write_r(hGEOSCtxt, writer, g[i].get(), &len);
		std::string s(reinterpret_cast<char const*>(wkb), len);
		out.push_back(s);
		GEOSFree_r(hGEOSCtxt, wkb);
	}
	GEOSWKBWriter_destroy_r(hGEOSCtxt, writer);
	geos_finish(hGEOSCtxt);
	return out;
}

std::vector<std::vector<unsigned char>> SpatVector::wkb_raw() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<std::vector<unsigned char>> out; 
	size_t len = 0;
	for (size_t i = 0; i < g.size(); i++) {
		unsigned char *hex = GEOSGeomToWKB_buf_r(hGEOSCtxt, g[i].get(), &len);
		std::vector<unsigned char> raw; 
		raw = std::vector<unsigned char>(hex, hex+len);
		out.push_back(raw);
		free(hex);
	}
	geos_finish(hGEOSCtxt);
	return out;
}	
	
std::vector<std::string> SpatVector::hex() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<std::string> out;
	out.reserve(g.size());
	size_t len = 0;
	for (size_t i = 0; i < g.size(); i++) {
		unsigned char *hex = GEOSGeomToHEX_buf_r(hGEOSCtxt, g[i].get(), &len);
		std::string s( reinterpret_cast<char const*>(hex), len) ;
		out.push_back(s);
		free(hex);
	}
	geos_finish(hGEOSCtxt);
	return out;
}


SpatVector SpatVector::from_hex(std::vector<std::string> x, std::string srs) {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	size_t n = x.size();
	std::vector<GeomPtr> p;
	p.resize(n);
	for (size_t i = 0; i < n; i++) {
		const char* cstr = x[i].c_str();
		size_t len = strlen(cstr);
		const unsigned char *hex = (const unsigned char *) cstr;
		GEOSGeometry* r = GEOSGeomFromHEX_buf_r(hGEOSCtxt, hex, len);
		p[i] = geos_ptr(r, hGEOSCtxt);
	}
	SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	SpatVector out = coll.get(0);
	if (coll.size() > 1) {
		out.addWarning("not all geometries were transferred, use svc for a geometry collection");
	}
	out.setSRS(srs);
	return out;
}



SpatVectorCollection SpatVectorCollection::from_hex_col(std::vector<std::string> x, std::string srs) {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	size_t n = x.size();
	std::vector<GeomPtr> p;
	p.resize(n);
	for (size_t i = 0; i < n; i++) {
		const char* cstr = x[i].c_str();
		size_t len = strlen(cstr);
		const unsigned char *hex = (const unsigned char *) cstr;
		GEOSGeometry* r = GEOSGeomFromHEX_buf_r(hGEOSCtxt, hex, len);
		p[i] = geos_ptr(r, hGEOSCtxt);
	}
	SpatVectorCollection out = coll_from_geos(p, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	for (size_t i = 0; i < out.size(); i++) {
		out.v[i].setSRS(srs);
	}
	return out;
}



std::vector<bool> SpatVector::geos_isvalid() {
	GEOSContextHandle_t hGEOSCtxt = geos_init2();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<bool> out;
	out.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		char v = GEOSisValid_r(hGEOSCtxt, g[i].get());
		out.push_back(v);
	}
	geos_finish(hGEOSCtxt);
	return {out};
}

std::vector<std::string> SpatVector::geos_isvalid_msg() {
	GEOSContextHandle_t hGEOSCtxt = geos_init2();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<std::string> out;
	out.reserve(2 * g.size());
	for (size_t i = 0; i < g.size(); i++) {
		char v = GEOSisValid_r(hGEOSCtxt, g[i].get());
		std::string valid = {v};
		out.push_back(valid);
		if (!v) {
			char *r = GEOSisValidReason_r(hGEOSCtxt, g[i].get());
			std::string reason = r;
			free(r);
			out.push_back(reason);
		} else {
			out.push_back("");
		}
	}
	geos_finish(hGEOSCtxt);
	return {out};
}


SpatVector SpatVector::make_valid2() {

	SpatVector out;
#ifndef GEOS380
	out.setError("make_valid is not available for GEOS < 3.8");
#else
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);

	size_t n = size();
	std::vector<long> ids;
	ids.reserve(n);

	for (size_t i=0; i<n; i++) {
		GEOSGeometry* r = GEOSMakeValid_r(hGEOSCtxt, x[i].get());
		if (r != NULL) {
			if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
				x[i] = geos_ptr(r, hGEOSCtxt);
				//ids.push_back(i);
			} else {
				GEOSGeom_destroy_r(hGEOSCtxt, r);
			}
		}
	}
	SpatVectorCollection coll = coll_from_geos(x, hGEOSCtxt, ids, false, false);
	out = coll.get(0);
	geos_finish(hGEOSCtxt);
	out.srs = srs;
//	if (ids.size() != n) {
//		out.df = df.subset_rows(out.df.iv[0]);
//	} else {
		out.df = df;
//	}
#endif
	return out;
}


/*
need to find the correct version required 
does not work with 3.4.2 (see #734)
#ifndef GEOS350
	out.setError("GEOS 3.5 required");
#else

SpatVector SpatVector::set_precision(double gridSize) {
	SpatVector out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r = GEOSGeom_setPrecision_r(hGEOSCtxt, g[i].get(), gridSize, 0);
		if (r == NULL) {
			out.setError("something bad happened");
			geos_finish(hGEOSCtxt);
			return out;
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
			p.push_back(geos_ptr(r, hGEOSCtxt));
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	if (p.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.df = df;
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	return out;
}
*/


std::vector<std::vector<size_t>> SpatVector::index_2d(SpatVector v) {
	std::vector<std::vector<size_t>> out(2);
	size_t n = std::max(size(), v.size()) * 2;
	out[0].reserve(n);
	out[1].reserve(n);
	size_t k = 0;
	for (size_t i=0; i<size(); i++) {
		for (size_t j=0; j<size(); j++) {
			if (geoms[i].extent.intersects(v.geoms[j].extent)) {
				out[0].push_back(i);
				out[1].push_back(j);
				k++;
				if (k > n) {
					n += std::max(size(), v.size());
					out[0].reserve(n);
					out[1].reserve(n);
				}
			}
		}
	}
	return out;
}


std::vector<std::vector<size_t>> SpatVector::index_sparse(SpatVector v) {
	std::vector<std::vector<size_t>> out(v.size());
	for (size_t i=0; i<size(); i++) {
		for (size_t j=0; j<size(); j++) {
			if (geoms[i].extent.intersects(v.geoms[j].extent)) {
				out[i].push_back(j);
			}
		}
	}
	return out;
}



SpatVector SpatVector::crop(SpatExtent e, bool wrap) {

	SpatVector out;
	if ((e.xmin == e.xmax) || (e.ymin == e.ymax)) {
		return *this;	
	}


#ifndef GEOS350
	out.setError("GEOS 3.5 required for crop");
	return out;
#else
	if (is_lonlat() & wrap) {
		SpatVector first = crop(e, false);
		if ((extent.xmin > -180.001) && (extent.xmax < 180.001)) {
			if ((e.xmin < -180) && (e.xmax < 180)) {
				double xmn = e.xmin + 360;
				if (xmn > -180) {
					if (xmn < e.xmax) xmn = e.xmax;
					SpatExtent e2 = e;
					e2.xmax = 180;
					e2.xmin = xmn;
					SpatVector second = crop(e2, false);
					first = first.append(second, true);
				}
			} else if ((e.xmax > 180) && (e.xmin > -180)) {
				double xmx = -360 + e.xmax;
				if (xmx < 180) {
					if (xmx > e.xmin) xmx = e.xmin;
					SpatExtent e2 = e;
					e2.xmin = -180;
					e2.xmax = xmx;
					SpatVector second = crop(e2, false);
					first = first.append(second, true);
				}
			} 
		}
		return(first);						
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	std::vector<long> id;
	id.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r = GEOSClipByRect_r(hGEOSCtxt, g[i].get(), e.xmin, e.ymin, e.xmax, e.ymax);
		if (r == NULL) {
			out.setError("something bad happened");
			geos_finish(hGEOSCtxt);
			return out;
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
			p.push_back(geos_ptr(r, hGEOSCtxt));
			id.push_back(i);
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	if (!p.empty()) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt, id);
		out = coll.get(0);
		out.df = df.subset_rows(out.df.iv[0]);
	}
	geos_finish(hGEOSCtxt);
	
	out.srs = srs;
	return out;
#endif
}




SpatVector SpatVector::make_nodes() {

	SpatVector out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r = GEOSNode_r(hGEOSCtxt, g[i].get());
		if (r == NULL) {
			out.setError("something bad happened");
			geos_finish(hGEOSCtxt);
			return out;
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
			p.push_back(geos_ptr(r, hGEOSCtxt));
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	if (!p.empty()) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.df = df;
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	return out;
}




SpatVector SpatVector::boundary() {

	SpatVector out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r = GEOSBoundary_r(hGEOSCtxt, g[i].get());
		if (r == NULL) {
			out.setError("something bad happened");
			geos_finish(hGEOSCtxt);
			return out;
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
			p.push_back(geos_ptr(r, hGEOSCtxt));
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	if (!p.empty()) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.df = df;
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	return out;
}


SpatVector SpatVector::normalize() {

	SpatVector out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r = g[i].get();
		if (GEOSNormalize_r(hGEOSCtxt, r)) {
			g[i] = geos_ptr(r, hGEOSCtxt);
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	out = vect_from_geos(g, hGEOSCtxt, type());
	geos_finish(hGEOSCtxt);
	out.df = df;
	out.srs = srs;
	return out;
}



SpatVector SpatVector::line_merge() {

	SpatVector out;
	if (type() != "lines") {
		out.setError("input must be lines");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r = GEOSLineMerge_r(hGEOSCtxt, g[i].get());
		if (r == NULL) {
			out.setError("something bad happened");
			geos_finish(hGEOSCtxt);
			return out;
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
			p.push_back(geos_ptr(r, hGEOSCtxt));
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	if (!p.empty()) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.df = df;
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	return out;
}



SpatVector SpatVector::simplify(double tolerance, bool preserveTopology) {
	SpatVector out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r;
		if (preserveTopology) {
			r = GEOSTopologyPreserveSimplify_r(hGEOSCtxt, g[i].get(), tolerance);
		} else {
			r = GEOSSimplify_r(hGEOSCtxt, g[i].get(), tolerance);
		}
		if (r == NULL) {
			out.setError("something bad happened");
			geos_finish(hGEOSCtxt);
			return out;
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
			p.push_back(geos_ptr(r, hGEOSCtxt));
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	if (!p.empty()) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.df = df;
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	return out;
}


SpatVector SpatVector::shared_paths(bool index) {

	if (type() == "polygons") {
		SpatVector x = as_lines();
		return x.shared_paths(index);
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);

	size_t s = size();
	std::vector<long> id1, id2;
	std::vector<GeomPtr> p;

	if (!index) {
		// calculate shared paths
		for (size_t i=0; i<(s-1); i++) {
			for (size_t j=(i+1); j<s; j++) {
				GEOSGeometry* r = GEOSSharedPaths_r(hGEOSCtxt, x[i].get(), x[j].get());
				if (r != NULL) {
					if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
						p.push_back(geos_ptr(r, hGEOSCtxt));
						id1.push_back(i+1);
						id2.push_back(j+1);
					} else {
						GEOSGeom_destroy_r(hGEOSCtxt, r);
					}
				}
			}
		}
	} else {
		// use spatial index
		std::vector<size_t> items(x.size());
		TreePtr tree1 = geos_ptr(GEOSSTRtree_create_r(hGEOSCtxt, 10), hGEOSCtxt);
		for (size_t i = 0; i < s; i++) {
			items[i] = i;
			if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get()))
				GEOSSTRtree_insert_r(hGEOSCtxt, tree1.get(), x[i].get(), &(items[i]));
		}
		for (size_t i = 0; i < s; i++) {
			// pre-select x's using tree:
			std::vector<size_t> tree_sel, sel;
			if (!GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
				GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), x[i].get(), callbck, &tree_sel);
			}
			if (! tree_sel.empty()) {
				for (size_t j = 0; j < tree_sel.size(); j++) {
					if (tree_sel[j] > i) {
						GEOSGeometry* r = GEOSSharedPaths_r(hGEOSCtxt, x[i].get(), x[tree_sel[j]].get());
						if (r != NULL) {
							if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
								p.push_back(geos_ptr(r, hGEOSCtxt));
								id1.push_back(i+1);
								id2.push_back(tree_sel[j]+1);
							} else {
								GEOSGeom_destroy_r(hGEOSCtxt, r);
							}
						}
					}
				}
			}
		}
	}

	SpatVector out;
	if (!p.empty()) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt, std::vector<long>(), false, false);
		out = coll.get(0);
		out = out.line_merge();
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	out.df.add_column(id1, "id1");
	out.df.add_column(id2, "id2");
	return out;
}


SpatVector SpatVector::shared_paths(SpatVector x, bool index) {

	if (x.type() == "polygons") {
		x = x.as_lines();
	}
	if (type() == "polygons") {
		SpatVector v = as_lines();
		return v.shared_paths(x, index);
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> a = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> b = geos_geoms(&x, hGEOSCtxt);

	size_t sa = size();
	size_t sb = b.size();
	std::vector<long> id1, id2;
	std::vector<GeomPtr> p;

	if (!index) {
		for (size_t i=0; i<sa; i++) {
			for (size_t j=0; j<sb; j++) {
				GEOSGeometry* r = GEOSSharedPaths_r(hGEOSCtxt, a[i].get(), b[j].get());
				if (r != NULL) {
					if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
						p.push_back(geos_ptr(r, hGEOSCtxt));
						id1.push_back(i+1);
						id2.push_back(j+1);
					} else {
						GEOSGeom_destroy_r(hGEOSCtxt, r);
					}
				}
			}
		}
	} else {
		// use spatial index
		std::vector<size_t> items(x.size());
		TreePtr tree1 = geos_ptr(GEOSSTRtree_create_r(hGEOSCtxt, 10), hGEOSCtxt);
		for (size_t i = 0; i < sb; i++) {
			items[i] = i;
			if (! GEOSisEmpty_r(hGEOSCtxt, b[i].get()))
				GEOSSTRtree_insert_r(hGEOSCtxt, tree1.get(), b[i].get(), &(items[i]));
		}
		for (size_t i = 0; i < sa; i++) {
			// pre-select x's using tree:
			std::vector<size_t> tree_sel, sel;
			if (!GEOSisEmpty_r(hGEOSCtxt, a[i].get())) {
				GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), a[i].get(), callbck, &tree_sel);
			}
			if (! tree_sel.empty()) {
				for (size_t j = 0; j < tree_sel.size(); j++) {
					GEOSGeometry* r = GEOSSharedPaths_r(hGEOSCtxt, a[i].get(), b[tree_sel[j]].get());
					if (r != NULL) {
						if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
							p.push_back(geos_ptr(r, hGEOSCtxt));
							id1.push_back(i+1);
							id2.push_back(tree_sel[j]+1);
						} else {
							GEOSGeom_destroy_r(hGEOSCtxt, r);
						}
					}
				}
			}
		}
	}

	SpatVector out;
	if (!p.empty()) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt, std::vector<long>(), false, false);
		out = coll.get(0);
		out = out.line_merge();
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	out.df.add_column(id1, "id1");
	out.df.add_column(id2, "id2");
	return out;
}


#ifndef GEOS370
bool find_segments(GEOSContextHandle_t hGEOSCtxt, std::vector<double> &x, std::vector<double> &y, std::vector<double> &cx, std::vector<double> &cy, std::vector<size_t> &si, std::vector<double> &sx, std::vector<double> &sy) {
	return false;
}	
#else
bool find_segments(GEOSContextHandle_t hGEOSCtxt, std::vector<double> &x, std::vector<double> &y, std::vector<double> &cx, std::vector<double> &cy, std::vector<size_t> &si, std::vector<double> &sx, std::vector<double> &sy) {
	size_t n = x.size() - 1;
	size_t m = cx.size() - 1;
	double ix, iy;
	si.resize(0);
	sx.resize(0);
	sy.resize(0);
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j<m; j++) {
			if (GEOSSegmentIntersection_r(hGEOSCtxt, x[i], y[i], x[i+1], cy[j+1], cx[j], cy[j], cx[j+1], cy[j+1], &ix, &iy) == 1) {
				si.push_back(i);
				sx.push_back(ix);
				sy.push_back(iy);
			}
		}
	}
	return si.size() > 0;
}
#endif

SpatVector SpatVector::split_lines(SpatVector v) {

// check GEOS version 

	SpatVector out;
	#ifndef GEOS370
	out.setError("not available with GEOS < 3.7");
	return out;
	#endif

/*
	std::vector<size_t> si;
	std::vector<double> sx, sy;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	for (size_t i=0; i<v.size(); i++) {
		SpatVector tmp = v.subset_rows(i);
		std::vector<int> x = relate(tmp, "intersects", true, true);
		std::vector<std::vector<double>> xy1 = tmp.coordinates();
		for (size_t j=0; j<x.size(); j++) {
			if (x[j] == 1) {
				std::vector<std::vector<double>> xy2 = tmp.coordinates();
				if (find_segments(hGEOSCtxt, xy1[0], xy1[1], xy2[0], xy2[1], si, sx, sy)) {
					
				}
			}
		}
	}
	
*/	
	return out;

}	


/*
SpatVector SpatVector::split_polygons(SpatVector lns) {
		SpatGeom glns;
		glns.gtype = lines;
		glns.setPart(SpatPart(x, y), 0);
		std::vector<double> xln = {180, 180};
		std::vector<double> yln = {-91, 91};
		glns.setPart(SpatPart(xln, yln), 1);
		SpatVector v;
		v.addGeom(glns);
		v = v.line_merge();
		v = v.aggregate(false);
		v = v.polygonize();
		g = v.geoms[0];
*/



SpatVector polygonize_one(const GEOSGeometry* gi, GEOSContextHandle_t hGEOSCtxt) {

	size_t ngeoms = 1;
	std::vector<GeomPtr> p(1);
	SpatVector out;

	GEOSGeometry* r = GEOSPolygonize_r(hGEOSCtxt, &gi, ngeoms);
	if (r == NULL) {
		out.setError("something bad happened");
		geos_finish(hGEOSCtxt);
		return out;
	}
	if (GEOSisEmpty_r(hGEOSCtxt, r)) {
		GEOSGeom_destroy_r(hGEOSCtxt, r);
	} else {
		p[0] = geos_ptr(r, hGEOSCtxt);
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.aggregate(false);
	}
	return out;
}

SpatVector SpatVector::polygonize() {

	if (type() == "polygons") {
		return *this;
	}

	SpatVector out;
	if (type() == "points") {
		out.setError("cannot make polygons from points (make lines first)");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);

	const GEOSGeometry* gi = g[0].get();
	out = polygonize_one(gi, hGEOSCtxt);

	for (size_t i = 1; i < g.size(); i++) {
		const GEOSGeometry* gi = g[i].get();
		SpatVector onegeom = polygonize_one(gi, hGEOSCtxt);
		out.addGeom(onegeom.getGeom(0));
	}
	geos_finish(hGEOSCtxt);

	out.srs = srs;
	if (df.nrow() == out.size()) {
		out.df = df;
	}
	return out;
}


SpatVector SpatVector::snap(double tolerance) {

	size_t s = size();
	SpatVector out;
	if (s == 0) {
		return out;
	}

	tolerance = std::max(0.0, tolerance);
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);

	std::vector<long> ids;
	ids.reserve(s);

	for (size_t i=0; i<(s-1); i++) {
		GEOSGeometry* r = x[i].get();
		for (size_t j=(i+1); j<s; j++) {
			r = GEOSSnap_r(hGEOSCtxt, r, x[j].get(), tolerance);
		}
		if (r != NULL) {
			if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
				x[i] = geos_ptr(r, hGEOSCtxt);
				ids.push_back(i);
			} else {
				GEOSGeom_destroy_r(hGEOSCtxt, r);
			}
		}
	}
	ids.push_back(s-1);
	SpatVectorCollection coll = coll_from_geos(x, hGEOSCtxt, ids, false, false);
	out = coll.get(0);
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	if (ids.size() != s) {
		out.df = df.subset_rows(out.df.iv[0]);
	} else {
		out.df = df;
	}
	return out;
}

SpatVector SpatVector::snapto(SpatVector y, double tolerance) {

	y = y.aggregate(false);
	size_t s = size();

	SpatVector out;
	if (s == 0) {
		return out;
	}

	tolerance = std::max(0.0, tolerance);
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> to = geos_geoms(&y, hGEOSCtxt);

	std::vector<long> ids;
	ids.reserve(s);

	GEOSGeometry* gto = to[0].get();
	for (size_t i=0; i<s; i++) {
		GEOSGeometry* r = GEOSSnap_r(hGEOSCtxt, x[i].get(), gto, tolerance);
		if (r != NULL) {
			if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
				x[i] = geos_ptr(r, hGEOSCtxt);
				ids.push_back(i);
			} else {
				GEOSGeom_destroy_r(hGEOSCtxt, r);
			}
		}
	}
	SpatVectorCollection coll = coll_from_geos(x, hGEOSCtxt, ids);
	out = coll.get(0);
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	if (ids.size() != s) {
		out.df = df.subset_rows(out.df.iv[0]);
	} else {
		out.df = df;
	}
	return out;
}


//GEOSPolygonizer_getCutEdges_r(GEOSContextHandle_t extHandle, const Geometry * const * g, unsigned int ngeoms)

//Geometry * GEOSPolygonize_full_r(GEOSContextHandle_t extHandle, const Geometry* g, Geometry** cuts, Geometry** dangles, Geometry** invalid)


SpatVector SpatVector::crop(SpatVector v) {

	SpatVector out;
	out.srs = srs;
	if (v.empty()) return *this;

	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
//	if ((type() != "polygons") & (type() != "mutlipolygons")) {
	if ((v.type() != "polygons")) {
		v = v.hull("convex");
	} else if (v.nrow() > 1) {
		v = v.aggregate(false);
	}
	
//	return intersect(v, false);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	std::vector<GeomPtr> result;
	std::vector<long> ids;
	size_t nx = size();
	ids.reserve(nx);

	for (size_t i = 0; i < nx; i++) {
		GEOSGeometry* geom = GEOSIntersection_r(hGEOSCtxt, x[i].get(), y[0].get());
		if (geom == NULL) {
			out.setError("GEOS exception");
			geos_finish(hGEOSCtxt);
			return(out);
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, geom)) {
			result.push_back(geos_ptr(geom, hGEOSCtxt));
			ids.push_back(i);
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, geom);
		}
	}

//	SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);

	if (!result.empty()) {
//		SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);
//		SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt, ids);
		SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt, ids, true, false);
		out = coll.get(0);	
//		std::vector<std::string> nms = out.get_names();
//		out = out.aggregate(nms[0], true);
		out.df = df.subset_rows(out.df.iv[0]);
		out.srs = srs;
	}
	geos_finish(hGEOSCtxt);
	return out;
}



SpatVector SpatVector::hull(std::string htype, std::string by, double param, bool allowHoles, bool tight) {

	SpatVector out;
	if (nrow() == 0) {
		out.srs = srs;
		return out;
	}

	std::vector<std::string> methods = {"convex", "rectangle", "circle", "concave_ratio", "concave_length"};
	if (std::find(methods.begin(), methods.end(), htype) == methods.end()) {
		out.setError("unknown hull type");
		return out;
	}

	if (!by.empty()) {
		SpatVector tmp = aggregate(by, false);
		if (tmp.hasError()) {
			return tmp;
		}
		for (size_t i=0; i<tmp.size(); i++) {
			SpatVector x = tmp.subset_rows(i);
			x = x.hull(htype, "");
			if (x.hasError()) {
				return x;
			}
			if (!x.geoms.empty() && (x.geoms[0].gtype == polygons)) {
				out.addGeom(x.geoms[0]);
			} else {
				SpatGeom g(polygons);
				out.addGeom(g);
			}
		}
		out.df = tmp.df;
		out.srs = srs;
		return out;
	}

	out.reserve(size());

/*
	if (is_lonlat()) {
		if ((extent.ymin > -85) && (extent.ymax < 85)) {
			SpatVector tmp = project("+proj=merc", false);
			tmp = tmp.hull(htype, "");
			tmp = tmp.project(srs.wkt, false);
			return tmp;
		}
	}
*/

	SpatVector a = aggregate(false);

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(&a, hGEOSCtxt);
	//std::string vt = type();
	GEOSGeometry* h;
	if (htype == "convex") {
		h = GEOSConvexHull_r(hGEOSCtxt, g[0].get());	
	} else if (htype == "circle") {
	#ifndef GEOS380
		geos_finish(hGEOSCtxt);
		out.setError("GEOS 3.8 required for bounding circle");
		return out;
	#else
		h = GEOSMinimumBoundingCircle_r(hGEOSCtxt, g[0].get(), NULL, NULL);
	#endif
	} else if (htype == "rectangle") {
	#ifndef GEOS361
		geos_finish(hGEOSCtxt);
		out.setError("GEOS 3.6.1 required for rotated rectangle");
		return out;
	#else
		h = GEOSMinimumRotatedRectangle_r(hGEOSCtxt, g[0].get());
	#endif
	} else if (htype == "concave_ratio") {
	#ifndef GEOS3110
		geos_finish(hGEOSCtxt);
		out.setError("GEOS 3.11 required for concave hull");
		return out;
	#else
		h = GEOSConcaveHull_r(hGEOSCtxt, g[0].get(), param, allowHoles);
	#endif
	} else if (htype == "concave_length") {
	#ifndef GEOS3110
		geos_finish(hGEOSCtxt);
		out.setError("GEOS 3.11 required for concave_length hull");
		return out;
	#else
		if (type() == "polygons") {
			h = GEOSConcaveHullOfPolygons_r(hGEOSCtxt, g[0].get(), param, tight, allowHoles);
		} else {
		#ifndef GEOS3120
				geos_finish(hGEOSCtxt);
				out.setError("GEOS 3.12 required for concave_length hull for points and lines");
				return out;
		#else 
				h = GEOSConcaveHullByLength_r(hGEOSCtxt, g[0].get(), param, allowHoles);
		#endif	
		}
	#endif
	} else {
		geos_finish(hGEOSCtxt);
		out.setError("unknown hull type");
		return out;
	}
	
	std::vector<GeomPtr> b(1);
	b[0] = geos_ptr(h, hGEOSCtxt);
	SpatVectorCollection coll = coll_from_geos(b, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	out = coll.get(0);
	out.srs = srs;
	return out;

}




SpatVector SpatVector::voronoi(SpatVector bnd, double tolerance, int onlyEdges) {
	SpatVector out;

	if (nrow() == 0) {
		out.addWarning("input SpatVector has no geometries");
		return out;
	}

#ifndef GEOS350
	out.setError("GEOS 3.5 required for voronoi");
	return out;
#else

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	SpatVector a = aggregate(false);
	std::vector<GeomPtr> g = geos_geoms(&a, hGEOSCtxt);
	GEOSGeometry* v;
	if (bnd.empty()) {
		v = GEOSVoronoiDiagram_r(hGEOSCtxt, g[0].get(), NULL, tolerance, onlyEdges);
	} else {
		if (bnd.type() != "polygons") {
			out.setError("boundary must have a polygon geometry");
			geos_finish(hGEOSCtxt);
			return out;
		}
		std::vector<GeomPtr> ge = geos_geoms(&bnd, hGEOSCtxt);
		v = GEOSVoronoiDiagram_r(hGEOSCtxt, g[0].get(), ge[0].get(), tolerance, onlyEdges);
	}
	if (v == NULL) {
		out.setError("GEOS exception");
		geos_finish(hGEOSCtxt);
		return(out);
	}
	std::vector<GeomPtr> b(1);
	b[0] = geos_ptr(v, hGEOSCtxt);
	SpatVectorCollection coll = coll_from_geos(b, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	out = coll.get(0);
	out.srs = srs;
	if (!out.hasError()) {
		out = out.disaggregate(false);
		if (!bnd.empty()) {
			SpatDataFrame empty;
			bnd.df = empty;
			out = out.intersect(bnd, true);
		}
		if ((type() == "points") && (!onlyEdges)) {
			std::vector<int> atts = out.relateFirst(*this, "intersects");
			std::vector<size_t> a;
			a.reserve(atts.size());
			for (size_t i=0; i<atts.size(); i++) {
				if (atts[i] >=0) a.push_back(atts[i]);
			}
			if (a.size() == out.size()) {
				out.df = df.subset_rows(a);
			}
		}
	}
	return out;
#endif
}




SpatVector SpatVector::delaunay(double tolerance, int onlyEdges, bool constrained) {
	SpatVector out;
	if (nrow() == 0) {
		out.addWarning("input SpatVector has no geometries");
		return out;
	}

#ifndef GEOS350
	out.setError("GEOS 3.5 required for delaunay");
	return out;
#endif 

#ifndef GEOS3100
	if (constrained) {
		out.setError("GEOS 3.10 required for constrained delaunay");
		return out;
	}
#else

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	SpatVector a = aggregate(false);
	std::vector<GeomPtr> g = geos_geoms(&a, hGEOSCtxt);

	GEOSGeometry* v;
	if (constrained) {
		v = GEOSConstrainedDelaunayTriangulation_r(hGEOSCtxt, g[0].get());
	} else {
		v = GEOSDelaunayTriangulation_r(hGEOSCtxt, g[0].get(), tolerance, onlyEdges);
	}	
	if (v == NULL) {
		out.setError("GEOS exception");
		geos_finish(hGEOSCtxt);
		return(out);
	}
	std::vector<GeomPtr> b(1);
	b[0] = geos_ptr(v, hGEOSCtxt);
	SpatVectorCollection coll = coll_from_geos(b, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	out = coll.get(0);
	out.srs = srs;
	if (!out.hasError()) {
		out = out.disaggregate(false);
		// associate with attributes
	}
	return out;
#endif
}



SpatVector SpatVector::buffer(std::vector<double> d, unsigned quadsegs, std::string capstyle, std::string joinstyle, double mitrelimit, bool singlesided) {

//	quadsegs = std::min(quadsegs, (unsigned) 180);
	SpatVector out;
	if (srs.is_empty()) {
		out.addWarning("unknown CRS. Results may be wrong");
	} 
	if (d.empty()) {
		out.setError("no buffer distance provided");
		return out;
	} 
	
	bool islonlat = is_lonlat();
	if (d.size() == 1 && d[0] == 0) {
		islonlat = false; //faster
	}
	std::string vt = type();
	if (vt == "points" || vt == "lines") {
		for (size_t i=0; i<d.size(); i++) {
			if (d[i] <= 0) {
				out.setError("a negative buffer is only meaningful with polygons");
				return out;
			}
		}
	}
	recycle(d, size());
	if (islonlat) {
		return buffer_lonlat(vt, d, quadsegs);
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();

	GEOSBufferParams* bufparms = GEOSBufferParams_create_r(hGEOSCtxt);	
	GEOSBufferParams_setQuadrantSegments_r(hGEOSCtxt, bufparms, quadsegs);
	if (capstyle == "flat") {
		GEOSBufferParams_setEndCapStyle_r(hGEOSCtxt, bufparms, GEOSBUF_CAP_FLAT);
	} else if (capstyle == "square") { 
		GEOSBufferParams_setEndCapStyle_r(hGEOSCtxt, bufparms, GEOSBUF_CAP_SQUARE);
	} else {
		GEOSBufferParams_setEndCapStyle_r(hGEOSCtxt, bufparms, GEOSBUF_CAP_ROUND);
	}
	if (joinstyle == "mitre") {
		GEOSBufferParams_setJoinStyle_r(hGEOSCtxt, bufparms, GEOSBUF_JOIN_MITRE);
	} else if (joinstyle == "bevel") { 
		GEOSBufferParams_setJoinStyle_r(hGEOSCtxt, bufparms, GEOSBUF_JOIN_BEVEL);
	} else {
		GEOSBufferParams_setJoinStyle_r(hGEOSCtxt, bufparms, GEOSBUF_JOIN_ROUND);
	}
	if (!std::isnan(mitrelimit)) { 
		GEOSBufferParams_setMitreLimit_r(hGEOSCtxt, bufparms, mitrelimit);
	}
	if (singlesided) { 
		GEOSBufferParams_setSingleSided_r(hGEOSCtxt, bufparms, 1);
	} 
	
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> b(size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* pt = GEOSBufferWithParams_r(hGEOSCtxt, g[i].get(), bufparms, d[i]);		
//		GEOSGeometry* pt = GEOSBuffer_r(hGEOSCtxt, g[i].get(), d[i], quadsegs);
		if (pt == NULL) {
			out.setError("GEOS exception");
			geos_finish(hGEOSCtxt);
			return(out);
		}
		b[i] = geos_ptr(pt, hGEOSCtxt);
	}
	const std::vector<long> ids = std::vector<long>();
	SpatVectorCollection coll = coll_from_geos(b, hGEOSCtxt, ids, false);

	GEOSBufferParams_destroy_r(hGEOSCtxt, bufparms);	
	geos_finish(hGEOSCtxt);
	out = coll.get(0);
	out.srs = srs;
	out.df = df;

// revdep problem with ShapePattern
	if (std::isnan(out.extent.xmin)) {
		SpatVector empty;
		empty.srs = srs;
		return empty;
	}
	
	return out;
}


// basic version of buffer, for debugging
SpatVector SpatVector::buffer2(std::vector<double> d, unsigned quadsegs) {

	SpatVector out;
	recycle(d, size());
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> b(size());
	for (size_t i = 0; i < g.size(); i++) {
		
//		Rcpp::Rcout << "buffer " << i;		
		GEOSGeometry* pt = GEOSBuffer_r(hGEOSCtxt, g[i].get(), d[i], (int) quadsegs);
//		Rcpp::Rcout << " done" << std::endl;

		if (pt == NULL) {
			out.setError("GEOS exception");
			geos_finish(hGEOSCtxt);
			return(out);
		}
		b[i] = geos_ptr(pt, hGEOSCtxt);
	}
	SpatVectorCollection coll = coll_from_geos(b, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	out = coll.get(0);
	out.srs = srs;
	out.df = df;
	return out;
}



SpatVector SpatVector::intersect(SpatVector v, bool values) {

	SpatVector out;
	out.srs = srs;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	//v = v.aggregate(false);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	std::vector<GeomPtr> result;
//	size_t nx = size();
//	size_t ny = v.size();
	std::vector<size_t> idx, idy;

	std::vector<std::vector<double>> r = which_relate(v, "intersects", true);
	size_t n = r[0].size();
	idx.reserve(n);
	idy.reserve(n);
	for (size_t i=0; i<n; i++) {
		idx.push_back( r[0][i] );
		idy.push_back( r[1][i] );
	}
	r.resize(0);
	std::vector<long> ids;
	ids.reserve(n);


	if (type() == "points") {
//		idx = wr[0];
//		idy = wr[1];

/*		for (size_t j = 0; j < ny; j++) {
			PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, y[j].get()), hGEOSCtxt);
			for (size_t i = 0; i < nx; i++) {
				if (GEOSPreparedIntersects_r(hGEOSCtxt, pr.get(), x[i].get())) {
					idx.push_back(i);
					idy.push_back(j);
				}
			}
		}
*/
		out = subset_rows(idx);

	} else {

		//long k = 0;
		for (size_t i = 0; i < n; i++) {
			GEOSGeometry* geom = GEOSIntersection_r(hGEOSCtxt, x[idx[i]].get(), y[idy[i]].get());
			if (geom == NULL) {
				out.setError("GEOS exception");
				geos_finish(hGEOSCtxt);
				return(out);
			}
			if (!GEOSisEmpty_r(hGEOSCtxt, geom)) {
				result.push_back(geos_ptr(geom, hGEOSCtxt));
				//idx.push_back(i);
				//idy.push_back(j);
				ids.push_back(i);
				//k++;
			} else {
				GEOSGeom_destroy_r(hGEOSCtxt, geom);
			}
		}

	//SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);
		if (!result.empty()) {
			SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt, ids, false, false);
			out = coll.get(0);
			out.srs = srs;
		}
	}
	geos_finish(hGEOSCtxt);

	if (!srs.is_same(v.srs, true)) {
		out.addWarning("different crs");
	}

	if ((type() == "polygons") && (v.type() == "polygons") && (out.type() != "polygons")) {
		// intersection is point or line, return empty
		out = SpatVector();
		out.addWarning("no intersection");
		out.srs = srs;
	}

	SpatDataFrame df1, df2;
	n = out.nrow();
	if (values) {
		if (n < idx.size()) {
			std::vector<size_t> idx2, idy2;
			idx2.reserve(n);
			idy2.reserve(n);
			for (size_t i=0; i<n; i++) {
				idx2.push_back( idx[ out.df.iv[0][i] ]);
				idy2.push_back( idy[ out.df.iv[0][i] ]);
			}
			df1 = df.subset_rows(idx2);
			df2 = v.df.subset_rows(idy2);
		} else {
			df1 = df.subset_rows(idx);
			df2 = v.df.subset_rows(idy);
		}
		if (!df1.cbind(df2)) {
			out.addWarning("could not combine attributes");
		}
	} else {
		if (n < idx.size()) {
			std::vector<size_t> idx2;
			idx2.reserve(n);
			for (size_t i=0; i<n; i++) {
				idx2.push_back( idx[ out.df.iv[0][i] ]);
			}
			df1 = df.subset_rows(idx2);
		} else {
			df1 = df.subset_rows(idx);
		}
	}
	std::vector<std::string> nms = df1.get_names();
	make_unique_names(nms);
	df1.set_names(nms);	
	out.df = df1;
	return out;
}




std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> getRelateFun(const std::string rel) {
	std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> rfun;
	if (rel == "equals") {
		rfun = GEOSEquals_r;
//	} else if (rel == "equalidentical") {
//		rfun = GEOSEqualsIdentical_r;
	} else if (rel == "intersects") {
		rfun = GEOSIntersects_r;
	} else if (rel == "disjoint") {
		rfun = GEOSDisjoint_r;
	} else if (rel == "touches") {
		rfun = GEOSTouches_r;
	} else if (rel == "crosses") {
		rfun = GEOSCrosses_r;
	} else if (rel == "within") {
		rfun = GEOSWithin_r;
	} else if (rel == "contains") {
		rfun = GEOSContains_r;
//	} else if (rel == "containsproperly") {
//		rfun = GEOSContainsProperly_r;
	} else if (rel == "overlaps") {
		rfun = GEOSOverlaps_r;
	} else if (rel == "covers") {
		rfun = GEOSCovers_r;
	} else if (rel == "coveredby") {
		rfun = GEOSCoveredBy_r;
	}
	return rfun;
}


std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> getPrepRelateFun(const std::string rel) {
	std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> rfun;
	if (rel == "intersects") {
		rfun = GEOSPreparedIntersects_r;
	} else if (rel == "disjoint") {
		rfun = GEOSPreparedDisjoint_r;
	} else if (rel == "touches") {
		rfun = GEOSPreparedTouches_r;
	} else if (rel == "crosses") {
		rfun = GEOSPreparedCrosses_r;
	} else if (rel == "within") {
		rfun = GEOSPreparedWithin_r;
	} else if (rel == "contains") {
		rfun = GEOSPreparedContains_r;
	} else if (rel == "containsproperly") {
		rfun = GEOSPreparedContainsProperly_r;
	} else if (rel == "overlaps") {
		rfun = GEOSPreparedOverlaps_r;
	} else if (rel == "covers") {
		rfun = GEOSPreparedCovers_r;
	} else if (rel == "coveredby") {
		rfun = GEOSPreparedCoveredBy_r;
	}
	return rfun;
}


int getRel(std::string &relation) {
	int pattern = 1;
	std::string rel = relation;
	std::transform(rel.begin(), rel.end(), rel.begin(), ::tolower);
	std::vector<std::string> f {"rook", "queen", "intersects", "touches", "crosses", "overlaps", "within", "contains", "covers", "coveredby", "disjoint", "equals"}; // tbd: "equalsexact", "equals"
	if (std::find(f.begin(), f.end(), rel) == f.end()) {
		if (relation.size() != 9) {
			pattern = 2;
		} else {
			std::string r = relation;
			for (size_t i=0; i<9; i++) {
				if (!(r.at(i) == 'T' || r.at(i) == 'F' || r.at(i) == '0' || r.at(i) == '1' || r.at(i) == '2' || r.at(i) == '*')) {
					pattern = 2;
					break;
				}
			}
		}
	} else if (rel == "rook") {
		relation = "F***1****";
	} else if (rel == "queen") {
		relation = "F***T****";
	} else {
		pattern = 0;
		relation = rel;
	}
	return pattern;
}


std::vector<int> SpatVector::pointInPolygon(std::vector<double> &x, std::vector<double> &y) {

	std::vector<int> out;

#ifdef GEOS3120
	size_t ng = size();
	size_t np = x.size();
	out.reserve(np);
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);

	for (size_t i = 0; i < ng; i++) {
		PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, g[i].get()), hGEOSCtxt);
		for (size_t j = 0; j < np; j++) {
			out.push_back( GEOSPreparedIntersectsXY_r(hGEOSCtxt, pr.get(), x[j], y[j]));
		}
	}

# else 
	SpatVector pnts;
	pnts.srs = srs;
	pnts.setPointsGeometry(x, y);
	out = relate(pnts, "intersects", true, true);

# endif

	return out;
}

std::vector<int> SpatVector::relate(SpatVector v, std::string relation, bool prepared, bool index) {

	// this method is redundant with "which_relate")
	std::vector<int> out;
	int pattern = getRel(relation);
	if (pattern == 2) {
		setError("'" + relation + "'" + " is not a valid relate name or pattern");
		return out;
	}
	if ((relation == "FF*FF****") || (relation == "disjoint")) index = false;
	if (relation.substr(0, 5) == "equal") {
		prepared = false;
	}
	
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();

	out.reserve(nx * ny);

	if (!index) {
		out.reserve(nx*ny);
		if (pattern == 1) {
			for (size_t i = 0; i < nx; i++) {
				for (size_t j = 0; j < ny; j++) {
					out.push_back( GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[j].get(), relation.c_str()));
				}
			}
		} else if (prepared) {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i = 0; i < nx; i++) {
				PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
				for (size_t j = 0; j < ny; j++) {
					out.push_back( relFun(hGEOSCtxt, pr.get(), y[j].get()));
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> relFun = getRelateFun(relation);
			for (size_t i = 0; i < nx; i++) {
				for (size_t j = 0; j < ny; j++) {
					out.push_back( relFun(hGEOSCtxt, x[i].get(), y[j].get()));
				}
			}
		}
	} else { // use spatial index
		out.resize(nx*ny);
		std::vector<size_t> items(y.size());
		TreePtr tree1 = geos_ptr(GEOSSTRtree_create_r(hGEOSCtxt, 10), hGEOSCtxt);
		for (size_t i = 0; i < y.size(); i++) {
			items[i] = i;
			if (! GEOSisEmpty_r(hGEOSCtxt, y[i].get()))
				GEOSSTRtree_insert_r(hGEOSCtxt, tree1.get(), y[i].get(), &(items[i]));
		}

		if (pattern == 1) {
			for (size_t i = 0; i < nx; i++) {
				// pre-select y's using tree:
				std::vector<size_t> tree_sel, sel;
				if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
					GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), x[i].get(), callbck, &tree_sel);
				}
				for (size_t j = 0; j < tree_sel.size(); j++) {
					if (GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[tree_sel[j]].get(), relation.c_str())) {
						out[i * nx + tree_sel[j]] = 1; //.push_back(tree_sel[j]);
					}
				}
			}
		} else if (prepared) {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i=0; i<nx; i++) {
				std::vector<size_t> tree_sel, sel;
				if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
					GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), x[i].get(), callbck, &tree_sel);
				}
				if (! tree_sel.empty()) {
					PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
					for (size_t j=0; j < tree_sel.size(); j++) {
						int r = relFun(hGEOSCtxt, pr.get(), y[tree_sel[j]].get());
						if (r == 2) {
							setError("an exception occurred");
							return out;
						}
						out[i*ny + tree_sel[j]] = r;
					}
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> relFun = getRelateFun(relation);
			for (size_t i=0; i<nx; i++) {
				std::vector<size_t> tree_sel, sel;
				if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
					GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), x[i].get(), callbck, &tree_sel);
				}
				if (! tree_sel.empty()) {
					for (size_t j=0; j < tree_sel.size(); j++) {
						int r = relFun(hGEOSCtxt, x[i].get(), y[tree_sel[j]].get());
						if (r == 2) {
							setError("an exception occurred");
							return out;
						}
						out[i * ny + tree_sel[j]] = r;
					}
				}
			}
		}
	}
	geos_finish(hGEOSCtxt);
	return out;
}


std::vector<std::vector<double>> SpatVector::which_relate(SpatVector v, std::string relation, bool narm) {

	bool index=true;
	if ((relation == "FF*FF****") || (relation == "disjoint")) index = false;

	std::vector<std::vector<double>> out(2);
	int pattern = getRel(relation);
	if (pattern == 2) {
		setError("'" + relation + "'" + " is not a valid relate name or pattern");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();
	out[0].reserve(nx * 1.5);
	out[1].reserve(nx * 1.5);


	if (relation.substr(0, 5) == "equal") {
		std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> relFun = getRelateFun(relation);
		for (size_t i = 0; i < nx; i++) {
			bool none = !narm;
			for (size_t j = 0; j < ny; j++) {
				if ( relFun(hGEOSCtxt, x[i].get(), y[j].get())) {
					out[0].push_back(i);
					out[1].push_back(j);
					none = false;
				}
			}
			if (none) {
				out[0].push_back(i);
				out[1].push_back(NAN);
			}
		}
	} else if (!index) {
		if (pattern == 1) {
			for (size_t i = 0; i < nx; i++) {
				bool none = !narm;
				for (size_t j = 0; j < ny; j++) {
					if ( GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[j].get(), relation.c_str())) {
						out[0].push_back(i);
						out[1].push_back(j);
						none = false;
					}
				}
				if (none) {
					out[0].push_back(i);
					out[1].push_back(NAN);
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i = 0; i < nx; i++) {
				bool none = !narm;
				PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
				for (size_t j = 0; j < ny; j++) {
					if (relFun(hGEOSCtxt, pr.get(), y[j].get())) {
						out[0].push_back(i);
						out[1].push_back(j);
						none = false;
					}
				}
				if (none) {
					out[0].push_back(i);
					out[1].push_back(NAN);
				}
			}
		} 
	} else {
		std::vector<size_t> items(y.size());
		TreePtr tree1 = geos_ptr(GEOSSTRtree_create_r(hGEOSCtxt, 10), hGEOSCtxt);
		for (size_t i = 0; i < y.size(); i++) {
			items[i] = i;
			if (! GEOSisEmpty_r(hGEOSCtxt, y[i].get()))
				GEOSSTRtree_insert_r(hGEOSCtxt, tree1.get(), y[i].get(), &(items[i]));
		}

		if (pattern == 1) {
			for (size_t i = 0; i < nx; i++) {
				// pre-select y's using tree:
				std::vector<size_t> tree_sel, sel;
				if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
					GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), x[i].get(), callbck, &tree_sel);
				}
				bool none = !narm;
				for (size_t j = 0; j < tree_sel.size(); j++) {
					if (GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[tree_sel[j]].get(), relation.c_str())) {
						out[0].push_back(i);
						out[1].push_back(tree_sel[j]);
						none = false;
					}
				}
				if (none) {
					out[0].push_back(i);
					out[1].push_back(NAN);
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i=0; i<nx; i++) {
				// pre-select y
				std::vector<size_t> tree_sel, sel;
				if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
					GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), x[i].get(), callbck, &tree_sel);
				}
				bool none = !narm;
				if (! tree_sel.empty()) {
					PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
					for (size_t j=0; j < tree_sel.size(); j++) {
						if (relFun(hGEOSCtxt, pr.get(), y[tree_sel[j]].get())) {
							out[0].push_back(i);
							out[1].push_back(tree_sel[j]);
							none = false;
						}
					}
				}
				if (none) {
					out[0].push_back(i);
					out[1].push_back(NAN);
				}
			}
		}
	}
	geos_finish(hGEOSCtxt);
	return out;
}



std::vector<std::vector<double>> SpatVector::which_relate(std::string relation, bool narm) {

	bool index=true;
	if ((relation == "FF*FF****") || (relation == "disjoint")) index = false;

	std::vector<std::vector<double>> out(2);
	int pattern = getRel(relation);
	if (pattern == 2) {
		setError("'" + relation + "'" + " is not a valid relate name or pattern");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	size_t nx = size();
	out[0].reserve(nx * 1.5);
	out[1].reserve(nx * 1.5);

	if (relation.substr(0, 5) == "equal") {
		std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> relFun = getRelateFun(relation);
		for (size_t i = 0; i < nx; i++) {
			bool none = !narm;
			for (size_t j = 0; j < nx; j++) {
				if ( relFun(hGEOSCtxt, x[i].get(), x[j].get())) {
					out[0].push_back(i);
					out[1].push_back(j);
					none = false;
				}
			}
			if (none) {
				out[0].push_back(i);
				out[1].push_back(NAN);
			}
		}
	} else if (!index) {
		if (pattern == 1) {
			for (size_t i=0; i<nx; i++) {
				bool none = !narm;
				for (size_t j=0; j<nx; j++) {
					if ( GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), x[j].get(), relation.c_str())) {
						out[0].push_back(i);
						out[1].push_back(j);
						none = false;
					}
				}
				if (none) {
					out[0].push_back(i);
					out[1].push_back(NAN);
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i = 0; i < nx; i++) {
				bool none = !narm;
				PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
				for (size_t j = 0; j < nx; j++) {
					if (relFun(hGEOSCtxt, pr.get(), x[j].get())) {
						out[0].push_back(i);
						out[1].push_back(j);
						none = false;
					}
				}
				if (none) {
					out[0].push_back(i);
					out[1].push_back(NAN);
				}
			}
		} 
	} else {
		std::vector<size_t> items(x.size());
		TreePtr tree1 = geos_ptr(GEOSSTRtree_create_r(hGEOSCtxt, 10), hGEOSCtxt);
		for (size_t i = 0; i < x.size(); i++) {
			items[i] = i;
			if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
				GEOSSTRtree_insert_r(hGEOSCtxt, tree1.get(), x[i].get(), &(items[i]));
			}
		}

		if (pattern == 1) {
			for (size_t i = 0; i < nx; i++) {
				std::vector<size_t> tree_sel, sel;
				if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
					GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), x[i].get(), callbck, &tree_sel);
				}
				bool none = !narm;
				for (size_t j = 0; j < tree_sel.size(); j++) {
					if (GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), x[tree_sel[j]].get(), relation.c_str())) {
						out[0].push_back(i);
						out[1].push_back(tree_sel[j]);
						none = false;
					}
				}
				if (none) {
					out[0].push_back(i);
					out[1].push_back(NAN);
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i=0; i<nx; i++) {
				std::vector<size_t> tree_sel, sel;
				if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
					GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), x[i].get(), callbck, &tree_sel);
				}
				bool none = !narm;
				if (! tree_sel.empty()) {
					PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
					for (size_t j=0; j < tree_sel.size(); j++) {
						if (relFun(hGEOSCtxt, pr.get(), x[tree_sel[j]].get())) {
							out[0].push_back(i);
							out[1].push_back(tree_sel[j]);
							none = false;
						}
					}
				}
				if (none) {
					out[0].push_back(i);
					out[1].push_back(NAN);
				}
			}
		}
	}
	geos_finish(hGEOSCtxt);
	return out;
}



std::vector<int> SpatVector::relate(std::string relation, bool symmetrical) {
	// this method is redundant with "which_relate")

	std::vector<int> out;
	int pattern = getRel(relation);
	if (pattern == 2) {
		setError("'" + relation + "'" + " is not a valid relate name or pattern");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);

	if (symmetrical) {
		size_t s = size();
		size_t n = ((s-1) * s)/2;
		out.reserve(n);
		if (pattern == 1) {
			for (size_t i=0; i<(s-1); i++) {
				for (size_t j=(i+1); j<s; j++) {
					out.push_back( GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), x[j].get(), relation.c_str()));
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i=0; i<(s-1); i++) {
				PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
				for (size_t j=(i+1); j<s; j++) {
					out.push_back( relFun(hGEOSCtxt, pr.get(), x[j].get()));
				}
			}
		}
	} else {
		size_t nx = size();
		out.reserve(nx*nx);
		if (pattern == 1) {
			for (size_t i = 0; i < nx; i++) {
				for (size_t j = 0; j < nx; j++) {
					out.push_back( GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), x[j].get(), relation.c_str()));
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i = 0; i < nx; i++) {
				PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
				for (size_t j = 0; j < nx; j++) {
					out.push_back( relFun(hGEOSCtxt, pr.get(), x[j].get()));
				}
			}
		}
	}
	geos_finish(hGEOSCtxt);
	return out;
}


/*
std::vector<std::vector<double>> SpatVector::which_relate(SpatVector v, std::string relation) {

	std::vector<std::vector<double>> out(2);
	int pattern = getRel(relation);
	if (pattern == 2) {
		setError("'" + relation + "'" + " is not a valid relate name or pattern");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();
	out[0].reserve(nx * 1.5);
	out[1].reserve(nx * 1.5);
	if (pattern == 1) {
		for (size_t i=0; i<nx; i++) {
			bool none = true;
			for (size_t j=0; j<ny; j++) {
				if (GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[j].get(), relation.c_str())) {
					out[0].push_back(i);
					out[1].push_back(j);
					none = false;
				}
			}
			if (none) {
				out[0].push_back(i);
				out[1].push_back(NAN);
			}
		}
	} else {
		std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
		for (size_t i=0; i<nx; i++) {
			PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
			bool none = true;
			for (size_t j=0; j<ny; j++) {
				if ( relFun(hGEOSCtxt, pr.get(), y[j].get())) {
					out[0].push_back(i);
					out[1].push_back(j);
					none = false;
				}
			}
			if (none) {
				out[0].push_back(i);
				out[1].push_back(NAN);
			}
		}
	}
	geos_finish(hGEOSCtxt);
	return out;
}
*/



/*
std::vector<int> SpatVector::relateFirst(SpatVector v, std::string relation) {

	int pattern = getRel(relation);
	if (pattern == 2) {
		setError("'" + relation + "'" + " is not a valid relate name or pattern");
		std::vector<int> out;
		return out;
	}
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();
	std::vector<int> out(nx, -1);
	if (pattern == 1) {
		for (size_t i = 0; i < nx; i++) {
			for (size_t j = 0; j < ny; j++) {
				if (GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[j].get(), relation.c_str())) {
					out[i] = j;
					continue;
				}
			}
		}
	} else {
		std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);

		for (size_t i = 0; i < nx; i++) {
			PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
			for (size_t j = 0; j < ny; j++) {
				if (relFun(hGEOSCtxt, pr.get(), y[j].get())) {
					out[i] = j;
					continue;
				}
			}
		}
	}
	geos_finish(hGEOSCtxt);
	return out;
}
*/


std::vector<int> SpatVector::relateFirst(SpatVector v, std::string relation) {

	bool index=true;
	if ((relation == "FF*FF****") || (relation == "disjoint")) index = false;

	std::vector<int> out;
	int pattern = getRel(relation);
	if (pattern == 2) {
		setError("'" + relation + "'" + " is not a valid relate name or pattern");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();
	out.resize(nx, -1);

	if (!index) {
		if (pattern == 1) {
			for (size_t i = 0; i < nx; i++) {
				for (size_t j = 0; j < ny; j++) {
					if ( GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[j].get(), relation.c_str())) {
						out[i] = j;
						continue;
					}
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i = 0; i < nx; i++) {
				PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
				for (size_t j = 0; j < ny; j++) {
					if (relFun(hGEOSCtxt, pr.get(), y[j].get())) {
						out[i] = j;
						continue;
					}
				}
			}
		} 
	} else {
		std::vector<size_t> items(y.size());
		TreePtr tree1 = geos_ptr(GEOSSTRtree_create_r(hGEOSCtxt, 10), hGEOSCtxt);
		for (size_t i = 0; i < y.size(); i++) {
			items[i] = i;
			if (! GEOSisEmpty_r(hGEOSCtxt, y[i].get()))
				GEOSSTRtree_insert_r(hGEOSCtxt, tree1.get(), y[i].get(), &(items[i]));
		}

		if (pattern == 1) {
			for (size_t i = 0; i < nx; i++) {
				std::vector<size_t> tree_sel, sel;
				if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
					GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), x[i].get(), callbck, &tree_sel);
				}
				for (size_t j = 0; j < tree_sel.size(); j++) {
					if (GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[tree_sel[j]].get(), relation.c_str())) {
						out[i] = tree_sel[j];
						continue;
					}
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i=0; i<nx; i++) {
				std::vector<size_t> tree_sel, sel;
				if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
					GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), x[i].get(), callbck, &tree_sel);
				}
				if (! tree_sel.empty()) {
					PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
					for (size_t j=0; j < tree_sel.size(); j++) {
						if (relFun(hGEOSCtxt, pr.get(), y[tree_sel[j]].get())) {
							out[i] = tree_sel[j];
							continue;
						}
					}
				}
			}
		}
	}
	geos_finish(hGEOSCtxt);
	return out;
}


/*
std::vector<bool> SpatVector::is_related(SpatVector v, std::string relation) {

	std::vector<bool> out;
	int pattern = getRel(relation);
	if (pattern == 2) {
		setError("'" + relation + "'" + " is not a valid relate name or pattern");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();
	out.resize(nx, false);
	if (pattern == 1) {
		for (size_t i = 0; i < nx; i++) {
			for (size_t j = 0; j < ny; j++) {
				bool isrel = GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[j].get(), relation.c_str());
				if (isrel) {
					out[i] = true;
					continue;
				}
			}
		}
	} else {
		std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
		for (size_t i = 0; i < nx; i++) {
			PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
			for (size_t j = 0; j < ny; j++) {
				bool isrel = relFun(hGEOSCtxt, pr.get(), y[j].get());
				if (isrel) {
					out[i] = true;
					continue;
				}
			}
		}
	}
	geos_finish(hGEOSCtxt);

	return out;
}
*/



std::vector<bool> SpatVector::is_related(SpatVector v, std::string relation) {

	bool index=true;
	if ((relation == "FF*FF****") || (relation == "disjoint")) index = false;

	std::vector<bool> out;
	int pattern = getRel(relation);
	if (pattern == 2) {
		setError("'" + relation + "'" + " is not a valid relate name or pattern");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();
	out.resize(nx, false);

	if (!index) {
		if (pattern == 1) {
			for (size_t i = 0; i < nx; i++) {
				for (size_t j = 0; j < ny; j++) {
					if ( GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[j].get(), relation.c_str())) {
						out[i] = true;
						continue;
					}
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i = 0; i < nx; i++) {
				PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
				for (size_t j = 0; j < ny; j++) {
					if (relFun(hGEOSCtxt, pr.get(), y[j].get())) {
						out[i] = true;
						continue;
					}
				}
			}
		} 
	} else {
		std::vector<size_t> items(y.size());
		TreePtr tree1 = geos_ptr(GEOSSTRtree_create_r(hGEOSCtxt, 10), hGEOSCtxt);
		for (size_t i = 0; i < y.size(); i++) {
			items[i] = i;
			if (! GEOSisEmpty_r(hGEOSCtxt, y[i].get()))
				GEOSSTRtree_insert_r(hGEOSCtxt, tree1.get(), y[i].get(), &(items[i]));
		}

		if (pattern == 1) {
			for (size_t i = 0; i < nx; i++) {
				std::vector<size_t> tree_sel, sel;
				if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
					GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), x[i].get(), callbck, &tree_sel);
				}
				for (size_t j = 0; j < tree_sel.size(); j++) {
					if (GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[tree_sel[j]].get(), relation.c_str())) {
						out[i] = true;
						continue;
					}
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i=0; i<nx; i++) {
				// pre-select y
				std::vector<size_t> tree_sel, sel;
				if (! GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
					GEOSSTRtree_query_r(hGEOSCtxt, tree1.get(), x[i].get(), callbck, &tree_sel);
				}
				if (! tree_sel.empty()) {
					PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
					for (size_t j=0; j < tree_sel.size(); j++) {
						if (relFun(hGEOSCtxt, pr.get(), y[tree_sel[j]].get())) {
							out[i] = true;
							continue;
						}
					}
				}
			}
		}
	}
	geos_finish(hGEOSCtxt);
	return out;
}


std::vector<size_t> SpatVector::equals_exact(SpatVector v, double tol) {
	std::vector<size_t> out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();
	out.reserve(nx*ny);
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			out.push_back( GEOSEqualsExact_r(hGEOSCtxt, x[i].get(), y[j].get(), tol));
		}
	}
	geos_finish(hGEOSCtxt);
	return out;
}


std::vector<size_t> SpatVector::equals_exact(bool symmetrical, double tol) {
	std::vector<size_t> out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);

	if (symmetrical) {
		size_t s = size();
		size_t n = ((s-1) * s)/2;
		out.reserve(n);
		for (size_t i=0; i<(s-1); i++) {
			for (size_t j=(i+1); j<s; j++) {
				out.push_back( GEOSEqualsExact_r(hGEOSCtxt, x[i].get(), x[j].get(), tol));
			}
		}
	} else {
		size_t nx = size();
		out.reserve(nx*nx);
		for (size_t i = 0; i < nx; i++) {
			for (size_t j = 0; j < nx; j++) {
				out.push_back( GEOSEqualsExact_r(hGEOSCtxt, x[i].get(), x[j].get(), tol));
			}
		}
	}
	geos_finish(hGEOSCtxt);
	return out;
}


SpatVector SpatVector::mask(SpatVector x, bool inverse) {
	std::vector<bool> b = is_related(x, "intersects");
	if (inverse) {
		for (size_t i=0; i<b.size(); i++) {
			b[i] = !b[i];
		}
	}
	std::vector<size_t> r;
	r.reserve(b.size());
	for (size_t i=0; i<b.size(); i++) {
		if (b[i]) r.push_back(i);
	}
	return subset_rows(r);
}


typedef int (* dist_fn)(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *, double *);

bool get_dist_fun(dist_fn &f, std::string s) {
	if ((s == "Euclidean") || s.empty())
		f = GEOSDistance_r;
	else if (s == "Hausdorff")
		f = GEOSHausdorffDistance_r;
#ifdef GEOS380
	else if (s == "Frechet")
		f = GEOSFrechetDistance_r;
		// GEOSFrechetDistanceDensify_r
#endif
	else {
		return false;
	}
	return true;
}


std::vector<double> SpatVector::geos_distance(SpatVector v, bool parallel, std::string fun, double m) {

	std::vector<double> out;

	dist_fn distfun;
	if (!get_dist_fun(distfun, fun)) {
		setError("invalid distance function");
		return out;
	}

	size_t nx = size();
	size_t ny = v.size();

	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);

	double d;

	if (parallel) {
		bool nyone = false;
		if (nx != ny) {
			if (ny == 1) {
				nyone = true;
			} else if ((nx == 1) && (ny > 1)) {
				std::swap(x, y);
				std::swap(nx, ny);
				nyone = true;
			} else {
				setError("SpatVectors have different lengths");
				return out;
			}
		}
		if (nyone) {
			out.reserve(nx);
			for (size_t i = 0; i < nx; i++) {
				if ( distfun(hGEOSCtxt, x[i].get(), y[0].get(), &d)) {
					out.push_back(d);
				} else {
					out.push_back(NAN);
				}
			}
		} else {
			out.reserve(nx);
			for (size_t i = 0; i < nx; i++) {
				if ( distfun(hGEOSCtxt, x[i].get(), y[i].get(), &d)) {
					out.push_back(d);
				} else {
					out.push_back(NAN);
				}
			}
		}
	} else {
		out.reserve(nx*ny);
		for (size_t i = 0; i < nx; i++) {
			for (size_t j = 0; j < ny; j++) {
				if ( distfun(hGEOSCtxt, x[i].get(), y[j].get(), &d)) {
					out.push_back(d);
				} else {
					out.push_back(NAN);
				}
			}
		}
	}
	geos_finish(hGEOSCtxt);
	if (m != 1) {
		for (double &d : out) d *= m;
	}
	return out;
 }

std::vector<double> SpatVector::geos_distance(bool sequential, std::string fun, double m) {

	std::vector<double> out;
	dist_fn distfun;
	if (!get_dist_fun(distfun, fun)) {
		setError("invalid distance function");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	size_t s = size();
	double d;
	if (sequential) {
		out.reserve(s);
		out.push_back(0);
		for (size_t i=0; i<(s-1); i++) {
			if ( distfun(hGEOSCtxt, x[i].get(), x[i+1].get(), &d)) {
				out.push_back(d);
			} else {
				out.push_back(NAN);
			}
		}
	} else {
		out.reserve((s-1) * s / 2);
		for (size_t i=0; i<(s-1); i++) {
			for (size_t j=(i+1); j<s; j++) {
				if ( distfun(hGEOSCtxt, x[i].get(), x[j].get(), &d)) {
					out.push_back(d);
				} else {
					out.push_back(NAN);
				}
			}
		}
	}
	if (s == 1) {
		out.push_back(0);
	}
	geos_finish(hGEOSCtxt);
	if (m != 1) {
		for (double &d : out) d *= m;
	}
	return out;
 }


SpatVector SpatVector::unite(SpatVector v) {

	SpatVector out;
	std::string gtp = type();
	if (gtp != v.type()) {
		out.setError("cannot unite different geom types");
		return out;
	}

	if (gtp != "polygons") {
		// should make records unique here. This is now done in R
		return append(v, true);
	}

	out = intersect(v, true);
	if (out.hasError()) {
		return out;
	}
	if (out.nrow() == 0) {
		return append(v, true);
	}

	if (gtp != out.type()) {
		SpatVector empty;
		out = empty;
	}
	SpatVector sdif = symdif(v);
	if (sdif.hasError()) {
		return sdif;
	}
	if ((sdif.nrow() > 0) && (gtp == sdif.type())) {	
		return sdif.append(out, true);
	}
	return out;
}


SpatVector SpatVector::unite() {
	int n = size();

	std::vector<long> x(1, 1);
	SpatDataFrame d;
	d.add_column(x, "id_1");
	SpatVector out = subset_rows(0);
	out.df = d;
	for (int i=1; i<n; i++) {
		std::string name = "id_" + std::to_string(i+1);
		SpatDataFrame d;
		d.add_column(x, name);
		SpatVector r = subset_rows(i);
		r.df = d;
		out = out.unite(r);
		if (out.hasError()) {
			return out;
		}
	}

	for (size_t i=0; i<out.df.iv.size(); i++) {
		for (size_t j=0; j<out.df.iv[i].size(); j++) {
			if (out.df.iv[i][j] != 1) {
				out.df.iv[i][j] = 0;
			}
		}
	}

	return out;
}



SpatVector SpatVector::symdif(SpatVector v) {

	if ((type() != "polygons") || (v.type() != "polygons")) {
		SpatVector out;
		out.setError("expected two polygon geometries");
		return out;
	}
	SpatVector out = erase(v);
	if (out.hasError()) {
		return out;
	}
	SpatVector ve = v.erase(*this);
	if (ve.hasError()) {
		return ve;
	}
	if ((out.type() == "polygons") && (ve.type() == "polygons")) {
		return out.append(ve, true);
	} else if (ve.type() == "polygons") {
		return ve;
	} else if (out.type() == "polygons") {
		return out;
	} else {
		SpatVector empty;
		return empty;
	} 

/*
	SpatVector out;
	out.srs = srs;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	std::vector<GeomPtr> result;
	std::vector<unsigned> ids;
	ids.reserve(size());
	size_t nx = size();
	size_t ny = v.size();


	for (size_t i = 0; i < nx; i++) {
		GEOSGeometry* geom = x[i].get();
		for (size_t j = 0; j < ny; j++) {
			geom = GEOSDifference_r(hGEOSCtxt, geom, y[j].get());
			if (geom == NULL) {
				out.setError("GEOS exception");
				geos_finish(hGEOSCtxt);
				return(out);
			}
			if (GEOSisEmpty_r(hGEOSCtxt, geom)) {
				break;
			}
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, geom)) {
			result.push_back(geos_ptr(geom, hGEOSCtxt));
			ids.push_back(i);
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, geom);
		}
	}

	if (result.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);
		out = coll.get(0);
		out.srs = srs;
		out.df = df.subset_rows(ids);
	}
	geos_finish(hGEOSCtxt);
	if (!srs.is_same(v.srs, true)) {
		out.addWarning("different crs");
	}

	return out.append(v, true);
	*/
}




SpatVector SpatVector::cover(SpatVector v, bool identity, bool expand) {
	if (v.srs.is_empty()) {
		v.srs = srs;
	}
	SpatVector out = erase(v);
	if (identity) {
		SpatVector insect = intersect(v, true);
		out = out.append(insect, true);
		if (expand) {
			v = v.erase(insect);
			out = out.append(v, true);
		}
	} else {
		if (!expand) {
			v = v.crop(*this);
		}
		out = out.append(v, true);
	}
	return out;
}


SpatVector SpatVector::erase_agg(SpatVector v) {

	if ((nrow()==0) || (v.nrow() == 0)) {
		return(*this);
	}


	if ((type() == "points") || (v.type() == "points")) {
		std::vector<bool> b = is_related(v, "intersects");
		std::vector<size_t> r;
		r.reserve(b.size());
		for (size_t i=0; i < b.size(); i++) {
			if (!b[i]) r.push_back(i);
		}
		return subset_rows(r);
	}

	SpatVector out;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);

// this approach is nicer than the below in ::erase
// but it fails if polys overlap
//	v = v.aggregate(false);

// so we do
	v = v.aggregate(true);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	std::vector<size_t> rids;
	size_t nx = size();
	std::vector<GeomPtr> result;

	for (size_t i = 0; i < nx; i++) {
		GEOSGeometry* geom = GEOSDifference_r(hGEOSCtxt, x[i].get(), y[0].get());
		if (geom == NULL) {
			out.setError("GEOS exception");
			geos_finish(hGEOSCtxt);
			return(out);
		}
		if (GEOSisEmpty_r(hGEOSCtxt, geom)) {
			GEOSGeom_destroy_r(hGEOSCtxt, geom);
		} else {
			result.push_back(geos_ptr(geom, hGEOSCtxt));
			rids.push_back(i);
		}
	}
	if (!result.empty()) {
		std::vector<long> ids;
		SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt, ids, true, false);
		out = coll.get(0);
		out.srs = srs;
		out.df = df.subset_rows(rids);
	} else {
		std::vector<long> none(1, -1);
		out = subset_rows(none);
	}
	geos_finish(hGEOSCtxt);
	if (!srs.is_same(v.srs, true)) {
		out.addWarning("different crs");
	}
	return out;
}



SpatVector SpatVector::erase(SpatVector v) {

	if ((nrow()==0) || (v.nrow() == 0)) {
		return(*this);
	}

	if ((type() == "points") || (v.type() == "points")) {
		std::vector<bool> b = is_related(v, "intersects");
		std::vector<size_t> r;
		r.reserve(b.size());
		for (size_t i=0; i<b.size(); i++) {
			if (!b[i]) r.push_back(i);
		}
		return subset_rows(r);
	}

	SpatVector out;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();
	std::vector<long> rids;
	rids.reserve(nx);

	for (size_t i = 0; i < nx; i++) {
		bool good=true;
		for (size_t j = 0; j < ny; j++) {
			GEOSGeometry* geom = GEOSDifference_r(hGEOSCtxt, x[i].get(), y[j].get());
			if (geom == NULL) {
				out.setError("GEOS exception");
				geos_finish(hGEOSCtxt);
				return(out);
			}
			if (GEOSisEmpty_r(hGEOSCtxt, geom)) {
				GEOSGeom_destroy_r(hGEOSCtxt, geom);
				good = false;
				break;
			}
			x[i] = geos_ptr(geom, hGEOSCtxt);
		}
		if (good) rids.push_back(i);
	}

	if (rids.empty()) {
		std::vector<long> none(1, -1);
		out = subset_rows(none);
	} else {
		SpatVectorCollection coll = coll_from_geos(x, hGEOSCtxt);
		out = coll.get(0);
		out.srs = srs;
		out.df = df;
		if (rids.size() != out.nrow()) {
			out = out.subset_rows(rids);
		}
	}
	geos_finish(hGEOSCtxt);
	if (!srs.is_same(v.srs, true)) {
		out.addWarning("different crs");
	}
	return out;
}



/*
SpatVector SpatVector::erase(SpatVector v) {

	if ((type() == "points") || (v.type() == "points")) {
		std::vector<int> b = relateFirst(v, "intersects");
		std::vector<unsigned> r;
		r.reserve(b.size());
		for (size_t i=0; i<b.size(); i++) {
			if (b[i] == -1) r.push_back(i);
		}
		return subset_rows(r);
	}

	SpatVector out;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	std::vector<unsigned> rids;
	size_t nx = size();
	size_t ny = v.size();

	for (size_t i = 0; i < nx; i++) {
		//GEOSGeometry* geom = x[i].get();
		for (size_t j = 0; j < ny; j++) {
			GEOSGeometry* geom = GEOSDifference_r(hGEOSCtxt, x[i].get(), y[j].get());
			if (geom == NULL) {
				out.setError("GEOS exception");
				geos_finish(hGEOSCtxt);
				return(out);
			}
			if (GEOSisEmpty_r(hGEOSCtxt, geom)) {
				GEOSGeom_destroy_r(hGEOSCtxt, geom);
				rids.push_back(i);
				break;
			}
			x[i] = geos_ptr(geom, hGEOSCtxt);
		}
	}

	if (rids.size() < nx) {
		SpatVectorCollection coll = coll_from_geos(x, hGEOSCtxt);
		out = coll.get(0);
		out.df = df;
		out.df.remove_rows(rids);
	}
	geos_finish(hGEOSCtxt);

	if (!srs.is_same(v.srs, true)) {
		out.addWarning("different crs");
	}
	out.srs = srs;

	return out;
}

*/

SpatVector SpatVector::erase(bool sequential) {

	if (nrow()==0) {
		return(*this);
	}

	SpatVector out;

	if (type() != "polygons") {
		out.setError("not polygons");
		return out;
	}
	size_t n = size();
	if (n < 2) {
		return *this;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<size_t> rids;

	if (sequential) {
		for (size_t i = 0; i < (n-1); i++) {
			for (size_t j = (i+1); j < n; j++) {
				GEOSGeometry* geom = GEOSDifference_r(hGEOSCtxt, x[i].get(), x[j].get());
				if (geom == NULL) {
					out.setError("GEOS exception");
					geos_finish(hGEOSCtxt);
					return(out);
				} else if (GEOSisEmpty_r(hGEOSCtxt, geom)) {
					GEOSGeom_destroy_r(hGEOSCtxt, geom);
					rids.push_back(i);
					break;
				} else {
					x[i] = geos_ptr(geom, hGEOSCtxt);
				}
			}
		}
	} else {
		std::vector<GeomPtr> y = geos_geoms(this, hGEOSCtxt);
		for (size_t i=0; i<n; i++) {
			for (size_t j=0; j<n; j++) {
				if (j == i) continue;
				GEOSGeometry* geom = GEOSDifference_r(hGEOSCtxt, x[i].get(), y[j].get());
				if (geom == NULL) {
					out.setError("GEOS exception");
					geos_finish(hGEOSCtxt);
					return(out);
				} else if (GEOSisEmpty_r(hGEOSCtxt, geom)) {
					GEOSGeom_destroy_r(hGEOSCtxt, geom);
					rids.push_back(i);
					//Rcpp::Rcout << i << std::endl;
					break;
				} else {
					x[i] = geos_ptr(geom, hGEOSCtxt);
				}
			}
		}
	}

	SpatVectorCollection coll = coll_from_geos(x, hGEOSCtxt);
	out = coll.get(0);
	out.srs = srs;
	out.df = df;
	out.df.remove_rows(rids);
	//SpatVector last = subset_rows(n-1);
	//out = out.append(last, true);

	geos_finish(hGEOSCtxt);
	return out;
}


SpatVector SpatVector::gaps() {
	SpatVector out;

	if (type() != "polygons") {
		out.setError("not polygons");
		return out;
	}

	size_t n = size();
	if (n < 2) {
		out.srs = srs;
		return out;
	}
	out = aggregate(true);
	return out.get_holes();
/*
	SpatExtent e = extent;
	e.xmin -= 11;
	e.xmax += 10;
	e.ymin -= 10;
	e.ymax += 10;
	SpatVector p(e, "");

	p = p.erase(*this);
	p = p.disaggregate(false);
	double exmin = e.xmin + 1;
	unsigned j;
	for (size_t i=0; i<p.size(); i++) {
		if (p.geoms[i].extent.xmin < exmin) {
			j = i;
			break;
		}
	}
	std::vector<unsigned> r(1, j);
	p.srs = srs;
	return p.remove_rows(r);
*/
}



// also use GEOSPreparedNearestPoints_r()

SpatVector SpatVector::nearest_point(SpatVector v, bool parallel, const std::string method) {
	SpatVector out;


	if ((size() == 0) || v.empty()) {
		out.setError("empty SpatVecor(s)");
		return out;
	}
	if (!srs.is_equal(v.srs)) {
		out.setError("CRSs do not match");
		return out;
	}
	bool lonlat = is_lonlat();
	if (lonlat) {
		std::vector<std::string> methods = {"geo", "cosine", "haversine"};
		if (!is_in_vector(method, methods)) {
			out.setError("invalid distance method");
		}
	}
	out.srs = srs;

	
	if (lonlat) {
		SpatVector x = *this;
		if (x.type() != "points") {
			x = x.densify(1000, false, false);
			x = x.as_points(false);
		}
		if (v.type() != "points") {
			v = v.densify(1000, false, false);
			v = v.as_points(false);
		}
		std::vector<double> nlon, nlat, dist;
		std::vector<long> id;
		std::vector<std::vector<double>> p = coordinates();
		std::vector<std::vector<double>> pv = v.coordinates();
		nearest_lonlat(id, dist, nlon, nlat, p[0], p[1], pv[0], pv[1], method);
		out.setPointsGeometry(nlon, nlat);
		std::vector<long> fromid(id.size());
		std::iota(fromid.begin(), fromid.end(), 0);
		out.df.add_column(fromid, "from_id");
		out.df.add_column(p[0], "from_x");
		out.df.add_column(p[1], "from_y");
		out.df.add_column(id, "to_id");
		out.df.add_column(nlon, "to_x");
		out.df.add_column(nlat, "to_y");
		out.df.add_column(dist, "distance");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	if (parallel) {
		if ((size() != v.size())) {
			out.setError("SpatVecors do not have the same size");
			return out;
		}
		std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
		std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
		std::vector<GeomPtr> b(size());
		for (size_t i=0; i < x.size(); i++) {
			GEOSCoordSequence* csq = GEOSNearestPoints_r(hGEOSCtxt, x[i].get(), y[i].get());
			GEOSGeometry* geom = GEOSGeom_createLineString_r(hGEOSCtxt, csq);
			b[i] = geos_ptr(geom, hGEOSCtxt);
		}
		out = vect_from_geos(b, hGEOSCtxt, "lines");

	} else {
		SpatVector mp = v.aggregate(false);
		std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
		std::vector<GeomPtr> y = geos_geoms(&mp, hGEOSCtxt);
		std::vector<GeomPtr> b(size());
		for (size_t i = 0; i < x.size(); i++) {
			GEOSCoordSequence* csq = GEOSNearestPoints_r(hGEOSCtxt, x[i].get(), y[0].get());
			GEOSGeometry* geom = GEOSGeom_createLineString_r(hGEOSCtxt, csq);
			b[i] = geos_ptr(geom, hGEOSCtxt);
		}
		out = vect_from_geos(b, hGEOSCtxt, "lines");
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	return out;
}

SpatVector SpatVector::nearest_point(const std::string method) {
	SpatVector out;
	if ((size() == 0)) {
		out.addWarning("empty SpatVecor");
		return out;
	}
	if ((size() == 1)) {
		out.setError("single geometry");
		return out;
	}
	size_t n = size();
	out.srs = srs;


	bool lonlat = is_lonlat();
	if (lonlat) {
		std::vector<std::string> methods = {"geo", "cosine", "haversine"};
		if (!is_in_vector(method, methods)) {
			out.setError("invalid distance method");
		}
	}


	if (lonlat) {
		if (type() == "points") {
			std::vector<double> nlon, nlat, dist;
			std::vector<long> id;
			std::vector<std::vector<double>> p = coordinates();
			nearest_lonlat_self(id, dist, nlon, nlat, p[0], p[1], method);
			out.setPointsGeometry(nlon, nlat);
			out.df.add_column(id, "id");
			out.df.add_column(dist, "distance");
			return out;
		} else {
			out.setError("not yet implement for non-point lonlat vector data");
			return out;
		}
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> b(n);
	for (size_t i = 0; i < n; i++) {
		SpatVector xa = remove_rows({i});
		xa = xa.aggregate(false);
		std::vector<GeomPtr> y = geos_geoms(&xa, hGEOSCtxt);
		GEOSCoordSequence* csq = GEOSNearestPoints_r(hGEOSCtxt, x[i].get(), y[0].get());
		GEOSGeometry* geom = GEOSGeom_createLineString_r(hGEOSCtxt, csq);
		b[i] = geos_ptr(geom, hGEOSCtxt);
	}
	out = vect_from_geos(b, hGEOSCtxt, "lines");
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	return out;
}

#ifdef GEOS361

// helper struct for STRtree:
typedef struct { GEOSGeom g; size_t id; } item_g;

int distance_fn(const void *item1, const void *item2, double *distance, void *userdata) {
	return GEOSDistance_r( (GEOSContextHandle_t) userdata, ((item_g *)item1)->g, ((item_g *)item2)->g, distance);
}

std::vector<int> SpatVector::nearest_geometry(SpatVector v) {

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);

	TreePtr tree = geos_ptr(GEOSSTRtree_create_r(hGEOSCtxt, 10), hGEOSCtxt);
	std::vector<item_g> items(y.size());
	bool tree_is_empty = true;
	for (size_t i = 0; i < y.size(); i++) {
		items[i].id = i;
		items[i].g = y[i].get();
		if (!GEOSisEmpty_r(hGEOSCtxt, y[i].get())) {
			GEOSSTRtree_insert_r(hGEOSCtxt, tree.get(), y[i].get(), &(items[i]));
			tree_is_empty = false;
		}
	}
	std::vector<int> out;
	if (tree_is_empty) {
		setError("cannot make spatial index");
		return out;
	}

	out.resize(nrow(), -1);
	for (size_t i = 0; i < x.size(); i++) {
		if (!GEOSisEmpty_r(hGEOSCtxt, x[i].get())) {
			item_g item, *ret_item;
			item.id = -99;
			item.g = x[i].get();
			ret_item = (item_g *) GEOSSTRtree_nearest_generic_r(hGEOSCtxt, tree.get(), &item,
					x[i].get(), distance_fn, hGEOSCtxt);
			if (ret_item != NULL) {
				out[i] = ret_item->id; 
			} else {
				setError("GEOS error");
				return out;
			}
		} 
	}
	geos_finish(hGEOSCtxt);

//	SpatVector out = v.subset_rows(sel);
	return out;
}
#else
std::vector<int> SpatVector::nearest_geometry(SpatVector v) {
	setError("you need GEOS 3.6.1 for this method");
	std::vector<int> out;
	return out;
}
#endif // GEOS361


SpatVector SpatVector::cross_dateline(bool &fixed) {
	SpatVector out;
	fixed = false;
	if (type() == "points") {
		return out;
	}
	out.reserve(size());

	for (size_t i=0; i<geoms.size(); i++) {
		if ((geoms[i].size() > 1) &&
			((geoms[i].extent.xmax - geoms[i].extent.xmin) > 180)) {
			SpatGeom g = geoms[i];
			for (size_t j=0; j<g.size(); j++) {
				if (g.parts[j].extent.xmax < 0) {
					for (size_t k=0; k<g.parts[j].x.size(); k++) {
						g.parts[j].x[k] += 360;
					}
					for (size_t k=0; k<g.parts[j].holes.size(); k++) {
						for (size_t m=0; m<g.parts[j].holes[k].x.size(); m++) {
							g.parts[j].holes[k].x[m] += 360;
						}
					}
					g.parts[j].extent.xmin += 360;
					g.parts[j].extent.xmax += 360;
					g.setPart(g.parts[j], j);
					fixed = true;
				}
			}
			out.addGeom(g);
		} else {
			out.addGeom(geoms[i]);
		}
	}
	out.srs = srs;
	out.df = df;
	return out;
}



SpatVector SpatVector::centroid(bool check_lonlat) {

	SpatVector out;
	if (nrow() == 0) {
		out.setError("input has no geometries");
		return out;
	}

	if (check_lonlat && could_be_lonlat()) {
		bool changed = false;
		SpatVector v = cross_dateline(changed);
		if (changed) {
			out = v.centroid(false);
			out.fix_lonlat_overflow();
			return out;
		}
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> b(size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* pt = GEOSGetCentroid_r(hGEOSCtxt, g[i].get());
		if (pt == NULL) {
			out.setError("NULL geom");
			geos_finish(hGEOSCtxt);
			return out;
		}
		b[i] = geos_ptr(pt, hGEOSCtxt);
	}
	out = vect_from_geos(b, hGEOSCtxt, "points");
	geos_finish(hGEOSCtxt);

	out.srs = srs;
	out.df = df;
	return out;
}


SpatVector SpatVector::point_on_surface(bool check_lonlat) {

	SpatVector out;
	if (nrow() == 0) {
		out.setError("input has no geometries");
		return out;
	}

	if (check_lonlat && could_be_lonlat()) {
		bool changed = false;
		SpatVector v = cross_dateline(changed);
		if (changed) {
			out = v.point_on_surface(false);
			out.fix_lonlat_overflow();
			return out;
		}
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> b(size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* pt = GEOSPointOnSurface_r(hGEOSCtxt, g[i].get());
		if (pt == NULL) {
			out.setError("NULL geom");
			geos_finish(hGEOSCtxt);
			return out;
		}
		b[i] = geos_ptr(pt, hGEOSCtxt);
	}
	out = vect_from_geos(b, hGEOSCtxt, "points");
	geos_finish(hGEOSCtxt);

	out.srs = srs;
	out.df = df;
	return out;
}


SpatVector SpatVector::unaryunion() {
	SpatVector out;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> gout(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* u = GEOSUnaryUnion_r(hGEOSCtxt, g[i].get());
		if (u == NULL) {
			out.setError("NULL geom");
			geos_finish(hGEOSCtxt);
			return out;
		}
		gout[i] = geos_ptr(u, hGEOSCtxt);
	}
	SpatVectorCollection coll = coll_from_geos(gout, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	out = coll.get(0);
	out.srs = srs;
	return out;
}


/*
bool geos_buffer(GEOSContextHandle_t hGEOSCtxt, std::vector<GeomPtr> &g, double dist, unsigned nQuadSegs) {
	std::vector<GeomPtr> g(size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* pt = GEOSBuffer_r(hGEOSCtxt, g[i].get(), dist, nQuadSegs);
		if (pt == NULL) {
			return false;
		}
		g[i] = geos_ptr(pt, hGEOSCtxt);
	}
	return true;
}
*/



SpatVector SpatVector::width() {

	SpatVector tmp;

#ifndef GEOS361
	tmp.setError("GEOS 3.6.1 required for width");
	return tmp;
#else

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> gout(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* w = GEOSMinimumWidth_r(hGEOSCtxt, g[i].get());
		if (w == NULL) {
			tmp.setError("found NULL geom");
			geos_finish(hGEOSCtxt);
			return tmp;
		}
		gout[i] = geos_ptr(w, hGEOSCtxt);
	}
	SpatVectorCollection coll = coll_from_geos(gout, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	tmp = coll.get(0);
	tmp.srs = srs;
	return tmp;

#endif
}

SpatVector SpatVector::clearance() {
	SpatVector tmp;
#ifndef GEOS361
	tmp.setError("GEOS 3.6 required for clearance");
	return tmp;
#else

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> gout(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* w = GEOSMinimumClearanceLine_r(hGEOSCtxt, g[i].get());
		if (w == NULL) {
			tmp.setError("NULL geom");
			geos_finish(hGEOSCtxt);
			return tmp;
		}
		gout[i] = geos_ptr(w, hGEOSCtxt);
	}
	SpatVectorCollection coll = coll_from_geos(gout, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	tmp = coll.get(0);
	tmp.srs = srs;
	return tmp;

#endif
}


bool SpatPart::is_CCW() {
#ifndef GEOS370
	return true;
#else
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	GEOSCoordSequence *pseq;
	size_t n = size();
	pseq = GEOSCoordSeq_create_r(hGEOSCtxt, n, 2);
	for (size_t i = 0; i < n; i++) {
		GEOSCoordSeq_setX_r(hGEOSCtxt, pseq, i, x[i]);
		GEOSCoordSeq_setY_r(hGEOSCtxt, pseq, i, y[i]);
	}
	char is_ccw;
	bool success = GEOSCoordSeq_isCCW_r(hGEOSCtxt, pseq, &is_ccw);
	GEOSCoordSeq_destroy_r(hGEOSCtxt, pseq);
	geos_finish(hGEOSCtxt);
	if (success) {
		return is_ccw != 0;
	} else {
		return true;
	}
#endif
}



void SpatVector::make_CCW() {
	#ifndef GEOS370
		setError("GEOS >= 3.7 needed for CCW");
		return;
	#else
	size_t n = size();
	if (n == 0) return;
	if (geoms[0].gtype != polygons) return;
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j<geoms[i].parts.size(); j++) {
			if (!geoms[i].parts[j].is_CCW()) {
				std::reverse(geoms[i].parts[j].x.begin(), geoms[i].parts[j].x.end());
				std::reverse(geoms[i].parts[j].y.begin(), geoms[i].parts[j].y.end());
				for (size_t k=0; k<geoms[i].parts[j].nHoles(); k++) {
					std::reverse(geoms[i].parts[j].holes[k].x.begin(), geoms[i].parts[j].holes[k].x.end());
					std::reverse(geoms[i].parts[j].holes[k].y.begin(), geoms[i].parts[j].holes[k].y.end());
				}
			}
		}
	}
	#endif
}


