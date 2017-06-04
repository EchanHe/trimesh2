//#pragma once
//
//#include <iostream>
//#include <fstream>
//#include "TriMesh.h"
//
//#include <boost\iterator\zip_iterator.hpp>
//#include <boost\tuple\tuple.hpp>
//#include <boost\tuple\tuple_io.hpp>
//#include <boost/foreach.hpp>
//#include <boost/range/combine.hpp>
//
//using namespace std;
//namespace trimesh {
//
//	void TriMesh::writeCurvature() {
//		
//		std::ofstream file;
//		file.open("curvature.csv");
//		float a, b;
//		BOOST_FOREACH(boost::tie(a, b),
//			boost::combine(gaus_curv, gaus_curv)) {
//			file << a << "  " << b;
//
//		}
//	}
//}
//#include <geodesic\geodesic_mesh_elements.h>
//#include <geodesic\geodesic_mesh.h>
//namespace geodesic {
//	void Mesh::initialize() {
//
//	}
//}