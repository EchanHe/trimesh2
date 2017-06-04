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
//		file << "gaussian,mean"<<"\n";
//		BOOST_FOREACH(boost::tie(a, b),
//			boost::combine(gaus_curv, mean_curv)) {
//			file << a << "," << b<<"\n";
//
//		}
//	}
//}