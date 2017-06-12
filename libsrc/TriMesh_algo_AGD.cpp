#include "TriMesh.h"
#include "TriMesh_algo.h"
using namespace geodesic;
namespace trimesh {


	void init_geod_mesh(TriMesh * mesh , geodesic::Mesh & geoMesh) {

		std::vector<double> points;
		std::vector<unsigned> faces;
		points.resize(mesh->vertices.size() * 3);
		faces.resize(mesh->faces.size() * 3);
		for (int i = 0; i < mesh->vertices.size(); i++) {
			int id = i * 3;
			points[id] = mesh->vertices[i][0];
			points[id + 1] = mesh->vertices[i][1];
			points[id + 2] = mesh->vertices[i][2];
			//points.push_back(mesh->vertices[i][0]);
			//points.push_back(mesh->vertices[i][1]);
			//points.push_back(mesh->vertices[i][2]);
		}

		for (int i = 0; i < mesh->faces.size(); i++) {
			int id = i * 3;
			faces[id] = mesh->faces[i][0];
			faces[id + 1] = mesh->faces[i][1];
			faces[id + 2] = mesh->faces[i][2];
		}

		clock_t begin = clock();

		geoMesh.initialize_mesh_data(points, faces);

		std::cout << "Geodesic Mesh constructing time: " << double(clock() - begin) / CLOCKS_PER_SEC << std::endl;
		begin = clock();

	};

	template<class Points, class Faces>
	std::vector<geodesic::SurfacePoint> cal_geo_dis(int start, int end, Points points, Faces faces) {
		clock_t begin = clock();

		geodesic::Mesh mesh;
		mesh.initialize_mesh_data(points, faces);

		cout << "Geodesic Mesh constructing time: " << double(clock() - begin) / CLOCKS_PER_SEC << std::endl;
		begin = clock();

		geodesic::GeodesicAlgorithmExact algorithm(&mesh);

		cout << "Geodesic Algorithm constructing time: " << double(clock() - begin) / CLOCKS_PER_SEC << std::endl;
		begin = clock();

		geodesic::SurfacePoint source(&mesh.vertices()[start]);
		std::vector<geodesic::SurfacePoint> all_sources(1, source);

		unsigned target_vertex_index = 2;
		geodesic::SurfacePoint target(&mesh.vertices()[end]);		//create target 

		std::vector<geodesic::SurfacePoint> path;	//geodesic path is a sequence of SurfacePoints

		bool const lazy_people_flag = false;		//there are two ways to do exactly the same
		if (lazy_people_flag)
		{
			algorithm.geodesic(source, target, path); //find a single source-target path
		}
		else		//doing the same thing explicitly for educational reasons
		{
			double const distance_limit = geodesic::GEODESIC_INF;			// no limit for propagation
			std::vector<geodesic::SurfacePoint> stop_points(1, target);	//stop propagation when the target is covered
			algorithm.propagate(all_sources, distance_limit, &stop_points);	//"propagate(all_sources)" is also fine, but take more time because covers the whole mesh

			algorithm.trace_back(target, path);		//trace back a single path 
		}
		return path;


		//print_info_about_path(path);
		//for (unsigned i = 0; i<path.size(); ++i)
		//{
		//	geodesic::SurfacePoint& s = path[i];
		//	std::cout << s.x() << "\t" << s.y() << "\t" << s.z() << std::endl;
		//}

	}

}