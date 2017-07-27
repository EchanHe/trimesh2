#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "kdTree_face.h"
#include "octree.h"
#include <numeric>
//#include "Vec.h"
using namespace std;
namespace trimesh {
	const int CONE_ANGLE = 60;
	const int RAYS_PER_ANGLE = 20;
	const int ANGLE_INTERVALS = 4;

	//rotation( vec dir, float angle, vec rotateAxis)
	//Goal: return the dirction vector after rotation *angle along *rotateaxis
	/*
	input:
		dir:		unit vector of direction
		angle:		Rotation angle in degree
		rotateAxis:	unit vector of the rotation axis
	*/
	inline vec rotation( vec dir, float angle, vec rotateAxis) {
		//vec oDir = dir - origin;
		float cos_angle = cos(angle*M_PIf/180);
		float sin_angle = sin(angle*M_PIf / 180);
		//stack overflow 
		//https://stackoverflow.com/questions/42421611/3d-vector-rotation-in-c
		vec result = (cos_angle* dir) + ((rotateAxis CROSS dir )*sin_angle) + (rotateAxis * (rotateAxis ^ dir) *(1- cos_angle));

		//wikipedia way
		vec Rr1 = vec(cos_angle + (rotateAxis[0] * rotateAxis[0])*(1 - cos_angle),
			(rotateAxis[0] * rotateAxis[1])*(1 - cos_angle) - rotateAxis[2] * sin_angle,
			(rotateAxis[0] * rotateAxis[2])*(1 - cos_angle) + rotateAxis[1] * sin_angle);

		vec Rr2 = vec((rotateAxis[0] * rotateAxis[1])*(1 - cos_angle) + rotateAxis[2] * sin_angle,
			cos_angle + (rotateAxis[1] * rotateAxis[1])*(1 - cos_angle),
			(rotateAxis[1] * rotateAxis[2])*(1 - cos_angle) - rotateAxis[0] * sin_angle);

		vec Rr3 = vec((rotateAxis[0] * rotateAxis[2])*(1 - cos_angle) - rotateAxis[1] * sin_angle,
			(rotateAxis[1] * rotateAxis[2])*(1 - cos_angle) + rotateAxis[0] * sin_angle,
			cos_angle + (rotateAxis[2] * rotateAxis[2])*(1 - cos_angle));

		result =vec( Rr1 ^ dir, Rr2 ^ dir, Rr3 ^ dir );
		//result = result + origin;
		normalize(result);
		return result;
	}
	//Goal: make a cone of rays
	/*input:
	halfAngle:	Half angle
	rings:		rays on
	v1,v2,v3:	points of the face
	*/
	void make_cone( vec origin, vec dir, vector<vec>& output) {
		output.clear();
		//float halfAngle = 60;

		//int totalRings = 20;
		//int intervals = 5;	

		output.resize(RAYS_PER_ANGLE*ANGLE_INTERVALS +1 );

		//find perpendicular  vector as rotation axis

		vec rotateAxis = vec(1, 1, (-dir[0] - dir[1]) / dir[2]);
		if (dir[2] == 0) {
			rotateAxis = vec(0, 0, 1);
		}
		normalize(rotateAxis);
		normalize(dir);
		
		for (int j = 0; j < ANGLE_INTERVALS; j++) {
			vec ray1 = rotation(dir, CONE_ANGLE-(15*j), rotateAxis);
			for (int i = 0; i < RAYS_PER_ANGLE; i++) {
				vec rayI = rotation(ray1, 360 * i / RAYS_PER_ANGLE, dir);
				output[(j*RAYS_PER_ANGLE)+i] = rayI;
			}
		}
		output[output.size()-1] = dir;
	}



	//Goal: Calculate the distance
	/*input:
		vectex:	vector start point
		dir:		vecter direction
		v1,v2,v3:	points of the face
	*/
	/*Moller¨CTrumbore intersection algorithm
		p + t * d = (1-u-v) * p0 + u * p1 + v * p2
		t , u v are parameters
		p0 p1 p2 are apex of triangles
		p are point on ray, d are dir
	
	*/
	 float triangle_inter(vec vertex ,vec dir , vec v1, vec v2, vec v3) {
		
		vec e[2] = { v2 - v1,
			v3 - v1 };

		vec pvec = e[1] CROSS dir;
		float det = e[0] ^ pvec;

		if (det < 1e-8 && det > -1e-8) {
			//std::cout << "parrel";
			return 0;
		}
		float inv_det = 1 / det;
		vec tvec = vertex - v1;
		float u = (tvec ^ pvec) * inv_det;
		if (u < 0 || u > 1) {
			return 0;
		}

		vec3 qvec =tvec CROSS e[0];
		float v = (dir ^ qvec) * inv_det;
		if (v < 0 || u + v > 1) {
			return 0;
		}
		float t = inv_det *(e[1] ^ qvec);
		
		//calculate the intersection points.
		vec inter_point = vertex + (t*dir);
		return len(t*dir);
	}

	template < class T>
	float sdf_stat_mean(vector<T> input) {
		float sum = std::accumulate(input.begin(), input.end(), 0.0f);
		float mean = sum / input.size();
		float sq_sum = std::inner_product(input.begin(), input.end(), input.begin(), 0.0);
		float stdev = std::sqrt(sq_sum / input.size() - mean * mean);

		return mean;
	}

	template < class T>
	float sdf_stat_mean(vector<T> input , float median , float stdev) {
		vector<T> new_input;
		int n_stdev = 1;
		float lower_threshold = median - stdev*n_stdev;
		float upper_threshold = median + stdev*n_stdev;
		for (int i = 0; i < input.size(); i++) {
			if (input[i] >= lower_threshold && input[i] <= upper_threshold) {
				new_input.push_back(input[i]);
			}
		}

		float sum = std::accumulate(new_input.begin(), new_input.end(), 0.0f);
		float mean = sum / new_input.size();
		return mean;
	}

	template < class T>
	float sdf_stat_stdev(vector<T> input) {
		float sum = std::accumulate(input.begin(), input.end(), 0.0f);
		float mean = sum / input.size();
		float sq_sum = std::inner_product(input.begin(), input.end(), input.begin(), 0.0);
		float stdev = std::sqrt(sq_sum / input.size() - mean * mean);

		return stdev;
	}

	template < class T>
	float sdf_stat_median(vector<T> vals) {
		int n = vals.size();
		if (n & 1) {
			nth_element(vals.begin(),
				vals.begin() + n / 2,
				vals.end());
			return vals[n / 2];
		}
		else {
			nth_element(vals.begin(),
				vals.begin() + n / 2 - 1,
				vals.end());
			float tmp = vals[n / 2 - 1];
			nth_element(vals.begin(),
				vals.begin() + n / 2,
				vals.end());
			return 0.5f * (tmp + vals[n / 2]);
		}
	}

	//---------function for sdf-------------------------
	void TriMesh::need_sdf() {
		std::clock_t begin = clock();
		need_normals();
		need_inwardNormals();
		need_faceNormals();
		need_faceMidPts();
		std::cout << std::endl << "Build kd tree: ";
		KDtreeFace *k = new KDtreeFace(faceMidPts);
		std::cout << double(clock() - begin) / CLOCKS_PER_SEC << " s" << std::endl;
		begin = clock();
		if (sdf.size() == vertices.size())
			return;
		sdf.resize(vertices.size());

		int nv = vertices.size();
		//nv = 1;
		int nf = faces.size();

		std::vector<vec> rays;

		std::cout << std::endl << "computing the shape diameters using kd tree: " ;
		
//#pragma omp parallel for

		for (int i = 0; i < nv; i++) {
			//initialize the vertex, in normal and rays for SDF
			vec inwardNormal = inwardNormals[i];
			vec vertex = vertices[i];
			make_cone(vertex, inwardNormal, rays);
		
			
			//vector<float>distances;
			//distances.resize(rays.size());
			//for (int j = 0;  j < rays.size(); j++) {
			//	const float* a = k->closest_to_ray(vertex, rays[j]);
			//	if (a != NULL) {
			//		vec intersect(a[0], a[1], a[2]);
			//		distances[j] = len(intersect - vertex);
			//	}
			//	if (i % 10000 == 0&&i/10000!=0)
			//	std::cout << "finish one "<<i << endl;
			//}

			const float* a = k->closest_to_ray(vertex, inwardNormal);
			vec intersect(a[0], a[1], a[2]);
			sdf[i] = len(intersect - vertex);
			//sdf[i] = sdf_stat_mean(distances);
			if (i% (nv / 10)==0)
				std::cout << "Finish "<< i / (nv / 10)<<" % of ht SDF" << std::endl;

		}

		std::cout << double(clock() - begin) / CLOCKS_PER_SEC << " s" << std::endl;
	}

	void TriMesh::need_sdf_brute() {
		need_normals();
		need_inwardNormals();

		need_faceNormals();
		if (sdf_brute.size() == vertices.size())
			return;
		sdf_brute.resize(vertices.size());

		

		int nv = vertices.size();
		//nv = 1;
		int nf = faces.size();

		std::vector<vec> rays;
		make_cone(vertices[8], inwardNormals[8], rays);
		std::cout << std::endl << "computing the shape diameters using brute: ";
		std::clock_t begin = clock();
		//#pragma omp parallel for

		//Octree * octreeFace = new Octree(faces, vertices, faceNormals);
		int diffCount = 0;
		
		for (int i = 82; i < nv; i++) {
			//initialize the vertex, in normal and rays for SDF
			vec inwardNormal = inwardNormals[i];
			vec normal = normals[i];
			vec vertex = vertices[i];
			make_cone(vertex, inwardNormal, rays);

			vector<float>distances;
	//		distances.resize(rays.size());
//#pragma omp parallel for		
			for (int j = 0; j < rays.size(); j++) {
				float minDist = 0;
				int faceId = -1;
				vec ray = rays[j];
//#pragma omp parallel for
				for (int k = 0; k < nf; k++) { //iterate through faces
					float angle = acos(normal ^ faceNormals[k])* 180.0 / M_PIf;
					if (angle > 90) {

						vec v1 = vertices[faces[k][0]]; vec v2 = vertices[faces[k][1]]; vec v3 = vertices[faces[k][2]];

						float distance = triangle_inter(vertex, ray, v1, v2, v3);
						if (distance != 0) {
							if (minDist == 0) {
								minDist = distance;
								faceId = k;
							}
							else {
								if (distance < minDist) {
									minDist = distance;
									faceId = k;
								}
							}
						}
					}
				}
				//float d_from_octree = octreeFace->intersect_face_from_raycast(vertex, ray, normal).closest_d;
				//if (minDist != d_from_octree) {
				//	diffCount++;
				//}

				if (minDist != 0 || faceId !=-1) {//--If there is a intersect between a face and the ray
					//--check the intersect is with the inward faces
						distances.push_back(minDist);
				//	minDist = minDist;
					//distances[j] = minDist;
						//int ring = (j) / RAYS_PER_ANGLE;
						//if (j == rays.size() - 1)
						//	distances.push_back(minDist);
						////distances[j] = minDist;
						//else
						//	distances.push_back(minDist);// ((ANGLE_INTERVALS - ring) * 2));
						//	//distances[j] = minDist / ((ANGLE_INTERVALS - ring) * 2);
				}


			}
			if (distances.size() <= 10) {
				cout << "The vertices " << i;
				cout << "Have " << distances.size() << " intersect"<<endl;
				
			}

			float median = sdf_stat_median(distances);
			float stdev = sdf_stat_stdev(distances);
			
			sdf_brute[i] = sdf_stat_mean(distances, median , stdev);
			if (i % (nv / 10) == 0)
				std::cout << "Finish " << i / (nv / 10) << "0% of ht SDF" << std::endl;
		}
		//cout << "The different percentatge:" << float(diffCount) / (nv * (RAYS_PER_ANGLE*ANGLE_INTERVALS+1))<<endl;
		std::cout << double(clock() - begin) / CLOCKS_PER_SEC << " s" << std::endl;
	}

	void TriMesh::need_sdf_from_simple(TriMesh *simple) {
		if (sdf.size() == vertices.size())
			return;
		int nv = vertices.size();
		int nf = faces.size();

		sdf.resize(vertices.size());
		std::clock_t begin = clock();
		need_normals();
		need_inwardNormals();
		need_faceNormals();
		need_faceMidPts();
		//using the Vertices to build the Octree
		//Octree *k = new Octree(vertices);
		Octree * octreeFace = new Octree(simple->faces, simple->vertices, simple->faceNormals);

		begin = clock();
		std::vector<vec> rays;
		std::cout << std::endl << "computing the shape diameters using octree tree: ";

		for (int i = 0; i < nv; i++) {
			//---initialize the vertex, in normal and rays for SDF
			vec inwardNormal = inwardNormals[i];
			vec vertex = vertices[i];
			vec normal = normals[i];
			make_cone(vertex, inwardNormal, rays);


			vector<float>distances;
			//distances.resize(rays.size());
			distances.clear();
			for (int j = 0; j < rays.size(); j++) {
				//Octree::Traversal_Info  info = k->find_cube_from_raycast(vertex, rays[j]);
				float d_from_octree = octreeFace->intersect_face_from_raycast(vertex, rays[j], normal).closest_d;
				if (d_from_octree != 0)
					distances.push_back(d_from_octree);
				//distances[j] = d_from_octree;
			}
			//const float* a = 
			//	vec intersect(a[0], a[1], a[2]);
			float median = sdf_stat_median(distances);
			float stdev = sdf_stat_stdev(distances);
			//sdf[i] = sdf_stat_mean(distances);
			sdf[i] = sdf_stat_mean(distances, median, stdev);
			if (i % (nv / 100) == 0)
				std::cout << "Finish " << i / (nv / 100) << "% of ht SDF" << std::endl;
		}

		std::cout << double(clock() - begin) / CLOCKS_PER_SEC << " s" << std::endl;
	}

	void TriMesh::need_sdf_octree() {
		if (sdf.size() == vertices.size())
			return;
		int nv = vertices.size();
		int nf = faces.size();

		sdf.resize(vertices.size());
		std::clock_t begin = clock();
		need_normals();
		need_inwardNormals();
		need_faceNormals();
		need_faceMidPts();
		//using the Vertices to build the Octree
		//Octree *k = new Octree(vertices);
		Octree * octreeFace = new Octree(faces, vertices, faceNormals);
		Octree::Traversal_Info info = octreeFace->find_all_leaves();
		cout << "branch, leaves, correct size leaves :" << info.branch_count << " " << info.leaves_count << " " << info.equal_max_faces_count << endl;
		begin = clock();
		std::vector<vec> rays;
		std::cout << std::endl << "computing the shape diameters using octree tree: ";
		int cal_count=0;
		for (int i = 0; i < nv; i++) {
			//---initialize the vertex, in normal and rays for SDF
			vec inwardNormal = inwardNormals[i];
			vec vertex = vertices[i];
			vec normal = normals[i];
			make_cone(vertex, inwardNormal, rays);


			vector<float>distances;
			//distances.resize(rays.size());
			distances.clear();
			for (int j = 0;  j < rays.size(); j++) {
				//Octree::Traversal_Info  info = k->find_cube_from_raycast(vertex, rays[j]);
				Octree::Traversal_Info info = octreeFace->intersect_face_from_raycast(vertex, rays[j], normal);
				float d_from_octree = info.closest_d;
				cal_count += info.cal_count;
				if (d_from_octree != 0)
					distances.push_back(d_from_octree);
					//distances[j] = d_from_octree;
			}
			//const float* a = 
		//	vec intersect(a[0], a[1], a[2]);
			float median = sdf_stat_median(distances);
			float stdev = sdf_stat_stdev(distances);
			//sdf[i] = sdf_stat_mean(distances);
			sdf[i] = sdf_stat_mean(distances, median, stdev);
			if (i % (nv / 10) == 0)
				std::cout << "Finish " << i / (nv / 10) << "0% of ht SDF" << std::endl;
		}
		std::cout << double(clock() - begin) / CLOCKS_PER_SEC << " s" << std::endl;
	//	std::cout << "The brute force calculation times: "<<nv*nf*81<<endl;
		std::cout << "The faces traverse per ray and vertices " << cal_count/(nv*81);

	}

	void TriMesh::compare_sdfs()
	{
		int decimal = 3;
		if (sdf.size() == 0 || sdf_brute.size() == 0) {
			need_sdf_octree();
			need_sdf_brute();
		}
		int count = 0;
		for (int i = 0; i < sdf.size(); i++) {
			float divider = 10 ^ decimal;
			float round_sdf = roundf(sdf[i] * divider) / divider;
			float round_sdf_brute = roundf(sdf_brute[i] * divider) / divider;
			if (round_sdf != round_sdf_brute)
				count++;
		}
		cout << "The different percent SDF between brute force and octree are: " << float(count )/ sdf.size()<<endl;


	}

	void TriMesh::writeSDF() {
		std::cout << "Writing curvature to sdf.csv" << std::endl;
		std::ofstream file;
		file.open("sdf.csv");
		file << "sdf" << "\n";
		for (int i = 0; i < sdf.size(); i++) {
			file << sdf[i] << "\n";
		}
		file.close();

		file.open("sdf_brute.csv");
		file << "sdf_brute" << "\n";
		for (int i = 0; i < sdf.size(); i++) {
			file << sdf_brute[i] << "\n";
		}
		file.close();
	}


}//end of namespace

/*
		need_normals();
		need_inwardNormals();

		need_faceNormals();
		if (colors.size() == vertices.size())
			return;
		colors.resize(vertices.size());

		int nv = vertices.size();
		nv = 1;
		int nf = faces.size();

		std::vector<vec> rays;
		for (int i = 0; i < nv; i++) {
			//initialize the vertex, in normal and rays for SDF
			vec inwardNormal = inwardNormals[i];
			vec vertex = vertices[i];
			make_cone(vertex, inwardNormal, rays);

			float minDist = 0;
			int faceId=0;
			for (int j = 0; j < nf; j++) {
				vec v1 = vertices[faces[j][0]]; vec v2 = vertices[faces[j][1]]; vec v3 = vertices[faces[j][2]];

				float distance = triangle_inter(vertex, inwardNormal, v1, v2, v3);
				if (distance != 0) {
					if (minDist == 0) {
						minDist = distance;
						faceId = j;
					}
					else {
						if (distance < minDist) {
							minDist = distance;
							faceId = j;
						}
					}
				}

			}
			// the distant and faces
			std::cout << "face number : " << faceId<< " distance: "<<minDist<<std::endl;
			colors[faces[faceId][0]] = Color(0.0f, 1.0f, 0.0f);
			colors[faces[faceId][1]] = Color(0.0f, 1.0f, 0.0f);
			colors[faces[faceId][2]] = Color(0.0f, 1.0f, 0.0f);

		}

*/