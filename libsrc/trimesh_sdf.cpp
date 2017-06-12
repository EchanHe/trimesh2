#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "kdTree_face.h"
#include <numeric>
//#include "Vec.h"
using namespace std;
namespace trimesh {
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
	inline void make_cone(float halfAngle, int rings, int intervals,vec origin,vec normal, vector<vec>& output) {
		output.resize(rings);

		//find perpendicular  vector as rotation axis
		vec rotateAxis = vec(1, 1, (-normal[0] - normal[1]) / normal[2]);
		normalize(rotateAxis);

		vec ray1 = rotation(normal, halfAngle, rotateAxis);

		for (int i = 0; i < rings; i++) {
			vec rayI = rotation(ray1, 360*i/rings, normal);
			output[i] = rayI;
		}
		
	}
	void make_cone( vec origin, vec normal, vector<vec>& output) {
		output.clear();
		float halfAngle = 60;
		int totalRings = 20;
		int intervals = 5;	
		output.resize(totalRings*intervals);

		//find perpendicular  vector as rotation axis
		vec rotateAxis = vec(1, 1, (-normal[0] - normal[1]) / normal[2]);
		normalize(rotateAxis);

		
		for (int j = 0; j < intervals; j++) {
			vec ray1 = rotation(normal, halfAngle-(20*j), rotateAxis);
			for (int i = 0; i < totalRings; i++) {
				vec rayI = rotation(ray1, 360 * i / totalRings, normal);
				output[(j*totalRings)+i] = rayI;
			}
		}
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
	float sdf_stat(vector<T> input) {
		float sum = std::accumulate(input.begin(), input.end(), 0.0f);
		float mean = sum / input.size();
		float sq_sum = std::inner_product(input.begin(), input.end(), input.begin(), 0.0);
		float stdev = std::sqrt(sq_sum / input.size() - mean * mean);

		return mean;
	}


	void TriMesh::need_sdf() {
		std::clock_t begin = clock();
		need_normals();
		need_inwardNormals();
		need_faceNormals();
		need_faceMidPts();
		std::cout << std::endl << "Budld kd tree: ";
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
		
			
			vector<float>distances;
			distances.resize(rays.size());
			for (int j = 0;  j < rays.size(); j++) {
				const float* a = k->closest_to_ray(vertex, rays[j]);
				if (a != NULL) {
					vec intersect(a[0], a[1], a[2]);
					distances[j] = len(intersect - vertex);
				}
				if (i % 10000 == 0&&i/10000!=0)
				std::cout << "finish one "<<i << endl;



			}
			sdf[i] = sdf_stat(distances);

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

		std::cout << std::endl << "computing the shape diameters using brute: ";
		std::clock_t begin = clock();
		//#pragma omp parallel for

		for (int i = 0; i < nv; i++) {
			//initialize the vertex, in normal and rays for SDF
			vec inwardNormal = inwardNormals[i];
			vec vertex = vertices[i];
			make_cone(vertex, inwardNormal, rays);

			vector<float>distances;
			distances.resize(rays.size());
			for (int j = 0; j < rays.size(); j++) {
				float minDist = 0;
				int faceId = 0;
				vec ray = rays[j];
				for (int j = 0; j < nf; j++) {
					vec v1 = vertices[faces[j][0]]; vec v2 = vertices[faces[j][1]]; vec v3 = vertices[faces[j][2]];

					float distance = triangle_inter(vertex, ray, v1, v2, v3);
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
				distances[j] = minDist;
			}
			sdf_brute[i] = sdf_stat(distances);



		}

		std::cout << double(clock() - begin) / CLOCKS_PER_SEC << " s" << std::endl;
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