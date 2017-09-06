#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "kdTree_face.h"
#include "octree.h"
#include <numeric>
#include "cuda_runtime.h"
#include<cuda.h>
#include<device_launch_parameters.h>
#include<conio.h>
//#include "Vec.h"
using namespace std;
namespace trimesh {
	//#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
	//
	//	inline void gpuAssert(cudaError_t code, char *file, int line, bool abort = true)
	//	{
	//		if (code != cudaSuccess)
	//		{
	//			fprintf(stderr, "GPUassert: %s %s %dn\n", cudaGetErrorString(code), file, line);
	//			//if (abort) { getch(); exit(code); }
	//		}
	//	}



	const int CONE_ANGLE = 60;
	const int RAYS_PER_ANGLE = 20;
	const int ANGLE_INTERVALS = 4;
	void vec_to_pointer(std::vector<vec>input, float***result) {
		int Ncol = 3;
		int Nrow = input.size();
		int nf = Nrow;
		*result = new float*[Nrow];
		for (int i = 0; i < Nrow; i++) {
			(*result)[i] = new float[Ncol];
			for (int j = 0; j < 3; j++) {
				(*result)[i][j] = input[i][j];
			}
		}
	};

	//rotation( vec dir, float angle, vec rotateAxis)
	//Goal: return the dirction vector after rotation *angle along *rotateaxis
	/*
	input:
	dir:		unit vector of direction
	angle:		Rotation angle in degree
	rotateAxis:	unit vector of the rotation axis
	*/
	inline vec rotation(vec dir, float angle, vec rotateAxis) {
		//vec oDir = dir - origin;
		float cos_angle = cos(angle*M_PIf / 180);
		float sin_angle = sin(angle*M_PIf / 180);
		//stack overflow 
		//https://stackoverflow.com/questions/42421611/3d-vector-rotation-in-c
		vec result = (cos_angle* dir) + ((rotateAxis CROSS dir)*sin_angle) + (rotateAxis * (rotateAxis ^ dir) *(1 - cos_angle));

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

		result = vec(Rr1 ^ dir, Rr2 ^ dir, Rr3 ^ dir);
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
	inline void make_cone(float halfAngle, int rings, int intervals, vec origin, vec normal, vector<vec>& output) {
		output.resize(rings);

		//find perpendicular  vector as rotation axis
		vec rotateAxis = vec(1, 1, (-normal[0] - normal[1]) / normal[2]);
		normalize(rotateAxis);

		vec ray1 = rotation(normal, halfAngle, rotateAxis);

		for (int i = 0; i < rings; i++) {
			vec rayI = rotation(ray1, 360 * i / rings, normal);
			output[i] = rayI;
		}

	}
	int iDivUp(int hostPtr, int b) { return ((hostPtr % b) != 0) ? (hostPtr / b + 1) : (hostPtr / b); }

	__device__ void  cross_device(float* v1, float * v2, float(&result)[3]) {
		//float result[3];
		result[0] = v1[1] * v2[2] - v1[2] * v2[1];
		result[1] = v1[2] * v2[0] - v1[0] * v2[2];
		result[2] = v1[0] * v2[1] - v1[1] * v2[0];
		//return result;
	}

	__device__ float dot_device(float * v1, float *v2) {
		float sum = v1[0] * v2[0];
		for (size_t i = 1; i < 3; i++)
			sum += v1[i] * v2[i];
		return sum;
	}

	__device__ void plus_device(float *v1, float* v2, float(&result)[3], bool plus = true) {
		//float * result2 = new float[3];
		if (plus) {
			for (int i = 0; i < 3; i++) {
				result[i] = v1[i] + v2[i];
			}
		}
		else {
			for (int i = 0; i < 3; i++) {
				result[i] = v1[i] - v2[i];
			}
		}

	}

	__device__ void scale_device(float s, float* v1, float(&result)[3]) {
		//float * result = new float[3];
		for (int i = 0; i < 3; i++) {
			result[i] = s*v1[i];
		}
		//return result;
	}
	__device__ float len_device(float* v1) {
		float result = dot_device(v1, v1);
		return sqrt(result);
	}


	__device__ float triangle_inter_device(float * v1, float * v2, float* v3,
		float *dir, float * vertex) {
		float e[2][3];
		//e[0] = new float[3];
		//e[1] = new float[3];
		for (int i = 0; i < 3; i++) {
			e[0][i] = v2[i] - v1[i];
			e[1][i] = v3[i] - v1[i];
		}

		//float * pvec;
		//pvec = cross_device(e[1], dir);
		float  pvec[3];

		cross_device(e[1], dir, pvec);

		float det = dot_device(e[0], pvec);
		if (det < 1e-8 && det > -1e-8) {
			//std::cout << "parrel";
			return 0;
		}

		float inv_det = 1 / det;
		float tvec[3];
		plus_device(vertex, v1, tvec, false);
		float u = dot_device(tvec, pvec) * inv_det;
		if (u < 0 || u > 1) {
			return 0;
		}
		float qvec[3];
		cross_device(tvec, e[0], qvec);
		float v = dot_device(dir, qvec) * inv_det;
		if (v < 0 || u + v > 1) {
			return 0;
		}
		float t = inv_det * dot_device(e[1], qvec);

		//calculate the intersection points.
		//float * inter_point =plus( vertex , scale(t,dir));
		float result[3];
		scale_device(t, dir, result);
		return len_device(result);

	}

	__device__ float cal_dist_on_node_device(float * p, float * dir, float * pNormal,
		float * p1X, float * p1Y, float * p1Z,
		float * p2X, float * p2Y, float * p2Z,
		float * p3X, float * p3Y, float * p3Z,
		float * faceNX, float * faceNY, float * faceNZ, int size) {
		float temp = 10000000;
		for (int faceid = 0; faceid < size; faceid++) {
			float sum = faceNX[faceid] * pNormal[0] + faceNY[faceid] * pNormal[1] + faceNZ[faceid] * pNormal[2];
			float angle = acos(sum)* 180.0 / M_PIf;

			float v1[3] = { p1X[faceid] , p1Y[faceid] , p1Z[faceid] };
			float v2[3] = { p2X[faceid] , p2Y[faceid] , p2Z[faceid] };
			float v3[3] = { p3X[faceid] , p3Y[faceid] , p3Z[faceid] };
			//--Check whether the ray and the faces normal is above 90 degree.
			if (angle>90) {
				//int distance = 0;
				float distance = triangle_inter_device(v1, v2, v3, dir, p);
				if (distance != 0 && distance < temp) {
					temp = distance;
				}
			}
		}
		if (temp == 10000000)
			return -1;
		return temp;
	}

	__device__ int kd_right_device(int index) {
		return ((index + 1) * 2);
	}
	__device__ int kd_left_device(int index) {
		return ((index + 1) * 2) - 1;
	}

	__global__ void intersect_gpu_greater(float *pt1x, float *pt1y, float *pt1z,
		float *pt2x, float *pt2y, float *pt2z,
		float *pt3x, float *pt3y, float *pt3z,
		float *fNx, float *fNy, float *fNz,
		float *vx, float*vy, float*vz,
		float * normal_x, float *normal_y, float *normal_z,
		float *dirX, float *dirY, float *dirZ,
		float *minDist,
		int nv, int nRay,
		int nf_max, int nf_min)
	{
		//int i = threadIdx.x;

		//int i_face = blockIdx.x * blockDim.x + threadIdx.x;
		//int i_rays = blockIdx.z * blockDim.z + threadIdx.z;
		//int i_vertices = blockIdx.y * blockDim.y + threadIdx.y;

		int i_face = blockIdx.x * blockDim.x + threadIdx.x;
		int i_rays = blockIdx.z * blockDim.z + threadIdx.z;
		int i_vertices = (blockIdx.y * blockDim.y + threadIdx.y) + nf_min;

		if (i_face < nv && i_vertices<nf_max && i_rays<nRay) {
			int j = i_vertices*nRay + i_rays;
			float v1[3] = { pt1x[i_face] , pt1y[i_face] , pt1z[i_face] };
			float v2[3] = { pt2x[i_face] , pt2y[i_face] , pt2z[i_face] };
			float v3[3] = { pt3x[i_face] , pt3y[i_face] , pt3z[i_face] };


			float dir[3] = { dirX[j] , dirY[j] , dirZ[j] };
			float vertex[3] = { vx[i_vertices] , vy[i_vertices] , vz[i_vertices] };

			float sum = fNx[i_face] * normal_x[i_vertices] + fNy[i_face] * normal_y[i_vertices] + fNz[i_face] * normal_z[i_vertices];
			//float sum = fNx[1] * normal_x[1] + fNy[1] * normal_y[1] + fNz[1] * normal_z[1];
			float angle = acos(sum)* 180.0 / M_PIf;
			//float result;
			if (angle>90) {
				//int distance = 0;
				float distance = triangle_inter_device(v1, v2, v3, dir, vertex);
				if (distance != 0) {
					if (minDist[j] == 0) {
						minDist[j] = distance;
					}
					//else {
					//	if (distance < minDist[j]) {
					//		minDist[j] = distance;
					//	}
					//}
				}
				//else if (distance < minDist[j]) {
				//	minDist[j] = distance;
				//}
			}

		}

	}
	

	__global__ void intersect_gpu_kdtree(
		float *pt1x, float *pt1y, float *pt1z,
		float *pt2x, float *pt2y, float *pt2z,
		float *pt3x, float *pt3y, float *pt3z,
		float *fNx, float *fNy, float *fNz,
		float *split, int * split_axis, int* tri_count, int * tri_index,
		float *aabb_max, float *aabb_min,
		float *vx, float*vy, float*vz,
		float * normal_x, float *normal_y, float *normal_z,
		float *dirX, float *dirY, float *dirZ,
		float *minDist,
		int nv, int nRay,
		int n_kd)
	{

		int i_vertices = blockIdx.x * blockDim.x + threadIdx.x;
		int i_rays = (blockIdx.y * blockDim.y + threadIdx.y);

		if (i_vertices<nv && i_rays<nRay) {
			int j = i_vertices*nRay + i_rays;
			float dir[3] = { dirX[j] , dirY[j] , dirZ[j] };
			float vertex[3] = { vx[i_vertices] , vy[i_vertices] , vz[i_vertices] };
			float pNormal[3] = { normal_x[i_vertices],normal_y[i_vertices] , normal_z[i_vertices] };
			float t_i0[3];
			float t_i1[3];

			for (int i = 0; i < 3; i++) {
				if (vertex[i] == 0 || dir[i] == 0) {
					vertex[i] = 0; dir[i] = 0;
				}
				if (dir[i] >= 0) {
					t_i0[i] = (aabb_min[i] - vertex[i]) / dir[i];
					t_i1[i] = (aabb_max[i] - vertex[i]) / dir[i];
				}
				else {
					t_i0[i] = (aabb_max[i] - vertex[i]) / dir[i];
					t_i1[i] = (aabb_min[i] - vertex[i]) / dir[i];
				}
			}
			float t_min = max(t_i0[0], max(t_i0[1], t_i0[2]));
			float t_max = min(t_i1[0], min(t_i1[1], t_i1[2]));
			float global_max = t_max;
			int index = 0; float distance = -1;
			if (t_max > t_min) {
				int id_stack[64];
				int * stackPtr = id_stack; *stackPtr++ = -1;
				float t_min_stack[64]; float t_max_stack[64];
				float * stackPtr_tmin = t_min_stack; *stackPtr_tmin++ = -1.0;
				float * stackPtr_tmax = t_max_stack; *stackPtr_tmax++ = -1.0;
				while (distance <= 0 && index <= n_kd) {
					if (tri_count[index] == 0) {
						//--KD restart:
						if (t_max == global_max)
							break;
						else {
							t_min = t_max;
							t_max = global_max;
							index = 0;
						}
						//NULL node;
						index = *--stackPtr; t_min = *--stackPtr_tmin; t_max = *--stackPtr_tmax;
						if (index == -1)
							break;
					}
					else if (tri_count[index] != -1) {
						//is leaf 
						int leaf_size = tri_count[index];
						int leaf_index = tri_index[index];
						distance = cal_dist_on_node_device(vertex, dir, pNormal,
							&pt1x[leaf_index], &pt1y[leaf_index], &pt1z[leaf_index],
							&pt2x[leaf_index], &pt2y[leaf_index], &pt2z[leaf_index],
							&pt3x[leaf_index], &pt3y[leaf_index], &pt3z[leaf_index],
							&fNx[leaf_index], &fNy[leaf_index], &fNz[leaf_index],
							leaf_size);

						//change the index by POPPing the value in the stack
						//if (*(stackPtr - 1) != -1) {
						//--kd restart search:
						if (t_max == global_max)
							break;
						else {
							t_min = t_max;
							t_max = global_max;
							index = 0;
						}


						//index = *--stackPtr; t_min = *--stackPtr_tmin; t_max = *--stackPtr_tmax;
						////}
						//if (index == -1)
						//	break;

					}
					else {
						//NOT leaf
						int first_index, second_index;
						int axis = split_axis[index];
						float thit = (split[index] - vertex[axis]) / dir[axis];
						if ((split[index] - vertex[axis]) >= 0) {
							//first_index = kd_left_device(index); second_index = kd_right_device(index);
							first_index = ((index + 1) * 2) - 1; second_index = ((index + 1) * 2);
						}
						else {
							//first_index = kd_right_device(index); second_index = kd_left_device(index);
							first_index = ((index + 1) * 2); second_index = ((index + 1) * 2) - 1;
						}

						if (thit >= t_max || thit < 0)
							index = first_index;
						else if (thit <= t_min)
							index = second_index;
						else {
							//--push the index of Right node into Stack
							//*stackPtr++ = second_index;
							//*stackPtr_tmin++ = thit; *stackPtr_tmax++ = t_max;
							index = first_index;
							t_max = thit;

						}

					}
				}

				minDist[j] = distance;
			}

		}

	}

	//__global__ void octree_ray_gpu(Octree * d_oct) {
	//	printf("%d ", d_oct->text);
	//}

	//void make_cone(vec origin, vec dir, vector<vec>& output) {
	//	output.clear();
	//	//float halfAngle = 60;

	//	//int totalRings = 20;
	//	//int intervals = 5;	

	//	output.resize(RAYS_PER_ANGLE*ANGLE_INTERVALS + 1);

	//	//find perpendicular  vector as rotation axis

	//	vec rotateAxis = vec(1, 1, (-dir[0] - dir[1]) / dir[2]);
	//	if (dir[2] == 0) {
	//		rotateAxis = vec(0, 0, 1);
	//	}
	//	normalize(rotateAxis);
	//	normalize(dir);

	//	for (int j = 0; j < ANGLE_INTERVALS; j++) {
	//		vec ray1 = rotation(dir, CONE_ANGLE - (15 * j), rotateAxis);
	//		for (int i = 0; i < RAYS_PER_ANGLE; i++) {
	//			vec rayI = rotation(ray1, 360 * i / RAYS_PER_ANGLE, dir);
	//			output[(j*RAYS_PER_ANGLE) + i] = rayI;
	//		}
	//	}
	//	output[output.size() - 1] = dir;
	//}

	void make_cone_total(vec origin, vec dir, vector<vec>& output) {
		//output.clear();
		//float halfAngle = 60;

		//int totalRings = 20;
		//int intervals = 5;	
		//vector<vec> output;
		//output.resize(RAYS_PER_ANGLE*ANGLE_INTERVALS + 1);

		//find perpendicular  vector as rotation axis

		vec rotateAxis = vec(1, 1, (-dir[0] - dir[1]) / dir[2]);
		if (dir[2] == 0) {
			rotateAxis = vec(0, 0, 1);
		}
		normalize(rotateAxis);
		normalize(dir);

		for (int j = 0; j < ANGLE_INTERVALS; j++) {
			vec ray1 = rotation(dir, CONE_ANGLE - (15 * j), rotateAxis);
			for (int i = 0; i < RAYS_PER_ANGLE; i++) {
				vec rayI = rotation(ray1, 360 * i / RAYS_PER_ANGLE, dir);
				output.push_back(rayI);
			}
		}
		output.push_back(dir);
	}

	//Goal: Calculate the distance
	/*input:
	vectex:	vector start point
	dir:		vecter direction
	v1,v2,v3:	points of the face
	*/
	/*Moller¡§CTrumbore intersection algorithm
	p + t * d = (1-u-v) * p0 + u * p1 + v * p2
	t , u v are parameters
	p0 p1 p2 are apex of triangles
	p are point on ray, d are dir

	*/


	template < class T>
	float sdf_stat_mean(vector<T> input) {
		float sum = std::accumulate(input.begin(), input.end(), 0.0f);
		float mean = sum / input.size();
		float sq_sum = std::inner_product(input.begin(), input.end(), input.begin(), 0.0);
		float stdev = std::sqrt(sq_sum / input.size() - mean * mean);

		return mean;
	}

	template < class T>
	float sdf_stat_mean(vector<T> input, float median, float stdev) {
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


	void TriMesh::need_sdf_gpu_1D() {
		std::clock_t begin = clock();
		std::cout << std::endl << "computing the SDF using gpu without loop: ";
		if (sdf_brute_gpu.size() == vertices.size())
			return;
		sdf_brute_gpu.resize(vertices.size());



		typedef long long llong;
		int nf = faces.size();
		int nv = vertices.size();
		int n_rays = 81;
		llong total_N = (llong)nf*(llong)nv *(llong)n_rays;
		std::cout << total_N << std::endl;
		llong nv_per_grid = ((llong)LONG_MAX) / (llong)(nf*n_rays);
		llong nf_per_grid = ((llong)LONG_MAX) / (llong)(nv*n_rays);
		if (nf>nf_per_grid) {
			std::cout << "greater than max block size" << endl;
		}
		else {
			nf_per_grid = nf;
		}
		if (nv>nv_per_grid) {
			std::cout << "greater than max block size" << endl;
		}
		else {
			nv_per_grid = nv;
		}


		//std::cout << "the max f " << nf_per_grid <<endl;
		//long n_blocks = LONG_MAX - 1;
		//	nf_per_grid = 2000;
		long n_grids = nf / nf_per_grid;

		long n_grids_v = nv / nv_per_grid;
		//std::cout << "the number of grids " << n_grids << endl;
		//	std::cout << n_grids << std::endl;


		need_faceNormals();

		//---allocate the points of faces and face normals on Memory
		float * p_faces_1pt[3], *p_faces_2pt[3], *p_faces_3pt[3], *p_face_normal[3];
		for (int j = 0; j < 3; j++) {
			p_faces_1pt[j] = new float[nf];
			p_faces_2pt[j] = new float[nf];
			p_faces_3pt[j] = new float[nf];
			p_face_normal[j] = new float[nf];

		}
		for (int j = 0; j < 3; j++) {
			for (int i = 0; i < nf; i++) {
				p_faces_1pt[j][i] = vertices[faces[i][0]][j];
				p_faces_2pt[j][i] = vertices[faces[i][1]][j];
				p_faces_3pt[j][i] = vertices[faces[i][2]][j];
				p_face_normal[j][i] = faceNormals[i][j];
			}
		}

		float * d_pt1[3], *d_pt2[3], *d_pt3[3], *d_face_normals[3];


		//allocate the points of faces and face normals on CUDA
		for (int i = 0; i < 3; i++) {
			cudaMalloc((void**)&(d_pt1[i]), nf * sizeof(float));
			cudaMemcpy(d_pt1[i], p_faces_1pt[i], nf * sizeof(float), cudaMemcpyHostToDevice);

			cudaMalloc((void**)&(d_pt2[i]), nf * sizeof(float));
			cudaMemcpy(d_pt2[i], p_faces_2pt[i], nf * sizeof(float), cudaMemcpyHostToDevice);

			cudaMalloc((void**)&(d_pt3[i]), nf * sizeof(float));
			cudaMemcpy(d_pt3[i], p_faces_3pt[i], nf * sizeof(float), cudaMemcpyHostToDevice);

			cudaMalloc((void**)&(d_face_normals[i]), nf * sizeof(float));
			cudaMemcpy(d_face_normals[i], p_face_normal[i], nf * sizeof(float), cudaMemcpyHostToDevice);
		}



		//make the total cones(nv * 81 )
		//	std::vector<vec> rays;
		std::vector<vec> vertices_rays;
		for (int i = 0; i < nv; i++) {
			vec inwardNormal = inwardNormals[i];
			vec vertex = vertices[i];
			make_cone_total(vertex, inwardNormal, vertices_rays);

		}
		int total_nRay = nv* n_rays;
		float * total_rayX = new float[total_nRay], *total_rayY = new float[total_nRay], *total_rayZ = new float[total_nRay];
		float *host_result = new float[total_nRay];
		//initialize the HOST pointers of total rays of vertices.
		for (int i = 0; i < total_nRay; i++) {
			total_rayX[i] = vertices_rays[i][0];
			total_rayY[i] = vertices_rays[i][1];
			total_rayZ[i] = vertices_rays[i][2];
			host_result[i] = 0;
		}
		//initialize the HOST pointers of vertices.
		float *host_v_x = new float[nv], *host_v_y = new float[nv], *host_v_z = new float[nv];
		float *host_normal_x = new float[nv], *host_normal_y = new float[nv], *host_normal_z = new float[nv];
		for (int i = 0; i < nv; i++) {
			host_v_x[i] = vertices[i][0];
			host_v_y[i] = vertices[i][1];
			host_v_z[i] = vertices[i][2];
			host_normal_x[i] = normals[i][0];
			host_normal_y[i] = normals[i][1];
			host_normal_z[i] = normals[i][2];
		}
		float * d_total_rayX, *d_total_rayY, *d_total_rayZ;
		float *d_v_x, *d_v_y, *d_v_z;
		float *d_normal_x, *d_normal_y, *d_normal_z;
		float * d_result;

		//locate and copy total rays
		gpuErrchk(cudaMalloc((void**)&(d_total_rayX), total_nRay * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_total_rayY), total_nRay * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_total_rayZ), total_nRay * sizeof(float)));
		gpuErrchk(cudaMemcpy(d_total_rayX, total_rayX, total_nRay * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_total_rayY, total_rayY, total_nRay * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_total_rayZ, total_rayZ, total_nRay * sizeof(float), cudaMemcpyHostToDevice));

		//locate and copy vertices
		gpuErrchk(cudaMalloc((void**)&(d_v_x), nv * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_v_y), nv * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_v_z), nv * sizeof(float)));
		gpuErrchk(cudaMemcpy(d_v_x, host_v_x, nv * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_v_y, host_v_y, nv * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_v_z, host_v_z, nv * sizeof(float), cudaMemcpyHostToDevice));

		//locate and copy normals
		gpuErrchk(cudaMalloc((void**)&(d_normal_x), nv * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_normal_y), nv * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_normal_z), nv * sizeof(float)));
		gpuErrchk(cudaMemcpy(d_normal_x, host_normal_x, nv * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_normal_y, host_normal_y, nv * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_normal_z, host_normal_z, nv * sizeof(float), cudaMemcpyHostToDevice));

		//locate result
		gpuErrchk(cudaMalloc((void**)&(d_result), total_nRay * sizeof(float)));
		gpuErrchk(cudaMemcpy(d_result, host_result, total_nRay * sizeof(float), cudaMemcpyHostToDevice));

		for (int i = 0; i <= (int)n_grids_v; i++) {
			//-----
			//initial the block size and block numbers.
			int blockSizeX = 128;
			int blockNumX = (nf + blockSizeX - 1) / blockSizeX;
			//blockNumX = nf / blockSizeX;

			int blockSizeY = 8;
			int blockNumY = (nv_per_grid + blockSizeY - 1) / blockSizeY;


			int blockSizeZ = 1;
			int blockNumZ = (n_rays + blockSizeZ - 1) / blockSizeZ;

			dim3 thread(blockSizeX, blockSizeY, blockSizeZ);
			dim3 blockNum(blockNumX, blockNumY, blockNumZ);
			//int i = 0;
			int nf_min = i*nv_per_grid;
			int nf_max = ((i + 1)*nv_per_grid) - 1;
			intersect_gpu_greater << <blockNum, thread >> > (d_pt1[0], d_pt1[1], d_pt1[2],
				d_pt2[0], d_pt2[1], d_pt2[2],
				d_pt3[0], d_pt3[1], d_pt3[2],
				d_face_normals[0], d_face_normals[1], d_face_normals[2],
				d_v_x, d_v_y, d_v_z,
				d_normal_x, d_normal_y, d_normal_z,
				d_total_rayX, d_total_rayY, d_total_rayZ,
				d_result,
				nv, n_rays, nf_max, nf_min);
			//if (i % (n_grids / 10) == 0)
			//	std::cout << "Finish " << i / (n_grids / 10) << "0% of ht SDF" << std::endl;
			//	gpuErrchk(cudaMemcpy(host_result, d_result, total_nRay * sizeof(float), cudaMemcpyDeviceToHost));
			//gpuErrchk(cudaMemcpy(d_result, host_result, total_nRay * sizeof(float), cudaMemcpyHostToDevice));
		}


		gpuErrchk(cudaMemcpy(host_result, d_result, total_nRay * sizeof(float), cudaMemcpyDeviceToHost));

		for (int i = 0; i < nv; i++) {
			vector<float>result;
			for (int j = 0; j < n_rays; j++) {
				if (host_result[i*n_rays + j] != 0)
					result.push_back(host_result[i*n_rays + j]);
			}
			if (result.size() == 0) {
				sdf_brute_gpu[i] = 0;// sdf_stat_mean(result);
			}
			else {
				sdf_brute_gpu[i] = sdf_stat_mean(result);
			}

		}

		int size0_dist_count = 0;


		//cuda free memory
		cudaFree(d_normal_x); cudaFree(d_normal_y); cudaFree(d_normal_z);
		cudaFree(d_v_x); cudaFree(d_v_y); cudaFree(d_v_z);
		cudaFree(d_total_rayX); cudaFree(d_total_rayY); cudaFree(d_total_rayZ);
		cudaFree(d_result);
		for (int i = 0; i < 3; i++) {
			cudaFree(d_pt1[i]);
			cudaFree(d_pt2[i]);
			cudaFree(d_pt3[i]);
			cudaFree(d_face_normals[i]);
		}
		std::cout << double(clock() - begin) / CLOCKS_PER_SEC << " s" << std::endl;
	}

	// use cuda to accelaterated sdf calculation.

	void TriMesh::need_sdf_kd_tree_gpu() {
		if (sdf.size() == vertices.size())
			return;
		sdf.resize(vertices.size());

		clock_t begin = clock();
		std::cout << std::endl << "computing the shape diameters using KD TREE with GPU: ";
		need_faceNormals();
		KD_tree kd_tree_pointer;
		trimesh::buildKDTree(kd_tree_pointer, faces, vertices, faceNormals);
		KD_tree_array * kd_array = KDTreeToArray(kd_tree_pointer);

		sdf_brute_gpu.resize(vertices.size());


		int n_tris = kd_array->pt1X_1d.size();
		typedef long long llong;
		int N_kd_array = kd_array->pt1X.size();
		int nv = vertices.size();
		int n_rays = 81;
		llong total_N = (llong)nv *(llong)n_rays;
		std::cout << total_N << std::endl;
		llong nv_per_grid = ((llong)LONG_MAX) / (llong)(n_rays);


		if (nv>nv_per_grid) {
			std::cout << "greater than max block size" << endl;
		}
		else {
			nv_per_grid = nv;
		}


		//std::cout << "the max f " << nf_per_grid <<endl;
		//long n_blocks = LONG_MAX - 1;
		//	nf_per_grid = 2000;


		long n_grids_v = nv / nv_per_grid;
		//std::cout << "the number of grids " << n_grids << endl;
		//	std::cout << n_grids << std::endl;





		float * d_pt1[3], *d_pt2[3], *d_pt3[3], *d_face_normals[3];


		//make the total cones(nv * 81 )
		//	std::vector<vec> rays;
		std::vector<vec> vertices_rays;
		for (int i = 0; i < nv; i++) {
			vec inwardNormal = inwardNormals[i];
			vec vertex = vertices[i];
			make_cone_total(vertex, inwardNormal, vertices_rays);

		}
		int total_nRay = nv* n_rays;
		float * total_rayX = new float[total_nRay], *total_rayY = new float[total_nRay], *total_rayZ = new float[total_nRay];
		float *host_result = new float[total_nRay];
		//initialize the HOST pointers of total rays of vertices.
		for (int i = 0; i < total_nRay; i++) {
			total_rayX[i] = vertices_rays[i][0];
			total_rayY[i] = vertices_rays[i][1];
			total_rayZ[i] = vertices_rays[i][2];
			host_result[i] = 0;
		}
		//initialize the HOST pointers of vertices.
		float *host_v_x = new float[nv], *host_v_y = new float[nv], *host_v_z = new float[nv];
		float *host_normal_x = new float[nv], *host_normal_y = new float[nv], *host_normal_z = new float[nv];
		for (int i = 0; i < nv; i++) {
			host_v_x[i] = vertices[i][0];
			host_v_y[i] = vertices[i][1];
			host_v_z[i] = vertices[i][2];
			host_normal_x[i] = normals[i][0];
			host_normal_y[i] = normals[i][1];
			host_normal_z[i] = normals[i][2];
		}
		float * d_total_rayX, *d_total_rayY, *d_total_rayZ;
		float *d_v_x, *d_v_y, *d_v_z;
		float *d_normal_x, *d_normal_y, *d_normal_z;
		float * d_result;

		//locate and copy total rays
		gpuErrchk(cudaMalloc((void**)&(d_total_rayX), total_nRay * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_total_rayY), total_nRay * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_total_rayZ), total_nRay * sizeof(float)));
		gpuErrchk(cudaMemcpy(d_total_rayX, total_rayX, total_nRay * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_total_rayY, total_rayY, total_nRay * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_total_rayZ, total_rayZ, total_nRay * sizeof(float), cudaMemcpyHostToDevice));

		//locate and copy vertices
		gpuErrchk(cudaMalloc((void**)&(d_v_x), nv * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_v_y), nv * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_v_z), nv * sizeof(float)));
		gpuErrchk(cudaMemcpy(d_v_x, host_v_x, nv * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_v_y, host_v_y, nv * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_v_z, host_v_z, nv * sizeof(float), cudaMemcpyHostToDevice));

		//locate and copy normals
		gpuErrchk(cudaMalloc((void**)&(d_normal_x), nv * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_normal_y), nv * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_normal_z), nv * sizeof(float)));
		gpuErrchk(cudaMemcpy(d_normal_x, host_normal_x, nv * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_normal_y, host_normal_y, nv * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_normal_z, host_normal_z, nv * sizeof(float), cudaMemcpyHostToDevice));

		//locate result
		gpuErrchk(cudaMalloc((void**)&(d_result), total_nRay * sizeof(float)));
		gpuErrchk(cudaMemcpy(d_result, host_result, total_nRay * sizeof(float), cudaMemcpyHostToDevice));

		//------------------------
		//--Locate the attributes in KD tree.
		float * d_p1x, *d_p1y, *d_p1z;
		float * d_p2x, *d_p2y, *d_p2z;
		float * d_p3x, *d_p3y, *d_p3z;
		float * d_fNx, *d_fNy, *d_fNz;



		gpuErrchk(cudaMalloc((void**)&(d_p1x), n_tris * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_p1y), n_tris * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_p1z), n_tris * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_p2x), n_tris * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_p2y), n_tris * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_p2z), n_tris * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_p3x), n_tris * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_p3y), n_tris * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_p3z), n_tris * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_fNx), n_tris * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_fNy), n_tris * sizeof(float)));
		gpuErrchk(cudaMalloc((void**)&(d_fNz), n_tris * sizeof(float)));

		gpuErrchk(cudaMemcpy(d_p1x, &(kd_array->pt1X_1d[0]), n_tris * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_p1y, &(kd_array->pt1Y_1d[0]), n_tris * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_p1z, &(kd_array->pt1Z_1d[0]), n_tris * sizeof(float), cudaMemcpyHostToDevice));

		gpuErrchk(cudaMemcpy(d_p2x, &(kd_array->pt2X_1d[0]), n_tris * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_p2y, &(kd_array->pt2Y_1d[0]), n_tris * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_p2z, &(kd_array->pt2Z_1d[0]), n_tris * sizeof(float), cudaMemcpyHostToDevice));

		gpuErrchk(cudaMemcpy(d_p3x, &(kd_array->pt3X_1d[0]), n_tris * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_p3y, &(kd_array->pt3Y_1d[0]), n_tris * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_p3z, &(kd_array->pt3Z_1d[0]), n_tris * sizeof(float), cudaMemcpyHostToDevice));

		gpuErrchk(cudaMemcpy(d_fNx, &(kd_array->fNomarlX_1d[0]), n_tris * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_fNy, &(kd_array->fNomarlY_1d[0]), n_tris * sizeof(float), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_fNz, &(kd_array->fNomarlZ_1d[0]), n_tris * sizeof(float), cudaMemcpyHostToDevice));


		float * d_split;
		int * d_triCount, *d_split_axis;
		int * d_tri_index;
		gpuErrchk(cudaMalloc((void**)  &(d_split), sizeof(float) *N_kd_array));
		gpuErrchk(cudaMalloc((void**)  &(d_triCount), sizeof(int) * N_kd_array));
		gpuErrchk(cudaMalloc((void**)  &(d_split_axis), sizeof(int) * N_kd_array));
		gpuErrchk(cudaMalloc((void**)  &(d_tri_index), sizeof(int) * N_kd_array));

		gpuErrchk(cudaMemcpy(d_split, &(kd_array->split[0]), sizeof(float) * N_kd_array, cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_triCount, &(kd_array->triCount[0]), sizeof(int)* N_kd_array, cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_split_axis, &(kd_array->split_axis[0]), sizeof(int)* N_kd_array, cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_tri_index, &(kd_array->triIndex[0]), sizeof(int)* N_kd_array, cudaMemcpyHostToDevice));

		float *d_max, *d_min;
		gpuErrchk(cudaMalloc((void**)  &(d_max), sizeof(float) * 3));
		gpuErrchk(cudaMalloc((void**)  &(d_min), sizeof(float) * 3));
		gpuErrchk(cudaMemcpy(d_max, &(kd_array->max), sizeof(float) * 3, cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_min, &(kd_array->min), sizeof(float) * 3, cudaMemcpyHostToDevice));

		int blockSizeX = 128;
		int blockNumX = (nv + blockSizeX - 1) / blockSizeX;
		//blockNumX = nf / blockSizeX;

		int blockSizeY = 8;
		int blockNumY = (n_rays + blockSizeY - 1) / blockSizeY;

		dim3 thread(blockSizeX, blockSizeY);
		dim3 blockNum(blockNumX, blockNumY);
		std::cout << double(clock() - begin) / CLOCKS_PER_SEC << " s" << std::endl;
		intersect_gpu_kdtree << <blockNum, thread >> > (d_p1x, d_p1y, d_p1z,
			d_p2x, d_p2y, d_p2z,
			d_p3x, d_p3y, d_p3z,
			d_fNx, d_fNy, d_fNz,
			d_split, d_split_axis, d_triCount, d_tri_index,
			d_max, d_min,
			d_v_x, d_v_y, d_v_z,
			d_normal_x, d_normal_y, d_normal_z,
			d_total_rayX, d_total_rayY, d_total_rayZ,
			d_result,
			nv, n_rays, N_kd_array

			);

		gpuErrchk(cudaMemcpy(host_result, d_result, total_nRay * sizeof(float), cudaMemcpyDeviceToHost));

		for (int i = 0; i < nv; i++) {
			vector<float>result;
			for (int j = 0; j < n_rays; j++) {
				if (host_result[i*n_rays + j] >= 0)
					result.push_back(host_result[i*n_rays + j]);
				else {
					result.push_back(bsphere.r * 2);
				}

			}
			if (result.size() == 0) {
				sdf[i] = 0;// sdf_stat_mean(result);
			}
			else {
				sdf[i] = sdf_stat_mean(result);
			}

		}
		std::cout << double(clock() - begin) / CLOCKS_PER_SEC << " s" << std::endl;
	}


}//end of namespace

 