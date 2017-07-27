/*
Szymon Rusinkiewicz
Princeton University

KDtree.cc
A K-D tree for points, with limited capabilities (find nearest point to
a given point, or to a ray).
*/

#include <cstring>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
#include "octree.h"
#include "TriMesh.h"
#include "mempool.h"
using namespace std;


namespace trimesh {

	// Small utility fcns - including them keeps this file independent of Vec.h
	static inline float sqr(float x)
	{
		return x*x;
	}

	static inline float dist2(const float *x, const float *y)
	{
		return sqr(x[0] - y[0]) + sqr(x[1] - y[1]) + sqr(x[2] - y[2]);
	}

	static inline float dot(const float *x, const float *y)
	{
		return (x[0] * y[0]) + (x[1] * y[1]) + (x[2] * y[2]);
	}

	static inline float dist2ray2(const float *x, const float *p, const float *d)
	{
		float xp0 = x[0] - p[0], xp1 = x[1] - p[1], xp2 = x[2] - p[2];
		return sqr(xp0) + sqr(xp1) + sqr(xp2) -
			sqr(xp0*d[0] + xp1*d[1] + xp2*d[2]);
	}
	// reference:https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
	static inline bool rayToPlane(vec nP , vec centerP , vec lDir , vec l0 , vec &intersect) {
		nP = normalize(nP);
		lDir = normalize(lDir);
		if ((nP ^ lDir) >1e-6 || (nP ^ lDir)<-0.000001) {
			float t = -((l0 - centerP) ^ nP) / (lDir ^ nP);
			if (t >= 0.0f) {
				intersect = l0 + t*(lDir); //interect points
			}
			return(t >= 0);
		}
		else return false;
	}

	static inline float distance_ray_face(vec p, vec dir, vec pNormal, vec faceNormal, vector<vec> pts)
	{
		float angle = acos(pNormal ^ faceNormal)* 180.0 / M_PIf;
		if (angle > 90) {
			vec v1 = pts[0]; vec v2 = pts[1]; vec v3 = pts[2];

			float distance = triangle_inter(p, dir, v1, v2, v3);
			return distance;
		}
		else
			return 0.0f;
	}
	//An Efficient Parametric Algorithm for Octree Traversal
	static int first_sub_node(vector<float> t_i0, vector<float> t_i_m) {
		//to do
		//--- check the which one of the 8 sub node
//		std::distance(t_i_m, std::max_element(t_i_m.begin(), t_i_m.end()));

		int index = std::max_element(t_i0.begin(), t_i0.end()) - t_i0.begin();

		string bits = "000";
		int x = stoi(bits.c_str(), nullptr, 2);

		switch (index) {
			//First entry plane is YZ plane, with t_x0 is the max
		case 0: {
			if (t_i_m[1] < t_i0[0])
				bits.replace(1, 1, "1"); 
			if (t_i_m[2] < t_i0[0])
				bits.replace(2, 1, "1");
			break;
		}
		case 1: {//First entry plane is XZ plane, with t_y0 is the max
			if (t_i_m[0] < t_i0[1])
				bits.replace(0, 1, "1");
			if (t_i_m[2] < t_i0[1])
				bits.replace(2, 1, "1");
			break;
		}
		case 2: {//First entry plane is XY plane, with t_z0 is the max
			if (t_i_m[0] < t_i0[2])
				bits.replace(0, 1, "1");
			if (t_i_m[1] < t_i0[2])
				bits.replace(1, 1, "1");
			break;
		}
		}
		//std::reverse(bits.begin(), bits.end());
		x = stoi(bits.c_str(), nullptr, 2);
		return x;
	}

	static int next_sub_node(int current_nodeID, vector<float> t_i1) {
		//return the next node ID using current Node ID and its t_i1
		//The index should be 1-7. The 8 means no more node is interacted.
		int index = std::min_element(t_i1.begin(), t_i1.end()) - t_i1.begin();
	
		switch (index) {
		case 0: {
			switch (current_nodeID) {
			case 0: return 4;
			case 1: return 5;
			case 2: return 6;
			case 3: return 7;
			case 4:
			case 5:
			case 6:
			case 7: return 8;
			}
			break;
		}
		case 1: {
			switch (current_nodeID) {
			case 0: return 2;
			case 1: return 3;
			case 2: 
			case 3: return 8;
			case 4: return 6;
			case 5: return 7;
			case 6:
			case 7: return 8;
			}
			break;
		}

		case 2: {
			switch (current_nodeID) {
			case 0: return 1;
			case 1: return 8;
			case 2: return 3;
			case 3: return 8;
			case 4: return 5;
			case 5: return 8;
			case 6: return 7;
			case 7: return 8;
			}
			break;
		}
		}
		return 8;
	}

	

	//---------------------


	class Octree::Node {
	public:
		enum Octants { R_UP_F = 0 ,R_DOWN_F=1, R_UP_B=2 , R_DOWN_B=3,
		L_UP_F=4, L_DOWN_F=5, L_UP_B=6, L_DOWN_B=7,
		ERROR = -1};

		std::vector<Node*> children;
		Node *father =NULL;
		
		bool isLeaf;
		bool isEmpty=false;
		bool has_faces_branch;
		int nPts;
		int nFaces;

		int totaln;

		vector<vec> points;
		vector<vector<vec>> node_face_pts;
		vector<vec> node_faceNormal;

		vector<int> node_face_id;

		struct {
			vec max;
			vec min;
			vec center;
			float r;
			float planeYZ;
			float planeXZ;
			float planeXY;
			//The max and min X,Y,Z of...
			//... front back,right,left,up,down
			//vector<vec> planesMax;
			//vector<vec> planesMin;
		} bBox;

		//member function
		Node() {  };

		Node(vector<vec> pts, int n , float xmax, float ymax , float zmax , float xmin , float ymin , float zmin);
		Node(vector<vector<vec>> facesPts , vec centroid, float edgeLen , std::vector<vec> fNormals , Node* father);
		Node(vector<vector<vec>> facesPts, vector<int> facesIDs, vec centroid, float edgeLen, std::vector<vec> fNormals);

		~Node();
		//void find_closest_to_ray(Traversal_Info &ti);
		Octants find_octant_point(const float *p);
		void distance_within_octants(const float *p , const float *dir ,  Octants octantP, Traversal_Info &info);

		bool is_pt_inside_cube(const float *p);


		void ray_intersect_cube(const float *p, const float *dir, Traversal_Info &info);

		void find_point_leaf(const float *p , Node *&result);
		void find_ray_intersect_originInside_bbox(const float *p, const float *dir, Traversal_Info &info);

		void find_ray_intersect_triangle_top_down(const float *p, const float *dir, const float *pNormal, Traversal_Info &info);
		void find_ray_intersect_triangle_bottom_up(const float *p, const float *dir, const float *pNormal, Traversal_Info &info);

		bool is_ray_intersect_with_octant(vec l0, vec lDir);

		void traverse_all_leaves(Traversal_Info &info){
			if (isEmpty)
				return;
			if (isLeaf) {
				info.bBox_maxs.push_back(bBox.max);
				info.bBox_mins.push_back(bBox.min);
				//cout << "reach the leaf";
				if (nFaces <= MAX_FACES_PER_NODE)
					info.equal_max_faces_count++;
				info.leaves_count++;
				return;
			}
			if (has_faces_branch) {
				info.branch_count++;
			}
			for (int i = 0; i < 8; i++) {
				children[i]->traverse_all_leaves(info);
			}
		}



		void ray_crossing_octants(vector<float> t_i0 , vector<float> t_i1,Traversal_Info &info) {
			//The first distance calculated is the closest distance
			if (info.closest_d != MIN_DIST_INIT)
				return;
			//--calculate the origin t1s and the min t1 for exit plane 
			vec originV = vec(info.origin_p[0], info.origin_p[1], info.origin_p[2]);
			vec originRay = vec(info.origin_dir[0], info.origin_dir[1], info.origin_dir[2]);
			vector<float> t_new_i1(3);
			for (int i = 0; i < 3; i++) {
				if(originRay[i]<0)
					t_new_i1[i] = (this->bBox.min[i] - originV[i]) / originRay[i];
				else
					t_new_i1[i] = (this->bBox.max[i] - originV[i]) / originRay[i];
			}
			float t1_new_min = *std::min_element(t_new_i1.begin(), t_new_i1.end());
			//Use minimum t1 to make sure no octant behind the origin points of the ray.
			if (t1_new_min <= 0) {
				return;
			}
			if (isEmpty)
				return;

			if (isLeaf) {
				//info.bBox_maxs.push_back(bBox.max);
				//info.bBox_mins.push_back(bBox.min);
				this->cal_sdf_on_node(info);
				return;
			}		
			//if has branch:


			//calculate the parameter t for hiting the half of X,Y,Z
			vector<float> t_i_m(3);
			for (int i = 0; i < 3; i++) {
				t_i_m[i] = 0.5*(t_i0[i] + t_i1[i]);
			}
			
			int currentNode_id = first_sub_node(t_i0, t_i_m);
			int first_node_id = currentNode_id;
			unsigned char negative_convert_bits = info.ray_neg_convert_bits;
		//	std::cout<<(currentNode_id ^ negative_convert_bits);
			//cout << "the Number of cubes intercted: " << (currentNode_id ^ negative_convert_bits) << endl;
			vector<float> t_i1_sub_node(3);
			vector<float> t_i0_sub_node(3);
			bool checked = false;
			do {
				switch (currentNode_id) {
				case 0: {
					t_i1_sub_node[0] = t_i_m[0]; t_i1_sub_node[1] = t_i_m[1]; t_i1_sub_node[2] = t_i_m[2];
					t_i0_sub_node[0] = t_i0[0]; t_i0_sub_node[1] = t_i0[1]; t_i0_sub_node[2] = t_i0[2];
					break;
				}
				case 1: {
					t_i1_sub_node[0] = t_i_m[0]; t_i1_sub_node[1] = t_i_m[1]; t_i1_sub_node[2] = t_i1[2];
					t_i0_sub_node[0] = t_i0[0]; t_i0_sub_node[1] = t_i0[1]; t_i0_sub_node[2] = t_i_m[2];
					break;
				}
				case 2: {
					t_i1_sub_node[0] = t_i_m[0]; t_i1_sub_node[1] = t_i1[1]; t_i1_sub_node[2] = t_i_m[2];
					t_i0_sub_node[0] = t_i0[0]; t_i0_sub_node[1] = t_i_m[1]; t_i0_sub_node[2] = t_i0[2];
					break;
				}
				case 3: {
					t_i1_sub_node[0] = t_i_m[0]; t_i1_sub_node[1] = t_i1[1]; t_i1_sub_node[2] = t_i1[2];
					t_i0_sub_node[0] = t_i0[0]; t_i0_sub_node[1] = t_i_m[1]; t_i0_sub_node[2] = t_i_m[2];
					break;
				}
				case 4: {
					t_i1_sub_node[0] = t_i1[0]; t_i1_sub_node[1] = t_i_m[1]; t_i1_sub_node[2] = t_i_m[2];
					t_i0_sub_node[0] = t_i_m[0]; t_i0_sub_node[1] = t_i0[1]; t_i0_sub_node[2] = t_i0[2];
					break;
				}
				case 5: {
					t_i1_sub_node[0] = t_i1[0]; t_i1_sub_node[1] = t_i_m[1]; t_i1_sub_node[2] = t_i1[2];
					t_i0_sub_node[0] = t_i_m[0]; t_i0_sub_node[1] = t_i0[1]; t_i0_sub_node[2] = t_i_m[2];
					break;
				}
				case 6: {
					t_i1_sub_node[0] = t_i1[0]; t_i1_sub_node[1] = t_i1[1]; t_i1_sub_node[2] = t_i_m[2];
					t_i0_sub_node[0] = t_i_m[0]; t_i0_sub_node[1] = t_i_m[1]; t_i0_sub_node[2] = t_i0[2];
					break;
				}
				case 7: {
					t_i1_sub_node[0] = t_i1[0]; t_i1_sub_node[1] = t_i1[1]; t_i1_sub_node[2] = t_i1[2];
					t_i0_sub_node[0] = t_i_m[0]; t_i0_sub_node[1] = t_i_m[1]; t_i0_sub_node[2] = t_i_m[2];
					break;
				}

				}
				children[currentNode_id ^ negative_convert_bits]->ray_crossing_octants(t_i0_sub_node, t_i1_sub_node, info);
				currentNode_id = next_sub_node(currentNode_id, t_i1_sub_node);
				if ( currentNode_id<8 && !checked) {
					bool a = check_which_first(info, std::min_element(t_i1_sub_node.begin(), t_i1_sub_node.end()) - t_i1_sub_node.begin());
				//	cout << "exit plant" << std::min_element(t_i1_sub_node.begin(), t_i1_sub_node.end()) - t_i1_sub_node.begin() << endl;
					//cout << "pass " << a << endl;;
					//
					if (a && has_faces_branch && info.closest_d == MIN_DIST_INIT) {
						this->cal_sdf_on_node(info);
						checked = true;
					}
				}
				//if (count == 1 && info.closest_d == MIN_DIST_INIT) {
				//	if (has_faces_branch) {
				//		//info.bBox_maxs.push_back(bBox.max);
				//		//info.bBox_mins.push_back(bBox.min);
				//		this->cal_sdf_on_nodes(info);
				//	}
				//}
		
			} while (currentNode_id < 8 && currentNode_id >= 0);

		}
		//static float distance_ray_face(vec p, vec dir, vec pNormal,vec node_faceNormal, vector<vec> pts );

		void cal_sdf_on_node(Traversal_Info &info) {
			vec vNormal = vec(info.ptNormal[0], info.ptNormal[1], info.ptNormal[2]);
			vec vertex = vec(info.origin_p[0], info.origin_p[1], info.origin_p[2]);
			vec ray = vec(info.origin_dir[0], info.origin_dir[1], info.origin_dir[2]);
			info.cal_count += nFaces;
		//	int faceid = -1;
			if (this->nFaces != 0) {
				#pragma omp parallel for
				for (int i = 0; i < nFaces; i++) { //Loop through all the faces in this Node
					float distance = distance_ray_face(vertex, ray, vNormal, node_faceNormal[i], this->node_face_pts[i]);
					if (distance != 0) {
						if (info.closest_d == MIN_DIST_INIT) {
							info.closest_d = distance;
							//info.faceID = k;
						}
						else if (distance < info.closest_d) {
							info.closest_d = distance;
							//faceId = k;
						}
					}
				}
			}
		}

		bool check_which_first(Traversal_Info &info, int index) {
			vec originV = vec(info.origin_p[0], info.origin_p[1], info.origin_p[2]);
			vec originRay = vec(info.origin_dir[0], info.origin_dir[1], info.origin_dir[2]);

			vector<float> t_new_i1(3);

			vec max = vec(this->bBox.max[0], this->bBox.max[1],this->bBox.max[2]);
			switch (index) {
			case 0: { // plane YZ
				max[0] = (bBox.max[0] + bBox.min[0]) / 2;
				break;
			}
			case 1: { // plane XZ
				max[1] = (bBox.max[1] + bBox.min[1]) / 2;
				break;
			}
			case 2: { //plane XY
				max[2] = (bBox.max[2] + bBox.min[2]) / 2;
				break;
			}
			}

			for (int i = 0; i < 3; i++) {
				if (originRay[i]<0)
					t_new_i1[i] = (this->bBox.min[i] - originV[i]) / originRay[i];
				else
					t_new_i1[i] = (max[i] - originV[i]) / originRay[i];
			}
			float t1_new_min = *std::min_element(t_new_i1.begin(), t_new_i1.end());
			return t1_new_min >= 0;

		}
		

	};// END the declaration of NODE

	Octree::Node::Node(vector<vec> pts, int n, 
		float xmax, float ymax, float zmax,
		float xmin, float ymin, float zmin) {
		//Find bounding box and initialize the bounding box attribute
		float planeYZ = 0.5f * (xmin + xmax);
		float planeXZ = 0.5f * (ymin + ymax);
		float planeXY = 0.5f * (zmin + zmax);
		//initialize the bounding box attribute
		bBox.max[0] = xmax; bBox.max[1] = ymax; bBox.max[2] = zmax;
		bBox.min[0] = xmin; bBox.min[1] = ymin; bBox.min[2] = zmin;
		bBox.planeYZ = planeYZ; bBox.planeXZ = planeXZ; bBox.planeXY = planeXY;
		bBox.center[0] = planeYZ; bBox.center[1] = planeXZ; bBox.center[2] = planeXY;

		float dx = xmax - xmin;
		float dy = ymax - ymin;
		float dz = zmax - zmin;
		float edgeLen = dx;
		// decide the longest
		if (dx >= dy && dx >= dz) {
			edgeLen = dx;
			bBox.r = dx;
		}
		if (dy >= dx && dy >= dz) {
			edgeLen = dy;
			bBox.r = dy;
		}
		if (dz >= dy && dz >= dx) {
			edgeLen = dz;
			bBox.r = dz;
		}
		isEmpty = false;
		//is leaf node
		if (n <= MAX_PTS_PER_NODE) {
			nPts = n;
			points.resize(n);
			points = pts;
			isLeaf = true;
			if (n == 0)
				isEmpty = true;
			return;
		}

		//is intermediate node
		nPts = 0;
		isLeaf = false;

		// Partition
		vector<vec> v[8];
		for (int i = 0; i < n; i++) {

			float pX = pts[i][0];
			float pY = pts[i][1];
			float pZ = pts[i][2];
			//right up front
			if (pX >= planeYZ && pY >= planeXZ && pZ >= planeXY) {
				v[0].push_back(pts[i]);
			}
			//right down front
			if (pX >= planeYZ && pY < planeXZ && pZ >= planeXY) {
				v[1].push_back(pts[i]);
			}

			//right up back
			if (pX >= planeYZ && pY >= planeXZ && pZ < planeXY) {
				v[2].push_back(pts[i]);
			}
			//right down back
			if (pX >= planeYZ && pY < planeXZ && pZ < planeXY) {
				v[3].push_back(pts[i]);
			}
			//----------------
			//left up front
			if (pX < planeYZ && pY >= planeXZ && pZ >= planeXY) {
				v[4].push_back(pts[i]);
			}
			//left down front
			if (pX < planeYZ && pY < planeXZ && pZ >= planeXY) {
				v[5].push_back(pts[i]);
			}

			//left up back
			if (pX < planeYZ && pY >= planeXZ && pZ < planeXY) {
				v[6].push_back(pts[i]);
			}
			//left down back
			if (pX < planeYZ && pY < planeXZ && pZ < planeXY) {
				v[7].push_back(pts[i]);
			}

		}
		//initialize the children bounding box
		float xMid = 0.5f * (xmin + xmax);float yMid = 0.5f * (ymin + ymax);float zMid = 0.5f * (zmin + zmax);
		float xMaxs[8], yMaxs[8],zMaxs[8];
		float xMins[8], yMins[8],zMins[8];
		//right up front
		xMaxs[0] = xmax; yMaxs[0] = ymax; zMaxs[0] = zmax;
		xMins[0] = xMid; yMins[0] = yMid; zMins[0] = zMid;
		//right down front
		xMaxs[1] = xmax; yMaxs[1] = yMid; zMaxs[1] = zmax;
		xMins[1] = xMid; yMins[1] = ymin; zMins[1] = zMid;
		//right up back
		xMaxs[2] = xmax; yMaxs[2] = ymax; zMaxs[2] = zMid;
		xMins[2] = xMid; yMins[2] = yMid; zMins[2] = zmin;
		//right down back
		xMaxs[3] = xmax; yMaxs[3] = yMid; zMaxs[3] = zMid;
		xMins[3] = xMid; yMins[3] = ymin; zMins[3] = zmin;
		//----------------
		//left up front
		xMaxs[4] = xMid; yMaxs[4] = ymax; zMaxs[4] = zmax;
		xMins[4] = xmin; yMins[4] = yMid; zMins[4] = zMid;
		//left down front
		xMaxs[5] = xMid; yMaxs[5] = yMid; zMaxs[5] = zmax;
		xMins[5] = xmin; yMins[5] = ymin; zMins[5] = zMid;
		//left up back
		xMaxs[6] = xMid; yMaxs[6] = ymax; zMaxs[6] = zMid;
		xMins[6] = xmin; yMins[6] = yMid; zMins[6] = zmin;
		//left down back
		xMaxs[7] = xMid; yMaxs[7] = yMid; zMaxs[7] = zMid;
		xMins[7] = xmin; yMins[7] = ymin; zMins[7] = zmin;

		children.resize(8);
		for (int i = 0; i < 8; i++) {
			//children[i] = new Node(v[i], v[i].size());
			children[i] = new Node(v[i], v[i].size(),
				xMaxs[i], yMaxs[i], zMaxs[i],
				xMins[i], yMins[i], zMins[i]);
			children[i]->father = this;

		}

	}
	

	//Build the Node using faces of the 3D model
	//No intermediate node.
	Octree::Node::Node(vector<vector<vec>> facesPts, vec centroid, float edgeLen ,vector<vec> fNormals, Node* Father)
	{
		vector<int> facesIds;

		//cout << "the size " << facesPts.size();
		int nf = facesPts.size();
		totaln = nf;
		float planeYZ = centroid[0];
		float planeXZ = centroid[1];
		float planeXY = centroid[2];

		for (int i = 0; i < 3; i++) {
			bBox.max[i] = centroid[i] + edgeLen / 2; 
			bBox.min[i] = centroid[i] - edgeLen / 2;
			bBox.center[i] = centroid[i];
		}
		//--is leaf node
		this->father = Father;
		bool isSame =this->father!=NULL && this->father->father != NULL && nf == Father->totaln && nf == Father->father->totaln;

		if (nf <= MAX_FACES_PER_NODE || isSame) {
			nFaces = nf;
			node_face_pts.resize(nf);
			node_faceNormal.resize(nf);
			//#pragma omp parallel for
			for (int i = 0; i < nf; i++) {
				for (int j = 0; j < 3; j++) {
					node_face_pts[i].push_back(facesPts[i][j]);
				}
				node_faceNormal[i] = fNormals[i];
			}
			has_faces_branch = false;
			isLeaf = true;
			if (nf == 0)
				isEmpty = true;
			return;
		}


		//--is branch node
		isEmpty = false;
		isLeaf = false;
		vector<vector<vec>> faces_in_octants[8];
		vector<vec> facesNormal_in_octants[8];

		vector<vector<vec>> faces_crossing_planes;
		vector<vec> faces_crossing_normals;
		//#pragma omp parallel for
		for (int i = 0; i < nf; i++) {
			//Decide whether the triangle is in the partition.
			bool tri_in_octants = true;
			bool all_pts_greater = true;
			bool all_pts_less = true;
			vector<vec> pts_per_tri = facesPts[i];
			vec faceNormal = fNormals[i];
			//bool ruf = true; bool rdf = true; bool rub = true; bool rdb = true;
			//bool luf = true; bool ldf = true; bool lub = true; bool ldb = true;
			bool ruf = false; bool rdf =false; bool rub = false; bool rdb = false;
			bool luf = false; bool ldf =false; bool lub = false; bool ldb = false;
			for (int j = 0; j < 3; j++) {
				// x,y,z of 3 points of the triangle
				float pX = pts_per_tri[j][0];
				float pY = pts_per_tri[j][1];
				float pZ = pts_per_tri[j][2];
				//right up front
				ruf = ruf||(pX >=planeYZ && pY >=planeXZ && pZ >=planeXY);
				rdf = rdf||(pX >=planeYZ && pY <=planeXZ && pZ >=planeXY);
				rub = rub||(pX >=planeYZ && pY >=planeXZ && pZ <=planeXY);
				rdb = rdb||(pX >=planeYZ && pY <=planeXZ && pZ <=planeXY);
   				   				   
				luf = luf||(pX <=planeYZ && pY >=planeXZ && pZ >=planeXY);
				ldf = ldf||(pX <=planeYZ && pY <=planeXZ && pZ >=planeXY);
				lub = lub|| (pX <=planeYZ && pY >=planeXZ && pZ <=planeXY);
				ldb = ldb||(pX <=planeYZ && pY <=planeXZ && pZ <=planeXY);
			}
			//right up front
			if (ruf) {
				faces_in_octants[7].push_back(pts_per_tri);
				facesNormal_in_octants[7].push_back(faceNormal);
			}
			//right down front
			  if (rdf) {
				faces_in_octants[5].push_back(pts_per_tri);
				facesNormal_in_octants[5].push_back(faceNormal);
			}

			//right up back
			  if (rub) {
				faces_in_octants[6].push_back(pts_per_tri);
				facesNormal_in_octants[6].push_back(faceNormal);
			}
			//right down back
			  if (rdb) {
				faces_in_octants[4].push_back(pts_per_tri);
				facesNormal_in_octants[4].push_back(faceNormal);
			}
			//----------------
			//left up front
			  if (luf) {
				faces_in_octants[3].push_back(pts_per_tri);
				facesNormal_in_octants[3].push_back(faceNormal);
			}
			//left down front
			  if (ldf) {
				faces_in_octants[1].push_back(pts_per_tri);
				facesNormal_in_octants[1].push_back(faceNormal);
			}

			//left up back
			  if (lub) {
				faces_in_octants[2].push_back(pts_per_tri);
				facesNormal_in_octants[2].push_back(faceNormal);
			}
			//left down back
			  if (ldb) {
				faces_in_octants[0].push_back(pts_per_tri);
				facesNormal_in_octants[0].push_back(faceNormal);
			}
			//else {
			//	faces_crossing_planes.push_back(pts_per_tri);
			//	faces_crossing_normals.push_back(faceNormal);
			//}
			//other wise:the triangle is crossing the dividing planes

		}

		//-- is branch Node with faces (faces which intersect with the partition planes)
		if (faces_crossing_planes.size() != 0) {
			has_faces_branch = true;
			int n_crossing = faces_crossing_planes.size();
			node_face_pts.resize(n_crossing);
			node_faceNormal.resize(n_crossing);
			//#pragma omp parallel for
			for (int i = 0; i < n_crossing; i++) {
				for (int j = 0; j < 3; j++) {
					node_face_pts[i].push_back(faces_crossing_planes[i][j]);
				}
				node_faceNormal[i] = faces_crossing_normals[i];
			}
			nFaces = n_crossing;
			if (n_crossing == nf)
				isLeaf = true;
			if(n_crossing>1000)
				cout << "the size " << n_crossing;
		}
		else {
			has_faces_branch = false;
		}

		//Calculate the parameters (centroid and edge length) for childern
		//Use the centroid and edge length is this Node.
		vec centroids[8];
		float cx = centroid[0]; float cy = centroid[1]; float cz = centroid[2]; float next_edge_len = edgeLen / 2; float c_diff = edgeLen / 4;
		//right up front
		centroids[7] = vec(cx + c_diff, cy + c_diff, cz + c_diff);
		//right down front
		centroids[5] = vec(cx + c_diff, cy - c_diff, cz + c_diff);
		//right up back
		centroids[6] = vec(cx + c_diff, cy + c_diff, cz - c_diff);
		//right down back
		centroids[4] = vec(cx + c_diff, cy - c_diff, cz - c_diff);
		//----------------
		//left up front
		centroids[3] = vec(cx - c_diff, cy + c_diff, cz + c_diff);
		//left down front
		centroids[1] = vec(cx - c_diff, cy - c_diff, cz + c_diff);
		//left up back
		centroids[2] = vec(cx - c_diff, cy + c_diff, cz - c_diff);
		//left down back
		centroids[0] = vec(cx - c_diff, cy - c_diff, cz - c_diff);

		children.resize(8);
		//for (int i = 0; i < 8; i++) {
		//	//if (faces_in_octants[i].size() != 0) {
		//	if (faces_in_octants[i].size() == nf) {
		//		if (nf > 100)
		//			cout << "size of child : " << nf;
		//		nFaces = nf;
		//		node_face_pts.resize(nf);
		//		node_faceNormal.resize(nf);
		//		//#pragma omp parallel for
		//		for (int i = 0; i < nf; i++) {
		//			for (int j = 0; j < 3; j++) {
		//				node_face_pts[i].push_back(facesPts[i][j]);
		//			}
		//			node_faceNormal[i] = fNormals[i];
		//		}
		//		has_faces_branch = false;
		//		isLeaf = true;
		//		if (nf == 0)
		//			isEmpty = true;
		//		return;
		//	}
		//	//cout << "no changes";
		//}
		for (int i = 0; i < 8; i++) {
			children[i] = new Node(faces_in_octants[i], centroids[i], next_edge_len, facesNormal_in_octants[i] , this);
			//children[i]->father = this;
			//}
		}
	}


	//Build the Node using faces of the 3D model 
	// With the face ID.
	Octree::Node::Node(vector<vector<vec>> facesPts, vector<int> facesIDs, vec centroid, float edgeLen, vector<vec> fNormals)
	{
	

		//cout << "the size " << facesPts.size();
		int nf = facesPts.size();
		totaln = nf;
		float planeYZ = centroid[0];
		float planeXZ = centroid[1];
		float planeXY = centroid[2];

		for (int i = 0; i < 3; i++) {
			bBox.max[i] = centroid[i] + edgeLen / 2;
			bBox.min[i] = centroid[i] - edgeLen / 2;
			bBox.center[i] = centroid[i];
		}
		//--is leaf node

		if (nf <= MAX_FACES_PER_NODE) {
			nFaces = nf;
			node_face_pts.resize(nf);
			node_faceNormal.resize(nf);
			node_face_id.resize(nf);
			//#pragma omp parallel for
			for (int i = 0; i < nf; i++) {
				for (int j = 0; j < 3; j++) {
					node_face_pts[i].push_back(facesPts[i][j]);
				}
				node_faceNormal[i] = fNormals[i];
				node_face_id[i] = facesIDs[i];
			}
			has_faces_branch = false;
			isLeaf = true;
			if (nf == 0)
				isEmpty = true;
			return;
		}


		//--is branch node
		isEmpty = false;
		isLeaf = false;
		vector<vector<vec>> faces_in_octants[8];
		vector<vec> facesNormal_in_octants[8];
		vector<int> faceID_in_octants[8];

		vector<vector<vec>> faces_crossing_planes;
		vector<vec> faces_crossing_normals;
		vector<int> faceID_crossing;
		
		
		//#pragma omp parallel for
		for (int i = 0; i < nf; i++) {
			//Decide whether the triangle is in the partition.
			bool tri_in_octants = true;
			bool all_pts_greater = true;
			bool all_pts_less = true;


			vector<vec> pts_per_tri = facesPts[i];
			vec face_normal = fNormals[i];
			int faceID = i;
			bool ruf = true; bool rdf = true; bool rub = true; bool rdb = true;
			bool luf = true; bool ldf = true; bool lub = true; bool ldb = true;
			for (int j = 0; j < 3; j++) {
				// x,y,z of 3 points of the triangle
				float pX = pts_per_tri[j][0];
				float pY = pts_per_tri[j][1];
				float pZ = pts_per_tri[j][2];
				//right up front
				ruf &= pX > planeYZ && pY > planeXZ && pZ > planeXY;
				rdf &= pX > planeYZ && pY < planeXZ && pZ > planeXY;
				rub &= pX > planeYZ && pY > planeXZ && pZ < planeXY;
				rdb &= pX > planeYZ && pY < planeXZ && pZ < planeXY;

				luf &= pX < planeYZ && pY > planeXZ && pZ > planeXY;
				ldf &= pX < planeYZ && pY < planeXZ && pZ > planeXY;
				lub &= pX < planeYZ && pY > planeXZ && pZ < planeXY;
				ldb &= pX < planeYZ && pY < planeXZ && pZ < planeXY;
			}
			//right up front
			if (ruf) {
				faces_in_octants[7].push_back(pts_per_tri);
				facesNormal_in_octants[7].push_back(face_normal);
				faceID_in_octants[7].push_back(faceID);
			}
			//right down front
			else if (rdf) {
				faces_in_octants[5].push_back(pts_per_tri);
				facesNormal_in_octants[5].push_back(face_normal);
				faceID_in_octants[5].push_back(faceID);
			}

			//right up back
			else if (rub) {
				faces_in_octants[6].push_back(pts_per_tri);
				facesNormal_in_octants[6].push_back(face_normal);
				faceID_in_octants[6].push_back(faceID);
			}
			//right down back
			else if (rdb) {
				faces_in_octants[4].push_back(pts_per_tri);
				facesNormal_in_octants[4].push_back(face_normal);
				faceID_in_octants[4].push_back(faceID);
			}
			//----------------
			//left up front
			else if (luf) {
				faces_in_octants[3].push_back(pts_per_tri);
				facesNormal_in_octants[3].push_back(face_normal);
				faceID_in_octants[3].push_back(faceID);
			}
			//left down front
			else if (ldf) {
				faces_in_octants[1].push_back(pts_per_tri);
				facesNormal_in_octants[1].push_back(face_normal);
				faceID_in_octants[1].push_back(faceID);
			}

			//left up back
			else if (lub) {
				faces_in_octants[2].push_back(pts_per_tri);
				facesNormal_in_octants[2].push_back(face_normal);
				faceID_in_octants[2].push_back(faceID);
			}
			//left down back
			else if (ldb) {
				faces_in_octants[0].push_back(pts_per_tri);
				facesNormal_in_octants[0].push_back(face_normal);
				faceID_in_octants[0].push_back(faceID);
			}
			else {
				faces_crossing_planes.push_back(pts_per_tri);
				faces_crossing_normals.push_back(face_normal);
				faceID_crossing.push_back(faceID);
			}
			//other wise:the triangle is crossing the dividing planes

		}

		//-- is branch Node with faces (faces which intersect with the partition planes)
		if (faces_crossing_planes.size() != 0) {
			has_faces_branch = true;
			int n_crossing = faces_crossing_planes.size();
			node_face_pts.resize(n_crossing);
			node_faceNormal.resize(n_crossing);
			node_face_id.resize(n_crossing);
			//#pragma omp parallel for
			for (int i = 0; i < n_crossing; i++) {
				for (int j = 0; j < 3; j++) {
					node_face_pts[i].push_back(faces_crossing_planes[i][j]);
				}
				node_faceNormal[i] = faces_crossing_normals[i];
				node_face_id[i] = faceID_crossing[i];
			}
			nFaces = n_crossing;
			if (n_crossing == nf)
				isLeaf = true;
			if (n_crossing>1000)
				cout << "the size " << n_crossing;
		}
		else {
			has_faces_branch = false;
		}

		//Calculate the parameters (centroid and edge length) for childern
		//Use the centroid and edge length is this Node.
		vec centroids[8];
		float cx = centroid[0]; float cy = centroid[1]; float cz = centroid[2]; float next_edge_len = edgeLen / 2; float c_diff = edgeLen / 4;
		//right up front
		centroids[7] = vec(cx + c_diff, cy + c_diff, cz + c_diff);
		//right down front
		centroids[5] = vec(cx + c_diff, cy - c_diff, cz + c_diff);
		//right up back
		centroids[6] = vec(cx + c_diff, cy + c_diff, cz - c_diff);
		//right down back
		centroids[4] = vec(cx + c_diff, cy - c_diff, cz - c_diff);
		//----------------
		//left up front
		centroids[3] = vec(cx - c_diff, cy + c_diff, cz + c_diff);
		//left down front
		centroids[1] = vec(cx - c_diff, cy - c_diff, cz + c_diff);
		//left up back
		centroids[2] = vec(cx - c_diff, cy + c_diff, cz - c_diff);
		//left down back
		centroids[0] = vec(cx - c_diff, cy - c_diff, cz - c_diff);

		children.resize(8);
		for (int i = 0; i < 8; i++) {
			//if (faces_in_octants[i].size() != 0) {

			children[i] = new Node(faces_in_octants[i] , faceID_in_octants[i], centroids[i], next_edge_len, facesNormal_in_octants[i]);
			children[i]->father = this;
			//}
		}
	}


	//void Octree::Node::box_of_point(Octree::Traversal_Info &ti) const {
	//	//check if inside the Bounding box of this node
	//	bool insideBoundBox = true;
	//	for (int i = 0; i < 3; i++) {
	//		insideBoundBox = insideBoundBox && (ti.origin_p[i] <= bBox.max[i]) && (ti.origin_p[i] >= bBox.min[i]);
	//	} 

	//	if (insideBoundBox&&isLeaf) {
	//		cout << "level: " << ti.level << endl;;
	//		cout << "The bound box is  ";
	//		cout <<endl<< "max  ";
	//		for (int i = 0; i < 3; i++) {
	//			cout << " " << bBox.max[i] << ",";
	//		}
	//		cout << endl << "min ";
	//		for (int i = 0; i < 3; i++) {
	//			cout << " " << bBox.min[i]<<",";
	//		}
	//		cout << endl;
	//		ti.level--;
	//		return;
	//	}
	//	if (insideBoundBox) {
	//		ti.level++;
	//		for (int i = 0; i < 8; i++)
	//			children[i]->find_closest_to_ray(ti);
	//	}

	//}

	Octree::Node::Octants Octree::Node::find_octant_point(const float *p) {
		float pX = p[0], pY = p[1], pZ = p[2];
		//float planeYZ = bBox.planeYZ, planeXZ = bBox.planeXZ, planeXY = bBox.planeXY;
		float planeYZ = bBox.center[0], planeXZ = bBox.center[1], planeXY = bBox.center[2];
		//right up front
		if (pX >= planeYZ && pY >= planeXZ && pZ >= planeXY) {
			return R_UP_F;
		}
		//right down front
		if (pX >= planeYZ && pY < planeXZ && pZ >= planeXY) {
			return R_DOWN_F;
		}

		//right up back
		if (pX >= planeYZ && pY >= planeXZ && pZ < planeXY) {
			return R_UP_B;
		}
		//right down back
		if (pX >= planeYZ && pY < planeXZ && pZ < planeXY) {
			return R_DOWN_B;
		}
		//----------------
		//left up front
		if (pX < planeYZ && pY >= planeXZ && pZ >= planeXY) {
			return L_UP_F;
		}
		//left down front
		if (pX < planeYZ && pY < planeXZ && pZ >= planeXY) {
			return L_DOWN_F;
		}

		//left up back
		if (pX < planeYZ && pY >= planeXZ && pZ < planeXY) {
			return L_UP_B;
		}
		//left down back
		if (pX < planeYZ && pY < planeXZ && pZ < planeXY) {
			return L_DOWN_B;
		}
		return ERROR;
	}

	void Octree::Node::distance_within_octants(const float *p, const float *dir, 
		Octree::Node::Octants octantP, Octree::Traversal_Info &info) {
		
		vec origin = vec(0, 0, 0);
		vec l0 = vec(p[0], p[1], p[2]);
		vec lDir = vec(dir[0], dir[1], dir[2]);
		// calculate line intersection with Plane Y Z
		vec nYZ = vec(bBox.center[0], 0, 0);// bBox.center - origin;
		vec pInterPlaneYZ;
		if (rayToPlane(nYZ, bBox.center, lDir, l0, pInterPlaneYZ)) {
			std::cout << "ray intersect with planeYZ";
		}

		// calculate line intersection with Plane X Z
		vec nXZ = vec(0, bBox.center[1], 0);// bBox.center - origin;
		vec pInterPlaneXZ;
		if (rayToPlane(nXZ, bBox.center, lDir, l0, pInterPlaneXZ)) {
			std::cout << "ray intersect with planeXZ";
		}

		// calculate line intersection with Plane X Y
		vec nXY = vec(0, 0, bBox.center[2]);// bBox.center - origin;
		vec pInterPlaneXY;
		if (rayToPlane(nXY, bBox.center, lDir, l0, pInterPlaneXY)) {
			std::cout << "ray intersect with planeXY";
		}

		for (int i = 0; i < 8; i++) {
			//calculate the distance between ray to no empty neighbours .
			if (children[i]->isEmpty != 0 && i != octantP) {
				switch (i) {
				case 0: {// right up front
				}
				}
			}
		}

	}
	//A function to check the point is inside this cube or not.
	bool Octree::Node::is_pt_inside_cube(const float * p)
	{
		bool insideBoundBox = true;
		for (int i = 0; i < 3; i++) {
			insideBoundBox = insideBoundBox && (p[i] <= bBox.max[i]) && (p[i] >= bBox.min[i]);
		}
		return insideBoundBox;
	}

	//check the whether the ray is intersect with the six faces of the Node region
	// p: origin point of the ray. Dir: direction vector
	void Octree::Node::ray_intersect_cube(const float *p, const float *dir, Traversal_Info &info) {
		if (isEmpty)
			return;

		vec l0 = vec(p[0], p[1], p[2]);
		vec lDir = vec(dir[0], dir[1], dir[2]);
		//The max and min X,Y,Z of...
		//... front back,right,left,up,down
		vector<vec> planesMax;
		vector<vec> planesMin;
		planesMax.push_back( vec(bBox.max[0], bBox.max[1], bBox.max[2]));
		planesMax.push_back( vec(bBox.max[0], bBox.max[1], bBox.min[2]));
		planesMax.push_back(  vec(bBox.max[0], bBox.max[1], bBox.max[2]));
		planesMax.push_back(  vec(bBox.min[0], bBox.max[1], bBox.max[2]));
		planesMax.push_back(  vec(bBox.max[0], bBox.max[1], bBox.max[2]));
		planesMax.push_back(  vec(bBox.max[0], bBox.min[1], bBox.max[2]));

		planesMin.push_back(  vec(bBox.min[0], bBox.min[1], bBox.max[2]));
		planesMin.push_back(  vec(bBox.min[0], bBox.min[1], bBox.min[2]));
		planesMin.push_back(  vec(bBox.max[0], bBox.min[1], bBox.min[2]));
		planesMin.push_back(  vec(bBox.min[0], bBox.min[1], bBox.min[2]));
		planesMin.push_back(  vec(bBox.min[0], bBox.max[1], bBox.min[2]));
		planesMin.push_back(  vec(bBox.min[0], bBox.min[1], bBox.min[2]));

		std::vector<vec> plane;
		//front and back planes
		plane.push_back(  vec(0, 0, bBox.max[2]));
		plane.push_back(  vec(0, 0, bBox.min[2]));
		//right and left planes
		plane.push_back(  vec(bBox.max[0], 0, 0));
		plane.push_back(  vec(bBox.min[0], 0, 0));
		//up and down planes;
		plane.push_back(  vec(0, bBox.max[1], 0));
		plane.push_back(  vec(0, bBox.min[1], 0));

		std::vector<vec> planeN;
		//front and back planes
		planeN.push_back(vec(0, 0, 1));
		planeN.push_back(vec(0, 0, -1));
		//right and left planes
		planeN.push_back(vec(1, 0, 0));
		planeN.push_back(vec(-1, 0, 0));
		//up and down planes;
		planeN.push_back(vec(0, 1, 0));
		planeN.push_back(vec(0, -1, 0));

		bool isIntersect = false;
		vector<vec> intersects;
		//Check the intersection between ray and 6 faces
		for (int i = 0; i < 6; i++) {
			vec intersect;
			bool temp = rayToPlane(planeN[i], plane[i], lDir, l0, intersect);
			bool inside = false;
			if (temp) {
				inside = intersect[0] <= planesMax[i][0] && intersect[1] <= planesMax[i][1] && intersect[2] <= planesMax[i][2];
				inside = inside && intersect[0] >= planesMin[i][0] && intersect[1] >= planesMin[i][1] && intersect[2] >= planesMin[i][2];
				if(inside)
					intersects.push_back(intersect);
			}
			isIntersect = isIntersect || temp || inside;
		}
		//IF it is leaf and intersect with the cube
		//calculate the distance from start of the ray to the intersect point
		if (isLeaf && isIntersect) {
			for (int i = 0; i < intersects.size(); i++) {
				vec d = intersects[i] - l0;
				float dis = len(intersects[i] - l0);
				if (dis < info.closest_d) {
					//Update the infomation of the intersecting octant
					info.closest_d = dis;
					info.closest_Pt = intersects[i];
					info.max = bBox.max;
					info.min = bBox.min;
				}
			}

		}
		if(isIntersect && !isLeaf){
			info.level += 1;
			for (int i = 0; i < 8; i++) {
				info.iChild = i;
				children[i]->ray_intersect_cube(p, dir, info);
			}
			info.level -= 1;
		}
		return;
	}

	//Find the cube(leaf) where the point are inside the bounding box of octree
	void Octree::Node::find_point_leaf(const float * p , Node *&result)
	{
		//check if inside the Bounding box of this node
		bool insideBoundBox = true;
		for (int i = 0; i < 3; i++) {
			insideBoundBox = insideBoundBox && (p[i] <= bBox.max[i]) && (p[i] >= bBox.min[i]);
		}
		Octants octantP;
		Octants direction;
		bool isOctantLeaf;
		bool isOctantEmpty;
		if (insideBoundBox && !isLeaf) {
			//check  in which octants in children
			octantP = this->find_octant_point(p);

			isOctantLeaf = this->children[octantP]->isLeaf == true;
			isOctantEmpty = this->children[octantP]->isEmpty == true;
			if (!isOctantEmpty) {
				//if result octant is not a leaf
				// traverse down to the Octant
				children[octantP]->find_point_leaf(p , result);
			}
			else {
				result = this;
				return;
			}
			//if (isOctantLeaf) {
			//	//check which other octant have the shortest distance
			//	//The octant must have node;

			//}
		}
		if (insideBoundBox && isLeaf) {
			result = this;
			return;
		}

	}
	// Find the closet bounding box with point in inside.
	void Octree::Node::find_ray_intersect_originInside_bbox(const float * p, const float * dir, Traversal_Info & info)
	{
		//if the node is ROOT, do nothing and return.
		if (this->father == NULL)
			return;
		Node *father = this->father;
		for (int i = 0; i < 8; i++) {
			//check whether ray is intersect with its siblings
			if (!(father->children[i]->isEmpty) && !(father->children[i]->is_pt_inside_cube(p))) {
				std::cout << i<<" ";
				father->children[i]->ray_intersect_cube(p,dir,info);
			}

		}
	}

	void Octree::Node::find_ray_intersect_triangle_top_down(const float * p, const float * dir, const float * pNormal, Traversal_Info & info)
	{
		vec vNormal = vec(pNormal[0], pNormal[1], pNormal[2]);
		vec vertex = vec(p[0], p[1], p[2]);
		vec ray = vec(dir[0], dir[1], dir[2]);
		int faceid = -1;
		if (this->nFaces != 0) {
			for (int i = 0; i < nFaces; i++) { //Loop through all the faces in this Node
				float distance = distance_ray_face(vertex,ray,vNormal, node_faceNormal[i], node_face_pts[i]);
				if (distance != 0) {
					if (info.closest_d == MIN_DIST_INIT) {
						info.closest_d = distance;
						faceid = i;
						//info.faceID = k;
					}
					else if (distance < info.closest_d) {
							info.closest_d = distance;
							faceid = i;
							//faceId = k;
					}
				}
			}
		}
		//Loop through all the intersected siblings.
		//for (int i = 0; i < 8; i++) {
		//	if (this->father->children[i] != this &&
		//		this->father->children[i]->is_ray_intersect_with_octant(vertex, ray))
		//		this->father->children[i]->find_ray_intersect_triangle_top_down(p, dir, pNormal, info);
		//}
		if (faceid != -1) {
			return;
		}

		if (this->isLeaf == false) {
			//If the node is not a leaf.
			//traverse itself and its childern.
			for (int i = 0; i < 8; i++) {
				if (children[i]->isEmpty == false)
					children[i]->find_ray_intersect_triangle_top_down(p, dir, pNormal, info);
			}
		}

	}

	void Octree::Node::find_ray_intersect_triangle_bottom_up(const float * p, const float * dir, const float * pNormal, Traversal_Info & info)
	{
		vec vNormal = vec(pNormal[0], pNormal[1], pNormal[2]);
		vec vertex = vec(p[0], p[1], p[2]);
		vec ray = vec(dir[0], dir[1], dir[2]);
		int faceid = -1;
		if (this->nFaces != 0) {
			for (int i = 0; i < nFaces; i++) { //Loop through all the faces in this Node
				float distance = distance_ray_face(vertex, ray, vNormal, node_faceNormal[i], node_face_pts[i]);
				//float angle = acos(vNormal ^ faceNormal[i])* 180.0 / M_PIf;
				//vector<vec> pts = node_face_pts[i];
				//if (angle > 90) {
				//	vec v1 = node_face_pts[i][0]; vec v2 = node_face_pts[i][1]; vec v3 = node_face_pts[i][2];

				//	float distance = triangle_inter(vertex, ray, v1, v2, v3);
				if (distance != 0) {
					if (info.closest_d == MIN_DIST_INIT) {
						info.closest_d = distance;
						faceid = i;
						//info.faceID = k;
					}
					else if (distance < info.closest_d) {
							info.closest_d = distance;
							faceid = i;
							//faceId = k;
					}
				}
				
			}
		}

		if (faceid != -1) {
			return;
		}

		if (this->father != NULL) {
			
			this->father->find_ray_intersect_triangle_bottom_up(p, dir, pNormal, info);
			//for (int i = 0; i < 8; i++) {
			//	bool is_node = this->father->children[i] == this;
			//	if (!is_node && this->father->children[i]->isEmpty == false) {
			//		this->father->children[i]->find_ray_intersect_triangle_bottom_up(p, dir, pNormal, info);
			//	}
			//}
			 
		}
	}
	//Check the ray whether is intersect with the Octant( 6 faces of the cube)
	bool Octree::Node::is_ray_intersect_with_octant(vec l0, vec lDir)
	{
		if (isEmpty)
			return false;

		//vec l0 = vec(p[0], p[1], p[2]);
		//vec lDir = vec(dir[0], dir[1], dir[2]);
		//The max and min X,Y,Z of...
		//... front back,right,left,up,down
		vector<vec> planesMax;
		vector<vec> planesMin;
		planesMax.push_back(vec(bBox.max[0], bBox.max[1], bBox.max[2]));
		planesMax.push_back(vec(bBox.max[0], bBox.max[1], bBox.min[2]));
		planesMax.push_back(vec(bBox.max[0], bBox.max[1], bBox.max[2]));
		planesMax.push_back(vec(bBox.min[0], bBox.max[1], bBox.max[2]));
		planesMax.push_back(vec(bBox.max[0], bBox.max[1], bBox.max[2]));
		planesMax.push_back(vec(bBox.max[0], bBox.min[1], bBox.max[2]));

		planesMin.push_back(vec(bBox.min[0], bBox.min[1], bBox.max[2]));
		planesMin.push_back(vec(bBox.min[0], bBox.min[1], bBox.min[2]));
		planesMin.push_back(vec(bBox.max[0], bBox.min[1], bBox.min[2]));
		planesMin.push_back(vec(bBox.min[0], bBox.min[1], bBox.min[2]));
		planesMin.push_back(vec(bBox.min[0], bBox.max[1], bBox.min[2]));
		planesMin.push_back(vec(bBox.min[0], bBox.min[1], bBox.min[2]));

		std::vector<vec> plane;
		//front and back planes
		plane.push_back(vec(0, 0, bBox.max[2]));
		plane.push_back(vec(0, 0, bBox.min[2]));
		//right and left planes
		plane.push_back(vec(bBox.max[0], 0, 0));
		plane.push_back(vec(bBox.min[0], 0, 0));
		//up and down planes;
		plane.push_back(vec(0, bBox.max[1], 0));
		plane.push_back(vec(0, bBox.min[1], 0));

		std::vector<vec> planeN;
		//front and back planes
		planeN.push_back(vec(0, 0, 1));
		planeN.push_back(vec(0, 0, -1));
		//right and left planes
		planeN.push_back(vec(1, 0, 0));
		planeN.push_back(vec(-1, 0, 0));
		//up and down planes;
		planeN.push_back(vec(0, 1, 0));
		planeN.push_back(vec(0, -1, 0));

		bool isIntersect = false;
		vector<vec> intersects;
		//Check the intersection between ray and 6 faces
		for (int i = 0; i < 6; i++) {
			vec intersect;
			bool temp = rayToPlane(planeN[i], plane[i], lDir, l0, intersect);
			bool inside = false;
			if (temp) {
				inside = intersect[0] <= planesMax[i][0] &&
					intersect[1] <= planesMax[i][1] &&
					intersect[2] <= planesMax[i][2];
				inside = inside &&
					intersect[0] >= planesMin[i][0] &&
					intersect[1] >= planesMin[i][1] &&
					intersect[2] >= planesMin[i][2];
			}
			isIntersect = isIntersect || temp || inside;
		}
		return isIntersect;
	}



	//
	/*float Octree::Node::distance_ray_face(vec p, vec dir, vec pNormal, vec faceNormal, vector<vec> pts)
	{
		float angle = acos(pNormal ^ faceNormal)* 180.0 / M_PIf;
		if (angle > 90) {
			vec v1 = pts[0]; vec v2 = pts[1]; vec v3 = pts[2];

			float distance = triangle_inter(p, dir, v1, v2, v3);
			return distance;
		}
		else
			return 0.0f;
	}*/


	// The the cube (leaf of Octree) that first intersect with ray.
	Octree::Traversal_Info Octree::find_cube_from_raycast(const float *p, const float *dir,
		float maxdist2 /* = 0.0f */) const
	{
		Traversal_Info ti;
		
		ti.origin_dir = dir;
		ti.origin_p = p;
		ti.level = 0;
		ti.closest_d = 1000000;

		bool originInsideBoundBox = true;
		for (int i = 0; i < 3; i++) {
			originInsideBoundBox = originInsideBoundBox && (p[i] <= root->bBox.max[i]) && (p[i] >= root->bBox.min[i]);
		}
		if (originInsideBoundBox) {
			//cout << "The emission point is inside of the bounding box" << endl;
			// to do find the shape diameters;
			Node *leaf =new Node();
			root->find_point_leaf(p, leaf);
			leaf->find_ray_intersect_originInside_bbox(p, dir, ti);
			return ti;
		}
		else {
			//cout << "The emission point is outside the bounding box" << endl;
			root->ray_intersect_cube(p, dir, ti);
			return ti;
		}

	}

	//Find the closest intersecting face for a ray.
	//Reference: An Efficient Parametric Algorithm for Octree Traversal

	Octree::Traversal_Info Octree::intersect_face_from_raycast( float * p, float * dir, 
		const float * pNormal, float maxdist2)
	{
		unsigned char a = 0;
		vec size = root->bBox.max - root->bBox.min;

		vec vertex = vec(p[0], p[1], p[2]); vec ray = vec(dir[0], dir[1], dir[2]);
		normalize(ray);
		vec ptNomral = vec(pNormal[0], pNormal[1], pNormal[2]);
		// calculate the converting bits for negative ray direction. 
		if (ray[0] < 0.0f) {
			//p[0] = size[0] - p[0];
			vertex[0] = 2 * root->bBox.center[0] - vertex[0];
			ray[0] = -ray[0];
			a = a | 4;
		}

		if (ray[1] < 0.0f) {
			vertex[1] = 2 * root->bBox.center[1] - vertex[1];
			ray[1] = -ray[1];
			a = a | 2;
		}

		if (ray[2] < 0.0f) {
			vertex[2] = 2 * root->bBox.center[2] - vertex[2];
			ray[2] = -ray[2];
			a = a | 1;
		}
		for (int i = 0; i < 3; i++) {
			if (ray[i] == 0) {
				ray[i] = 0;
			}
		}

		Traversal_Info info;
		info.origin_dir = dir;
		info.origin_p = p;
		

		info.level = 0;
		info.closest_d = MIN_DIST_INIT;
		info.ptNormal = ptNomral;
		info.ray_neg_convert_bits = a;
		bool originInsideBoundBox = true;
		for (int i = 0; i < 3; i++) {
			originInsideBoundBox = originInsideBoundBox && (p[i] <= root->bBox.max[i]) && (p[i] >= root->bBox.min[i]);
		}

		
		if (originInsideBoundBox) {
			//vector<float> t_i0(3);	
			vector<float>  t_i0(3);
			vector<float> t_i1(3);
			for (int i = 0; i < 3; i++) {
				t_i0[i] = (root->bBox.min[i] - vertex[i]) / ray[i];
				t_i1[i] = (root->bBox.max[i] - vertex[i]) / ray[i];
			}
			float t0_max = *std::max_element(t_i0.begin(), t_i0.end());
			vec outO = vertex + ((t0_max) * ray);

			//info.origin_t = t0_max;
			vector<float> new_t_i0(3);
			vector<float> new_t_i1(3);
			for (int i = 0; i < 3; i++) {
				new_t_i0[i] = (root->bBox.min[i] - outO[i]) / ray[i];
				new_t_i1[i] = (root->bBox.max[i] - outO[i]) / ray[i];
			}
			 //*std::min_element(t_i0, t_i0.end())
			float t1_min = *std::min_element(new_t_i1.begin(), new_t_i1.end());
			t0_max =  *std::min_element(new_t_i0.begin(), new_t_i0.end());
			if (t1_min > t0_max)
				root->ray_crossing_octants(new_t_i0, new_t_i1, info);
			else
				cout << "Vertex inside the bounding box. BUT not intersect with the Bounding" << endl;
			info.vertex = outO;
			info.ray = ray;
		}
		else {// the origin point is outside of the bounding box
			
			float tx0 = (root->bBox.min[0] - vertex[0]) / ray[0];
		//	float tx0 = (root->bBox.max[0] - vertex[0]) / ray[0];
			vector<float> t_i0(3);	vector<float> t_i1(3);
			for (int i = 0; i < 3; i++) {
				t_i0[i] = (root->bBox.min[i] - vertex[i]) / ray[i];
				t_i1[i] = (root->bBox.max[i] - vertex[i]) / ray[i];
			}
			float t0_max = *std::max_element(t_i0.begin(), t_i0.end()); //*std::min_element(t_i0, t_i0.end())
			float t1_min = *std::min_element(t_i1.begin(), t_i1.end());
			if (t1_min > t0_max)
				root->ray_crossing_octants(t_i0, t_i1, info);
				//cout << "inside" << endl;
			else
				cout << "not intersect with the Bounding" << endl;
			//root->ray_crossing_octants(t_i0 , t_i1, info);
		}

		if (info.closest_d == MIN_DIST_INIT)
			info.closest_d = 0;
		return info;
	}

	Octree:: Octree::Node * Octree::getRoot()
	{
		float a = root->bBox.center[0];
		return root;
	}

	Octree::Traversal_Info Octree::find_all_leaves()
	{
		Traversal_Info info;
		root->traverse_all_leaves(info);

		return info;
	}

	// A point together with a distance - default comparison is by "first",
	// i.e., distance
	typedef pair<float, const float *> pt_with_d;
	// Destroy a Octree tree node
	Octree::Node::~Node()
	{
		if (!nPts) {
			for (int i = 0; i < 8; i++) {
				delete children[i];
			}
		}
	}

	// Delete a Octree
	Octree::~Octree()
	{
		delete root;
	}


	//Build the octree with vertices as input
	void Octree::build(std::vector<vec> pts, size_t n)
	{
		// Find bbox of this node
		float xmin = pts[0][0], xmax = pts[0][0];
		float ymin = pts[0][1], ymax = pts[0][1];
		float zmin = pts[0][2], zmax = pts[0][2];
		for (size_t i = 1; i < n; i++) {
			if (pts[i][0] < xmin)  xmin = pts[i][0];
			if (pts[i][0] > xmax)  xmax = pts[i][0];
			if (pts[i][1] < ymin)  ymin = pts[i][1];
			if (pts[i][1] > ymax)  ymax = pts[i][1];
			if (pts[i][2] < zmin)  zmin = pts[i][2];
			if (pts[i][2] > zmax)  zmax = pts[i][2];
		}
		float dx = xmax - xmin;
		float dy = ymax - ymin;
		float dz = zmax - zmin;
		float edgeLen = dx;
		// decide the longest
		root = new Node(pts, n, xmax, ymax, zmax, xmin, ymin, zmin);
		root->father = NULL;
	}
	//build the octree with faces
	void Octree::build_with_faces(std::vector<TriMesh::Face> faces, std::vector<vec> pts , std::vector<vec> fNormals)
	{
		std::clock_t time = std::clock();
		cout << "constructing the octree using model faces... ";
		// Find bbox of this node
		int n = pts.size();
		float xmin = pts[0][0], xmax = pts[0][0];
		float ymin = pts[0][1], ymax = pts[0][1];
		float zmin = pts[0][2], zmax = pts[0][2];
		for (size_t i = 1; i < n; i++) {
			if (pts[i][0] < xmin)  xmin = pts[i][0];
			if (pts[i][0] > xmax)  xmax = pts[i][0];
			if (pts[i][1] < ymin)  ymin = pts[i][1];
			if (pts[i][1] > ymax)  ymax = pts[i][1];
			if (pts[i][2] < zmin)  zmin = pts[i][2];
			if (pts[i][2] > zmax)  zmax = pts[i][2];
		}
		float dx = xmax - xmin;
		float dy = ymax - ymin;
		float dz = zmax - zmin;
		float edgeLen = max(dx,dy);
		edgeLen = max(edgeLen, dz);
		vec centroid = vec(0.5f * (xmin + xmax), 0.5f * (ymin + ymax), 0.5f * (zmin + zmax));

		//-----initial face points
		int nf = faces.size();
		vector<vector<vec>> facesPts(nf);
		vector<int> total_faces_index(nf);
		#pragma omp parallel for
		for (int i = 0; i < nf; i++) {
			facesPts[i].resize(3);
			for (int j = 0; j < 3; j++) {
				facesPts[i][j]= pts[faces[i][j]];
			}

			total_faces_index[i] = i;
		}
		
		if (interNode) {
			root = new Node(facesPts, total_faces_index, centroid, edgeLen, fNormals);
		}
		else {
			root = new Node(facesPts,  centroid, edgeLen, fNormals,NULL);
		}
		
		//Node(facesPt1, facesPt2, facesPt3, centroid , edgeLen);
		std::cout << endl << "constructing Octree " << double(clock() - time) / CLOCKS_PER_SEC << " s" << std::endl;
	}

	//// Class for nodes in the K-D tree
	//class Octree::Node {
	//private:
	//	static PoolAlloc memPool;

	//public:
	//	// A place to put all the stuff required while traversing the K-D
	//	// tree, so we don't have to pass tons of variables at each fcn call
	//	struct Traversal_Info {
	//		const float *p, *dir;
	//		const float *closest;
	//		float closest_d, closest_d2;
	//		const Octree::CompatFunc *iscompat;
	//		size_t k;
	//		vector<pt_with_d> knn;
	//	};

	//	enum { MAX_PTS_PER_NODE = 1 };


	//	// The node itself

	//	int npts; // If this is 0, intermediate node.  If nonzero, leaf.

	//	union {
	//		struct {
	//			float center[3];
	//			float r;
	//			int splitaxis;
	//			Node *child1, *child2;
	//		} node;
	//		struct {
	//			const float *p[MAX_PTS_PER_NODE];
	//		} leaf;
	//	};

	//	Node(const float **pts, size_t n);
	//	~Node();

	//	void find_closest_to_pt(Traversal_Info &ti) const;
	//	void find_k_closest_to_pt(Traversal_Info &ti) const;
	//	void find_closest_to_ray(Traversal_Info &ti) const;

	//	void *operator new(size_t n) { return memPool.alloc(n); }
	//	void operator delete(void *p, size_t n) { memPool.free(p, n); }
	//};


	//// Class static variable
	//PoolAlloc Octree::Node::memPool(sizeof(Octree::Node));


	//// Create a KD tree from the points pointed to by the array pts
	//Octree::Node::Node(const float **pts, size_t n)
	//{
	//	// Leaf nodes
	//	if (n <= MAX_PTS_PER_NODE) {
	//		npts = n;
	//		memcpy(leaf.p, pts, n * sizeof(float *));
	//		return;
	//	}


	//	// Else, interior nodes
	//	npts = 0;

	//	// Find bbox
	//	float xmin = pts[0][0], xmax = pts[0][0];
	//	float ymin = pts[0][1], ymax = pts[0][1];
	//	float zmin = pts[0][2], zmax = pts[0][2];
	//	for (size_t i = 1; i < n; i++) {
	//		if (pts[i][0] < xmin)  xmin = pts[i][0];
	//		if (pts[i][0] > xmax)  xmax = pts[i][0];
	//		if (pts[i][1] < ymin)  ymin = pts[i][1];
	//		if (pts[i][1] > ymax)  ymax = pts[i][1];
	//		if (pts[i][2] < zmin)  zmin = pts[i][2];
	//		if (pts[i][2] > zmax)  zmax = pts[i][2];
	//	}

	//	// Find node center and size
	//	node.center[0] = 0.5f * (xmin + xmax);
	//	node.center[1] = 0.5f * (ymin + ymax);
	//	node.center[2] = 0.5f * (zmin + zmax);
	//	float dx = xmax - xmin;
	//	float dy = ymax - ymin;
	//	float dz = zmax - zmin;
	//	node.r = 0.5f * sqrt(sqr(dx) + sqr(dy) + sqr(dz));

	//	// Find longest axis
	//	node.splitaxis = 2;
	//	if (dx > dy) {
	//		if (dx > dz)
	//			node.splitaxis = 0;
	//	}
	//	else {
	//		if (dy > dz)
	//			node.splitaxis = 1;
	//	}

	//	// Partition
	//	const float splitval = node.center[node.splitaxis];
	//	const float **left = pts, **right = pts + n - 1;
	//	while (1) {
	//		while ((*left)[node.splitaxis] < splitval)
	//			left++;
	//		while ((*right)[node.splitaxis] > splitval)
	//			right--;
	//		if (right <= left)
	//			break;
	//		swap(*left, *right);
	//		left++; right--;
	//	}

	//	// Build subtrees
	//	node.child1 = new Node(pts, left - pts);
	//	node.child2 = new Node(left, n - (left - pts));
	//}




	//// Create a KDtree from a list of points (i.e., ptlist is a list of 3*n floats)
	//void Octree::build(const float *ptlist, size_t n)
	//{
	//	vector<const float *> pts(n);
	//	for (size_t i = 0; i < n; i++)
	//		pts[i] = ptlist + i * 3;

	//	root = new Node(&(pts[0]), n);
	//}








}; // namespace trimesh
