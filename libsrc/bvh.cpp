#include "bvh.h"

using namespace std;
namespace trimesh {
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


	class BVH::Node {
	public:
		Node* left,*right;

		bool isRoot = false;
		bool isLeaf;
		bool isEmpty = false;
		bool has_faces_branch;
		int nPts;
		int nFaces;

		int totaln;

	
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
	

		Node(vector<vector<vec>> facesPts, vector<vec>pts, vector<int> facesIDs,std::vector<vec> fNormals);
		Node(vector<vector<vec>> facesPts, vector<vec>pts, std::vector<vec> fNormals);
		~Node();

		void ray_crossing_octants(Traversal_Info & info) {
			info.bBox_maxs.push_back(bBox.max);
			info.bBox_mins.push_back(bBox.min);
			if (info.closest_d != MIN_DIST_INIT)
				return;
			if (isLeaf) {
				
				//intersect
				this->cal_sdf_on_node(info);
				return;
			}

			vector<float> left_t0(3); vector<float> right_t0(3);
			vector<float> left_t1(3); vector<float> right_t1(3);
			for (int i = 0; i < 3; i++) {
				if (info.ray[i] >= 0) {
					left_t0[i] = (this->left->bBox.min[i] - info.vertex[i]) / info.ray[i];
					right_t0[i] = (this->right->bBox.min[i] - info.vertex[i]) / info.ray[i];

					left_t1[i] = (this->left->bBox.max[i] - info.vertex[i] )/ info.ray[i];
					right_t1[i] =( this->right->bBox.max[i] - info.vertex[i] )/ info.ray[i];
				}
				else {
					left_t0[i] = (this->left->bBox.max[i] - info.vertex[i]) / info.ray[i];
					right_t0[i] = (this->right->bBox.max[i] - info.vertex[i]) / info.ray[i];

					left_t1[i] = (this->left->bBox.min[i] - info.vertex[i]) / info.ray[i];
					right_t1[i] = (this->right->bBox.min[i] - info.vertex[i]) / info.ray[i];
				}
			}
			float left_t0_max = *std::max_element(left_t0.begin(), left_t0.end());
			float right_t0_max = *std::max_element(right_t0.begin(), right_t0.end());
			float left_t1_min = *std::min_element(left_t1.begin(), left_t1.end());
			float right_t1_min = *std::min_element(right_t1.begin(), right_t1.end());

			if (left_t1_min > left_t0_max && right_t1_min > right_t0_max) {
				if (left_t0_max < right_t0_max) {
					this->left->ray_crossing_octants(info);
					this->right->ray_crossing_octants(info);
				}
				else {
					this->right->ray_crossing_octants(info);
					this->left->ray_crossing_octants(info);
				}
			}
			else if(left_t1_min > left_t0_max)
				this->left->ray_crossing_octants(info);
			else if(right_t1_min > right_t0_max)
				this->right->ray_crossing_octants(info);
		};
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
		};

		void traverse_all_info(Tree_Info & info) {
			if (isLeaf) {

				if (nFaces > info.max_faces_count)
					info.max_faces_count = nFaces;

				//cout << "reach the leaf";
				if (nFaces <= MAX_FACES_PER_NODE)
					info.equal_max_faces_count++;
				info.leaves_count++;
				if(nFaces>=0)
					info.faces_count += nFaces;
				return;
			}
			left->traverse_all_info(info);
			right->traverse_all_info(info);
		}

	};// END the declaration of NODE

	struct s_plane {
		float value;
		int index;
		bool isBegin;
	};
	enum Axis{X=0 ,  Y=1 , Z=2 ,NONE=-1};
	bool s_plane_sort(s_plane a, s_plane b) { return a.value < b.value; }
	BVH::Node::Node(vector<vector<vec>> facesPts, vector<vec>pts, vector<int> facesIDs,  std::vector<vec> fNormals) {
		int nf = facesPts.size();
		//calculate the bounding box

		//SURFACE AREA HEURISTIC
		s_plane sp; vector<s_plane> listX, listY, listZ;
		for (int i = 0; i < nf; i++) {
			vector<vec>  facePts = facesPts[i];
			sp.index = i;
			//x
			sp.value = min(facePts[0][0], min(facePts[1][0], facePts[2][0])); sp.isBegin = true; listX.push_back(sp);
			sp.value = max(facePts[0][0], max(facePts[1][0], facePts[2][0])); sp.isBegin = false; listX.push_back(sp);
			//y
			sp.value = min(facePts[0][1], min(facePts[1][1], facePts[2][1])); sp.isBegin = true; listY.push_back(sp);
			sp.value = max(facePts[0][1], max(facePts[1][1], facePts[2][1])); sp.isBegin = false; listY.push_back(sp);
			//Z
			sp.value = min(facePts[0][2], min(facePts[1][2], facePts[2][2])); sp.isBegin = true; listZ.push_back(sp);
			sp.value = max(facePts[0][2], max(facePts[1][2], facePts[2][2])); sp.isBegin = false; listZ.push_back(sp);
		}

		//find the best plane
		std::sort(listZ.begin(), listZ.end(), s_plane_sort);
		std::sort(listY.begin(), listY.end(), s_plane_sort);
		std::sort(listX.begin(), listX.end(), s_plane_sort);

		float xmax = listX[listX.size() - 1].value; float xmin = listX[0].value;
		float ymax = listY[listY.size() - 1].value; float ymin = listY[0].value;
		float zmax = listZ[listZ.size() - 1].value; float zmin = listZ[0].value;
		float dx = xmax - xmin;
		float dy = ymax - ymin;
		float dz = zmax - zmin;
		bBox.max[0] = xmax; bBox.min[0] = xmin;
		bBox.max[1] = ymax; bBox.min[1] = ymin;
		bBox.max[2] = zmax; bBox.min[2] = zmin;

		//--leaf node
		if (nf <= MAX_FACES_PER_NODE) {
			isLeaf = true;
			nFaces = nf;
			for (int i = 0; i < nf; i++) {
				node_face_pts.push_back(facesPts[i]);
				node_faceNormal.push_back(fNormals[i]);
			}
			return;
		}
		isLeaf = false;



		double traversal_cost = 1.0;
		double intersection_cost = 1.0;

		double best_cost = 1e8, cost, nosplit_cost = nf* intersection_cost;
		// best cost, current cost, and the cost of not splitting the node
		float split_position = 100000;

		double _width = dx;
		double _height = dy;
		double _depth = dz;
		double surface_area = 2.0 * (_width * _height + _height * _depth + _depth * _width), // surface area of node
			templ, tempr, surface_area_l, surface_area_r; // temp values for dimension split and surface area of child nodes
		int nl, nr; // number of primitives in child nodes
		char flag; // we are including primitives lying on a split plane in both children.. this should prevent a decrease in nr on the first pass

		nl = 0; nr = nf; flag = 0;
		Axis axis = NONE;
		if (dx >= dy && dx >= dz)axis = X; if (dy >= dx && dy >= dz)axis = Y; if (dz >= dx && dz >= dy)axis = Z;


		for (int i = 0; i < listZ.size(); i++) {
			//if (listZ[i].isBegin == true) { if (flag == 1) { nr--; flag = 0; } nl++; }
			//if (listZ[i].isBegin == false) { if (flag == 1) { nr--; } flag = 1; }
			if (axis != Z)continue;
		
			if (listZ[i].isBegin == false) { nl++; nr--; }


			if (listZ[i].value <= zmin || listZ[i].value >= zmax) continue; // ensure the split position is a proper one

			templ = listZ[i].value - zmin;
			tempr = zmax - listZ[i].value;
			surface_area_l = 2.0 * (_width * _height + _height * templ + templ * _width);
			surface_area_r = 2.0 * (_width * _height + _height * tempr + tempr * _width);

			cost = traversal_cost +  nl/nf  +  nr/nf;
			if (nl == nf / 2) {
			//	best_cost = cost;
				axis = Z;
				split_position = listZ[i].value;
			}
		}
		nl = 0; nr = nf; flag = 0;
		for (int i = 0; i < listY.size(); i++) {
			if (axis != Y)continue;
			
			if (listY[i].isBegin == false) { nl++; nr--; }

			if (listY[i].value <= ymin || listY[i].value >= ymax) continue; // ensure the split position is a proper one

			templ = listY[i].value - ymin;
			tempr = ymax - listY[i].value;
			surface_area_l = 2.0 * (_width * templ + templ * _depth + _depth * _width);
			surface_area_r = 2.0 * (_width * tempr + tempr * _depth + _depth * _width);

			cost = traversal_cost +  nl/nf +  nr/nf ;
			if (nl == nf / 2) {
			//	best_cost = cost;
				axis = Y;
				split_position = listY[i].value;
			}
		}
		nl = 0; nr = nf; flag = 0;
		for (int i = 0; i < listX.size(); i++) {
			if (axis != X)continue;


			if (listX[i].isBegin == false) { nl++; nr--; }
			if (listX[i].value <= xmin || listX[i].value >= xmax) continue; // ensure the split position is a proper one

			templ = listX[i].value - xmin;
			tempr = xmax - listX[i].value;
			surface_area_l = 2.0 * (templ * _height + _height * _depth + _depth * templ);
			surface_area_r = 2.0 * (tempr * _height + _height * _depth + _depth * tempr);

			
			if (nl == nf / 2) {
				//best_cost = cost;
				axis = X;
				split_position = listX[i].value;
			}
		}

		if (split_position == 100000) {
			isLeaf = true;
			nFaces = nf;
			for (int i = 0; i < nf; i++) {
				node_face_pts.push_back(facesPts[i]);
				node_faceNormal.push_back(fNormals[i]);
			}
			return;
		}

		vector<vector<vec>> faces_left;
		vector<vector<vec>> faces_right;
		vector<vec> fNormals_left;
		vector<vec> fNormals_right;
		//split the node
		bool left_right_flag = true;
		for (int i = 0; i <  listX.size(); i++) {
			//int list_id = 2 * i;
			bool not_crossing = true;
			if (axis == X) {
				if (listX[i].value <= split_position  && listX[i].isBegin == false) {
					faces_left.push_back(facesPts[listX[i].index]);
					fNormals_left.push_back(fNormals[listX[i].index]);
				}
				else if (listX[i].value> split_position&& listX[i].isBegin == false) {
					faces_right.push_back(facesPts[listX[i].index]);
					fNormals_right.push_back(fNormals[listX[i].index]);
				}
			}
			else if (axis == Y) {
				if (listY[i].value <= split_position  && listY[i].isBegin == false) {
					faces_left.push_back(facesPts[listY[i].index]);
					fNormals_left.push_back(fNormals[listY[i].index]);
				}
				else if (listY[i].value> split_position&& listY[i].isBegin == false) {
					faces_right.push_back(facesPts[listY[i].index]);
					fNormals_right.push_back(fNormals[listY[i].index]);
				}
			}
			else if (axis == Z) {
				if (listZ[i].isBegin == false) {
					if (listZ[i].value <= split_position) {
						faces_left.push_back(facesPts[listZ[i].index]);
						fNormals_left.push_back(fNormals[listZ[i].index]);
					}
					else if (listZ[i].value > split_position) {
						faces_right.push_back(facesPts[listZ[i].index]);
						fNormals_right.push_back(fNormals[listZ[i].index]);
					}
				}
			}
		
			//if (left_right_flag && !not_crossing) {
			//	faces_right.push_back(facesPts[i]);
			//	fNormals_right.push_back(fNormals[i]);
			//}
			//else if(!left_right_flag && !not_crossing){
			//	faces_left.push_back(facesPts[i]);
			//	fNormals_left.push_back(fNormals[i]);
			//}
			//left_right_flag = !left_right_flag;

		}

		if (faces_left.size() == nf || faces_right.size() == nf) {
			isLeaf = true;
			nFaces = nf;
			for (int i = 0; i < nf; i++) {
				node_face_pts.push_back(facesPts[i]);
				node_faceNormal.push_back(fNormals[i]);
			}
			return;
		}

		int breakpoint = faces_left.size() ;
		left = new Node(faces_left, pts, facesIDs, fNormals_left);
		right = new Node(faces_right, pts, facesIDs, fNormals_right);
	}


	BVH::Node::~Node()
	{
		if (!nPts) {
			delete left;
			delete right;
		}
	}


	// Delete a BVH
	BVH::~BVH()
	{
		delete root;
	}
	void BVH::build_with_faces(std::vector<TriMesh::Face> faces, std::vector<vec> pts, std::vector<vec> fNormals) {
		std::clock_t time = std::clock();
		cout << "constructing the BVH using model faces... ";
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
		float edgeLen = max(dx, dy);
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
				facesPts[i][j] = pts[faces[i][j]];
			}

			total_faces_index[i] = i;
		}

		root =new Node(facesPts, pts, total_faces_index, fNormals);
		std::cout << endl << "constructing Octree " << double(clock() - time) / CLOCKS_PER_SEC << " s" << std::endl;
	}

	BVH::Tree_Info BVH::find_tree_info() {
		Tree_Info info;
		root->traverse_all_info(info);

		cout << "faces count: " << info.faces_count<<endl;
		cout << "The leaves and the equal MAX_FACES leaves: " << info.leaves_count << " " << info.equal_max_faces_count << endl;

		return info;
	}


	BVH::Traversal_Info BVH::intersect_face_from_raycast(float * p, float * dir,
		const float * pNormal)
	{
		unsigned char a = 0;
		vec size = root->bBox.max - root->bBox.min;

		vec vertex = vec(p[0], p[1], p[2]); vec ray = vec(dir[0], dir[1], dir[2]);
		normalize(ray);
		vec ptNomral = vec(pNormal[0], pNormal[1], pNormal[2]);
		// calculate the converting bits for negative ray direction. 
		//if (ray[0] < 0.0f) {
		//	//p[0] = size[0] - p[0];
		//	vertex[0] = 2 * root->bBox.center[0] - vertex[0];
		//	ray[0] = -ray[0];
		//	a = a | 4;
		//}

		//if (ray[1] < 0.0f) {
		//	vertex[1] = 2 * root->bBox.center[1] - vertex[1];
		//	ray[1] = -ray[1];
		//	a = a | 2;
		//}

		//if (ray[2] < 0.0f) {
		//	vertex[2] = 2 * root->bBox.center[2] - vertex[2];
		//	ray[2] = -ray[2];
		//	a = a | 1;
		//}
		
		//adjust -0 to 0.
		for (int i = 0; i < 3; i++) {
			if (ray[i] == 0) {
				ray[i] = 0;
			}
		}

		Traversal_Info info;
		info.origin_dir = dir;
		info.origin_p = p;

		info.ptNormal = ptNomral;
		info.closest_d = MIN_DIST_INIT;
		
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
				if (ray[i] >= 0) {
					t_i0[i] = (root->bBox.min[i] - vertex[i]) / ray[i];
					t_i1[i] = (root->bBox.max[i] - vertex[i]) / ray[i];
				}
				else {
					t_i0[i] = (root->bBox.max[i] - vertex[i]) / ray[i];
					t_i1[i] = (root->bBox.min[i] - vertex[i]) / ray[i];
				}

			}
			float t0_max = *std::max_element(t_i0.begin(), t_i0.end());
			vec outO = vertex + ((t0_max)* ray);

			//info.origin_t = t0_max;
			vector<float> new_t_i0(3);
			vector<float> new_t_i1(3);
			for (int i = 0; i < 3; i++) {
				if (ray[i] >= 0) {
					new_t_i0[i] = (root->bBox.min[i] - outO[i]) / ray[i];
					new_t_i1[i] = (root->bBox.max[i] - outO[i]) / ray[i];
				}
				else {
					new_t_i0[i] = (root->bBox.max[i] - outO[i]) / ray[i];
					new_t_i1[i] = (root->bBox.min[i] - outO[i]) / ray[i];
				}
			}
			//*std::min_element(t_i0, t_i0.end())
			float t1_min = *std::min_element(new_t_i1.begin(), new_t_i1.end());
			t0_max = *std::min_element(new_t_i0.begin(), new_t_i0.end());
			info.vertex = outO;
			info.ray = ray;
			if (t1_min > t0_max)
				//cout << "in side inside";
				root->ray_crossing_octants(info);
			else
				cout << "Vertex inside the bounding box. BUT not intersect with the Bounding" << endl;

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
				//root->ray_crossing_octants(t_i0, t_i1, info);
				cout << "inside" << endl;
			else
				cout << "not intersect with the Bounding" << endl;
			//root->ray_crossing_octants(t_i0 , t_i1, info);
		}

		if (info.closest_d == MIN_DIST_INIT)
			info.closest_d = 0;
		return info;
	}


};