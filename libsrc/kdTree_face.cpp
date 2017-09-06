#include "kdTree_face.h"

using namespace std;
namespace trimesh {
	inline int kd_left(int index) {
		return ((index + 1) * 2)-1;

	};
	inline int kd_right(int index) {
		return ((index + 1) * 2) ;
	};
	 void  cross(float* v1, float * v2, float(&result)[3]) {
		//float result[3];
		result[0] = v1[1] * v2[2] - v1[2] * v2[1];
		result[1] = v1[2] * v2[0] - v1[0] * v2[2];
		result[2] = v1[0] * v2[1] - v1[1] * v2[0];
		//return result;
	}

	 float dot(float * v1, float *v2) {
		float sum = v1[0] * v2[0];
		for (size_t i = 1; i < 3; i++)
			sum += v1[i] * v2[i];
		return sum;
	}

	 void plus(float *v1, float* v2, float(&result)[3], bool plus = true) {
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

	 void scale(float s, float* v1, float(&result)[3]) {
		//float * result = new float[3];
		for (int i = 0; i < 3; i++) {
			result[i] = s*v1[i];
		}
		//return result;
	}
	 float len(float* v1) {
		float result = dot(v1, v1);
		return sqrt(result);
	}
	float triangle_inter(float * v1, float * v2, float* v3,
		float *dir, float * vertex) {
		float e[2][3];
		//e[0] = new float[3];
		//e[1] = new float[3];
		for (int i = 0; i < 3; i++) {
			e[0][i] = v2[i] - v1[i];
			e[1][i] = v3[i] - v1[i];
		}

		//float * pvec;
		//pvec = cross(e[1], dir);
		float  pvec[3];

		cross(e[1], dir, pvec);

		float det = dot(e[0], pvec);
		if (det < 1e-8 && det > -1e-8) {
			//std::cout << "parrel";
			return 0;
		}

		float inv_det = 1 / det;
		float tvec[3];
		plus(vertex, v1, tvec, false);
		float u = dot(tvec, pvec) * inv_det;
		if (u < 0 || u > 1) {
			return 0;
		}
		float qvec[3];
		cross(tvec, e[0], qvec);
		float v = dot(dir, qvec) * inv_det;
		if (v < 0 || u + v > 1) {
			return 0;
		}
		float t = inv_det * dot(e[1], qvec);

		//calculate the intersection points.
		//float * inter_point =plus( vertex , scale(t,dir));
		float result[3];
		scale(t, dir, result);
		return len(result);
	}

	//-end of inline 
	void buildKDTree(KD_tree & tree , vector<TriMesh::Face> faces, vector<vec>pts,vector<vec> fNormals) {
		int nv = pts.size();
		float xmin = pts[0][0], xmax = pts[0][0];
		float ymin = pts[0][1], ymax = pts[0][1];
		float zmin = pts[0][2], zmax = pts[0][2];
		for (size_t i = 1; i < nv; i++) {
			if (pts[i][0] < xmin)  xmin = pts[i][0];
			if (pts[i][0] > xmax)  xmax = pts[i][0];
			if (pts[i][1] < ymin)  ymin = pts[i][1];
			if (pts[i][1] > ymax)  ymax = pts[i][1];
			if (pts[i][2] < zmin)  zmin = pts[i][2];
			if (pts[i][2] > zmax)  zmax = pts[i][2];
		}
		tree.max[0] = xmax; tree.min[0] = xmin;
		tree.max[1] = ymax; tree.min[1] = ymin;
		tree.max[2] = zmax; tree.min[2] = zmin;

		int nf = faces.size();

		vector<float> fPt1[3]; vector<float> fPt2[3]; vector<float> fPt3[3];
		vector<float> faceNormals[3];
		for (int i = 0; i < 3; i++) {
			fPt1[i].resize(nf); fPt2[i].resize(nf); fPt3[i].resize(nf); faceNormals[i].resize(nf);
		}

		for (int i = 0; i < nf; i++) {
			
			for (int j = 0; j < 3; j++) {
				fPt1[j][i] = pts[faces[i][0]][j];
				fPt2[j][i] = pts[faces[i][1]][j];
				fPt3[j][i] = pts[faces[i][2]][j];
				faceNormals[j][i] = fNormals[i][j];
			}
		}

		tree.root = buildKDNode(tree.root, fPt1,fPt2,fPt3,faceNormals);

	};

	struct s_plane {
		float value;
		int index;
		bool isBegin;
	};
	enum Axis { X = 0, Y = 1, Z = 2, NONE = -1 };
	bool split_p_sort(s_plane a, s_plane b) { return a.value < b.value; }
	//fPti [nf][3]
	KD_Node * buildKDNode(KD_Node * node, vector<float> fPt1[], vector<float> fPt2[], vector<float> fPt3[], vector<float> fNormals[]) {
		int nf = fPt1[0].size();
		node = new KD_Node;
		s_plane sp; vector<s_plane> listX, listY, listZ;
		for (int i = 0; i < nf; i++) {

			sp.index = i;
			//x
			sp.value = min(fPt1[0][i], min(fPt2[0][i], fPt3[0][i])); sp.isBegin = true; listX.push_back(sp);
			sp.value = max(fPt1[0][i], max(fPt2[0][i], fPt3[0][i])); sp.isBegin = false; listX.push_back(sp);
			//y
			sp.value = min(fPt1[1][i], min(fPt2[1][i], fPt3[1][i])); sp.isBegin = true; listY.push_back(sp);
			sp.value = max(fPt1[1][i], max(fPt2[1][i], fPt3[1][i])); sp.isBegin = false; listY.push_back(sp);
			//Z
			sp.value = min(fPt1[2][i], min(fPt2[2][i], fPt3[2][i])); sp.isBegin = true; listZ.push_back(sp);
			sp.value = max(fPt1[2][i], max(fPt2[2][i], fPt3[2][i])); sp.isBegin = false; listZ.push_back(sp);
		}

		//find the best plane
		std::sort(listZ.begin(), listZ.end(), split_p_sort);
		std::sort(listY.begin(), listY.end(), split_p_sort);
		std::sort(listX.begin(), listX.end(), split_p_sort);

		float xmax = listX[listX.size() - 1].value; float xmin = listX[0].value;
		float ymax = listY[listY.size() - 1].value; float ymin = listY[0].value;
		float zmax = listZ[listZ.size() - 1].value; float zmin = listZ[0].value;
		float dx = xmax - xmin;
		float dy = ymax - ymin;
		float dz = zmax - zmin;

		float  traversal_cost = 1.0;
		float  intersection_cost = 1.0;
		float  best_cost = 1e8, cost, nosplit_cost = nf* intersection_cost;
		// best cost, current cost, and the cost of not splitting the node
		float split_position;

		float _width = dx;
		float _height = dy;
		float _depth = dz;
		float surface_area = 2.0 * (_width * _height + _height * _depth + _depth * _width), // surface area of node
			templ, tempr, surface_area_l, surface_area_r; // temp values for dimension split and surface area of child nodes
		int nl, nr; // number of primitives in child nodes
		int final_nl, final_nr;
		char flag; // we are including primitives lying on a split plane in both children.. this should prevent a decrease in nr on the first pass
		nl = 0; nr = nf; flag = 0; final_nl = 0; nr = nf;
		Axis axis = NONE;
		for (int i = 0; i < listZ.size(); i++) {
			if (listZ[i].isBegin == true) { if (flag == 1) { nr--; flag = 0; } nl++; }
			if (listZ[i].isBegin == false) { if (flag == 1) { nr--; } flag = 1; }

			//if (listZ[i].isBegin == true) { if (flag == 1) { nr--; flag = 0; } nl++; }
			//if (listZ[i].isBegin == false) { nl++; nr--; }


			if (listZ[i].value <= zmin || listZ[i].value >= zmax) continue; // ensure the split position is a proper one

			templ = listZ[i].value - zmin;
			tempr = zmax - listZ[i].value;
			surface_area_l = 2.0 * (_width * _height + _height * templ + templ * _width);
			surface_area_r = 2.0 * (_width * _height + _height * tempr + tempr * _width);

			cost = traversal_cost + surface_area_l / surface_area * nl * intersection_cost + surface_area_r / surface_area * nr * intersection_cost;
			if (cost < best_cost) {
				best_cost = cost;
				axis = Z;
				split_position = listZ[i].value;
				final_nl = nl;
				final_nr = nr;
			}
		}
		nl = 0; nr = nf; flag = 0;
		for (int i = 0; i < listY.size(); i++) {
			if (listY[i].isBegin == true) { if (flag == 1) { nr--; flag = 0; } nl++; }
			if (listY[i].isBegin == false) { if (flag == 1) { nr--; } flag = 1; }
			//if (listY[i].isBegin == false) { nl++; nr--; }

			if (listY[i].value <= ymin || listY[i].value >= ymax) continue; // ensure the split position is a proper one

			templ = listY[i].value - ymin;
			tempr = ymax - listY[i].value;
			surface_area_l = 2.0 * (_width * templ + templ * _depth + _depth * _width);
			surface_area_r = 2.0 * (_width * tempr + tempr * _depth + _depth * _width);

			cost = traversal_cost + surface_area_l / surface_area * nl * intersection_cost + surface_area_r / surface_area * nr * intersection_cost;
			if (cost < best_cost) {
				best_cost = cost;
				axis = Y;
				split_position = listY[i].value;
				final_nl = nl;
				final_nr = nr;
			}
		}
		nl = 0; nr = nf; flag = 0;
		for (int i = 0; i < listX.size(); i++) {
			if (listX[i].isBegin == true) { if (flag == 1) { nr--; flag = 0; } nl++; }
			if (listX[i].isBegin == false) { if (flag == 1) { nr--; } flag = 1; }

			if (listX[i].value <= xmin || listX[i].value >= xmax) continue; // ensure the split position is a proper one

			templ = listX[i].value - xmin;
			tempr = xmax - listX[i].value;

			surface_area_l = 2.0 * (templ * _height + _height * _depth + _depth * templ);
			surface_area_r = 2.0 * (tempr * _height + _height * _depth + _depth * tempr);

			cost = traversal_cost + surface_area_l / surface_area * nl * intersection_cost + surface_area_r / surface_area * nr * intersection_cost;
			if (cost < best_cost) {
				best_cost = cost;
				axis = X;
				split_position = listX[i].value;
				final_nl = nl;
				final_nr = nr;
			}
		}

		listX.clear();
		listY.clear();
		listZ.clear();
		for (int i = 0; i < nf; i++) {
			sp.index = i;
			//x
			sp.value = min(fPt1[0][i], min(fPt2[0][i], fPt3[0][i])); sp.isBegin = true; listX.push_back(sp);
			sp.value = max(fPt1[0][i], max(fPt2[0][i], fPt3[0][i])); sp.isBegin = false; listX.push_back(sp);
			//y
			sp.value = min(fPt1[1][i], min(fPt2[1][i], fPt3[1][i])); sp.isBegin = true; listY.push_back(sp);
			sp.value = max(fPt1[1][i], max(fPt2[1][i], fPt3[1][i])); sp.isBegin = false; listY.push_back(sp);
			//Z
			sp.value = min(fPt1[2][i], min(fPt2[2][i], fPt3[2][i])); sp.isBegin = true; listZ.push_back(sp);
			sp.value = max(fPt1[2][i], max(fPt2[2][i], fPt3[2][i])); sp.isBegin = false; listZ.push_back(sp);
		}

		vector<float> left_fPt1[3]; vector<float> left_fPt2[3]; vector<float> left_fPt3[3];
		vector<float> left_faceNormals[3];
		vector<float> right_fPt1[3]; vector<float> right_fPt2[3]; vector<float> right_fPt3[3];
		vector<float> right_faceNormals[3];
		//for (int i = 0; i < 3; i++) {
		//	left_fPt1[i].resize(nf); left_fPt2[i].resize(nf); left_fPt3[i].resize(nf); left_faceNormals[i].resize(nf);
		//	right_fPt1[i].resize(nf);
		//}
	
		if (best_cost < nosplit_cost && nf != final_nl && nf != final_nr) { // termination criteria
			node->axis = axis;
			node->split_value = split_position;
			for (int i = 0; i < nf; i++) {
				int list_id = 2 * i;
				bool not_crossing = true;
				bool left_condition = 
					(axis == X && (listX[list_id].value <= split_position && listX[list_id+1].value <= split_position) )||
					(axis == Y && (listY[list_id].value <= split_position && listY[list_id + 1].value <= split_position) )||
					(axis == Z && (listZ[list_id].value <= split_position && listZ[list_id + 1].value <= split_position));
				bool right_condition = 
					(axis == X && (listX[list_id].value >= split_position && listX[list_id + 1].value >= split_position)) ||
					(axis == Y && (listY[list_id].value >= split_position && listY[list_id + 1].value >= split_position))||
					(axis == Z && (listZ[list_id].value >= split_position && listZ[list_id + 1].value >= split_position));
				if (left_condition) {
					//--less
					left_fPt1[0].push_back(fPt1[0][i]); left_fPt1[1].push_back(fPt1[1][i]); left_fPt1[2].push_back(fPt1[2][i]);
					left_fPt2[0].push_back(fPt2[0][i]); left_fPt2[1].push_back(fPt2[1][i]); left_fPt2[2].push_back(fPt2[2][i]);
					left_fPt3[0].push_back(fPt3[0][i]); left_fPt3[1].push_back(fPt3[1][i]); left_fPt3[2].push_back(fPt3[2][i]);
					left_faceNormals[0].push_back(fNormals[0][i]);	left_faceNormals[1].push_back(fNormals[1][i]);	left_faceNormals[2].push_back(fNormals[2][i]);
				}
				else if (right_condition) {
					right_fPt1[0].push_back(fPt1[0][i]); right_fPt1[1].push_back(fPt1[1][i]); right_fPt1[2].push_back(fPt1[2][i]);
					right_fPt2[0].push_back(fPt2[0][i]); right_fPt2[1].push_back(fPt2[1][i]); right_fPt2[2].push_back(fPt2[2][i]);
					right_fPt3[0].push_back(fPt3[0][i]); right_fPt3[1].push_back(fPt3[1][i]); right_fPt3[2].push_back(fPt3[2][i]);
					right_faceNormals[0].push_back(fNormals[0][i]);	right_faceNormals[1].push_back(fNormals[1][i]);	right_faceNormals[2].push_back(fNormals[2][i]);
				}
			}
			node->left=buildKDNode(node->left, left_fPt1, left_fPt2, left_fPt3, left_faceNormals);
			node->right = buildKDNode(node->right, right_fPt1, right_fPt2, right_fPt3, right_faceNormals);
		}
		else {
			node->isLeaf = true;
			for (int i = 0; i < nf; i++) {
				node->pt1X.push_back(fPt1[0][i]); node->pt1Y.push_back(fPt1[1][i]); node->pt1Z.push_back(fPt1[2][i]);
				node->pt2X.push_back(fPt2[0][i]); node->pt2Y.push_back(fPt2[1][i]); node->pt2Z.push_back(fPt2[2][i]);
				node->pt3X.push_back(fPt3[0][i]); node->pt3Y.push_back(fPt3[1][i]); node->pt3Z.push_back(fPt3[2][i]);
				node->fNomarlX.push_back(fNormals[0][i]); node->fNomarlY.push_back(fNormals[1][i]); node->fNomarlZ.push_back(fNormals[2][i]);
			}
		
		}
		return node;


	};

	int heightKD_Tree(KD_Node* node)
	{
		if (node == NULL)
			return 0;
		else
		{
			/* compute the height of each subtree */
			int lheight = heightKD_Tree(node->left);
			int rheight = heightKD_Tree(node->right);

			/* use the larger one */
			if (lheight > rheight)
				return(lheight + 1);
			else return(rheight + 1);
		}
	}

	void printLevelOrder(KD_Node* root)
	{
		int h = heightKD_Tree(root);
		int i;
		
		for (i = 1; i <= h; i++)
			printGivenLevel(root, i);
	}
	//void find_all_bbox(KD_tree tree) {
	//	
	//	
	//	traverse_all_bbox(tree.root);
	//}
	//void traverse_all_bbox(KD_Node* node) {

	//}
	/* Print nodes at a given level */
	void printGivenLevel(KD_Node* node, int level)
	{
		if (node == NULL)
			return;
		if (level == 1)
			cout<< node->split_value <<" ";
		else if (level > 1)
		{
			printGivenLevel(node->left, level - 1);
			printGivenLevel(node->right, level - 1);
		}
	}
	
	KD_tree_array* KDTreeToArray(KD_tree tree) {
		
		int h = heightKD_Tree(tree.root);
		int count = 0;
		for (int i = 1; i <= h; i++) {
			count += pow(2, i - 1);
		}
		KD_tree_array* kdArray = new KD_tree_array(count);
		for (int i = 0; i < 3; i++) {
			kdArray->max[i] = tree.max[i];
			kdArray->min[i] = tree.min[i];
		}
		int i;
		
		kd_to_array(kdArray, tree.root, 1);
		//for (i = 1; i <= h; i++)
		return kdArray;
	}
	void kd_to_array(KD_tree_array * kdArray, KD_Node* node, int index) {
		if (node == NULL)
			return;
		int id = index - 1;
		if (node->isLeaf == false) {
			kdArray->split[id] = node->split_value;
			kdArray->split_axis[id] = node->axis;
			kdArray->triCount[id] = -1;
		}
		if (node->isLeaf == true) {
			{
				kdArray->triIndex[id] = kdArray->pt1X_1d.size();

				for (int i = 0; i < node->pt1X.size(); i++) {
					kdArray->pt1X_1d.push_back(node->pt1X[i]);
					kdArray->pt1Y_1d.push_back(node->pt1Y[i]);
					kdArray->pt1Z_1d.push_back(node->pt1Z[i]);

					kdArray->pt2X_1d.push_back(node->pt2X[i]);
					kdArray->pt2Y_1d.push_back(node->pt2Y[i]);
					kdArray->pt2Z_1d.push_back(node->pt2Z[i]);

					kdArray->pt3X_1d.push_back(node->pt3X[i]);
					kdArray->pt3Y_1d.push_back(node->pt3Y[i]);
					kdArray->pt3Z_1d.push_back(node->pt3Z[i]);

					kdArray->fNomarlX_1d.push_back(node->fNomarlX[i]);
					kdArray->fNomarlY_1d.push_back(node->fNomarlY[i]);
					kdArray->fNomarlZ_1d.push_back(node->fNomarlZ[i]);

					kdArray->pt1X[id].push_back(node->pt1X[i]); kdArray->pt1Y[id].push_back(node->pt1Y[i]); kdArray->pt1Z[id].push_back(node->pt1Z[i]);
					kdArray->pt2X[id].push_back(node->pt2X[i]); kdArray->pt2Y[id].push_back(node->pt2Y[i]); kdArray->pt2Z[id].push_back(node->pt2Z[i]);
					kdArray->pt3X[id].push_back(node->pt3X[i]); kdArray->pt3Y[id].push_back(node->pt3Y[i]); kdArray->pt3Z[id].push_back(node->pt3Z[i]);
					kdArray->fNomarlX[id].push_back(node->fNomarlX[i]); kdArray->fNomarlY[id].push_back(node->fNomarlY[i]); kdArray->fNomarlZ[id].push_back(node->fNomarlZ[i]);

				}
				kdArray->triCount[id] = node->pt1X.size();
				
			}
		}
		kd_to_array(kdArray, node->left, index * 2);
		kd_to_array(kdArray, node->right, index * 2+1);
		
	}
	float rayToFaces(KD_tree_array * kdArray, float * p, float * dir, float * pNormal) {
		vector<float> t_i0(3), t_i1(3);
		for (int i = 0; i < 3; i++) {
			if (p[i] == 0 || dir[i]==0) {
				p[i] = 0; dir[i] = 0;
			}
			if (dir[i] >= 0) {
				t_i0[i] = (kdArray->min[i] - p[i]) / dir[i];
				t_i1[i] = (kdArray->max[i] - p[i]) / dir[i];
			}
			else {
				t_i0[i] = (kdArray->max[i] - p[i]) / dir[i];
				t_i1[i] = (kdArray->min[i] - p[i]) / dir[i];
			}
		}
		float t_min = *std::max_element(t_i0.begin(), t_i0.end());
		float t_max = *std::min_element(t_i1.begin(), t_i1.end());
		
		int totalN = kdArray->pt1X.size();
		int index = 0; float distance = 0;

		if (t_max > t_min) {
			//cout << "inter sect";
			int id_stack[64];
			int * stackPtr = id_stack; *stackPtr++ = -1;
			float t_min_stack[64];float t_max_stack[64];
			float * stackPtr_tmin = t_min_stack; *stackPtr_tmin++ = -1.0;
			float * stackPtr_tmax = t_max_stack; *stackPtr_tmax++ = -1.0;
			int track_id = 0;
			while (distance <= 0 && index <= totalN) {
				if (kdArray->triCount[index] == 0) {
					//NULL node;
					index = *--stackPtr; t_min = *--stackPtr_tmin; t_max = *--stackPtr_tmax;
					if (index == -1)
						break;
				}
				else if (kdArray->triCount[index] != -1) {
					//is leaf 
					int leaf_size = kdArray->triCount[index];
					float * p_p1x = &(kdArray->pt1X[index][0]); float * p_p1y = &(kdArray->pt1Y[index][0]); float * p_p1z = &(kdArray->pt1Z[index][0]);
					float * p_p2x = &(kdArray->pt2X[index][0]); float * p_p2y = &(kdArray->pt2Y[index][0]); float * p_p2z = &(kdArray->pt2Z[index][0]);
					float * p_p3x = &(kdArray->pt3X[index][0]); float * p_p3y = &(kdArray->pt3Y[index][0]); float * p_p3z = &(kdArray->pt3Z[index][0]);
					float * p_fnx = &(kdArray->fNomarlX[index][0]); float * p_fny = &(kdArray->fNomarlY[index][0]); float * p_fnz = &(kdArray->fNomarlZ[index][0]);
					distance = cal_dist_on_node(p, dir, pNormal,
						p_p1x, p_p1y, p_p1z,
						p_p2x, p_p2y, p_p2z,
						p_p3x, p_p3y, p_p3z,
						p_fnx, p_fny, p_fnz, leaf_size);

					//change the index by POPPing the value in the stack
					//if (*(stackPtr - 1) != -1) {
					index = *--stackPtr; t_min = *--stackPtr_tmin; t_max = *--stackPtr_tmax;
					//}
					if (index == -1)
						break;
					
				}
				else {
					//NOT leaf
					int first_index, second_index;
					int axis = kdArray->split_axis[index];
					float thit = (kdArray->split[index] - p[axis]) / dir[axis];
					if (kdArray->split[index] - p[axis]>=0) {
						first_index = kd_left(index); second_index  = kd_right(index);
					}
					else {
						first_index = kd_right(index); second_index  = kd_left(index);
					}
					
					if (thit >= t_max || thit < 0)
						index = first_index;
					else if(thit <= t_min)
						index = second_index;
					else {
						//push the index of Right node into Stack
						*stackPtr++ = second_index;
						*stackPtr_tmin++ = thit; *stackPtr_tmax++ = t_max;
						index = first_index;
						t_max = thit;
						
					}

				}
				track_id++;
			}

		}
		return distance;

	}
	float cal_dist_on_node(float * p, float * dir, float * pNormal,
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
				float distance = triangle_inter(v1, v2, v3, dir, p);
				if (distance != 0 && distance < temp) {
					temp = distance;
				}
			}
		}
		if (temp == 10000000)
			return 0;
		return temp;

	}
	void len_tri(KD_tree_array * kdArray) {
		int sum = 0;
		for (int i = 0; i < kdArray->triCount.size(); i++) {
			if (kdArray->triCount[i] > 0)
				sum += kdArray->triCount[i];
			else
				sum += 1;
		}
		cout << "the length kd array is: " << kdArray->triCount.size()<<endl;
		cout << "the length of triangles per node ARRAY is: "<<sum;
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


	class KDtree_face::Node {
	public:
		Node* left, *right;

		bool isRoot = false;
		bool isLeaf;
		bool isEmpty = false;
		bool has_faces_branch;
		int nPts;
		int nFaces;
		int Axis;
		float split_value;
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


		Node(vector<vector<vec>> facesPts, vector<vec>pts, vector<int> facesIDs, std::vector<vec> fNormals);
		Node(vector<vector<vec>> facesPts, vector<vec>pts, std::vector<vec> fNormals);
		~Node();
		void traverse_all_bbox(Traversal_Info &info);
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

					left_t1[i] = (this->left->bBox.max[i] - info.vertex[i]) / info.ray[i];
					right_t1[i] = (this->right->bBox.max[i] - info.vertex[i]) / info.ray[i];
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
			else if (left_t1_min > left_t0_max)
				this->left->ray_crossing_octants(info);
			else if (right_t1_min > right_t0_max)
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
				//cout << "reach the leaf";
				if (nFaces <= MAX_FACES_PER_NODE)
					info.equal_max_faces_count++;
				info.leaves_count++;
				if (nFaces >= 0)
					info.faces_count += nFaces;
				return;
			}
			left->traverse_all_info(info);
			right->traverse_all_info(info);
		}

	};// END the declaration of NODE

	void KDtree_face::Node::traverse_all_bbox(Traversal_Info &info) {
		if (isEmpty)
			return;
		info.bBox_maxs.push_back(bBox.max);
		info.bBox_mins.push_back(bBox.min);
		if (isLeaf) {
			//cout << "reach the leaf";
			return;
		}
		left->traverse_all_bbox(info);
		right->traverse_all_bbox(info);
	}
	KDtree_face::Node::Node(vector<vector<vec>> facesPts, vector<vec>pts, vector<int> facesIDs, std::vector<vec> fNormals) {
		int nf = facesPts.size();
		//calculate the bounding box
		int nv = pts.size();
		//--leaf node
		if (nf < MAX_FACES_PER_NODE) {
			isLeaf = true;
			nFaces = nf;
			for (int i = 0; i < nf; i++) {
				node_face_pts.push_back(facesPts[i]);
				node_faceNormal.push_back(fNormals[i]);
			}
			return;
		}
		isLeaf = false;





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
		std::sort(listZ.begin(), listZ.end(), split_p_sort);
		std::sort(listY.begin(), listY.end(), split_p_sort);
		std::sort(listX.begin(), listX.end(), split_p_sort);

		float xmax = listX[listX.size() - 1].value; float xmin = listX[0].value;
		float ymax = listY[listY.size() - 1].value; float ymin = listY[0].value;
		float zmax = listZ[listZ.size() - 1].value; float zmin = listZ[0].value;
		float dx = xmax - xmin;
		float dy = ymax - ymin;
		float dz = zmax - zmin;
		bBox.max[0] = xmax; bBox.min[0] = xmin;
		bBox.max[1] = ymax; bBox.min[1] = ymin;
		bBox.max[2] = zmax; bBox.min[2] = zmin;

		double traversal_cost = 1.0;
		double intersection_cost = 1.0;

		double best_cost = 1e8, cost, nosplit_cost = nf* intersection_cost;
		// best cost, current cost, and the cost of not splitting the node
		double split_position;

		double _width = dx;
		double _height = dy;
		double _depth = dz;
		double surface_area = 2.0 * (_width * _height + _height * _depth + _depth * _width), // surface area of node
			templ, tempr, surface_area_l, surface_area_r; // temp values for dimension split and surface area of child nodes
		int nl, nr; // number of primitives in child nodes
		char flag; // we are including primitives lying on a split plane in both children.. this should prevent a decrease in nr on the first pass

		nl = 0; nr = nf; flag = 0;
		int final_nl = nl, final_nr = nr;
		trimesh::Axis axis = NONE;
		for (int i = 0; i < listZ.size(); i++) {
			if (listZ[i].isBegin == true) { if (flag == 1) { nr--; flag = 0; } nl++; }
			if (listZ[i].isBegin == false) { if (flag == 1) { nr--; } flag = 1; }

			//if (listZ[i].isBegin == true) { if (flag == 1) { nr--; flag = 0; } nl++; }
			//if (listZ[i].isBegin == false) { nl++; nr--; }


			if (listZ[i].value <= zmin || listZ[i].value >= zmax) continue; // ensure the split position is a proper one

			templ = listZ[i].value - zmin;
			tempr = zmax - listZ[i].value;
			surface_area_l = 2.0 * (_width * _height + _height * templ + templ * _width);
			surface_area_r = 2.0 * (_width * _height + _height * tempr + tempr * _width);

			cost = traversal_cost + surface_area_l / surface_area * nl * intersection_cost + surface_area_r / surface_area * nr * intersection_cost;
			if (cost < best_cost) {
				best_cost = cost;
				axis = Z;
				split_position = listZ[i].value;
				final_nl = nl;
				final_nr = nr;
			}
		}
		nl = 0; nr = nf; flag = 0;
		for (int i = 0; i < listY.size(); i++) {
			if (listY[i].isBegin == true) { if (flag == 1) { nr--; flag = 0; } nl++; }
			if (listY[i].isBegin == false) { if (flag == 1) { nr--; } flag = 1; }
			//if (listY[i].isBegin == false) { nl++; nr--; }

			if (listY[i].value <= ymin || listY[i].value >= ymax) continue; // ensure the split position is a proper one

			templ = listY[i].value - ymin;
			tempr = ymax - listY[i].value;
			surface_area_l = 2.0 * (_width * templ + templ * _depth + _depth * _width);
			surface_area_r = 2.0 * (_width * tempr + tempr * _depth + _depth * _width);

			cost = traversal_cost + surface_area_l / surface_area * nl * intersection_cost + surface_area_r / surface_area * nr * intersection_cost;
			if (cost < best_cost) {
				best_cost = cost;
				axis = Y;
				split_position = listY[i].value;
				final_nl = nl;
				final_nr = nr;
			}
		}
		nl = 0; nr = nf; flag = 0;
		for (int i = 0; i < listX.size(); i++) {
			if (listX[i].isBegin == true) { if (flag == 1) { nr--; flag = 0; } nl++; }
			if (listX[i].isBegin == false) { if (flag == 1) { nr--; } flag = 1; }

			if (listX[i].value <= xmin || listX[i].value >= xmax) continue; // ensure the split position is a proper one

			templ = listX[i].value - xmin;
			tempr = xmax - listX[i].value;
			surface_area_l = 2.0 * (templ * _height + _height * _depth + _depth * templ);
			surface_area_r = 2.0 * (tempr * _height + _height * _depth + _depth * tempr);

			cost = traversal_cost + surface_area_l / surface_area * nl * intersection_cost + surface_area_r / surface_area * nr * intersection_cost;
			if (cost < best_cost) {
				best_cost = cost;
				axis = X;
				split_position = listX[i].value;
				final_nl = nl;
				final_nr = nr;
			}
		}

		listX.clear();
		listY.clear();
		listZ.clear();
		for (int i = 0; i < nf; i++) {
			vector<vec>  facePts = facesPts[i];
			sp.index = i;
			//x
			sp.value = min(facePts[0][0], min(facePts[1][0], facePts[2][0])); sp.isBegin = true; listX.push_back(sp);
			sp.value = max(facePts[0][0], max(facePts[1][0], facePts[2][0])); sp.isBegin = false; listX.push_back(sp);
			//y
			sp.value = min(facePts[0][1], min(facePts[1][1], facePts[2][1])); sp.isBegin = true; listY.push_back(sp);
			sp.value = max(facePts[0][1], max(facePts[1][1], facePts[2][1])); sp.isBegin = false; listY.push_back(sp);
			//x
			sp.value = min(facePts[0][2], min(facePts[1][2], facePts[2][2])); sp.isBegin = true; listZ.push_back(sp);
			sp.value = max(facePts[0][2], max(facePts[1][2], facePts[2][2])); sp.isBegin = false; listZ.push_back(sp);
		}

		vector<vector<vec>> faces_left;
		vector<vector<vec>> faces_right;
		vector<vec> fNormals_left;
		vector<vec> fNormals_right;
		//split the node
		bool left_right_flag = true;
		if (best_cost < nosplit_cost && nf != final_nl && nf != final_nr) { // termination criteria
			Axis = axis;
			split_value = split_position;
			for (int i = 0; i < nf; i++) {
				int list_id = 2 * i;
				bool not_crossing = true;
				bool left_condition =
					(axis == X && (listX[list_id].value <= split_position && listX[list_id + 1].value <= split_position)) ||
					(axis == Y && (listY[list_id].value <= split_position && listY[list_id + 1].value <= split_position)) ||
					(axis == Z && (listZ[list_id].value <= split_position && listZ[list_id + 1].value <= split_position));
				bool right_condition =
					(axis == X && (listX[list_id].value >= split_position && listX[list_id + 1].value >= split_position)) ||
					(axis == Y && (listY[list_id].value >= split_position && listY[list_id + 1].value >= split_position)) ||
					(axis == Z && (listZ[list_id].value >= split_position && listZ[list_id + 1].value >= split_position));
				if (left_condition) {
					//--less
					faces_left.push_back(facesPts[i]); 
					fNormals_left.push_back(fNormals[i]);
				}
				else if (right_condition) {
					faces_right.push_back(facesPts[i]);
					fNormals_right.push_back(fNormals[i]);
				}
			}
			this->left = new Node(faces_left, pts,facesIDs, fNormals_left);
			right = new Node(faces_right, pts, facesIDs, fNormals_right);
			//node->left = buildKDNode(node->left, left_fPt1, left_fPt2, left_fPt3, left_faceNormals);
			//node->right = buildKDNode(node->right, right_fPt1, right_fPt2, right_fPt3, right_faceNormals);
		}
		else {
			isLeaf = true;
			nFaces = nf;
			for (int i = 0; i < nf; i++) {
				node_face_pts.push_back(facesPts[i]);
				node_faceNormal.push_back(fNormals[i]);
			}
			return;
		}
	

	/*	if (faces_left.size() == nf || faces_right.size() == nf) {
			isLeaf = true;
			nFaces = nf;
			for (int i = 0; i < nf; i++) {
				node_face_pts.push_back(facesPts[i]);
				node_faceNormal.push_back(fNormals[i]);
			}
			return;
		}

		int breakpoint = faces_left.size();
		left = new Node(faces_left, pts, facesIDs, fNormals_left);
		right = new Node(faces_right, pts, facesIDs, fNormals_right);*/
	}


	KDtree_face::Node::~Node()
	{
		if (!nPts) {
			delete left;
			delete right;
		}
	}


	// Delete a KDtree_face
	KDtree_face::~KDtree_face()
	{
		delete root;
	}
	void KDtree_face::build_with_faces(std::vector<TriMesh::Face> faces, std::vector<vec> pts, std::vector<vec> fNormals) {
		std::clock_t time = std::clock();
		cout << "constructing the KDtree_face using model faces... ";
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

		root = new Node(facesPts, pts, total_faces_index, fNormals);
		std::cout << endl << "constructing Octree " << double(clock() - time) / CLOCKS_PER_SEC << " s" << std::endl;
	}

	KDtree_face::Tree_Info KDtree_face::find_tree_info() {
		Tree_Info info;
		root->traverse_all_info(info);

		cout << "faces count: " << info.faces_count << endl;
		cout << "The leaves and the equal MAX_FACES leaves: " << info.leaves_count << " " << info.equal_max_faces_count << endl;

		return info;
	}

	KDtree_face::Traversal_Info KDtree_face::find_all_bbox() {
		Traversal_Info info;
		root->traverse_all_bbox(info);

		return info;
	}
	KDtree_face::Traversal_Info KDtree_face::intersect_face_from_raycast(float * p, float * dir,
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