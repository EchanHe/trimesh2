#pragma once

#include <vector>
#include "Vec.h"
#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <algorithm>
//#include <GL/glut.h>
using namespace std;
namespace trimesh {



	struct KD_Node {
		KD_Node *left = NULL;
		KD_Node * right = NULL;
		float split_value=0;
		int axis=-1;
		vector<float> pt1X; vector<float> pt1Y; vector<float> pt1Z;
		vector<float> pt2X; vector<float> pt2Y; vector<float> pt2Z;
		vector<float> pt3X; vector<float> pt3Y;	vector<float> pt3Z;
		vector<float> fNomarlX; vector<float> fNomarlY; vector<float> fNomarlZ;
		bool isLeaf = false;
	};


	struct KD_tree {
		KD_Node * root= NULL;
		float max[3];
		float min[3];
	};

	struct KD_tree_array {
		float max[3];
		float min[3];
		vector<float> split;
		vector<int> split_axis;
		vector<int> triCount;
		vector<int> triIndex;
		vector<vector<float>> pt1X;
		vector<vector<float>> pt1Y;
		vector<vector<float>> pt1Z;
			
		vector<vector<float>> pt2X;
		vector<vector<float>> pt2Y;
		vector<vector<float>> pt2Z;
			 
		vector<vector<float>> pt3X;
		vector<vector<float>> pt3Y;
		vector<vector<float>> pt3Z;
	
		vector<vector<float>> fNomarlX;
		vector<vector<float>> fNomarlY;
		vector<vector<float>> fNomarlZ;

		vector<float> pt1X_1d;
		KD_tree_array(int n) {
			pt1X.resize(n); pt1Y.resize(n);pt1Z.resize(n);
			pt2X.resize(n);pt2Y.resize(n);pt2Z.resize(n);
			pt3X.resize(n);pt3Y.resize(n);pt3Z.resize(n);
			fNomarlX.resize(n); fNomarlY.resize(n); fNomarlZ.resize(n);
			split.resize(n);
			split_axis.resize(n);
			triCount.resize(n); triIndex.resize(n);
		}
	};

	void buildKDTree(KD_tree & tree, vector<TriMesh::Face> faces, vector<vec>pts, vector<vec> fNormals);
	void buildKDNode(KD_Node & node, float ** fPt1, float ** fPt2, float ** fPt3, float **fNormals);
	KD_Node * buildKDNode(KD_Node * node, vector<float> fPt1[], vector<float> fPt2[], vector<float> fPt3[], vector<float> fNormals[]);
	int heightKD_Tree(KD_Node* node);
	void printLevelOrder(KD_Node* root);
	void printGivenLevel(KD_Node* node, int level);

	KD_tree_array* KDTreeToArray(KD_tree tree);
	void KDTreeToArrayPerLevel(KD_Node* node, int level  );
	//void cal_sdf_on_node();
	void kd_to_array(KD_tree_array * kdArray, KD_Node* node , int index);

	float rayToFaces(KD_tree_array * kdArray, float * p, float * dir,float * pNormal);
	float cal_dist_on_node(float * p, float * dir, float * pNormal,
		float * p1X, float * p1Y, float * p1Z,
		float * p2X, float * p2Y, float * p2Z, 
		float * p3X, float * p3Y, float * p3Z,
		float * faceNX, float * faceNY, float * faceNZ, int size
		);
	void len_tri(KD_tree_array * kdArray);


	//Class of kd tree.
	class KDtree_face {
	private:
		//class Node;
		//Node *root;
		//Build the octree with faces as input
		void build_with_faces(std::vector<TriMesh::Face> faces, std::vector<vec> pts, std::vector<vec> fNormals);
	public:
		class Node;
		Node *root;
		struct Traversal_Info {
			const float *origin_p, *origin_dir, *ptNormal;
			//const float *vertex, *ray;
			vec origin_vec_p, origin_vec_dir;
			float closest_d;
			int level;
			int iChild;
			int faceID = -1;
			size_t k;
			unsigned char ray_neg_convert_bits;
			std::vector<vec> bBox_maxs;
			std::vector<vec> bBox_mins;
			vec vertex;
			vec ray;
			//float origin_t;
			int cal_count = 0;// the times for counting.
			int level_count = 0;
			//vector<pt_with_d> knn;
		};

		struct Tree_Info {
			int highest_layer;
			int lowest_layer;
			int equal_max_faces_count = 0;
			int leaves_count = 0;
			int branch_count = 0;
			int faces_count = 0;
		};
		const static bool INTER_NODE = true;

		const static int MAX_FACES_PER_NODE = 2;

		const static int MIN_DIST_INIT = 1000000;


		// Constructor from a vector of faces -------------------
		template <class T> KDtree_face(const ::std::vector<T> &f, const ::std::vector<vec> &v, const ::std::vector<vec> &fNormals)
		{
			//MAX_PTS_PER_NODE = pts_per_node;
			build_with_faces(f, v, fNormals);
		}


		// Destructor - recursively frees the tree
		~KDtree_face();

		Tree_Info find_tree_info();
		Traversal_Info intersect_face_from_raycast(float * p, float * dir,
			const float * pNormal);
	};
}; // namespace trimesh