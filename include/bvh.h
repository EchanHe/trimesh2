#pragma once

#include <vector>
#include "Vec.h"
#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <algorithm>
//#include <GL/glut.h>
namespace trimesh {
	//bounding volume hierachy
	class BVH {
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

			int max_faces_count = 0;
		};
		const static bool INTER_NODE = true;

		const static int MAX_FACES_PER_NODE = 3;

		const static int MIN_DIST_INIT = 1000000;


		// Constructor from a vector of faces -------------------
		template <class T> BVH(const ::std::vector<T> &f, const ::std::vector<vec> &v, const ::std::vector<vec> &fNormals)
		{
			//MAX_PTS_PER_NODE = pts_per_node;
			build_with_faces(f, v, fNormals);
		}


		// Destructor - recursively frees the tree
		~BVH();

		Tree_Info find_tree_info();
		Traversal_Info intersect_face_from_raycast(float * p, float * dir,
			const float * pNormal);
	};
}; // namespace trimesh