#pragma once


#include <vector>
#include "Vec.h"
#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <algorithm>
//#include <GL/glut.h>
namespace trimesh {

	class Octree {
	private:
		//class Node;
		//Node *root;
		void build(const float *ptlist, size_t n);

		//Build the octree with vertices as input
		void build(std::vector<vec> pts, size_t n);
		//Build the octree with faces as input
		void build_with_faces(std::vector<TriMesh::Face> faces, std::vector<vec> pts , std::vector<vec> fNormals);
	public:
		class Node;
		Node *root;
		struct Traversal_Info {
			const float *origin_p, *origin_dir;
			vec closest_Pt;
			vec max;
			vec min;
			float closest_d;
			int level;
			int iChild;
			int faceID=-1;
			size_t k;
			unsigned char ray_neg_convert_bits;
			vec ptNormal;
			std::vector<vec> bBox_maxs;
			std::vector<vec> bBox_mins;
			vec vertex;
			vec ray;
			float origin_t;
			//vector<pt_with_d> knn;
		};

		const static int MAX_PTS_PER_NODE=10;

		const static int MAX_FACES_PER_NODE = 10;

		const static int MIN_DIST_INIT = 1000000;
		// Compatibility function for closest-compatible-point searches
		struct CompatFunc
		{
			virtual bool operator () (const float *p) const = 0;
			virtual ~CompatFunc() {}  // To make the compiler shut up
		};

		// Constructor from an array of points
		Octree(const float *ptlist, size_t n)
		{
			build(ptlist, n);
		}

		// Constructor from a vector of points---------------
		template <class T> Octree(const ::std::vector<T> &v ,int pts_per_node)
		{
			//MAX_PTS_PER_NODE = pts_per_node;
			build(v, v.size());
		}

		template <class T> Octree(const ::std::vector<T> &v)
		{
			//MAX_PTS_PER_NODE = pts_per_node;
			build(v, v.size());
		}

		// Constructor from a vector of faces -------------------
		template <class T> Octree(const ::std::vector<T> &f, const ::std::vector<vec> &v , const ::std::vector<vec> &fNormals)
		{
			//MAX_PTS_PER_NODE = pts_per_node;
			build_with_faces(f, v, fNormals);
		}


		// Destructor - recursively frees the tree
		~Octree();

		// The queries: returns closest point to a point or a ray,
		// provided it's within sqrt(maxdist2) and is compatible
		const float *closest_to_pt(const float *p,
			float maxdist2 = 0.0f,
			const CompatFunc *iscompat = NULL) const;
		const float *closest_to_ray(const float *p, const float *dir,
			float maxdist2 = 0.0f,
			const CompatFunc *iscompat = NULL) const;

		// Find the k nearest neighbors
		void find_k_closest_to_pt(::std::vector<const float *> &knn,
			int k,
			const float *p,
			float maxdist2 = 0.0f,
			const CompatFunc *iscompat = NULL) const;

		// Find the k nearest neighbors
		void find_leaf_closest_ray(const float *p, const float *dir,
			float maxdist2 = 0.0f /* = 0.0f */) const;

		Traversal_Info find_cube_from_raycast(const float * p, const float * dir, float maxdist2 = 0.0f) const;

		Traversal_Info intersect_face_from_raycast( float * p,  float * dir,const float *pNormal, float maxdist2 = 0.0f);

		Node * getRoot();

		Traversal_Info find_all_leaves();

	};

}; // namespace trimesh