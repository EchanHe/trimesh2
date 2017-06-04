// mesh_smooth.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <cstdio>
#include <cstdlib>

using namespace std;
using namespace trimesh;

void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s in.ply out.ply\n", myname);
	exit(1);
}

int main(int argc, char *argv[])
{
	if (argc < 3)
		usage(argv[0]);
	TriMesh *in = TriMesh::read(argv[1]);
	if (!in)
		usage(argv[0]);

	bool had_tstrips = !in->tstrips.empty();

	//TriMesh *out = new TriMesh;
	//float voxelsize = 2.0f * in->feature_size();
	//crunch(in, out, voxelsize);

	lmsmooth(in, 20);

	if (had_tstrips)
		in->need_tstrips();
	in->write(argv[2]);
}

