# TriMesh2 
with GPU KD tree and Shape diameter funciton.

gpu branch


**TriMesh2:** C++ library and set of utilities for input, output, and basic manipulation of 3D triangle meshes


This is a fork of the [TriMesh2 library](http://gfx.cs.princeton.edu/proj/trimesh2/) (originally by [Szymon Rusinkiewicz](https://www.cs.princeton.edu/~smr/)), which I use a lot in my other graphics projects. I like TriMesh2 because of the low setup costs required to do model loading and manipulation. I'm staying close to the original, and adding my fixes for MSVC and 64-bit compilation. See further down for details.

 * For the original TriMesh2 project, see [the Trimesh2 homepage](http://gfx.cs.princeton.edu/proj/trimesh2/).
 * Original MSVC 2012 project by Bengt Rosenberger.
 
## News
 * mesh_view tool now builds in Win64 without external dependencies!

## Info
Features: 

 * Support for reading/writing PLY, OFF, OBJ files. Read-only: 3DS, SM, RAY.
 * Vec: a templated C++ class for constant-length vectors, with support for the usual arithmetic operations and XForm: a class for rigid-body transformations.
 * An OpenGL trackball/arcball implementation, with automatic selection of rotation center.
 * Algorithms for subdivision, smoothing, curvature estimation, triangle stripping, and various other simple mesh manipulations.

The following utility programs are included:

 * **mesh_view**: a simple 3D mesh viewer. And shows the landmarks using mean curvature, gaussian curvature and shape diameter function.
 * **mesh_make**: create arbitrarily-tessellated meshes of various simple shapes
 * **mesh_filter**: applies a variety of simple transformations to a mesh, such as converting formats, flipping faces, subdivision, smoothing, rigid-body transformations, etc.
 * **mesh_cc**: list and/or extract connected components from a mesh
 * **mesh_cat**: combine several meshes into a single file
 * **mesh_align**: align 2 meshes using ICP
 * **mesh_shade**: a few procedural shaders for adding per-vertex color
 * **mesh_check**: check for some kinds of topological oddities (e.g., more than 2 faces at an edge) in a mesh file.
 * **mesh_crunch**: quick-n-dirty mesh decimation using the Rossignac-Borrel method of vertex collapse
 * **mesh_info**: print out some information about a mesh
 * **xf**: create or compose transformations in .xf files

## Fork Details

The original TriMesh2 project was very Linux/GCC-oriented, and only had limited Win32 support. The primary aim of this fork is to provide a stable MSVC solution, for both x86 and x64 targets. For several tools, TriMesh2 depends on a custom-cherrypicked freeglut version - dragging this kicking and screaming into the new century is also part of the challenge.


For todo/planned features, see todo.txt.
