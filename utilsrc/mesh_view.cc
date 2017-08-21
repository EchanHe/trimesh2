/*
Szymon Rusinkiewicz
Princeton University

mesh_view.cc
Simple viewer
*/

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#define FREEGLUT_STATIC

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "XForm.h"
#include "GLCamera.h"
#include "ICP.h"
#include "strutil.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <GL/glut.h>

//#include<omp.h>

#include <fstream>


#include "kdTree_face.h"
#include "octree.h"
#include "bvh.h"
//#include<geodesic\geodesic_mesh_elements.h>
//#include<geodesic\geodesic_constants_and_simple_functions.h>

using namespace std;
using namespace trimesh;


// Globals
vector<TriMesh *> meshes;
vector<xform> xforms;
vector<bool> visible;
vector<string> filenames;

//global geodist
std::vector<geodesic::SurfacePoint> path;

TriMesh::BSphere global_bsph;
xform global_xf;
GLCamera camera;

int current_mesh = -1;

bool draw_edges = false;
bool draw_points = false;
bool draw_2side = false;
bool draw_shiny = true;
bool draw_lit = false;
bool draw_falsecolor = true;
bool draw_index = false;
bool white_bg = true;
bool grab_only = false;
int point_size = 1, line_width = 1;

//--------yichen Draw
bool draw_normals = false;
bool draw_inwardNormals = true;
bool draw_geo_distance = true;
bool draw_feature = true;
bool draw_RayCasting = false;
bool draw_cone = true;
bool draw_Octree = false;

vec octRayO = vec(2,-0.5,-2);
vec octRayDir = vec(-1,1,1);
int octRay_index = 1;// 7999;// 49999;// 50000;


Octree* octree;
BVH * bvh;
BVH::Traversal_Info bvh_INFO;
trimesh::Octree::Traversal_Info octRayInfo;


template<class Points, class Faces>
void init_Geodesic_face_point(Points& points, Faces& faces, TriMesh *mesh)
{
	points.resize(mesh->vertices.size() * 3);
	faces.resize(mesh->faces.size() * 3);
	for (int i = 0; i < mesh->vertices.size(); i++) {
		int id = i * 3;
		points[id] = mesh->vertices[i][0];
		points[id + 1] = mesh->vertices[i][1];
		points[id + 2] = mesh->vertices[i][2];
		//points.push_back(mesh->vertices[i][0]);
		//points.push_back(mesh->vertices[i][1]);
		//points.push_back(mesh->vertices[i][2]);
	}

	for (int i = 0; i < mesh->faces.size(); i++) {
		int id = i * 3;
		faces[id] = mesh->faces[i][0];
		faces[id + 1] = mesh->faces[i][1];
		faces[id + 2] = mesh->faces[i][2];
	}

}
void make_geodesic_mesh(TriMesh *themesh , geodesic::Mesh &mesh) {
	std::vector<double> points;
	std::vector<unsigned> faces;
	init_Geodesic_face_point(points, faces, themesh);

	clock_t begin = clock();


	mesh.initialize_mesh_data(points, faces);

	cout << "Geodesic Mesh constructing time: " << double(clock() - begin) / CLOCKS_PER_SEC << std::endl;
	begin = clock();

}

// Signal a redraw
void need_redraw()
{
	glutPostRedisplay();
}


// Clear the screen
void cls()
{
	glDisable(GL_DITHER);
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_NORMALIZE);
	glDisable(GL_LIGHTING);
	glDisable(GL_NORMALIZE);
	glDisable(GL_COLOR_MATERIAL);
	if (white_bg)
		glClearColor(1, 1, 1, 0);
	else
		glClearColor(0.08f, 0.08f, 0.08f, 0);
	glClearDepth(1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}


// Set up lights and materials
void setup_lighting(int id)
{
	Color c(1.0f);
	if (draw_falsecolor)
		c = Color::hsv(-3.88f * id, 0.6f + 0.2f * sin(0.42f * id), 1);
	glColor3fv(c);

	if (!draw_lit || meshes[id]->normals.empty()) {
		glDisable(GL_LIGHTING);
		return;
	}

	GLfloat mat_specular[4] = { 0.18f, 0.18f, 0.18f, 0.18f };
	if (!draw_shiny) {
		mat_specular[0] = mat_specular[1] =
		mat_specular[2] = mat_specular[3] = 0.0f;
	}
	GLfloat mat_shininess[] = { 64 };
	GLfloat global_ambient[] = { 0.02f, 0.02f, 0.05f, 0.05f };
	GLfloat light0_ambient[] = { 0, 0, 0, 0 };
	GLfloat light0_diffuse[] = { 0.85f, 0.85f, 0.8f, 0.85f };
	if (current_mesh >= 0 && id != current_mesh) {
		light0_diffuse[0] *= 0.5f;
		light0_diffuse[1] *= 0.5f;
		light0_diffuse[2] *= 0.5f;
	}
	GLfloat light1_diffuse[] = { -0.01f, -0.01f, -0.03f, -0.03f };
	GLfloat light0_specular[] = { 0.85f, 0.85f, 0.85f, 0.85f };
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, draw_2side);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);
}


// Draw triangle strips.  They are stored as length followed by values.
void draw_tstrips(const TriMesh *themesh)
{
	static bool use_glArrayElement = false;
	static bool tested_renderer = false;
	if (!tested_renderer) {
		use_glArrayElement = !!strstr(
			(const char *) glGetString(GL_RENDERER), "Intel");
		tested_renderer = true;
	}

	const int *t = &themesh->tstrips[0];
	const int *end = t + themesh->tstrips.size();
	if (use_glArrayElement) {
		while (likely(t < end)) {
			glBegin(GL_TRIANGLE_STRIP);
			int striplen = *t++;
			for (int i = 0; i < striplen; i++)
				glArrayElement(*t++);
			glEnd();
		}
	} else {
		while (likely(t < end)) {
			int striplen = *t++;
			glDrawElements(GL_TRIANGLE_STRIP, striplen, GL_UNSIGNED_INT, t);
			t += striplen;
		}
	}
}


// Draw the mesh
void draw_mesh(int i)
{
	const TriMesh *themesh = meshes[i];

	glPushMatrix();
	glMultMatrixd(xforms[i]);

	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);

	if (draw_2side) {
		glDisable(GL_CULL_FACE);
	} else {
		glCullFace(GL_BACK);
		glEnable(GL_CULL_FACE);
	}

	// Vertices
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT,
			sizeof(themesh->vertices[0]),
			&themesh->vertices[0][0]);

	// Normals
	if (!themesh->normals.empty() && !draw_index) {
		glEnableClientState(GL_NORMAL_ARRAY);
		glNormalPointer(GL_FLOAT,
				sizeof(themesh->normals[0]),
				&themesh->normals[0][0]);
	} else {
		glDisableClientState(GL_NORMAL_ARRAY);
	}

	// Colors
	if (!themesh->colors.empty() && !draw_falsecolor && !draw_index) {
		glEnableClientState(GL_COLOR_ARRAY);
		glColorPointer(3, GL_FLOAT,
			       sizeof(themesh->colors[0]),
			       &themesh->colors[0][0]);
	} else {
		glDisableClientState(GL_COLOR_ARRAY);
	}


	// Main drawing pass
	if (draw_points || themesh->tstrips.empty()) {
		// No triangles - draw as points
		glPointSize(float(point_size));
		glDrawArrays(GL_POINTS, 0, themesh->vertices.size());
		glPopMatrix();
		return;
	}

	if (draw_edges) {
		glPolygonOffset(10.0f, 10.0f);
		glEnable(GL_POLYGON_OFFSET_FILL);
	}

	draw_tstrips(themesh);
	glDisable(GL_POLYGON_OFFSET_FILL);

	// Edge drawing pass
	if (draw_edges) {
		glPolygonMode(GL_FRONT, GL_LINE);
		glLineWidth(float(line_width));
		glDisableClientState(GL_COLOR_ARRAY);
		glDisable(GL_COLOR_MATERIAL);
		GLfloat global_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
		GLfloat light0_diffuse[] = { 0.8f, 0.8f, 0.8f, 0.0f };
		GLfloat light1_diffuse[] = { -0.2f, -0.2f, -0.2f, 0.0f };
		GLfloat light0_specular[] = { 0.0f, 0.0f, 0.0f, 0.0f };
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
		glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
		glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
		GLfloat mat_diffuse[4] = { 0.0f, 0.0f, 1.0f, 1.0f };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
		glColor3f(0, 0, 1); // Used iff unlit
		draw_tstrips(themesh);
		glPolygonMode(GL_FRONT, GL_FILL);
	}

	glPopMatrix();
}

void draw_path(int i ) {
	if (draw_geo_distance == false)
		return;
	TriMesh * m = meshes[i];

	glPushMatrix();

	for (int i = 0; i < path.size(); i++) {
		geodesic::SurfacePoint& s = path[i];
		//std::cout << s.x() << "\t" << s.y() << "\t" << s.z() << std::endl;
		double x = s.x();
		double y = s.y();
		double z = s.z();
		glPushMatrix();
		if(i==0 )
			glColor3f(0.0, 1.0, 0.0);
		else if(i + 1 == path.size())
			glColor3f(0.0, 0.0, 1.0);
		else
			glColor3f(1.0, 1.0, 1.0);
		glTranslatef(x, y, z);
		glutSolidSphere(0.1, 50, 50);
		glPopMatrix();

		if (i + 2 <= path.size()) {
			double x2 = path[i+1].x();
			double y2 = path[i + 1].y();
			double z2 = path[i + 1].z();

			glPushMatrix();
			glLineWidth(5.0);
			glColor3f(1.0, 0.0, 0.0);
			glBegin(GL_LINES);
			glVertex3f(x, y, z);
			glVertex3f(x2, y2, z2);
			glEnd();
			glPopMatrix();
			
		}

	}

	//----------draw all the verticies
	//for (int i = 0; i < m->vertices.size(); i++)
	//{
	//	glPushMatrix();
	//	glTranslatef(m->vertices[i][0], m->vertices[i][1], m->vertices[i][2]);
	//	glutSolidSphere(0.01, 50, 50);
	//	glPopMatrix();
	//}
	glPopMatrix();
}

//--draw coord
void drawCoord() {
	glPushMatrix();

	glLineWidth(5.0);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(50, 0, 0);
	glEnd();
	glPopMatrix();

	glPushMatrix();

	glLineWidth(5.0);
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 50, 0);
	glEnd();
	glPopMatrix();

	glPushMatrix();

	glLineWidth(5.0);
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 0, 50);
	glEnd();
	glPopMatrix();
}

//draw normals
void drawNormals(int i) {
	TriMesh * m = meshes[i];
	if (draw_normals) {
		//TriMesh * m = meshes[i];
		for (int i = 0; i < m->normals.size(); i++) {
			double x = m->vertices[i][0];
			double y = m->vertices[i][1];
			double z = m->vertices[i][2];

			vec endPoint = m->normals[i] + m->vertices[i];

			double x2 = endPoint[0];
			double y2 = endPoint[1];
			double z2 = endPoint[2];

			glPushMatrix();

			glLineWidth(2.0);
			glColor3f(0.0, 0.0, 1.0);
			glBegin(GL_LINES);
			glVertex3f(x, y, z);
			glVertex3f(x2, y2, z2);
			glEnd();
			glPopMatrix();

			if (draw_inwardNormals) {
				vec endPoint = m->vertices[i] - m->normals[i];
				double x2 = endPoint[0];
				double y2 = endPoint[1];
				double z2 = endPoint[2];

				glPushMatrix();

				glLineWidth(2.0);
				glColor3f(1.0, 0.0, 0.0);
				glBegin(GL_LINES);
				glVertex3f(x, y, z);
				glVertex3f(x2, y2, z2);
				glEnd();
				glPopMatrix();
			}
		}
	}
	bool draw_faceNormals=false;
	if (draw_faceNormals) {
		m->need_faceNormals();
		for (int i = 0; i < m->faceNormals.size(); i++) {
			vec origin = m->centroid(i);
			double x = origin[0];
			double y = origin[1];
			double z = origin[2];

			vec endPoint = (0.1f*m->faceNormals[i]) + origin;

			double x2 = endPoint[0];
			double y2 = endPoint[1];
			double z2 = endPoint[2];
			glPushMatrix();

			glLineWidth(1.0);
			glColor3f(1.0, 0.0, 0.0);
			//glScalef(0.1, 0.1, 0.1);
			glBegin(GL_LINES);
			glVertex3f(x, y, z);
			glVertex3f(x2, y2, z2);
			glEnd();
			glPopMatrix();
		}
	}
}
// draw cone
void drawCone(int i) {

	if (draw_cone) {
		TriMesh * m = meshes[i];
		int nv = m->vertices.size();
		int id = 0;// rand() % nv;

		double x = m->vertices[id][0];
		double y = m->vertices[id][1];
		double z = m->vertices[id][2];

		vec endPoint = m->normals[id] + m->vertices[id];
		//normalize(endPoint);
		double x2 = endPoint[0];
		double y2 = endPoint[1];
		double z2 = endPoint[2];

		double x3 = m->normals[id][0];
		double y3 = m->normals[id][1];
		double z3 = m->normals[id][2];
		glPushMatrix();
		//glTranslatef(x,y,z);
		glLineWidth(5.0);
		glColor3f(0.1, 0.5, 0.1);
		glBegin(GL_LINES);
		glVertex3f(x, y, z);
		glVertex3f(x2, y2, z2);
		glEnd();
		glPopMatrix();

		//---------draw inward normal
		endPoint = m->inwardNormals[id] + m->vertices[id];
		//normalize(endPoint);
		x2 = endPoint[0];
		y2 = endPoint[1];
		z2 = endPoint[2];

		glPushMatrix();
		//glTranslatef(x,y,z);
		glLineWidth(5.0);
		glColor3f(0.1, 0.1, 0.5);
		glBegin(GL_LINES);
		glVertex3f(x, y, z);
		glVertex3f(x2, y2, z2);
		glEnd();
		glPopMatrix();

		vector<vec> rays;
		//make_cone(45, 10, 0, m->vertices[id], m->normals[id], rays);
		make_cone(m->vertices[id], m->inwardNormals[id], rays);
		for (int j = 0; j < rays.size(); j++) {
			glPushMatrix();
			glLineWidth(2.0);
			glColor3f(0.0, 1.0, 0.0);
			glBegin(GL_LINES);

			//glVertex3f(0, 0, 0);
			//glVertex3f(v[0], v[1], v[2]);

			glVertex3f(x, y, z);
			glVertex3f(x + rays[j][0], y + rays[j][1], z + rays[j][2]);
			glEnd();
			glPopMatrix();
		}
	}
}

void drawFeaturePts(int i) {
	if(draw_feature){
		TriMesh * m = meshes[i];
		float lmSize = m->bsphere.r / 100;
		for (int i = 0; i < m->landmarks.size(); i++) {
			float x = m->landmarks[i][0]; float y = m->landmarks[i][1]; float z = m->landmarks[i][2];
			glPushMatrix();
				
			glColor3f(0.0, 0.0, 1.0);
			glTranslatef(x, y, z);
			glutSolidSphere(lmSize, 10, 10);
			glPopMatrix();
		}
		for (int i = 0; i < m->landmarks_2.size(); i++) {
			float x = m->landmarks_2[i][0]; float y = m->landmarks_2[i][1]; float z = m->landmarks_2[i][2];
			glPushMatrix();

			glColor3f(1.0, 0.0, 0.0);
			glTranslatef(x, y, z);
			glutSolidSphere(lmSize, 10, 10);
			glPopMatrix();
		}
		//int nv = m->vertices.size();
		//glPushMatrix();

		//glColor3f(0.0, 1.0, 0.0);
		//glTranslatef(m->vertices[nv-1][0], m->vertices[nv-1][1], m->vertices[nv-1][2]);
		//glutSolidSphere(10, 10, 10);
		//glPopMatrix();
		
	}
	
}

void drawRayCasting(int i) {
	if (draw_RayCasting) {

		TriMesh * m = meshes[i];

		octRayO = m->vertices[octRay_index];
		octRayDir = m->inwardNormals[octRay_index];

		bvh_INFO = bvh->intersect_face_from_raycast(octRayO, octRayDir, octRayO);
		float point_scale = 100.0f;

		//draw ray
		float scale = (m->bsphere.r * 3);
		vec end = (octRayO + scale *octRayDir);

		glPushMatrix();

		glColor3f(0.8, 0.0, 0.0);
		glTranslatef(octRayO[0], octRayO[1], octRayO[2]);
		//glutSolidSphere(scale/ point_scale, 10, 10);
		glPopMatrix();

		glPushMatrix();
		glLineWidth(5.0);
		glColor3f(0.1, 0.5, 0.1);
		glBegin(GL_LINES);
		glVertex3f(octRayO[0], octRayO[1], octRayO[2]);
		glVertex3f(end[0], end[1], end[2]);
		glEnd();
		glPopMatrix();

		// draw out bounding box ray

		end = (bvh_INFO.vertex + scale * bvh_INFO.ray);
		glPushMatrix();

		glColor3f(0.8, 0.0, 0.0);
		glTranslatef(bvh_INFO.vertex[0], bvh_INFO.vertex[1], bvh_INFO.vertex[2]);
		//glutSolidSphere(scale / point_scale, 10, 10);
		glPopMatrix();

		glPushMatrix();
		glLineWidth(5.0);
		glColor3f(0.8, 0.0, 0.0);
		glBegin(GL_LINES);
		glVertex3f(bvh_INFO.vertex[0], bvh_INFO.vertex[1], bvh_INFO.vertex[2]);
		glVertex3f(end[0], end[1], end[2]);
		glEnd();
		glPopMatrix();


		//--- draw nodes with ray crossing
		for (int i = 0; i < bvh_INFO.bBox_maxs.size(); i++) {
			float dx = bvh_INFO.bBox_maxs[i][0] - bvh_INFO.bBox_mins[i][0];
			float dy = bvh_INFO.bBox_maxs[i][1] - bvh_INFO.bBox_mins[i][1];
			float dz = bvh_INFO.bBox_maxs[i][2] - bvh_INFO.bBox_mins[i][2];

			float cx = (bvh_INFO.bBox_maxs[i][0] + bvh_INFO.bBox_mins[i][0]) / 2;
			float cy = (bvh_INFO.bBox_maxs[i][1] + bvh_INFO.bBox_mins[i][1]) / 2;
			float cz = (bvh_INFO.bBox_maxs[i][2] + bvh_INFO.bBox_mins[i][2]) / 2;
			glPushMatrix();
			glColor3d(0, 0, 0);
			glTranslatef(cx, cy, cz);
			glScalef(dx, dy, dz);
			glutWireCube(1);
			glPopMatrix();
		}

	}
	//----draw raycast of octree

	//if (draw_RayCasting) {

	//	TriMesh * m = meshes[i];
	//	
	//	octRayO = m->vertices[octRay_index];
	//	octRayDir = m->inwardNormals[octRay_index];

	//	octRayInfo = octree->intersect_face_from_raycast(octRayO, octRayDir, octRayO);
	//	float point_scale = 100.0f;
	//	
	//	//draw ray
	//	float scale = (m->bsphere.r*3);
	//	vec end =  (octRayO + scale *octRayDir);

	//	glPushMatrix();

	//	glColor3f(0.8, 0.0, 0.0);
	//	glTranslatef(octRayO[0], octRayO[1], octRayO[2]);
	//	//glutSolidSphere(scale/ point_scale, 10, 10);
	//	glPopMatrix();

	//	glPushMatrix();
	//	glLineWidth(5.0);
	//	glColor3f(0.1, 0.5, 0.1);
	//	glBegin(GL_LINES);
	//	glVertex3f(octRayO[0], octRayO[1], octRayO[2]);
	//	glVertex3f(end[0], end[1], end[2]);
	//	glEnd();
	//	glPopMatrix();
	//	
	//	//---draw reverse ray, if neccessary
	//	vec oRay = vec(octRayO[0], octRayO[1] , octRayO[2]); vec dirRay =  vec(octRayDir[0], octRayDir[1], octRayDir[2]);
	//	if (octRayDir[0] < 0.0f) {
	//		oRay[0] = m->bsphere.center[0]*2 - octRayO[0];
	//		dirRay[0] = -octRayDir[0];
	//	}
	//	if (octRayDir[1] < 0.0f) {
	//		oRay[1] = m->bsphere.center[1] * 2 - octRayO[1];
	//		dirRay[1] = -octRayDir[1];
	//	}
	//	if (octRayDir[2] < 0.0f) {
	//		oRay[2] = m->bsphere.center[2] * 2 - octRayO[2];
	//		dirRay[2] = -octRayDir[2];
	//	}
	//	if (octRayDir[0] < 0.0f || octRayDir[1] < 0.0f || octRayDir[2] < 0.0f) {
	//		end = (oRay + scale *dirRay);
	//		glPushMatrix();

	//		glColor3f(0.8, 0.0, 0.0);
	//		glTranslatef(oRay[0], oRay[1], oRay[2]);
	//		//glutSolidSphere(scale / point_scale, 10, 10);
	//		glPopMatrix();

	//		glPushMatrix();
	//		glLineWidth(5.0);
	//		glColor3f(0.1, 0.1, 0.5);
	//		glBegin(GL_LINES);
	//		glVertex3f(oRay[0], oRay[1], oRay[2]);
	//		glVertex3f(end[0], end[1], end[2]);
	//		glEnd();
	//		glPopMatrix();
	//	}

	//	// draw out bounding box ray

	//	end = (octRayInfo.vertex + scale * octRayInfo.ray);
	//	glPushMatrix();

	//	glColor3f(0.8, 0.0, 0.0);
	//	glTranslatef(octRayInfo.vertex[0], octRayInfo.vertex[1], octRayInfo.vertex[2]);
	//	//glutSolidSphere(scale / point_scale, 10, 10);
	//	glPopMatrix();

	//	glPushMatrix();
	//	glLineWidth(5.0);
	//	glColor3f(0.8, 0.0, 0.0);
	//	glBegin(GL_LINES);
	//	glVertex3f(octRayInfo.vertex[0], octRayInfo.vertex[1], octRayInfo.vertex[2]);
	//	glVertex3f(end[0], end[1], end[2]);
	//	glEnd();
	//	glPopMatrix();


	//	////----draw bounding box
	//	//float r = m->bsphere.r;
	//	//vec center = m->bsphere.center;
	//	////glPushMatrix();
	//	////glTranslatef(center[0]+r/2, center[1], center[2]);
	//	////glutWireCube(r);
	//	////glPopMatrix();
	//	//float x[8], y[8], z[8];
	//	//for (int i = 0; i < 8; i++) {
	//	//	if(i<=3)
	//	//		x[i] = center[0] - r / 2;
	//	//	else
	//	//		x[i] = center[0] + r / 2;
	//	//}
	//	//y[0] = center[1] - r / 2; y[1] = center[1] - r / 2; y[4] = center[1] - r / 2; y[5] = center[1] - r / 2;
	//	//y[2] = center[1] + r / 2; y[3] = center[1] + r / 2; y[6] = center[1] + r / 2; y[7] = center[1] + r / 2;

	//	//z[0] = center[2] - r / 2; z[2] = center[2] - r / 2; z[4] = center[2] - r / 2; z[6] = center[2] - r / 2;
	//	//z[1] = center[2] + r / 2; z[3] = center[1] + r / 2; z[5] = center[2] + r / 2; z[7] = center[2] + r / 2;
	//	//double red[4], green[4];
	//	//red[0] = 1.0; red[1] = 1.0; red[2] = 1.0; red[3] = 0.0;
	//	//green[0] = 0.0; green[1] = 0.5; green[2] = 1.0; green[3] = 1.0;
	//	//for (int i = 0; i < 8; i++) {

	//	//	glPushMatrix();
	//	//	if (i < 4) 
	//	//		glColor3d(red[i], green[i], 0);
	//	//	else
	//	//		glColor3d(0, 0, 0);
	//	//	glTranslatef(x[i], y[i], z[i]);
	//	//	glutWireCube(r);
	//	//	glPopMatrix();
	//	//}


	//	//--- draw nodes with ray crossing
	//	for (int i = 0; i < octRayInfo.bBox_maxs.size(); i++) {
	//		float dx = octRayInfo.bBox_maxs[i][0] - octRayInfo.bBox_mins[i][0];
	//		float dy = octRayInfo.bBox_maxs[i][1] - octRayInfo.bBox_mins[i][1];
	//		float dz = octRayInfo.bBox_maxs[i][2] - octRayInfo.bBox_mins[i][2];

	//		float cx = (octRayInfo.bBox_maxs[i][0] + octRayInfo.bBox_mins[i][0]) / 2;
	//		float cy = (octRayInfo.bBox_maxs[i][1] + octRayInfo.bBox_mins[i][1]) / 2;
	//		float cz = (octRayInfo.bBox_maxs[i][2] + octRayInfo.bBox_mins[i][2]) / 2;
	//		glPushMatrix();
	//		glColor3d(0,0, 0);
	//		glTranslatef(cx, cy, cz);
	//		glScalef(dx, dy, dz);
	//		glutWireCube(1);
	//		glPopMatrix();
	//	}

	//}
}

void drawOctree(int i) {
	TriMesh * m = meshes[i];
	
	if (draw_Octree) {
		octRayInfo = octree->find_all_leaves();
		for (int i = 0; i < octRayInfo.bBox_maxs.size(); i++) {
			float dx = octRayInfo.bBox_maxs[i][0] - octRayInfo.bBox_mins[i][0];
			float dy = octRayInfo.bBox_maxs[i][1] - octRayInfo.bBox_mins[i][1];
			float dz = octRayInfo.bBox_maxs[i][2] - octRayInfo.bBox_mins[i][2];

			float cx = (octRayInfo.bBox_maxs[i][0] + octRayInfo.bBox_mins[i][0]) / 2;
			float cy = (octRayInfo.bBox_maxs[i][1] + octRayInfo.bBox_mins[i][1]) / 2;
			float cz = (octRayInfo.bBox_maxs[i][2] + octRayInfo.bBox_mins[i][2]) / 2;
			glPushMatrix();
			glTranslatef(cx, cy, cz);
			glScalef(dx, dy, dz);
			glutWireCube(1);
			glPopMatrix();
		}

	}
}



// Draw the scene
void redraw()
{
	timestamp t = now();
	camera.setupGL(global_xf * global_bsph.center, global_bsph.r);
	glPushMatrix();
	glMultMatrixd(global_xf);
	cls();
	for (size_t i = 0; i < meshes.size(); i++) {
		if (!visible[i])
			continue;
		setup_lighting(i);
		draw_mesh(i);
		draw_path(i);
		drawCoord();
		drawNormals(i);
		drawCone(i);
		drawFeaturePts(i);
		drawRayCasting(i);
		drawOctree(i);
	}

	glPopMatrix();
	glutSwapBuffers();
	if (grab_only) {
		void dump_image();
		dump_image();
		exit(0);
	}
	printf("\r                        \r%.1f msec.", 1000.0f * (now() - t));
	fflush(stdout);
}


// Update global bounding sphere.
void update_bsph()
{
	point boxmin(1e38f, 1e38f, 1e38f);
	point boxmax(-1e38f, -1e38f, -1e38f);
	bool some_vis = false;
	for (size_t i = 0; i < meshes.size(); i++) {
		if (!visible[i])	
			continue;
		some_vis = true;
		point c = xforms[i] * meshes[i]->bsphere.center;
		float r = meshes[i]->bsphere.r;
		for (int j = 0; j < 3; j++) {
			boxmin[j] = min(boxmin[j], c[j]-r);
			boxmax[j] = max(boxmax[j], c[j]+r);
		}
	}
	if (!some_vis)
		return;
	point &gc = global_bsph.center;
	float &gr = global_bsph.r;
	gc = 0.5f * (boxmin + boxmax);
	gr = 0.0f;
	for (size_t i = 0; i < meshes.size(); i++) {
		if (!visible[i])	
			continue;
		point c = xforms[i] * meshes[i]->bsphere.center;
		float r = meshes[i]->bsphere.r;
		gr = max(gr, dist(c, gc) + r);
	}
}


// Set the view...
void resetview()
{
	camera.stopspin();

	// Reloc mesh xforms
	for (size_t i = 0; i < meshes.size(); i++)
		if (!xforms[i].read(xfname(filenames[i])))
			xforms[i] = xform();

	update_bsph();

	// Set camera to first ".camxf" if we have it...
	for (size_t i = 0; i < filenames.size(); i++) {
		if (global_xf.read(replace_ext(filenames[i], "camxf")))
			return;
	}

	// else default view
	global_xf = xform::trans(0, 0, -5.0f * global_bsph.r) *
		    xform::trans(-global_bsph.center);
}


// Make some mesh current
void set_current(int i)
{
	camera.stopspin();
	if (i >= 0 && i < (int)meshes.size() && visible[i])
		current_mesh = i;
	else
		current_mesh = -1;
	need_redraw();
}


// Change visiblility of a mesh
void toggle_vis(int i)
{
	if (i >= 0 && i < (int)meshes.size())
		visible[i] = !visible[i];
	if (current_mesh == i && !visible[i])
		set_current(-1);
	update_bsph();
	need_redraw();
}


// Save the current image to a PPM file.
// Uses the next available filename matching filenamepattern
void dump_image()
{
	// Find first non-used filename
	const char filenamepattern[] = "img%d.ppm";
	int imgnum = 0;
	FILE *f;
	while (1) {
		char filename[1024];
		sprintf(filename, filenamepattern, imgnum++);
		f = fopen(filename, "rb");
		if (!f) {
			f = fopen(filename, "wb");
			printf("\n\nSaving image %s... ", filename);
			fflush(stdout);
			break;
		}
		fclose(f);
	}

	// Read pixels
	GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	GLint width = V[2], height = V[3];
	char *buf = new char[width*height*3];
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(V[0], V[1], width, height, GL_RGB, GL_UNSIGNED_BYTE, buf);

	// Flip top-to-bottom
	for (int i = 0; i < height/2; i++) {
		char *row1 = buf + 3 * width * i;
		char *row2 = buf + 3 * width * (height - 1 - i);
		for (int j = 0; j < 3 * width; j++)
			swap(row1[j], row2[j]);
	}

	// Write out file
	fprintf(f, "P6\n#\n%d %d\n255\n", width, height);
	fwrite(buf, width*height*3, 1, f);
	fclose(f);
	delete [] buf;

	printf("Done.\n\n");
}


// Save scan transforms
void save_xforms()
{
	for (size_t i = 0; i < xforms.size(); i++) {
		string xffile = xfname(filenames[i]);
		printf("Writing %s\n", xffile.c_str());
		xforms[i].write(xffile);
	}
}


// Save camera xform
void save_cam_xform()
{
	std::string camfile = replace_ext(filenames[0], "camxf");
	printf("Writing %s\n", camfile.c_str());
	global_xf.write(camfile);
}


// ICP
void do_icp(int n)
{
	camera.stopspin();

	if (current_mesh < 0 || current_mesh >= (int)meshes.size())
		return;
	if (n < 0 || n >= (int)meshes.size())
		return;
	if (!visible[n] || !visible[current_mesh] || n == current_mesh)
		return;
	ICP(meshes[n], meshes[current_mesh], xforms[n], xforms[current_mesh], 2);
	update_bsph();
	need_redraw();
}


// Handle mouse button and motion events
static unsigned buttonstate = 0;

void doubleclick(int button, int x, int y)
{
	// Render and read back ID reference image
	camera.setupGL(global_xf * global_bsph.center, global_bsph.r);
	glDisable(GL_BLEND);
	glDisable(GL_LIGHTING);
	glClearColor(1,1,1,1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	draw_index = true;
	glPushMatrix();
	glMultMatrixd(global_xf);
	for (size_t i = 0; i < meshes.size(); i++) {
		if (!visible[i])
			continue;
		glColor3ub((i >> 16) & 0xff,
			   (i >> 8)  & 0xff,
			    i        & 0xff);
		draw_mesh(i);
	}
	glPopMatrix();
	draw_index = false;
	GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	y = int(V[1] + V[3]) - 1 - y;
	unsigned char pix[3];
	glReadPixels(x, y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, pix);
	int n = (pix[0] << 16) + (pix[1] << 8) + pix[2];

	if (button == 0 || buttonstate == (1 << 0)) {
		// Double left click - select a mesh
		set_current(n);
	} else if (button == 2 || buttonstate == (1 << 2)) {
		// Double right click - ICP current to clicked-on
		do_icp(n);
	}
}

void mousehelperfunc(int x, int y)
{
	static const Mouse::button physical_to_logical_map[] = {
		Mouse::NONE, Mouse::ROTATE, Mouse::MOVEXY, Mouse::MOVEZ,
		Mouse::MOVEZ, Mouse::MOVEXY, Mouse::MOVEXY, Mouse::MOVEXY,
	};

	Mouse::button b = Mouse::NONE;
	if (buttonstate & (1 << 3))
		b = Mouse::WHEELUP;
	else if (buttonstate & (1 << 4))
		b = Mouse::WHEELDOWN;
	else if (buttonstate & (1 << 30))
		b = Mouse::LIGHT;
	else
		b = physical_to_logical_map[buttonstate & 7];

	if (current_mesh < 0) {
		camera.mouse(x, y, b,
			     global_xf * global_bsph.center, global_bsph.r,
			     global_xf);
	} else {
		xform tmp_xf = global_xf * xforms[current_mesh];
		camera.mouse(x, y, b,
			     tmp_xf * meshes[current_mesh]->bsphere.center,
			     meshes[current_mesh]->bsphere.r,
			     tmp_xf);
		xforms[current_mesh] = inv(global_xf) * tmp_xf;
		update_bsph();
	}
}

void mousemotionfunc(int x, int y)
{
	mousehelperfunc(x,y);
	if (buttonstate)
		need_redraw();
}

void mousebuttonfunc(int button, int state, int x, int y)
{
	static timestamp last_click_time;
	static unsigned last_click_buttonstate = 0;
	static float doubleclick_threshold = 0.4f;

	if (glutGetModifiers() & GLUT_ACTIVE_CTRL)
		buttonstate |= (1 << 30);
	else
		buttonstate &= ~(1 << 30);

	if (state == GLUT_DOWN) {
		buttonstate |= (1 << button);
		if (buttonstate == last_click_buttonstate &&
		    now() - last_click_time < doubleclick_threshold) {
			doubleclick(button, x, y);
			last_click_buttonstate = 0;
		} else {
			last_click_time = now();
			last_click_buttonstate = buttonstate;
		}
	} else {
		buttonstate &= ~(1 << button);
	}

	mousehelperfunc(x, y);
	if (buttonstate & ((1 << 3) | (1 << 4))) // Wheel
		need_redraw();
	if ((buttonstate & 7) && (buttonstate & (1 << 30))) // Light
		need_redraw();
}


// Idle callback
void idle()
{
	xform tmp_xf = global_xf;
	if (current_mesh >= 0)
		tmp_xf = global_xf * xforms[current_mesh];

	if (camera.autospin(tmp_xf))
		need_redraw();
	else
		usleep(10000);

	if (current_mesh >= 0) {
		xforms[current_mesh] = inv(global_xf) * tmp_xf;
		update_bsph();
	} else {
		global_xf = tmp_xf;
	}
}
int file_index = 0;
const int fileN = 5;
char * filenames_input[fileN] = { "..\\..\\..\\data\\beak_sdf.ply" ,"..\\..\\..\\data\\beak_sdf2.ply" ,"..\\..\\..\\data\\frog.ply",
"..\\..\\..\\data\\2934_smooth_sdf.ply" ,"..\\..\\..\\data\\bird_one_component_sdf.ply" };
//char * f = filenames[0];
//read mesh()
void readMesh(int file_index) {
	TriMesh *themesh = TriMesh::read(filenames_input[file_index]);
	themesh->need_bsphere();
	for (int i = 0; i < themesh->vertices.size(); i++) {
		for (int j = 0; j < 3; j++) {
			themesh->vertices[i][j] = themesh->vertices[i][j] - themesh->bsphere.center[j];
		}

	}

	themesh->bsphere = {};
	themesh->normals.clear();
	themesh->tstrips.clear();


	themesh->colors.clear();
	themesh->need_faceMidPts();
	themesh->need_normals();
	themesh->need_inwardNormals();
	themesh->need_tstrips();
	themesh->need_bsphere();
	themesh->need_faceNormals();


	themesh->need_curvatures();

	//------FEATURES EXTRACTION
	//Assign the Vertex quality to sdf
	themesh->sdf = themesh->quality;

	if (file_index < 2)
		themesh->need_feature_points_curv_sdf_beak();
	else
		themesh->need_feature_points_curv_sdf();
	//std::cout << "The size of feature points :" << themesh->landmarkID.size() << endl;
	//std::cout << "colors? :" << themesh->colors.size() << endl;

	//------COLORING PART--------------

	//themesh->color_vertex(themesh->landmarkID);
	//themesh->color_vertex(themesh->quality);
	//themesh->color_vertex(themesh->sdf_brute,false);
	themesh->color_vertex(themesh->sdf, false);
	//cout << themesh->curv1[0];
	//calculate the color


	//cout << compSizes[0];
	meshes.clear();
	meshes.push_back(themesh);
	xforms.push_back(xform());
	visible.push_back(true);
	filenames.push_back(filenames_input[file_index]);
}

int mesh_attri_index = 0;
const int attriN = 3;

// Keyboard
#define Ctrl (1-'a')
void keyboardfunc(unsigned char key, int, int)
{
	switch (key) {
		case ' ':
			if (current_mesh < 0)
				resetview();
			else
				set_current(-1);
			break;
		case '@': // Shift-2
			draw_2side = !draw_2side; break;
		case 'e':
			draw_edges = !draw_edges; break;
		case 'f':
			draw_falsecolor = !draw_falsecolor; break;
		case Ctrl + 'f':
			draw_feature = !draw_feature; break;
		case 'l':
			draw_lit = !draw_lit; break;
		case 'p':
			point_size++; break;
		case 'P':
			point_size = max(1, point_size - 1); break;
		case Ctrl+'p':
			draw_points = !draw_points; break;
		case 's':
			draw_shiny = !draw_shiny; break;
		case 't':
			line_width++; break;
		case 'T':
			line_width = max(1, line_width - 1); break;
		case 'w':
			white_bg = !white_bg; break;
		case 'I':
			dump_image(); break;
		case Ctrl+'x':
			save_xforms();
			break;
		case Ctrl+'c':
			save_cam_xform();
			break;
		case '\033': // Esc
		case Ctrl+'q':
		case 'Q':
		case 'q':
			exit(0);
		case '0':
			toggle_vis(9); break;
		case '-':
			toggle_vis(10); break;
		case '=':
			toggle_vis(11); break;
		case '[': {
			//octRayDir=rotation(octRayDir, 10, vec(0,0,1));; break;
			octRayO[0] += 0.1;
			break;
		}
		case ']': {
			//octRayDir = rotation(octRayDir, -10, vec(0, 0, 1));; break;
			octRayO[0] -= 0.1;
			break;
		}
		case '{': {
			//octRayDir = rotation(octRayDir, 10, vec(0, 1, 0));; break;
			octRayO[1] += 0.1;
			break;
		}
		case '}': {
			octRayO[1] -= 0.1;
			break;
			//octRayDir = rotation(octRayDir, -10, vec(0, 1, 0));; break;
		}
		case '>': {
			octRay_index = (octRay_index + 1) % meshes[0]->vertices.size();
			//readMesh((file_index++)%fileN);
			break;
			//octRayDir = rotation(octRayDir, -10, vec(0, 1, 0));; break;
		}
		case '<': {
			//octRay_index = (octRay_index - 1) % meshes[0]->vertices.size();
		//	readMesh((file_index--) % fileN);
			break;
			//octRayDir = rotation(octRayDir, -10, vec(0, 1, 0));; break;
		}
		case'c': {
			mesh_attri_index = (mesh_attri_index +1) % attriN;
			meshes[0]->colors.clear();
			switch (mesh_attri_index) {
			case 0:		
				cout << "showing gaussian curv ..." << endl;
				meshes[0]->color_vertex(meshes[0]->gaus_curv, true);
				break;
			case 1:
				cout << "showing mean curv ..." << endl;
				meshes[0]->color_vertex(meshes[0]->mean_curv, true);
				break;
			case 2:
				cout << "showing shape diameter function ..." << endl;
				meshes[0]->color_vertex(meshes[0]->sdf, true);
				break;
			}
		}
		default:
			if (key >= '1' && key <= '9') {
				int m = key - '1';
				toggle_vis(m);
			}
	}
	need_redraw();
}


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s [-grab] infile...\n", myname);
	exit(1);
}



int main(int argc, char *argv[])
{
	glutInitWindowSize(512, 512);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInit(&argc, argv);

	if (argc < 2)
		usage(argv[0]);

	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-grab")) {
			grab_only = true;
			continue;
		}
		char *filename = argv[i];
		//filename = "..\\..\\..\\data\\torus_sdf.ply";
		//filename = "..\\..\\..\\data\\external_surface_only.ply";
	//	filename = "..\\..\\..\\data\\bird_one_component_sdf.ply";
	//	filename = "..\\..\\..\\data\\beak_sdf3.ply";
		filename = "..\\..\\..\\data\\torus2.ply";
		//filename = "..\\..\\..\\data\\cow2.ply";
	//filename = "..\\..\\..\\data\\elephant.off";
	//	filename = "..\\..\\..\\data\\cube.ply";
	//filename = "..\\..\\..\\data\\Armadillo.ply";
	//filename = "..\\..\\..\\data\\frog.ply";
//		filename = "..\\..\\..\\data\\2934_smooth_sdf.ply";
	//	filename = "..\\..\\..\\data\\Armadillo.ply";
	//	const char *filename_simple = "..\\..\\..\\data\\bird_one_component_simple_1.ply";
		TriMesh *themesh = TriMesh::read(filename);
		//TriMesh *simple = TriMesh::read(filename);
		//simple->need_faceNormals();

		if (!themesh)
			usage(argv[0]);
		//trimesh::lmsmooth(themesh, 40);
		//themesh->need_bsphere();
		//for (int i = 0; i < themesh->vertices.size(); i++) {
		//	for (int j = 0; j < 3; j++) {
		//		themesh->vertices[i][j] = themesh->vertices[i][j] -themesh->bsphere.center[j];
		//	}
		//	
		//}

		themesh->bsphere = {};
		themesh->normals.clear();
		themesh->tstrips.clear();
		themesh->colors.clear();
		themesh->need_faceMidPts();
		themesh->need_normals();
		themesh->need_inwardNormals();
		themesh->need_tstrips();
		themesh->need_bsphere();
		themesh->need_faceNormals();
		themesh->need_neighbors();

		themesh->need_curvatures();


		//KDtree_face * kdface_tree = new KDtree_face(themesh->faces, themesh->vertices, themesh->faceNormals);
		//kdface_tree->intersect_face_from_raycast(themesh->vertices[1], themesh->inwardNormals[1], themesh->normals[1]);


		//KD_tree aaa;
		//trimesh::buildKDTree(aaa, themesh->faces, themesh->vertices, themesh->faceNormals);
		//int height = heightKD_Tree(aaa.root);
		//KD_tree_array * kd_array =  KDTreeToArray(aaa);
		//rayToFaces(kd_array , themesh->vertices[1] , themesh->inwardNormals[1],themesh->normals[1]);
		
	
		themesh->need_sdf_bvh();
		//test octree

		//octree = new Octree(themesh->vertices, 1);
	//	Octree * octreeFace = new Octree(themesh->faces, themesh->vertices, themesh->faceNormals);
	/*	octRayInfo = octreeFace->intersect_face_from_raycast(themesh->vertices[0], themesh->inwardNormals[0], themesh->normals[0]);
*/
		//octree = new Octree(themesh->faces, themesh->vertices, themesh->faceNormals);
		//Octree::Node * a =octree->getRoot();
		
		//--ray is outside of the shape
		//octRayInfo = octree->find_cube_from_raycast(octRayO, octRayDir);

		//--ray is inside the shape
		//octRayInfo = octree->find_cube_from_raycast(themesh->vertices[0], themesh->inwardNormals[0]);
	//	octRayInfo = octree->find_cube_from_raycast(octRayO, octRayDir);
		//themesh->need_sdf_from_simple(simple);
	//themesh->need_sdf_octree();
	//	themesh->compare_sdfs();
		//themesh->writeSDF();
		
		//themesh->need_sdf_brute();
		//cout << "number of sdf " << themesh->sdf_brute.size()<<endl;
		//themesh->colors.resize(themesh->vertices.size());
		//themesh->colors[themesh->faces[5500][0]] =Color(1.0, 0.0, 0.0);
		//themesh->colors[themesh->faces[5500][1]] = Color(1.0, 0.0, 0.0);
		//themesh->colors[themesh->faces[5500][2]] = Color(1.0, 0.0, 0.0);

		//bvh = new BVH(themesh->faces, themesh->vertices, themesh->faceNormals);
		//BVH::Tree_Info b_info = bvh->find_tree_info();

		//themesh->need_sdf_bvh();
		//themesh->need_sdf_kd();
		//themesh->sdf.clear();
		//themesh->need_sdf_octree();
		//b->intersect_face_from_raycast(themesh->vertices[0], themesh->inwardNormals[0] , themesh->normals[0]);

		//------FEATURES EXTRACTION
		//Assign the Vertex quality to sdf
		//if (themesh->quality.size() != 0) {
		//	themesh->sdf = themesh->quality;
		//}
		
		themesh->need_feature_points_curv_sdf_beak();
		//std::cout << "The size of feature points :" << themesh->landmarkID.size() << endl;
		//std::cout << "colors? :" << themesh->colors.size() << endl;

		//------COLORING PART--------------

		//themesh->color_vertex(themesh->landmarkID);
		//themesh->color_vertex(themesh->quality);
		//themesh->color_vertex(themesh->sdf_brute,false);
	/*	themesh->color_vertex(themesh->gaus_curv, true);
		themesh->color_vertex(themesh->mean_curv, true);*/
		themesh->color_vertex(themesh->sdf, false);
		//cout << themesh->curv1[0];
		//calculate the color

		themesh->writeAttri();

		//cout << compSizes[0];
		
		meshes.push_back(themesh);
		xforms.push_back(xform());
		visible.push_back(true);
		filenames.push_back(filename);
	

		


		//---------geodesic distance

		
		//std::vector<double> points;
		//std::vector<unsigned> faces;
		//init_Geodesic_face_point(points, faces, themesh);
		
		//geodesic::Mesh geoMesh;

		//trimesh::init_geod_mesh(themesh, geoMesh);

		//KDtree *kd = new KDtree(themesh->vertices);
		//const float *a=kd->closest_to_pt(themesh->bsphere.center);
	//	//path = trimesh::cal_geo_dis(5, 6, points, faces);
	//	//trimesh::cal_geo_dis(points, faces);
		
		
	//	print_info_about_path(path);
		//for (unsigned i = 0; i<path.size(); ++i)
		//{
		//	geodesic::SurfacePoint& s = path[i];
		//	std::cout << s.x() << "\t" << s.y() << "\t" << s.z() << std::endl;
		//	
		//}
		//std::cout << "length " << geodesic::length(path) << endl;
		//themesh->write("a.ply");
	}

	glutCreateWindow(argv[1]);
	glutDisplayFunc(redraw);
	glutMouseFunc(mousebuttonfunc);
	glutMotionFunc(mousemotionfunc);
	glutKeyboardFunc(keyboardfunc);
	glutIdleFunc(idle);

	resetview();
	glutMainLoop();
}


