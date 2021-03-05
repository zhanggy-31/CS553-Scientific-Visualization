#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>


#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "polyhedron.h"
#include "trackball.h"
#include "tmatrix.h"

Polyhedron* poly;

/*scene related variables*/
const float zoomspeed = 0.9;
const int win_width = 1024;
const int win_height = 1024;
const int view_mode = 0;		// 0 = othogonal, 1=perspective
const double radius_factor = 0.9;

/*
Use keys 1 to 0 to switch among different display modes.
Each display mode can be designed to show one type 
visualization result.

Predefined ones: 
display mode 1: solid rendering
display mode 2: show wireframes
display mode 3: render each quad with colors of vertices
*/
int display_mode = 1;

/*User Interaction related variabes*/
float s_old, t_old;
float rotmat[4][4];
double zoom = 1.0;
double translation[2] = { 0, 0 };
int mouse_mode = -2;	// -1 = no action, 1 = tranlate y, 2 = rotate

/******************************************************************************
Forward declaration of functions
******************************************************************************/

void init(void);

/*glut attaching functions*/
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void mousewheel(int wheel, int direction, int x, int y);

/*functions for element picking*/
void display_vertices(GLenum mode, Polyhedron* poly);
void display_quads(GLenum mode, Polyhedron* poly);
void display_selected_vertex(Polyhedron* poly);
void display_selected_quad(Polyhedron* poly);

/*display vis results*/
void display_polyhedron(Polyhedron* poly);

/******************************************************************************
Main program.
******************************************************************************/
int main(int argc, char* argv[])
{
	/*load mesh from ply file*/
	

	/*read filename in folder "scalar_data"*/

	FILE* this_file = fopen("../quadmesh_2D/scalar_data/2x_square_plus_y_square.ply", "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);

	
	/*initialize the mesh*/
	poly->initialize(); // initialize the mesh
	poly->write_info();


	/*init glut and create window*/
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Scientific Visualization");


	/*initialize openGL*/
	init();
	
	/*the render function and callback registration*/
	glutKeyboardFunc(keyboard);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMouseWheelFunc(mousewheel);
	
	/*event processing loop*/
	glutMainLoop();
	
	/*clear memory before exit*/
	poly->finalize();	// finalize everything
	
	return 0;
}


/******************************************************************************
Set projection mode
******************************************************************************/

void set_view(GLenum mode)
{
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };

	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (view_mode == 0)
		glOrtho(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom, radius_factor * zoom, 0.0, 40.0);
	else
		glFrustum(-radius_factor * zoom, radius_factor * zoom, -radius_factor, radius_factor, -1000, 1000);

	GLfloat light_position[3];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

/******************************************************************************
Update the scene
******************************************************************************/

void set_scene(GLenum mode, Polyhedron* poly)
{
	glTranslatef(translation[0], translation[1], -3.0);

	/*multiply rotmat to current mat*/
	{
		int i, j, index = 0;

		GLfloat mat[16];

		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
				mat[index++] = rotmat[i][j];

		glMultMatrixf(mat);
	}

	glScalef(0.9 / poly->radius, 0.9 / poly->radius, 0.9 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}


/******************************************************************************
Init scene
******************************************************************************/

void init(void) {

	mat_ident(rotmat);

	/* select clearing color */
	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	
	//set pixel storage modes
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	
	glEnable(GL_NORMALIZE);
	if (poly->orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}


/******************************************************************************
Pick objects from the scene
******************************************************************************/

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, * ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	ptr = (GLuint*)buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	return seed_id;
}

/******************************************************************************
Diaplay all quads for selection
******************************************************************************/

void display_quads(GLenum mode, Polyhedron* this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	//glDisable(GL_LIGHTING);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	for (i = 0; i < this_poly->nquads; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Quad* temp_q = this_poly->qlist[i];
		{
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;
		}
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		
		glBegin(GL_POLYGON);
		for (j = 0; j < 4; j++) {
			Vertex* temp_v = temp_q->verts[j];
			//glColor3f(0, 0, 0);
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

/******************************************************************************
Diaplay all vertices for selection
******************************************************************************/

void display_vertices(GLenum mode, Polyhedron* this_poly)
{

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	CHECK_GL_ERROR();

	for (int i = 0; i < this_poly->nverts; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		CHECK_GL_ERROR();

		Vertex* temp_v = this_poly->vlist[i];

		CHECK_GL_ERROR();

		{
			GLUquadric* quad = gluNewQuadric();
			CHECK_GL_ERROR();

			glPushMatrix();
			CHECK_GL_ERROR();
			glTranslatef(temp_v->x, temp_v->y, temp_v->z);
			glColor4f(0, 0, 1, 1.0);
			gluSphere(quad, this_poly->radius * 0.01, 50, 50);
			glPopMatrix();
			CHECK_GL_ERROR();

			gluDeleteQuadric(quad);
			CHECK_GL_ERROR();

		}
		CHECK_GL_ERROR();
	}
	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay selected quad
******************************************************************************/

void display_selected_quad(Polyhedron* this_poly)
{
	if (this_poly->selected_quad == -1)
	{
		return;
	}

	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	Quad* temp_q = this_poly->qlist[this_poly->selected_quad];

	glBegin(GL_POLYGON);
	for (j = 0; j < 4; j++) {
		Vertex* temp_v = temp_q->verts[j];
		glColor3f(1.0, 0.0, 1.0);
		glVertex3d(temp_v->x, temp_v->y, 0.001);
	}
	glEnd();
}

/******************************************************************************
Diaplay selected vertex
******************************************************************************/

void display_selected_vertex(Polyhedron* this_poly)
{
	if (this_poly->selected_vertex == -1)
	{
		return;
	}

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	CHECK_GL_ERROR();

	Vertex* temp_v = this_poly->vlist[this_poly->selected_vertex];

	GLfloat mat_diffuse[4];

	{
		mat_diffuse[0] = 1.0;
		mat_diffuse[1] = 0.0;
		mat_diffuse[2] = 0.0;
		mat_diffuse[3] = 1.0;
	}
	CHECK_GL_ERROR();

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	CHECK_GL_ERROR();

	{
		GLUquadric* quad = gluNewQuadric();
		glPushMatrix();
		CHECK_GL_ERROR();
		glTranslatef(temp_v->x, temp_v->y, temp_v->z);
		gluSphere(quad, this_poly->radius * 0.01, 50, 50);
		glPopMatrix();
		CHECK_GL_ERROR();
		gluDeleteQuadric(quad);
	}

	CHECK_GL_ERROR();
}

/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

	/* set escape key to exit */
	switch (key) {
	case 27:
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '1':
		display_mode = 1;
		glutPostRedisplay();
		break;

	case '2':
		display_mode = 2;
		glutPostRedisplay();
		break;

	case '3':
	{
		display_mode = 3;

		double L = (poly->radius * 2) / 30;
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			for (int j = 0; j < 4; j++) {

				Vertex* temp_v = temp_q->verts[j];

				temp_v->R = int(temp_v->x / L) % 2 == 0 ? 1 : 0;
				temp_v->G = int(temp_v->y / L) % 2 == 0 ? 1 : 0;
				temp_v->B = 0.0;
			}
		}
		glutPostRedisplay();
		break;
	}
	
	case '4': /*P.1.a grey scale map*/
		display_mode = 4;
		{
			//codes to modify the color of each vertex
			double M = poly->vlist[0]->scalar, m = poly->vlist[0]->scalar;
			for (int i = 0; i < poly->nverts; i++) {
				Vertex* temp_v = poly->vlist[i];
				double scalar = temp_v->scalar;
				//printf("scalar:%f\n", scalar);
				if (scalar > M) {
					M = scalar;
				}
				if (scalar < m) {
					m = scalar;
				}
			}
			/*get all scalar values of vertices, and get the maxmum M and minmum m*/
			//printf("M:%f\n", M);
			//printf("m:%f\n", m);
			for (int i = 0; i < poly->nquads; i++) {
				Quad* temp_q = poly->qlist[i];
				for (int j = 0; j < 4; j++) {
					Vertex* temp_v = temp_q->verts[j];
					double scalar = temp_v->scalar;
					temp_v->R = temp_v->G = temp_v->B =(scalar - m) / (M-m);
				}
			}
			/*based on scalar values, set R = G = B = (scalar - m) / (M-m)*/
		}
		display();
		break;

	case '5':/*P.1.b bi-color map*/
		display_mode = 4;
		{
			//codes to modify the color of each vertex
			double M = poly->vlist[0]->scalar, m = poly->vlist[0]->scalar;
			for (int i = 0; i < poly->nverts; i++) {
				Vertex* temp_v = poly->vlist[i];
				double scalar = temp_v->scalar;
				//printf("scalar:%f\n", scalar);
				if (scalar > M) {
					M = scalar;
				}
				if (scalar < m) {
					m = scalar;
				}
			}
			/*get all scalar values of vertices, and get the maxmum M and minmum m*/
			//printf("M:%f\n", M);
			//printf("m:%f\n", m);
			for (int i = 0; i < poly->nquads; i++) {
				Quad* temp_q = poly->qlist[i];
				for (int j = 0; j < 4; j++) {
					Vertex* temp_v = temp_q->verts[j];
					double scalar = temp_v->scalar;
					temp_v->R = 1*(scalar - m) / (M - m) + 0*(M - scalar) / (M - m);
					temp_v->G = 1*(scalar - m) / (M - m) + 0*(M - scalar) / (M - m);
					temp_v->B = 0*(scalar - m) / (M - m) + 1*(M - scalar) / (M - m);
				}
			}
			/*using yellow and blue as c1 and c2 to set RGB with scalar*/
		}
		display();
		break;

	case '6':/*P.1.c rainbow map*/
		display_mode = 4;
		{
			//codes to modify the color of each vertex
			double M = poly->vlist[0]->scalar, m = poly->vlist[0]->scalar;
			for (int i = 0; i < poly->nverts; i++) {
				Vertex* temp_v = poly->vlist[i];
				double scalar = temp_v->scalar;
				//printf("scalar:%f\n", scalar);
				if (scalar > M) {
					M = scalar;
				}
				if (scalar < m) {
					m = scalar;
				}
			}
			/*get all scalar values of vertices, and get the maxmum M and minmum m*/
			//printf("M:%f\n", M);
			//printf("m:%f\n", m);
			for (int i = 0; i < poly->nquads; i++) {
				Quad* temp_q = poly->qlist[i];
				for (int j = 0; j < 4; j++) {
					Vertex* temp_v = temp_q->verts[j];
					double scalar = temp_v->scalar;
					double c = (scalar - m) / (M - m);
					if (c < 0.25) {
						temp_v->R = 0.f;
						temp_v->G = c * 4.f;
						temp_v->B = 1.f;
					}
					else if (c < 0.5) {
						temp_v->R = 0.f;
						temp_v->G = 1.f;
						temp_v->B = (0.5 - c) * 4.f;
					}
					else if (c < 0.75) {
						temp_v->R = (c - 0.5) * 4.f;
						temp_v->G = 1.f;
						temp_v->B = 0.f;
					}
					else{
						temp_v->R = 1.f;
						temp_v->G = (1 - c) * 4.f;
						temp_v->B = 0.f;
						//printf("R:%f,G:%f, B:%f\n", temp_v->R, temp_v->G, temp_v->B);
					}
					
				}
			}
			/*based on scalar values, for color between blue and cyan, increases G with scalar.
			for color between cyan and green, decreases B with scalar. for color between green and yellow
			, increases R with scalar. for color between yellow and red, decreases G with scalar*/
		}
		display();
		break;

	case '7':/*P.2 height field*/
		display_mode = 5;
		{
			//codes to modify the color of each vertex
			double M = poly->vlist[0]->scalar, m = poly->vlist[0]->scalar;
			for (int i = 0; i < poly->nverts; i++) {
				Vertex* temp_v = poly->vlist[i];
				double scalar = temp_v->scalar;
				//printf("scalar:%f\n", scalar);
				if (scalar > M) {
					M = scalar;
				}
				if (scalar < m) {
					m = scalar;
				}
			}
			/*get all scalar values of vertices, and get the maxmum M and minmum m*/
			//printf("M:%f\n", M);
			//printf("m:%f\n", m);
			for (int i = 0; i < poly->nquads; i++) {
				Quad* temp_q = poly->qlist[i];
				for (int j = 0; j < 4; j++) {
					Vertex* temp_v = temp_q->verts[j];
					double scalar = temp_v->scalar;
					temp_v->z = 10 *(scalar - m) / (M - m);
				}
			}
			/*height changes with scalar.
			  multiple by 10 to make curve clearly*/
		}
		display();
		break;

	case '8':/*P.3 combine color and height*/
		display_mode = 5;
		{
			//codes to modify the color of each vertex
			double M = poly->vlist[0]->scalar, m = poly->vlist[0]->scalar;
			for (int i = 0; i < poly->nverts; i++) {
				Vertex* temp_v = poly->vlist[i];
				double scalar = temp_v->scalar;
				//printf("scalar:%f\n", scalar);
				if (scalar > M) {
					M = scalar;
				}
				if (scalar < m) {
					m = scalar;
				}
			}
			/*get all scalar values of vertices, and get the maxmum M and minmum m*/
			//printf("M:%f\n", M);
			//printf("m:%f\n", m);
			for (int i = 0; i < poly->nquads; i++) {
				Quad* temp_q = poly->qlist[i];
				for (int j = 0; j < 4; j++) {
					Vertex* temp_v = temp_q->verts[j];
					double scalar = temp_v->scalar;
					double c = (scalar - m) / (M - m);
					if (c < 0.25) {
						temp_v->R = 0.f;
						temp_v->G = c * 4.f;
						temp_v->B = 1.f;
						temp_v->z = 10 * (scalar - m) / (M - m);
					}
					else if (c < 0.5) {
						temp_v->R = 0.f;
						temp_v->G = 1.f;
						temp_v->B = (0.5 - c) * 4.f;
						temp_v->z = 10 * (scalar - m) / (M - m);
					}
					else if (c < 0.75) {
						temp_v->R = (c - 0.5) * 4.f;
						temp_v->G = 1.f;
						temp_v->B = 0.f;
						temp_v->z = 10 * (scalar - m) / (M - m);
					}
					else {
						temp_v->R = 1.f;
						temp_v->G = (1 - c) * 4.f;
						temp_v->B = 0.f;
						temp_v->z = 10 * (scalar - m) / (M - m);
						//printf("R:%f,G:%f, B:%f\n", temp_v->R, temp_v->G, temp_v->B);
					}
				}
			}
			/*set the color as wheel between yellow and red.
			  height changes with scalar.
			  multiple by 10 to make curve clearly*/
		}
		display();
		break;
	
	case 'r':
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		glutPostRedisplay();
		break;
	}
}

/******************************************************************************
Callback function for dragging mouse
******************************************************************************/

void motion(int x, int y) {
	float r[4];
	float s, t;

	s = (2.0 * x - win_width) / win_width;
	t = (2.0 * (win_height - y) - win_height) / win_height;

	if ((s == s_old) && (t == t_old))
		return;

	switch (mouse_mode) {
	case 2:

		Quaternion rvec;

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;

	case 1:

		translation[0] += (s - s_old);
		translation[1] += (t - t_old);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

/******************************************************************************
Callback function for mouse clicks
******************************************************************************/

void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		
		if (state == GLUT_DOWN) {
			float xsize = (float)win_width;
			float ysize = (float)win_height;

			float s = (2.0 * x - win_width) / win_width;
			float t = (2.0 * (win_height - y) - win_height) / win_height;

			s_old = s;
			t_old = t;

			/*translate*/
			if (button == GLUT_LEFT_BUTTON)
			{
				mouse_mode = 1;
			}

			/*rotate*/
			if (button == GLUT_RIGHT_BUTTON)
			{
				mouse_mode = 2;
			}
		}
		else if (state == GLUT_UP) {

			if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_SHIFT) {  // build up the selection feedback mode

				/*select face*/

				GLuint selectBuf[win_width];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_quads(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_quad = processHits(hits, selectBuf);
				printf("Selected quad id = %d\n", poly->selected_quad);
				glutPostRedisplay();

				CHECK_GL_ERROR();

			}
			else if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_CTRL)
			{
				/*select vertex*/

				GLuint selectBuf[win_width];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*  create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_vertices(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_vertex = processHits(hits, selectBuf);
				printf("Selected vert id = %d\n", poly->selected_vertex);
				glutPostRedisplay();

			}

			mouse_mode = -1;
		}
	}
}

/******************************************************************************
Callback function for mouse wheel scroll
******************************************************************************/

void mousewheel(int wheel, int direction, int x, int y) {
	if (direction == 1) {
		zoom *= zoomspeed;
		glutPostRedisplay();
	}
	else if (direction == -1) {
		zoom /= zoomspeed;
		glutPostRedisplay();
	}
}

/******************************************************************************
Callback function for scene display
******************************************************************************/

void display(void)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting

	glClearDepth(1.0);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();

	set_view(GL_RENDER);
	set_scene(GL_RENDER, poly);

	/*display the mesh*/
	display_polyhedron(poly);

	/*display selected elements*/
	display_selected_vertex(poly);
	display_selected_quad(poly);


	glFlush();
	glutSwapBuffers();
	glFinish();

	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay the polygon with visualization results
******************************************************************************/

void display_polyhedron(Polyhedron* poly)
{
	unsigned int i, j;

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);

	for (i = 0; i < poly->nquads; i++) {

		Quad* temp_q = poly->qlist[i];

		switch (display_mode) {
		case 1:
		{
			glEnable(GL_LIGHTING);
			glEnable(GL_LIGHT0);
			glEnable(GL_LIGHT1);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			GLfloat mat_diffuse[4] = { 1.0, 1.0, 0.0, 0.0 };
			GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
			glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		}
		case 2:
			glDisable(GL_LIGHTING);
			glEnable(GL_LINE_SMOOTH);
			glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
				glColor3f(0.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();

			glDisable(GL_BLEND);
			break;

		case 3:
			glDisable(GL_LIGHTING);

			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 4:
			glDisable(GL_LIGHTING);
			glShadeModel(GL_SMOOTH);

			for (j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				temp_v->z = 0;
			}/*clear the height values for P.1*/
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 5:
			glDisable(GL_LIGHTING);
			glShadeModel(GL_SMOOTH);

			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		}

	}
}
