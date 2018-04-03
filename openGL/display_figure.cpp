#include <GL/glut.h>
#include<iostream>
void renderScene(void) {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glBegin(GL_TRIANGLES);
		glVertex3f(0,1,-7);
		glVertex3f(1,1,-7);
		glVertex3f(.5,2,-7);
	glEnd();

	glBegin(GL_TRIANGLES);
		glVertex3f(0,0,-7);
		glVertex3f(0,1,-7);
		glVertex3f(-1,0.5,-7);
	glEnd();

	glBegin(GL_TRIANGLES);
		glVertex3f(0,0,-7);
		glVertex3f(1,0,-7);
		glVertex3f(.5,-1,-7);
	glEnd();

	glBegin(GL_TRIANGLES);
		glVertex3f(1,1,-7);
		glVertex3f(1,0,-7);
		glVertex3f(2,.5,-7);
	glEnd();


	glBegin(GL_QUADS);
		glVertex3f(0,0,-7);
		glVertex3f(1,0.0,-7);
		glVertex3f(1,1,-7);
		glVertex3f(0,1,-7);
	glEnd();

	glutSwapBuffers();


	glutSwapBuffers();
}
void changeSize(int w, int h) {
	if(h == 0)
		{h = 1;}


	float ratio = 1.0* w / h;

	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);

        // Reset Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective.
	gluPerspective(45,ratio,1,1000);

	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
}


int main(int argc, char **argv) {

	// init GLUT and create Window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA); //depth buffer, smooth, colored
	glutInitWindowPosition(100,100);
	glutInitWindowSize(320,320);
	glutCreateWindow("The 1976 Annual World's Best SF - C++, reinvention is its own reward. ");
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);	
	glutMainLoop();
	return 1;

}
