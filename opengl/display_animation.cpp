#include <GL/glut.h>
#include<iostream>
float angle = 0;
float z = 0;
int x = 10;
float y = 0;
void renderScene(void) {

	// Clear Color and Depth 

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();


	gluLookAt(	0.0f, 0.0f, 10.0f,
			0.0f, 0.0f,  0.0f,
			0.0f, 1.0f,  0.0f);



    glRotatef(angle, 1.0f, 1.0f, 1.0f);
	//glTranslatef(0,0,angle*0.1);
	//glTranslatef(-2,0,0);
	/*
	glBegin(GL_QUADS); // of the color cube

	// Top-face
	glColor3f(0.0f, 1.0f, 0.0f); // green
	glVertex3f(1.0f, 1.0f, -1.0f);
	glVertex3f(-1.0f, 1.0f, -1.0f);
	glVertex3f(-1.0f, 1.0f, 1.0f);
	glVertex3f(1.0f, 1.0f, 1.0f);

	// Bottom-face
	glColor3f(1.0f, 0.5f, 0.0f); // orange
	glVertex3f(1.0f, -1.0f, 1.0f);
	glVertex3f(-1.0f, -1.0f, 1.0f);
	glVertex3f(-1.0f, -1.0f, -1.0f);
	glVertex3f(1.0f, -1.0f, -1.0f);

	// Front-face
	glColor3f(1.0f, 0.0f, 0.0f); // red
	glVertex3f(1.0f, 1.0f, 1.0f);
	glVertex3f(-1.0f, 1.0f, 1.0f);
	glVertex3f(-1.0f, -1.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, 1.0f);

	// Back-face
	glColor3f(1.0f, 1.0f, 0.0f); // yellow
	glVertex3f(1.0f, -1.0f, -1.0f);
	glVertex3f(-1.0f, -1.0f, -1.0f);
	glVertex3f(-1.0f, 1.0f, -1.0f);
	glVertex3f(1.0f, 1.0f, -1.0f);

	// Left-face
	glColor3f(0.0f, 0.0f, 1.0f); // blue
	glVertex3f(-1.0f, 1.0f, 1.0f);
	glVertex3f(-1.0f, 1.0f, -1.0f);
	glVertex3f(-1.0f, -1.0f, -1.0f);
	glVertex3f(-1.0f, -1.0f, 1.0f);

	// Right-face
	glColor3f(1.0f, 0.0f, 1.0f); // magenta
	glVertex3f(1.0f, 1.0f, -1.0f);
	glVertex3f(1.0f, 1.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, -1.0f);

	glEnd(); // of the color cube
	*/
	glBegin(GL_TRIANGLES);
		glColor3f(0.6f, int(angle*10)%x/100., int(angle*10)%x/10.); // green
		glVertex3f(0,1,z);
		glVertex3f(1,1,z);
		glVertex3f(.5,2,z);
	glEnd();

	glBegin(GL_TRIANGLES);
		glColor3f(0.0f, 1.0f, int(angle*10)%x/10.); // green
		glVertex3f(0,0,z);
		glVertex3f(0,1,z);
		glVertex3f(-1,0.5,z);
	glEnd();

	glBegin(GL_TRIANGLES);
		glColor3f(1.0f, int(angle*5)%x/100., int(angle*10)%x/10.); // green
		glVertex3f(0,0,z);
		glVertex3f(1,0,z);
		glVertex3f(.5,-1,z);
	glEnd();

	glBegin(GL_TRIANGLES);
		glColor3f(int(angle*20)%x/100., int(angle*10)%x/10., 1.0f); // green
		glVertex3f(1,1,z);
		glVertex3f(1,0,z);
		glVertex3f(2,.5,z);
	glEnd();
	
    //firkant
	glBegin(GL_QUADS);
		glColor3f(int(angle*20)%x/100., int(angle*10)%x/100., 0.0f); // green
		glVertex3f(0,0,z);
		glVertex3f(1,0.0,z);
		glVertex3f(1,1,z);
		glVertex3f(0,1,z);
	glEnd();
	
	angle+=0.8f;
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
	gluPerspective(45.0f,ratio,1.0f,100.0f);

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
	glutDisplayFunc(renderScene); //display
	glutReshapeFunc(changeSize);  //perspective / zoom
    glutIdleFunc(renderScene);	  //animation
	glutMainLoop();
	return 1;

}
