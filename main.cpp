#include "mpm.hpp"
#include "opengl.hpp"
#include <iostream>
#include <cstdio>

nt w = 800, h = 800;                           // size of window in pixels
double xmin = 0, xmax = 1, ymin = 0, ymax = 1; // range of coordinates drawn
double dt = 0.01;                              // time step
int lastTime = 0, prevTime = 0, frame = 0;

void reshape(int w, int h) {
    ::w = w;
    ::h = h;
    glViewport(0, 0, w, h);
}

void onKey(unsigned char key, int x, int y) {
	switch(key) {
		case 27: // escape
			exit(0);
            break;
	}
}

void display() {
    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(::xmin, ::xmax, ::ymin, ::ymax);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    /*
     * Do material drawing stuff here
     */
    glutSwapBuffers();
}

void idle() {
    //FPS
    frame++;
    int time = glutGet(GLUT_ELAPSED_TIME);
    if(time - lastTime > 1000) {
        double fps = frame * 1000.0 / (time-lastTime);
        lastTime = time;
        frame = 0;
        //printf("FPS: %f\n", fps);
        char fpsString[15];
        snprintf(fpsString, 15, "%f", fps);
        glutSetWindowTitle(fpsString);
    }
    double timeStep = (time-prevTime) / 1000.0;
    /*
     * Do a time step
     */
    prevTime = time;

    glutPostRedisplay();
}

int main(int argc, char** argv) {
    
    /******************************
     * Algorithm Process from Stomakhin et al. 2013
     * 
     * while not ended
     *      Rasterize particle data to the grid
     *      Compute Grid Forces
     *      Update grid velocities
     *      Solve semi-implicit integration system
     *      Update Deformation Gradient
     *      Update particle velocities
     *      Update particle positions
     * 
     * 
     * 
     * Rasterize particle data to the grid
     *      Transfer mass       m^n_i = Ep m_p * w^n_ip
     *      Transfer velocities v^n_i = Ep (v^n_p * m_p * w^n_ip) / m^n_i
     * 
     * Compute grid forces
     *      
     * 
     * Update grid velocities
     *      Update to v^*_i = v^n_i + dt * m^-1_i * f^n_i
     * 
     * Solve semi-implicit integration system
     * 
     * 
     * Update deformation gradient
     * 
     * 
     * Update particle velocities
     *      v^n+1_p = (1-a)*v_pic + a*v_flip
     * 
     * Update particle positions
     *      x^n+1_p = x^n_p + dt * v^n+1_p
     * 
     ******************************/
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_MULTISAMPLE);
    glutInitWindowSize(::w, ::h);
    glutCreateWindow("Animation");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
    /// glutMouseFunc(mouse);
    /// glutMotionFunc(motion);
    glutKeyboardFunc(onKey);
    glutMainLoop();
    
    return 0;
}
