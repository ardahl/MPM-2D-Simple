/// #include "mpm.hpp"
#include "material.hpp"
#include "opengl.hpp"
#include <iostream>
#include <cstdio>
//OpenCV Includes for saving frames
#include <cv.h>
#include <highgui.h>

int w = 800, h = 800;                           // size of window in pixels
double xmin = 0, xmax = 1, ymin = 0, ymax = 1; // range of coordinates drawn
int lastTime = 0, prevTime = 0, frame = 0;
double seconds = 3, curr = 0;

/// cv::Mat img(h, w, CV_8UC3);
/// cv::Mat flipped(h, w, CV_8UC3);
/// cv::VideoWriter output;

Material* m;

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
    if(curr >= seconds) {
        printf("\nRender Done\n");
        exit(0);
    }
    for(int i = 0; i < 10000; i++) {    //Hardcode for 30fps with dt of (1/3)e-5 
        m->step((1.0/3.0)*0.00001);
    }
    curr += (1.0/30.0);
    
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
     * Draw the actual stuff
     */
    //just draw the particles 
    std::vector<Particle*> ps = m->particles;
    glPointSize(5);
    glColor3f(0,0,0);
    glBegin(GL_POINTS);
    for(int i = 0; i < ps.size(); i++) {
        Particle* p = ps[i];
        glColor3f(p->color(0), p->color(1), p->color(2));
        glVertex3f(p->x(0), p->x(1), p->x(2));
    }
    glEnd();
    
    /// glReadPixels(0, 0, img.cols, img.rows, GL_BGR, GL_UNSIGNED_BYTE, img.data);
    /// cv::flip(img, flipped, 0);
    /// output << flipped;
    
    glutSwapBuffers();
    
    printf("Progress: %f%%\r", curr/seconds);
}

/// void idle() {
    /// //FPS
    /// frame++;
    /// int time = glutGet(GLUT_ELAPSED_TIME);
    /// if(time - lastTime > 1000) {
        /// double fps = frame * 1000.0 / (time-lastTime);
        /// lastTime = time;
        /// frame = 0;
        /// //printf("FPS: %f\n", fps);
        /// char fpsString[15];
        /// snprintf(fpsString, 15, "%f", fps);
        /// glutSetWindowTitle(fpsString);
    /// }
    /// ///double timeStep = (time-prevTime) / 1000.0;
    /// double timeStep = 1e-5;                     //time in seconds
    /// /*
     /// * Do a time step
     /// */
    /// for(int i = 0; i < 30.0/timeStep; i++) {    //render 30 frames for second
        
    /// }
    /// prevTime = time;

    /// glutPostRedisplay();
/// }

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_MULTISAMPLE);
    glutInitWindowSize(::w, ::h);
    glutCreateWindow("Animation");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    /// glutIdleFunc(idle);
    glutKeyboardFunc(onKey);
    
    ///initialize
    ///while not done
    ///    for i=0 : fps / timeStep
    ///        step
    ///    render opengl frame
    ///    save frame using opencv
    
    std::string config = std::string(argv[1]);
    /// int fps = 30;
    /// double dt = (1.0/3.0)e-5;
    /// int itersPerFrame = (int)(1.0/(fps*dt));
    
    
    
    /// //use fast 4-byte alignment (default anyway) if possible
    /// glPixelStorei(GL_PACK_ALIGNMENT, (img.step & 3) ? 1 : 4);
    /// //set length of one complete row in destination data (doesn't need to equal img.cols)
    /// glPixelStorei(GL_PACK_ROW_LENGTH, img.step/img.elemSize());
    /// output = cv::VideoWriter("test.avi", CV_FOURCC('P', 'I', 'M', '1'), 30, cv::Size(w, h));
    
    m = new Material(config);
    m->init();
    
    glutMainLoop();
    return 0;
}
