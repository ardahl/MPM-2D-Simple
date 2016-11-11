//Uncomment this to remove debug file output
/// #define NDEBUG

#include "mpm.hpp"
#include "opengl.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>
//OpenCV Includes for saving frames
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/videoio/videoio.hpp>

int w = 800, h = 800;                           // size of window in pixels
double xmin = 0, xmax = 1, ymin = 0, ymax = 1; // range of coordinates drawn
int lastTime = 0, prevTime = 0, frame = 0;
int seconds = 3*30, curr = 0;
bool next = true;

cv::Mat img(h, w, CV_8UC3);
cv::Mat flipped(h, w, CV_8UC3);
cv::VideoWriter output;

#ifndef NDEBUG
std::ofstream debug;
#endif
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
        case 'n':
            next = true;
            break;
	}
}

void display() {
    /// if(!next) {
        /// return;
    /// }
    /// else {
        /// next = false;
    /// }
    if(curr >= seconds) {
        printf("\nRender Done\n");
        #ifndef NDEBUG
        debug.close();
        #endif
        exit(0);
    }
    //Perform step
    int iters = 1000;
    for(int i = 0; i < iters; i++) {    //Hardcode for 30fps with dt of (1/3)e-4 
        m->step((1.0/3.0)*0.0001);
        printf("Step: %d/%d\r", i, iters);
    }
    printf("\n");
    #ifndef NDEBUG
    //output gradient for debug
    for(int i = 0; i < m->particles.size(); i++) {
        Particle* p = m->particles[i];
        debug << "particle " << i << "\n";
        debug << p->gradient << "\n";
    }
    debug << "\n\n";
    #endif
    curr++;
    
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
    //Draw grid structure
    glColor3f(0.8, 0.8, 0.8);
    glPointSize(10);
    for (int i = 0; i < m->m; i++) {
        for (int j = 0; j < m->n; j++) {
            Vector2d x = m->x0 + Vector2d(i+0.5, j+0.5)*m->h;
            glBegin(GL_POINTS);
            glVertex2f(x(0), x(1));
            glEnd();
        }
    }
    //just draw the particles 
    std::vector<Particle*> ps = m->particles;
    glPointSize(5);
    glColor3f(0,0,0);
    glBegin(GL_POINTS);
    for(int i = 0; i < ps.size(); i++) {
        Particle* p = ps[i];
        glColor3f(p->color(0), p->color(1), p->color(2));
        glVertex2f(p->x(0), p->x(1));
    }
    glEnd();
    
    glReadPixels(0, 0, img.cols, img.rows, GL_BGR, GL_UNSIGNED_BYTE, img.data);
    cv::flip(img, flipped, 0);
    output << flipped;
    
    glutSwapBuffers();
    printf("Frame %d/%d\n", curr, seconds);
    /// printf("Progress: %.2f%%\r", 100.0*(double)curr/seconds);
}

void idle() {
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    debug.open("grad.txt");
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_MULTISAMPLE);
    glutInitWindowSize(::w, ::h);
    glutCreateWindow("Animation");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
    glutKeyboardFunc(onKey);
    
    ///initialize
    ///while not done
    ///    for i=0 : fps / timeStep
    ///        step
    ///    render opengl frame
    ///    save frame using opencv
    
    std::string config = std::string(argv[1]);
    
    //use fast 4-byte alignment (default anyway) if possible
    glPixelStorei(GL_PACK_ALIGNMENT, (img.step & 3) ? 1 : 4);
    //set length of one complete row in destination data (doesn't need to equal img.cols)
    glPixelStorei(GL_PACK_ROW_LENGTH, img.step/img.elemSize());
    output = cv::VideoWriter("test.avi", CV_FOURCC('P', 'I', 'M', '1'), 30, cv::Size(w, h));
    
    m = new Material(config);
    m->init();
    
    glutMainLoop();
    return 0;
}
