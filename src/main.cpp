#include "defines.hpp"
#include "mpm.hpp"
#include "opengl.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>
//OpenCV Includes for saving frames
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/videoio/videoio.hpp>

using namespace Eigen;

int w = 800, h = 800;                           // size of window in pixels
double xmin = 0, xmax = 1, ymin = 0, ymax = 1; // range of coordinates drawn
int lastTime = 0, prevTime = 0, frame = 0;
int seconds = 7*30, curr = 0;
bool next = true;

cv::Mat img(h, w, CV_8UC3);
cv::Mat flipped(h, w, CV_8UC3);
cv::VideoWriter output;

#ifndef NDEBUG
extern std::ofstream debug;
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
    curr++;
    if(curr > seconds) {
        printf("\nRender Done\n");
        #ifndef NDEBUG
        debug.close();
        #endif
        exit(0);
    }
    //Perform step
    int iters = 2000;
    double itersInv = 1.0/iters;
    for(int i = 0; i < iters; i++) {    //Hardcode for 30fps with dt of (1/3)e-5
        /// printf("Step %d\n", i);
        m->step((1.0/30.0)*itersInv);
        printf("Frame %d/%d Step: %d/%d\r", curr, seconds, i+1, iters);
    }
    printf("\n");
    
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
    Vector2d x0 = m->x0 - Vector2d(m->h/2.0, m->h/2.0);
    Vector2d x1 = m->x1 - Vector2d(m->h/2.0, m->h/2.0);
    Vector2d x10 = x1 - x0;
    glColor3f(0.8, 0.8, 0.8);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < m->m+1; i++) {
        for (int j = 0; j < m->n+1; j++) {
            Vector2d gridPos = Vector2d(i, j);
            if(i % 2 == 1) {
                gridPos = Vector2d(i, m->n-j);
            }
            Vector2d x = x0 + gridPos*m->h;
            x = x - x0;
            glVertex2f(x(0)/x10(0), x(1)/x1(1));
        }
    }
    glEnd();
    glBegin(GL_LINE_STRIP);
    for (int j = 0; j < m->n+1; j++) {
        for (int i = 0; i < m->m+1; i++) {
            Vector2d gridPos = Vector2d(i, j);
            if(j % 2 == 1) {
                gridPos = Vector2d(m->m-i, j);
            }
            Vector2d x = x0 + gridPos*m->h;
            x = x - x0;
            glVertex2f(x(0)/x10(0), x(1)/x10(1));
        }
    }
    glEnd();
    //just draw the particles 
    std::vector<Particle*> ps = m->particles;
    glPointSize(5);
    glColor3f(0,0,0);
    glBegin(GL_POINTS);
    for(size_t i = 0; i < ps.size(); i++) {
        Particle* p = ps[i];
        glColor3f(p->color(0), p->color(1), p->color(2));
        glVertex2f((p->x(0)-x0(0))/x10(0), (p->x(1)-x0(1))/x10(1));
    }
    glEnd();
    
    glReadPixels(0, 0, img.cols, img.rows, GL_BGR, GL_UNSIGNED_BYTE, img.data);
    cv::flip(img, flipped, 0);
    output << flipped;
    
    glutSwapBuffers();
    /// printf("Frame %d/%d\n", curr, seconds);
}

void idle() {
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    if(argc < 3) {
        std::cout << "Usage: ./mpm <configfile> <debug and video outputs name>\n";
        std::exit(0);
    }
    std::string outfile = std::string(argv[2]);
    #ifndef NDEBUG
    std::string dbout = outfile + std::string(".txt");
    debug.open(dbout.c_str());
    #endif
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
    std::string videoout = outfile + std::string(".avi");
    output = cv::VideoWriter(videoout, CV_FOURCC('P', 'I', 'M', '1'), 30, cv::Size(w, h));
    
    m = new Material(config);
    m->init();
    
    glutMainLoop();
    return 0;
}
