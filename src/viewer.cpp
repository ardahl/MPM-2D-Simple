#include "defines.hpp"
#include "mpm.hpp"
#include "opengl.hpp"
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <iomanip>
//OpenCV Includes for saving frames
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/videoio/videoio.hpp>

using namespace Eigen;

// Global Variables
int width = 800, height = 800;                           // size of window in pixels
double xmin = 0, xmax = 1, ymin = 0, ymax = 1; // range of coordinates drawn
int frame = 0;
bool next = true;

World *world;
std::vector<std::vector<std::vector<Particle> > > frames;

cv::Mat img(height, width, CV_8UC3);
cv::Mat flipped(height, width, CV_8UC3);
cv::VideoWriter output;

void reshape(int w, int h) {
    ::width = w;
    ::height = h;
    glViewport(0, 0, width, height);
}

//TODO: Make this better. Currently requires being run in the directory of the files.
//Also requires the configuration file to be there as well
//First, remove world. We don't need to set up the world.
//Second, Update this so we can run the viewer from anywhere
bool readAnimation(const char *fname, std::vector<std::vector<std::vector<Particle> > > &frames) {
	for (int frame = 0; true; frame++) { 
	  std::vector<std::vector<Particle> > objects;
	  for (int obj = 0; true; obj++) {
   	    std::vector<Particle> parts;
   	    std::ostringstream ss;
		ss << std::setw(2) << std::setfill('0') << obj << "." << std::setw(6)<< frame;
		std::string pframe(ss.str());
		std::string parIn = std::string(fname) + "-" + pframe + ".bgeo";
        if(!std::ifstream(parIn)) break;
		if (!readParticles(parIn.c_str(), parts)) break;
		if (parts.size() == 0) break;
		objects.push_back(parts);
	  }
	  if (objects.size() == 0) break;
	  frames.push_back(objects);
	}
    
    return true;
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
     * Draw the actual stuff
     */
    //Draw grid structure
    Vector2d x0 = world->origin - Vector2d(world->h/2.0, world->h/2.0);
    Vector2d x1 = x0 + Vector2d(world->res[0]*world->h, world->res[1]*world->h);
    Vector2d x10 = x1 - x0;
    glColor3f(0.8, 0.8, 0.8);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < world->res[0]+1; i++) {
        for (int j = 0; j < world->res[1]+1; j++) {
            Vector2d gridPos = Vector2d(i, j);
            if(i % 2 == 1) {
                gridPos = Vector2d(i, world->res[1]-j);
            }
            Vector2d x = x0 + gridPos*world->h;
            x = x - x0;
            glVertex2f(x(0)/x10(0), x(1)/x10(1));
        }
    }
    glEnd();
    glBegin(GL_LINE_STRIP);
    for (int j = 0; j < world->res[1]+1; j++) {
        for (int i = 0; i < world->res[0]+1; i++) {
            Vector2d gridPos = Vector2d(i, j);
            if(j % 2 == 1) {
                gridPos = Vector2d(world->res[0]-i, j);
            }
            Vector2d x = x0 + gridPos*world->h;
            x = x - x0;
            glVertex2f(x(0)/x10(0), x(1)/x10(1));
        }
    }
    glEnd();
    //just draw the particles 
    glPointSize(5);
    glColor3f(0,0,0);
    glBegin(GL_POINTS);
	
    for(size_t i = 0; i < frames[frame].size(); i++) {
	  for(size_t j = 0; j < frames[frame][i].size(); j++) {
        Particle &p = frames[frame][i][j];
        glColor3f(p.color(0), p.color(1), p.color(2));
        glVertex2f((p.x(0)-x0(0))/x10(0), (p.x(1)-x0(1))/x10(1));
	  }
	}
    glEnd();
    glReadPixels(0, 0, img.cols, img.rows, GL_BGR, GL_UNSIGNED_BYTE, img.data);
    cv::flip(img, flipped, 0);
    output << flipped;
    
    glutSwapBuffers();
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

void myTimerFunc(int id) {
    frame++;
    if(frame > (int)frames.size()-1) {
        frame = 0;
    }
    glutPostRedisplay();
    glutTimerFunc(33, myTimerFunc, 0);
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
  
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_MULTISAMPLE);
    glutInitWindowSize(::width, ::height);
    glutCreateWindow("Animation");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(onKey);
    glutTimerFunc(1000, myTimerFunc, 0);
  
    //use fast 4-byte alignment (default anyway) if possible
    glPixelStorei(GL_PACK_ALIGNMENT, (img.step & 3) ? 1 : 4);
    //set length of one complete row in destination data (doesn't need to equal img.cols)
    glPixelStorei(GL_PACK_ROW_LENGTH, img.step/img.elemSize());

	std::string outfile, conffile;
    if(argc == 1 || argc > 3) {
        printf("Usage: ./viewer path/to/bgeo/<base_name> [/path/to/config]");
        std::exit(0);
    }
    else if(argc == 2){
        outfile = std::string(argv[1]);
        conffile = outfile + std::string(".json");
    }
    else {
        outfile = std::string(argv[1]);
        conffile = std::string(argv[2]);
    }

    /// std::string videoout = outfile.substr(outfile.find_last_of('/')+1,std::string::npos) + std::string(".avi");
    std::string videoout = outfile + std::string(".mp4");
    //Old AVI format
    /// output = cv::VideoWriter(videoout, CV_FOURCC('P', 'I', 'M', '1'), 30, cv::Size(width, height)); 
    //Theoretically this should be the H.264 codec with mp4 files, but ffmpeg uses different tags.
    //From docs: http://docs.opencv.org/3.2.0/dd/d9e/classcv_1_1VideoWriter.html#ac3478f6257454209fa99249cc03a5c59
        //FFMPEG backend with MP4 container natively uses other values as fourcc code: see ObjectType, so you may receive a warning message from OpenCV about fourcc code conversion. 
    /// output = cv::VideoWriter(videoout, CV_FOURCC('H', '2', '6', '4'), 30, cv::Size(width, height));
    //0x20 is the code for the mp4v codec (MPEG-4 Video)
    //0x21 should be the code for the H.264 codec
    output = cv::VideoWriter(videoout, 0x21, 30, cv::Size(width, height));
    if(!output.isOpened()) {
        printf("Video Failed to open\n");
        std::exit(1);
    }
    
    std::string config = std::string(conffile);
    world = new World(config);
    world->init();
    
    readAnimation(outfile.c_str(), frames);
    printf("Frames: %d\n", (int)frames.size());
    if(frames.size() == 0) { 
        std::exit(0);
    }
    glutMainLoop();
  
    return 0;
}

