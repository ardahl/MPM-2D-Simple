// ./bin/sampler <output filename>

#include "world.hpp"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <random>
#include <ctime>
#include <Partio.h>
#include <Eigen/Geometry>
#include "json/json.h"
#include "range.hpp"

#define RAND false

#define CIRCLE true //could change to SQUARE, LINE


using namespace Eigen;
using benlib::range;
enum VelType {
    constant,
    linear,
    rotation,
    none
};

int main(int argc, char* argv[]) {
    std::vector<Particle> parts;
    std::srand(std::time(0));

    // parameters of the object
    double size[2] = {0.15, 0.15};  //radius of object in x, y directions
    Vector2d object(0.0,0.5);       //position of the object center
    int ores[2] = {75,75};          //particle resolution for sampling
    double pmass = 0.001;             //mass of each particle
    double scale = 0.75;            //velocity scale
    #if RAND
    int numPart = 648;              //number of random particles in object
    #endif
    VelType vt = rotation;          //How to initialize velocity

    //Jitter
    bool jitter = false;            //Whether to jitter the particles
    Vector2d cellSize(0.025, 0.025);//Size of each grid cell
    double jitterFract = 0.25;      //How much of the cell size it should be jitter by at most

    //Constant veloctiy
    double xvel = 1.0;              //x and y values for constant veloctiy
    double yvel = 0.0;
    //Linear veloctiy v=(x,y)
    //Rotational velocity v=(-y,x)
    bool centered = true;           //whether linear and rotational velocity uses the center for (0,0)

    //Initialize B
    Matrix2d B, C;
    Matrix2d D = ((cellSize(0)*cellSize(1))/3.0)*Matrix2d::Identity();
    if(vt == constant) {
        C << 0, 0, 0, 0;
    }
    else if(vt == linear) {
        C << 1, 0, 0, 1;
    }
    else if(vt == rotation) {
        C << 0, -1, 1, 0;
    }
    else {
        C << 0, 0, 0, 0;
    }
    C = scale*C;
    B = C * D;

#if SQUARE
    //Vector2d center = object + (Vector2d(size[0],size[1]) * 0.5);
    Vector2d center = object;
    object = object - Vector2d(size[0], size[1]);
    //Set up particles at each object vertex
    double diffx, diffy;
    if(ores[0] < 2) {
        diffx = size[0];
    }
    else {
        diffx = 2*size[1] / (ores[1]-1);
    }
    if(ores[1] < 2) {
        diffy = size[1];
    }
    else {
        diffy = 2*size[1] / (ores[1]-1);
    }

    for(int i = 0; i < ores[0]; i++) {
        for(int j = 0; j < ores[1]; j++) {
            Vector2d pos = object + Vector2d(diffx*i, diffy*j);
            if(jitter) {
                //Jitter position by quarter grid cell size
                Vector2d jsize = jitterFract * cellSize;
                pos = pos + Vector2d(jsize(0)*(2*((double)rand()/RAND_MAX)-1), jsize(1)*(2*((double)rand()/RAND_MAX)-1));
            }
            Vector2d quad = pos - center;
            Vector3d col(1.0, 0.0, 0.0);
            if((quad(0) < 0 && quad(1) < 0) || (quad(0) >= 0 && quad(1) >= 0)) {
                col = Vector3d(0.0, 1.0, 0.0);
            }

            Vector2d d = pos-center;
            Vector2d vel;
            if(vt == constant) {
                vel = Vector2d(xvel, yvel);
            }
            else if(vt == linear) {
                if(centered) {
                    vel = Vector2d(d(0), d(1));
                }
                else {
                    vel = Vector2d(pos(0), pos(1));
                }
            }
            else if(vt == rotation) {
                if(centered) {
                    vel = Vector2d(-d(1), d(0));
                }
                else {
                    vel = Vector2d(-pos(1), pos(0));
                }
            }

            Particle par(pos, scale*vel, col, pmass);
            par.B = B;
            parts.push_back(par);
        }
    }
#elif CIRCLE
    Vector2d center = object;
        //non-randomly make a circle
        //technically this is an ellipse because I'm using the same data as the
        //square and just using the size vector as radius of the semi-major and semi-minor axes
        //just make a square and reject those outside the ellipse.
        //~78.5% of the resx*resy particles are accepted - Pi/4 * (l*w) particles
    double diffx, diffy;
    if(ores[0] < 2) {
        diffx = size[0];
    }
    else {
        diffx = 2*size[1] / (ores[1]-1);
    }
    if(ores[1] < 2) {
        diffy = size[1];
    }
    else {
        diffy = 2*size[1] / (ores[1]-1);
    }

    #if RAND
    while(parts.size() < numPart) {
        Vector2d pos = Matrix<double, 2, 1>::Random();
        pos(0) = pos(0)*size[0];
        pos(1) = pos(1)*size[1];
        pos = pos + center;
    #else
    for(int i = 0; i < ores[0]; i++) {
      for(int j = 0; j < ores[1]; j++) {
        Vector2d pos = center - Vector2d(size[0], size[1]) + Vector2d(diffx*i, diffy*j);
    #endif
        if(jitter) {
            //Jitter position by quarter grid cell size
            Vector2d jsize = jitterFract * cellSize;
            pos = pos + Vector2d(jsize(0)*(2*((double)rand()/RAND_MAX)-1), jsize(1)*(2*((double)rand()/RAND_MAX)-1));
        }
        Vector3d col(1.0, 0.0, 0.0);
        Vector2d quad = pos - center;
        if((quad(0) < 0 && quad(1) < 0) || (quad(0) > 0 && quad(1) > 0)) {
            col = Vector3d(0.0, 1.0, 0.0);
        }
        Vector2d ph = pos - object;
        if( ((ph(0)*ph(0))/(size[0]*size[0])) + ((ph(1)*ph(1))/(size[1]*size[1])) < 1+EPS) {
            Vector2d vel;
            if(vt == none) {
                vel = Vector2d::Zero();
            }
            else if(vt == constant) {
                vel = Vector2d(xvel, yvel);
            }
            else if(vt == linear) {
                if(centered) {
                    vel = Vector2d(ph(0), ph(1));
                }
                else {
                    vel = Vector2d(pos(0), pos(1));
                }
            }
            else if(vt == rotation) {
                if(centered) {
                    vel = Vector2d(-ph(1), ph(0));
                }
                else {
                    vel = Vector2d(-pos(1), pos(0));
                }
            }
            Particle par(pos, scale*vel, col, pmass);
            par.B = B;
            parts.push_back(par);
        }
      }
    #if !RAND
    }
    #endif
#elif LINE
    Vector2d center = object;
    double diffx, diffy;
    if(ores[0] < 2) {
        diffx = size[0];
    }
    else {
        diffx = 2*size[1] / (ores[1]-1);
    }
    if(ores[1] < 2) {
        diffy = size[1];
    }
    else {
        diffy = 2*size[1] / (ores[1]-1);
    }

    for(int i = 0; i < ores[0]; i++) {
        Vector2d pos = center - Vector2d(size[0], 0) + Vector2d(diffx*i, 0);
        Vector3d col(1.0, 0.0, 0.0);
        Vector2d ph = pos - center;
        Vector2d vel;
        if(vt == constant) {
            vel = Vector2d(xvel, yvel);
        }
        else if(vt == linear) {
            if(centered) {
                vel = Vector2d(d(0), d(1));
            }
            else {
                vel = Vector2d(pos(0), pos(1));
            }
        }
        else if(vt == rotation) {
            if(centered) {
                vel = Vector2d(-d(1), d(0));
            }
            else {
                vel = Vector2d(-pos(1), pos(0));
            }
        }
        Particle par(pos, scale*vel, col, pmass);
        par.B = B;
        parts.push_back(par);
    }
#endif
    //Single particle special case
    if(ores[0] == 0 && ores[1] == 0) {
        Vector2d vel = object;
        Vector3d col(1.0, 0.0, 0.0);
        Particle par(object, vel, col, pmass);
        parts.push_back(par);
    }
    if(jitter) {
        Particle par(center, scale*Vector2d(0,0), Vector3d(1,1,0), pmass);
        par.B = B;
        parts.push_back(par);
    }

    writeParticles(argv[1], parts);
    return 0;
}
