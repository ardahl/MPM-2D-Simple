/**************
 * This main will dump the particles to a file for a seperate viewer.
 * 
 * Output format:
 * num_frames num_particles particle_color
 * grid_x0 grid_y0 spacing dx dy
 * 
 * p0
 * p1
 * ...
 * pn
 **************/

#include "defines.hpp"
#include "mpm.hpp"
#include "configparser.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>

using namespace Eigen;

#ifndef NDEBUG
extern std::ofstream debug;
#endif

int main(int argc, char** argv) {
    if(argc < 3) {
        std::cout << "Usage: ./mpm <configfile> <output file>\n";
        std::exit(0);
    }
    std::string outfile = std::string(argv[2]);
    std::ofstream dumpOut;
    dumpOut.open(outfile.c_str());
    #ifndef NDEBUG
    std::string dbout = outfile + std::string(".txt");
    debug.open(dbout.c_str());
    #endif
    
    std::string config = std::string(argv[1]);
    
    ConfigParser conf(config);
    Material* m = new Material(config);
    m->init();
    
    int frame = 0, totalFrames = conf.getInt("time")*30;
    int iters = conf.getInt("iterations");
    double itersInv = 1.0/iters;
    //setup
    dumpOut << totalFrames << " " << m->particles.size() << "\n";
    dumpOut << m->x0(0) << " " << m->x0(1) << " " << m->h << " " << m->m << " " << m->n << "\n\n";
    while(frame < totalFrames) {
        frame++;
        //Perform step
        for(int i = 0; i < iters; i++) {
            /// printf("Step %d\n", i);
            m->step((1.0/30.0)*itersInv);
            printf("Frame %d/%d Step: %d/%d\r", frame, totalFrames, i+1, iters);
        }
        printf("\n");
        
        //Dump data
        for(size_t i = 0; i < m->particles.size(); i++) {
            Particle* p = m->particles[i];
            dumpOut << p->x(0) << " " << p->x(1) << "\n";
        }
        dumpOut << "\n";
    }
    
    #ifndef NDEBUG
    debug.close();
    #endif
    dumpOut.close();
    return 0;
}
