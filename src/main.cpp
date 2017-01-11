#include "defines.hpp"
#include "mpm.hpp"
#include "opengl.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>

using namespace Eigen;

#ifndef NDEBUG
extern std::ofstream debug;
#endif


int main(int argc, char** argv) {
	double timeSinceLastFrame = 1.0;
	int frame = 0;
	int iters = 0;

    if(argc < 2) {
        std::cout << "Usage: ./mpm <configfile> [debug and video outputs name]\n";
        std::exit(0);
    }
	std::string outfile;
	if (argc < 3) {
	  std::string inputfname = std::string(argv[1]);
	  auto const start = inputfname.find_last_of('/');
	  auto const end = inputfname.find_last_of('.');
	  outfile = inputfname.substr(start+1,end-start-1);
	} else {
	  outfile = std::string(argv[2]);
	}
#ifndef NDEBUG
    std::string dbout = outfile + std::string(".txt");
    debug.open(dbout.c_str());
#endif
  
    std::string config = std::string(argv[1]);
    World world(config);
    world.init();
  
	while (world.elapsedTime < world.totalTime) {
	  if (timeSinceLastFrame > 1.0/30.0) {
   	    std::ostringstream ss;
		ss << std::setw(5) << std::setfill('0') << frame;
		std::string pframe(ss.str());
		std::string parOut = outfile + "." + pframe + ".bgeo";
        writeParticles(parOut.c_str(), world.particles);
		frame++;
		iters = 0;
		timeSinceLastFrame = 0.0;
		printf("\n");
	  }
	  
	  world.step();
	  timeSinceLastFrame += world.dt;
	  world.elapsedTime += world.dt;
	  iters++;

	  //std::cout<<world.dt<<std::endl;
	  printf("Frame %d/%d Step: %d/%d\r", frame, (int)(30.0*world.totalTime), iters, (int)(1.0/(30.0*world.dt)));
	  std::cout<<std::flush;
	}

    return 0;
}


