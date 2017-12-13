#include "defines.hpp"
#include "eulerworld.hpp"
#include "opengl.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>

using namespace Eigen;
/// #define SLOMO

#ifndef NDEBUG
extern std::ofstream debug;
#endif

#ifdef SLOMO
//For slow motion, set start/end frame
int startFrame = 0;
int endFrame = 1;
int slowmoSteps = 100;
#endif

int main(int argc, char** argv) {
	double timeSinceLastFrame = 1.0;
	int frame = 0;
    int step = 0;
	int iters = 0;

    if(argc < 2) {
        std::cout << "Usage: ./euler <configfile> [debug and video outputs name]\n";
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
      #ifndef SLOMO
      if (timeSinceLastFrame > 1.0/30.0) {
      #else
      if(frame >= startFrame && frame <= endFrame && iters % slowmoSteps == 0) {
      #endif
        for (unsigned int obj = 0; obj<world.objects.size(); obj++) {
		  /// std::ostringstream ss;
          #ifndef SLOMO
		  /// ss << std::setw(2) << std::setfill('0') << obj << "." << std::setw(6)<< frame;
          #else
          /// ss << std::setw(2) << std::setfill('0') << obj << "." << std::setw(6)<< step;
		  #endif
          /// std::string pframe(ss.str());
		  /// std::string parOut = outfile + "-" + pframe + ".bgeo";
		  /// writeParticles(parOut.c_str(), world.objects[obj].particles);
          step++;
          #ifndef SLOMO
		  frame++;
		  iters = 0;
		  timeSinceLastFrame = world.dt;
          #endif
          printf("                                                                 \r");//clear out terminal line
		}
	  }
      #ifdef SLOMO
      if(timeSinceLastFrame > 1.0/30.0) {
          frame++;
          iters = 0;
          timeSinceLastFrame = world.dt;
      }
	  #endif
	  world.step();
	  timeSinceLastFrame += world.dt;
	  iters++;

	  printf("Frame: %d/%d\tStep: %d/%d\r", frame, (int)(30.0*world.totalTime), iters, (int)(1.0/(30.0*world.dt)));
	  std::cout<<std::flush;
	}
	std::cout<<std::endl<<std::endl;

    return 0;
}

