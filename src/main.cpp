#include "defines.hpp"
#include "mpm.hpp"
#include "opengl.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>

using namespace Eigen;

#ifndef NDEBUG
extern std::ofstream debug;
#endif

int main(int argc, char** argv) {
  int w = 800, h = 800;                           // size of window in pixels
  double xmin = 0, xmax = 1, ymin = 0, ymax = 1; // range of coordinates drawn
  int lastTime = 0, prevTime = 0, frame = 0;
  int seconds = 10*30, curr = 0;
  bool next = true;
  

  if(argc < 3) {
	std::cout << "Usage: ./mpm <configfile> <debug and video outputs name>\n";
	std::exit(0);
  }
  std::string outfile = std::string(argv[2]);
#ifndef NDEBUG
  std::string dbout = outfile + std::string(".txt");
  debug.open(dbout.c_str());
#endif
  
  std::string config = std::string(argv[1]);
  World world(config);
  world.init();
  
  int iters = 4000;
  double itersInv = 1.0/iters;
  for(int i = 0; i < iters; i++) {    //Hardcode for 30fps with dt of (1/3)e-5
	/// printf("Step %d\n", i);
	world.step((1.0/30.0)*itersInv);
	printf("Frame %d/%d Step: %d/%d\r", curr, seconds, i+1, iters);
  }
  printf("\n");
  
  return 0;
}
