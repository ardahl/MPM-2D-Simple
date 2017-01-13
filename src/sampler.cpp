#include "world.hpp"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <Partio.h>
#include <Eigen/Geometry>
#include "json/json.h"
#include "range.hpp"

using namespace Eigen;
using benlib::range;

int main(int argc, char* argv[]) {
  std::vector<Particle> parts;

  // parameters of the object
  double size[2] = {1.0, 1.0};
  Vector2d object(0.0,1.5);
  int ores[2] = {101,101};
  double pmass = 1.0;
  double rotation = 100.0;
  
#if 0
  Vector2d center = object + (Vector2d(size[0],size[1]) * 0.5);
  //Set up particles at each object vertex
  double diffx = size[0] / (ores[0]-1);
  double diffy = size[1] / (ores[1]-1);
  for(int i = 0; i < ores[0]; i++) {
	for(int j = 0; j < ores[1]; j++) {
	  Vector2d pos = object + Vector2d(diffx*i, diffy*j);
	  Vector3d col = ((double)j/(ores[1]-1))*Vector3d(1, 0, 0);

	  Vector2d d = pos-center;
	  Vector2d vel = rotation*Vector2d(-d(1), d(0));

	  Particle par(pos, vel, col, pmass, 0.0, 0.0);
	  parts.push_back(par);
	}
  }
#else
  printf("Making Circle\n");
  Vector2d center = object;
        //non-randomly make a circle
        //technically this is an ellipse because I'm using the same data as the
        //square and just using the size vector as radius of the semi-major and semi-minor axes
        //just make a square and reject those outside the ellipse.
        //~78.5% of the resx*resy particles are accepted - Pi/4 * (l*w) particles 
  double diffx = 2*size[0] / (ores[0]-1);
  double diffy = 2*size[1] / (ores[1]-1);
  for(int i = 0; i < ores[0]; i++) {
	for(int j = 0; j < ores[1]; j++) {
	  Vector2d pos = center - Vector2d(size[0], size[1]) + Vector2d(diffx*i, diffy*j);
	  Vector3d col = ((double)j/(ores[1]-1))*Vector3d(1, 0, 0);
	  Vector2d ph = pos - object;
	  if( ((ph(0)*ph(0))/(size[0]*size[0])) + ((ph(1)*ph(1))/(size[1]*size[1])) < 1+EPS) {
		Vector2d vel = rotation*Vector2d(-ph(1), ph(0));
		Particle par(pos, vel, col, pmass);
		parts.push_back(par);
	  }
	}
  }
#endif

  writeParticles(argv[1], parts);
  return 0;
}
