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

/// #define RAND

#define SQUARE 1 //could change to SQUARE, LINE

using namespace Eigen;
using benlib::range;

int main(int argc, char* argv[]) {
  std::vector<Particle> parts;

  // parameters of the object
  double size[2] = {0.25, 0.25};
  Vector2d object(0.0,0.5);
  int ores[2] = {50,50};
  double pmass = 1.0;
  double rotation = 0.333333333;
  int numPart = 648;
  
#if SQUARE
  //Vector2d center = object + (Vector2d(size[0],size[1]) * 0.5);
  Vector2d center = object;
  object = object - Vector2d(size[0], size[1]);
  //Set up particles at each object vertex
  double diffx = 2*size[0] / (ores[0]-1);
  double diffy = 2*size[1] / (ores[1]-1);
  for(int i = 0; i < ores[0]; i++) {
	for(int j = 0; j < ores[1]; j++) {
	  Vector2d pos = object + Vector2d(diffx*i, diffy*j);
	  Vector2d quad = pos - center;
	  Vector3d col(1.0, 0.0, 0.0);
	  if((quad(0) < 0 && quad(1) < 0) || (quad(0) >= 0 && quad(1) >= 0)) {
	    col = Vector3d(0.0, 1.0, 0.0);
	  }

	  Vector2d d = pos-center;
	  Vector2d vel = rotation*Vector2d(-d(1), d(0));

	  Particle par(pos, vel, col, pmass);
	  Matrix2d B;
	  B << 0, -rotation, rotation, 0;
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
  double diffx = 2*size[0] / (ores[0]-1);
  double diffy = 2*size[1] / (ores[1]-1);
#ifdef RAND
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
	  /// Vector3d col = ((double)j/(ores[1]-1))*Vector3d(1, 0, 0);
      Vector3d col(1.0, 0.0, 0.0);
      Vector2d quad = pos - center;
      if((quad(0) < 0 && quad(1) < 0) || (quad(0) > 0 && quad(1) > 0)) {
          col = Vector3d(0.0, 1.0, 0.0);
      }
	  Vector2d ph = pos - object;
	  if( ((ph(0)*ph(0))/(size[0]*size[0])) + ((ph(1)*ph(1))/(size[1]*size[1])) < 1+EPS) {
		Vector2d vel = rotation*Vector2d(-ph(1), ph(0));
        /// Vector2d vel(pos(0), pos(1));
        Particle par(pos, 1.0*vel, col, pmass);
        Matrix2d B;
        B << 0, -rotation, rotation, 0;
        par.B = B;
		parts.push_back(par);
	  }
	}
#ifndef RAND
  }
#endif
#elif LINE
  Vector2d center = object;
  double diffx = 2*size[0] / (ores[0]-1);
  double diffy = 2*size[1] / (ores[1]-1);
  for(int i = 0; i < ores[0]; i++) {
      Vector2d pos = center - Vector2d(size[0], 0) + Vector2d(diffx*i, 0);
      Vector3d col(1.0, 0.0, 0.0);
      Vector2d ph = pos - center;
      Vector2d vel = rotation*Vector2d(-ph(1), ph(0));
      Particle par(pos, vel, col, pmass);
      Matrix2d B;
      B << 0, -rotation, rotation, 0;
      par.B = B;
      parts.push_back(par);
  }
#endif
  if(ores[0] == 0 && ores[1] == 0) {
      Vector2d vel = object;
      Vector3d col(1.0, 0.0, 0.0);
      Particle par(object, vel, col, pmass);
      parts.push_back(par);
  }

  writeParticles(argv[1], parts);
  return 0;
}
