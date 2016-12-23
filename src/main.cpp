#include "defines.hpp"
#include "mpm.hpp"
#include "opengl.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <Partio.h>

using namespace Eigen;

#ifndef NDEBUG
extern std::ofstream debug;
#endif

void writeParticles(const char *fname, const std::vector<Particle> &particles);
void readParticles(const char *fname, std::vector<Particle> &particles);

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


void writeParticles(const char *fname, const std::vector<Particle> &particles) {
	Partio::ParticlesDataMutable *data = Partio::create();
	Partio::ParticleAttribute xattr;
	Partio::ParticleAttribute uattr;
	Partio::ParticleAttribute sattr;
	Partio::ParticleAttribute gattr;
	Partio::ParticleAttribute cattr;
	Partio::ParticleAttribute mattr;
	Partio::ParticleAttribute rattr;
	Partio::ParticleAttribute vattr;
	data->addParticles(particles.size());

	data->addAttribute("position", Partio::VECTOR, 2);
	data->addAttribute("velocity", Partio::VECTOR, 2);
	data->addAttribute("stress", Partio::VECTOR, 4);
	data->addAttribute("gradient", Partio::VECTOR, 4);
	data->addAttribute("color", Partio::VECTOR, 3);
	data->addAttribute("mass", Partio::FLOAT, 1);
	data->addAttribute("rho", Partio::FLOAT, 1);
	data->addAttribute("vol", Partio::FLOAT, 1);
	
	data->attributeInfo("position", xattr);
	data->attributeInfo("velocity", uattr);
	data->attributeInfo("stress", sattr);
	data->attributeInfo("gradient", gattr);
	data->attributeInfo("color", cattr);
	data->attributeInfo("mass", mattr);
	data->attributeInfo("rho", rattr);
	data->attributeInfo("vol", vattr);

	for (int i=0; i < particles.size(); i++) {
		const Particle &p = particles[i];
		float *x = data->dataWrite<float>(xattr, i);
		float *u = data->dataWrite<float>(uattr, i);
		float *s = data->dataWrite<float>(sattr, i);
		float *g = data->dataWrite<float>(gattr, i);
		float *c = data->dataWrite<float>(cattr, i);
		float *m = data->dataWrite<float>(mattr, i);
		float *r = data->dataWrite<float>(rattr, i);
		float *v = data->dataWrite<float>(vattr, i);

		x[0] = p.x(0), x[1] = p.x(1);
		v[0] = p.v(0), v[1] = p.v(1);
		s[0] = p.stress(0,0), s[1] = p.stress(0,1), s[2] = p.stress(1,0), s[3] = p.stress(1,1);
		g[0] = p.gradient(0,0), g[1] = p.gradient(0,1), g[2] = p.gradient(1,0), g[3] = p.gradient(1,1);
		c[0] = p.color(0), c[1] = p.color(1), s[2] = p.color(2);
		m[0] = p.m;
		r[0] = p.rho;
		v[0] = p.vol;
	}

	Partio::write(fname, *data);
	data->release();
}


void readParticles(const char *fname, std::vector<Particle> &particles) {
	Partio::ParticlesDataMutable *data = Partio::read(fname);
	Partio::ParticleAttribute xattr;
	Partio::ParticleAttribute uattr;
	Partio::ParticleAttribute sattr;
	Partio::ParticleAttribute gattr;
	Partio::ParticleAttribute cattr;
	Partio::ParticleAttribute mattr;
	Partio::ParticleAttribute rattr;
	Partio::ParticleAttribute vattr;

	bool position = data->attributeInfo("position", xattr);
	bool velocity = data->attributeInfo("velocity", uattr);
	bool stress = data->attributeInfo("stress", sattr);
	bool gradient = data->attributeInfo("gradient", gattr);
	bool color = data->attributeInfo("color", cattr);
	bool mass = data->attributeInfo("mass", mattr);
	bool rho = data->attributeInfo("rho", rattr);
	bool vol = data->attributeInfo("vol", vattr);

	particles.resize(data->numParticles());

	for (int i=0; i < particles.size(); i++) {
	  Particle &p = particles[i];
	  if (position) {
		float *x = data->dataWrite<float>(xattr, i);
		p.x[0] = x[0], p.x[1] = x[1];
	  } else {
		p.x = Vector2d(0.0, 0.0);
	  }
	  if (velocity) {
		float *v = data->dataWrite<float>(uattr, i);
		p.v[0] = v[0], p.v[1] = v[1];
	  } else {
		p.v = Vector2d(0.0, 0.0);
		}
	  if (stress) {
		float *s = data->dataWrite<float>(sattr, i);
		p.stress(0,0) = s[0], p.stress(0,1) = s[1], p.stress(1,0) = s[2], p.stress(1,1) = s[3];
	  } else {
		p.stress = Matrix2d::Zero();
	  }
	  if (gradient) {
		float *g = data->dataWrite<float>(gattr, i);
		p.gradient(0,0) = g[0], p.gradient(0,1) = g[1], p.gradient(1,0) = g[2], p.gradient(1,1) = g[3];
	  } else {
		p.gradient = Matrix2d::Identity();
	  }
	  if (color) {
		float *c = data->dataWrite<float>(cattr, i);
		p.color(0) = c[0], p.color(1) = c[1], p.color(2) = c[2];
	  } else {
		p.color = Vector3d(1.0,1.0,1.0);
	  }
	  if (mass) {
		float *m = data->dataWrite<float>(mattr, i);
		p.m = m[0];
	  } else {
		p.m = 1.0;
	  }
	  if (rho) {
		float *r = data->dataWrite<float>(rattr, i);
			p.rho = r[0];
	  } else {
		p.rho = 1.0;
	  }
	  if (vol) {
		float *v = data->dataWrite<float>(vattr, i);
		p.vol = v[0];
	  } else {
		p.vol = 1.0;
	  }
	}
}
