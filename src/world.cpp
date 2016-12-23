#include "world.hpp"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include "json/json.h"
#include "range.hpp"

using namespace Eigen;
using benlib::range;

#ifndef NDEBUG
std::ofstream debug;
int count = 0;
#endif

//TODO: Different colored particles
//TODO: Make circle

//TODO: Check if forces to cells are what you should be getting 
//applying it to particles
//Plot 2norm of (F-R)
//Rotation = (speed*dt)/(mass of body) * pi/2

World::World(std::string config) {
  stepNum = 0;
  elapsedTime = 0.0;
  filename = config;

  std::ifstream ins(filename.c_str());
  
  Json::Value root;
  Json::Reader jReader;

  if(!jReader.parse(ins, root)){
	std::cout << "couldn't read input file: " << filename << '\n'
			  << jReader.getFormattedErrorMessages() << std::endl;
	exit(1);
  }
  double size[2];
  Vector2d object;
  Vector3d c;
  Vector2d object, grav;
  std::string str, objType = "square";

  x0 = Vector2d(0, 0);
  Force* frc = NULL;
  int ores[2];

  auto objectsIn = root["objects"];
  for (auto i : range(objectsIn.size())) {
    auto typeIn = objectsIn[i]["type"];
	if (typeIn.size() != 1) {
	  std::cout<< "bad object type, skipping" << std::endl;
      std::cout << typeIn.size() << std::endl;
	  continue;
	}
    objType = typeIn[0].asString();
    auto locationIn = objectsIn[i]["location"];
	if (locationIn.size() != 2) {
	  std::cout<< "bad object location, skipping" << std::endl;
	  continue;
	}
	object(0) = locationIn[0].asDouble();
	object(1) = locationIn[1].asDouble();
	auto sizeIn = objectsIn[i]["size"];
	if (locationIn.size() != 2) {
	  std::cout<< "bad object size, skipping" << std::endl;
	  continue;
	}
	size[0] = sizeIn[0].asDouble();
	size[1] = sizeIn[1].asDouble();
	auto resIn = objectsIn[i]["resolution"];
	if (locationIn.size() != 2) {
	  std::cout<< "bad object resolution, skipping" << std::endl;
	  continue;
	}
	ores[0] = resIn[0].asInt();
	ores[1] = resIn[1].asInt();
  }

  auto gridIn = root["grid"];
  {
	auto originIn = gridIn["origin"];
	if (originIn.size() != 2) {
	  std::cout<< "bad grid origin, skipping" << std::endl;
	} else {
	  origin = Vector2d(originIn[0].asDouble(),originIn[1].asDouble());
	}
	auto sizeIn = gridIn["size"];
	if (sizeIn.size() != 2) {
	  std::cout<< "bad grid size, skipping" << std::endl;
	} else {
	  res[0] = sizeIn[0].asInt();
	  res[1] = sizeIn[1].asInt();
	}

	h = gridIn["h"].asDouble();
  }
  double pmass = root["mass"].asDouble();
  auto colorIn = root["color"];
  if (colorIn.size() != 3) {
	std::cout<< "bad color, skipping" << std::endl;
  } else {
	c = Eigen::Vector3d(colorIn[0].asDouble(), colorIn[1].asDouble(), colorIn[2].asDouble());
  }

  auto lameIn = root["lame"];
  if (lameIn.size() != 2) {
	std::cout<< "bad lame, skipping" << std::endl;
  } else {
	lambda = lameIn[0].asDouble();
	mu = lameIn[1].asDouble();
  }
  auto gravityIn = root["gravity"];
  if (gravityIn.isNull() || gravityIn.size() != 2) {
	std::cout<< "no gravity" << std::endl;
	gravityEnabled = false;
  } else {
    gravity(0)= gravityIn[0].asDouble();
    gravity(1)= gravityIn[1].asDouble();
	gravityEnabled = true;
  }

  auto rotationIn = root["rotate"];
  if (rotationIn.isNull()) {
	rotationEnabled = false;
  } else {
	rotationEnabled = true;
	rotation = rotationIn.asDouble();
  }
  rotspeed = root["rotate"].asDouble();
  if(objType == "square") {
    center = object + (Vector2d(size[0],size[1]) * 0.5);
    origin(0) += h/2.0; origin(1) += h/2.0;
  }
  else if(objType == "circle") {
      center = origin;
  }
  frc = new Rotate(center, rotspeed);
  forces.push_back(frc);

  Vector2d xGrid = x0 + Vector2d(h/2.0,h/2.0);
  x0 = xGrid;
  x1 = x0 + h*Vector2d(m, n);
  if(objType == "square") {
    //Set up particles at each object vertex
    double diffx = l / (ores[0]-1);
    double diffy = w / (ores[1]-1);
    for(int i = 0; i < ores[0]; i++) {
        for(int j = 0; j < ores[1]; j++) {
            Vector2d pos = object + Vector2d(diffx*i, diffy*j);
            Vector3d col = ((double)j/(ores[1]-1))*Vector3d(1, 0, 0);
            Particle* par = new Particle(pos, Vector2d(0,0), col, pmass);
            particles.push_back(par);
        }
    }
  }
  if(objType == "circle") {
      //non-randomly make a circle
      //technically this is an ellipse because I'm using the same data as the
      //square and just using the size vector as radius of the semi-major and semi-minor axes
      //just make a square and reject those outside the ellipse.
      //~78.5% of the resx*resy particles are accepted - Pi/4 * (l*w) particles 
      double diffx = 2*l / (ores[0]-1);
      double diffy = 2*w / (ores[1]-1);
      for(int i = 0; i < ores[0]; i++) {
        for(int j = 0; j < ores[1]; j++) {
            Vector2d pos = center - Vector2d(l, w) + Vector2d(diffx*i, diffy*j);
            Vector3d col = ((double)j/(ores[1]-1))*Vector3d(1, 0, 0);
            Vector2d ph = pos - object;
            if(i == resx/2.0 && j == resy/2.0) {
                printf("Pos: (%f, %f)\n", pos(0), pos(1));
                printf("Obj: (%f, %f)\n", object(0), object(1));
                printf("ph: (%f, %f)\n", ph(0), ph(1));
                printf("l, w: %f, %f\n", l, w);
            }
            /// double rx = l/2.0, ry = w/2.0;
            if( ((ph(0)*ph(0))/(l*l)) + ((ph(1)*ph(1))/(w*w)) < 1+EPS) {
                Particle* par = new Particle(pos, Vector2d(0,0), col, pmass);
                particles.push_back(par);
            }
        }
      }
  }
  
  //Average position for center of mass
  Vector2d avePos = Vector2d::Zero();
  for(size_t i = 0; i < particles.size(); i++) {
	avePos += particles[i].x;
  }
  avePos /= particles.size();
  
  printf("Number of particles: %d\n", (int)particles.size());
  printf("Dimentions: %dx%d\n", res[0], res[1]);
  printf("Grid Spacing: %f\n", h);
  printf("Lame Constants: %f, %f\n", lambda, mu);
  printf("Gravity: (%f, %f)\n", gravity(0), gravity(1));
  printf("Center of Mass: (%f, %f), (%f, %f)\n", center(0), center(1), avePos(0), avePos(1));
  printf("X0: (%f, %f)\n", origin(0), origin(1));
  printf("X1: (%f, %f)\n", origin(0)+h*(res[0]-1), origin(1)+h*(res[1]-1));
  mass = new double[res[0]*res[1]];
  vel = new Vector2d[res[0]*res[1]];
  velStar = new Vector2d[res[0]*res[1]];
  frc = new Vector2d[res[0]*res[1]];
  //TODO: Check gradient
}

void World::init() {
    particleVolumesDensities();
}

inline double weight(double x) {
    if(x < 1.0) {
        return 0.5*x*x*x - x*x + (2.0/3.0);
    } else if(x < 2.0) {
        return (-1.0/6.0)*x*x*x + x*x - 2.0*x + (4.0/3.0);
    }
    return 0;
}

inline double gradweight1d(double x) {
    double ax = std::abs(x);
    if(x < 1.0) {
        return 1.5*x*ax - 2.0*x;
    } else if(ax < 2.0) {
        return (-0.5)*x*ax + 2.0*x - 2.0*(x/ax);
    }
    return 0;
}

inline Vector2d gradweight(const Vector2d &offset, double h) {
	return Vector2d(gradweight1d(offset(0))*weight(offset(1))/h, 
		weight(offset(0))*gradweight1d(offset(1))/h);
}


inline void bounds(const Vector2d &offset, const int res[2], int *xbounds, int *ybounds) {
  xbounds[0] = std::max(0, ((int)std::ceil(-2.0 + offset(0))) - 1);
  xbounds[1] = std::min(res[0]-1, ((int)std::floor(2.0 + offset(0))) + 1);
  ybounds[0] = std::max(0, ((int)std::ceil(-2.0 + offset(1))) - 1);
  ybounds[1] = std::min(res[1]-1, ((int)std::floor(2.0 + offset(1))) + 1);
}



/******************************
 * Compute_Particle_Volumes_And_Densities
 *      Raserise_Particle_Data_To_Grid
 *      for each particle p do
 *          rho_p = 0
 *          for each grid node i s.t. w_ip > 0
 *              rho_p += w_ip * m_i
 *          end for
 *          rho_p /= h^3
 *          V_p = m_p / rho_p
 *      end for
 *****************************/
void World::particleVolumesDensities() {
    particlesToGrid();
    for(size_t i = 0; i < particles.size(); i++) {
        Particle &p = particles[i];
		p.rho = 0.0;
		Vector2d offset = (p.x - origin) / h;
		int xbounds[2], ybounds[2];
		bounds(offset, res, xbounds, ybounds);
        for(int j = xbounds[0]; j < xbounds[1]; j++) {
		  double w1 = weight(offset(0) - j);
            for(int k = ybounds[0]; k < ybounds[1]; k++) {
			  double r = mass[j*res[1] + k] * w1 * weight(offset(1) - k);
                p.rho += r;
            }
        }
        p.rho /= (h*h);            //3D: h*h*h
        #ifndef NDEBUG
        if(std::isnan(p.rho) || std::isinf(p.rho)) {
            printf("Paricle %d rho has NaN: %f\n", (int)i, p.rho);
            exit(0);
        }
        #endif
        p.vol = p.m / p.rho;
        #ifndef NDEBUG
        if(std::isnan(p.vol) || std::isinf(p.vol)) {
            printf("Paricle %d volume has NaN: %f\nMass: %f\nDensity: %f\n", (int)i, p.vol, p.m, p.rho);
            exit(0);
        }
        #endif
    }
}

/******************************
 * Algorithm Process from Stomakhin et al. 2013 and Yue et al. 2015
 *
 * Compute_Particle_Volumes_And_Densities
 * while true
 *      Rasterize_Particle_Data_to_Grid
 *      Compute_Grid_Forces
 *      Update_Grid_Velocities
 *      Update_Deformation_Gradient
 *      Update_Particle_Velocities
 *      Update_Particle_Positions
 ******************************/
void World::step(double dt) {
    particlesToGrid();
    computeGridForces();
    updateGridVelocities(dt);
    updateGradient(dt);
    gridToParticles(dt);
    stepNum++;
}

/******************************
 * Rasterize_Particle_Data_To_Grid
 *      for each grid node i do
 *          m_i = 0
 *          v_i = 0
 *      end for
 *      for each particle p do
 *          for each grid node i s.t. w_ip > 0 do
 *              m_i += w_ip * m_p
 *              v_i += w_ip * m_p * v_p
 *          end for
 *      end for
 *      for each grid cell i do
 *          v_i /= m_i
 *      end for
 *
 *****************************/
void World::particlesToGrid() {
    {Vector2d *v = vel; for (unsigned int i=0; i<res[0]*res[1]; i++, v++) (*v) = Vector2d(0.0,0.0);}  
    {double *m = mass; for (unsigned int i=0; i<res[0]*res[1]; i++, m++) (*m) = 0.0;}

    for(size_t i = 0; i < particles.size(); i++) {
		Particle &p = particles[i];
		Vector2d offset = (p.x - origin) / h;
		int xbounds[2], ybounds[2];
		bounds(offset, res, xbounds, ybounds);
        for(int j = xbounds[0]; j < xbounds[1]; j++) {
   		    double w1 = weight(offset(0) - j);
            for(int k = ybounds[0]; k < ybounds[1]; k++) {
    			double w = w1*weight(offset(1) - k);
                mass[j*res[1] + k] += w * p.m;
                vel [j*res[1] + k] += w * p.m * p.v;
            }
        }
	}
    /// #pragma omp parallel for collapse(2)
	for(int i = 0; i < res[0]; i++) {
		for(int j = 0; j < res[1]; j++) {
            if(mass[i*res[1] + j] < EPS) {
                vel[i*res[1] + j] = Vector2d(0.0, 0.0);
            }
            else {
                vel[i*res[1] + j] /= mass[i*res[1] + j];
            }
            #ifndef NDEBUG
            if(vel[i*res[1] + j].hasNaN()) {
                printf("interpolated vel NaN at (%d, %d)\n", i, j);
                std::cout << vel[i*res[1] + j] << std::endl;
                exit(0);
            }
            #endif
		}
	}
}

/******************************
 * Compute_Grid_Forces
 *      for each grid node i do
 *          f_i = 0
 *      end for
 *      for each particle p do
 *          J_p = det(F_p)
 *          epsilon_p = (1/2)*(F^T * F - I)
 *          stress_p = lambda * trace(epsilon_p) * I + 2 * mu * epsilon_p
 *          for each grid cell i s.t. w_ip > 0
 *              f_i -= V_p * J_p * stress_p * grad(w_ip)
 *          end for
 *      end for
 *****************************/
void World::computeGridForces() {
    {Vector2d *f = frc; for (unsigned int i=0; i<res[0]*res[1]; i++, f++) (*f) = Vector2d(0.0,0.0);}

    for(size_t i = 0; i < particles.size(); i++) {
        Particle &p = particles[i];
        double J = p.gradient.determinant();
        Matrix2d eps = 0.5 * (p.gradient.transpose() * p.gradient - Matrix2d::Identity());
        double trace = eps.trace();
        Matrix2d stress = lambda*trace*Matrix2d::Identity() + 2.0*mu*eps;

		Vector2d offset = (p.x - origin) / h;
		int xbounds[2], ybounds[2];
		bounds(offset, res, xbounds, ybounds);
        for(int j = xbounds[0]; j < xbounds[1]; j++) {
            for(int k = ybounds[0]; k < ybounds[1]; k++) {
			  Vector2d accumF = p.vol * J * stress * gradweight(Vector2d(offset(0)-i,offset(1)-j),h);
                frc[j*res[1] + k] -= accumF;
                #ifndef NDEBUG
                if(frc[j*res[1] + k].hasNaN()) {
                    printf("\nf NaN at (%d, %d)\n", j, k);
                    std::cout << "Force:\n" << frc[j*res[1] + k] << std::endl;
                    std::cout << "Volume: " << p.vol << std::endl;
                    std::cout << "Determinant: " << p.gradient.determinant() << std::endl;
                    std::cout << "Stress:\n" << p.stress << std::endl;
                    std::cout << "Gradient:\n" << gradweight(Vector2d(offset(0)-i, offset(1)-j),h) << std::endl;
                    exit(0);
                }
                #endif
            }
        }
    }
}

/******************************
 * Update_Grid_Velocities
 *      for each grid node i do
 *          v^*_i = v_i^n + (dt * f_i)/m_i + dt * g
 *      end for
 *****************************/
void World::updateGridVelocities(double dt) {
    /// #pragma omp parallel for collapse(2)
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
		  int index = i*res[1] + j;
		  if(mass[index] < EPS) {
                velStar[index] = Vector2d(0, 0);
            }
            else {
			    Vector2d extfrc(0.0,0.0);
				if (gravityEnabled) {
				  extfrc += mass[index]*gravity;
				}
				if (rotationEnabled) {
				  Vector2d d = origin+Vector2d(h*i,h*j)-center;
				  extfrc += rotation*Vector2d(-d(1), d(0));
				}
                velStar[index] = vel[index] + dt * (1.0/mass[index]) * (frc[index] + extfrc); //dt*g
            }
            #ifndef NDEBUG
		    if(velStar[i*res[1] + j].hasNaN()) {
                printf("velStar NaN at (%d, %d)\n", i, j);
                std::cout << velStar[index] << std::endl;
                exit(0);
            }
            #endif
        }
    }
}

/******************************
 * Update_Deformation_Gradient
 *      for each particle p do
 *          grad(v_p) = 0
 *          for each grid node i s.t. w_ip > 0
 *              grad(v_p) += v^*_i * (grad(w_ip))^T
 *          end for
 *          f_p,n+1 = I + dt*grad(v_p)
 *          F_p,n+1 *= f_p,n+1
 *      end for
 *****************************/
void World::updateGradient(double dt) {
    /// #pragma omp parallel for 
    for(size_t i = 0; i < particles.size(); i++) {
        Particle &p = particles[i];
        Matrix2d gradV = Matrix2d::Zero();
		Vector2d offset = (p.x - origin) / h;
		int xbounds[2], ybounds[2];
		bounds(offset, res, xbounds, ybounds);
        for(int j = xbounds[0]; j < xbounds[1]; j++) {
            for(int k = ybounds[0]; k < ybounds[1]; k++) {
			    int index = j*res[1]+k;
                #ifndef NDEBUG
                if(velStar[index].hasNaN()) {
                    printf("gradV velStar has NaN at (%d, %d)\n", j, k);
                    std::cout << velStar[index] << std::endl;
                    exit(0);
                }
                if(gradweight(Vector2d(offset(0)-i,offset(1)-j),h).transpose().hasNaN()) {
                    printf("gradV gradW has NaN at (%d, %d)\n", j, k);
                    std::cout << gradweight(Vector2d(offset(0)-i,offset(1)-j),h).transpose() << std::endl;
                    exit(0);
                }
                #endif
                Matrix2d accumGrad = velStar[index] * gradweight(Vector2d(offset(0)-i,offset(1)-j),h).transpose();
                gradV += accumGrad;
            }
        }
        #ifndef NDEBUG
        if(gradV.hasNaN()) {
            printf("gradV has NaN\n");
            std::cout << gradV << std::endl;
            exit(0);
        }
        #endif
        Matrix2d fp = Matrix2d::Identity() + dt*gradV;
        p.gradient = fp*p.gradient;
    }
}

/******************************
 * Update_Particle_Velocities
 *      for each particle p do
 *          v_pic = 0
 *          v_flip = v_p
 *          for all grid cells i s.t. w_ip > 0
 *              v_pic += w_ip * v^*_i
 *              v_flip += w_ip * (v^*_i - v_i)
 *          end for
 *          v_p = (1-alpha)*v_pic + alpha*v_flip
 *      end for
 *
 * Update_Particle_Positions
 *      for each particle p do
 *          x_p += dt * v_p
 *      end for
 *****************************/
void World::gridToParticles(double dt) {
    double alpha = 0.95;
    /// #pragma omp parallel for
    for(size_t i = 0; i < particles.size(); i++) {
		Particle &p = particles[i];
		//Update velocities
        Vector2d pic = Vector2d::Zero();
        Vector2d flip = p.v;
        Vector2d tmpPic, tmpFlip;
		Vector2d offset = (p.x - origin) / h;
		int xbounds[2], ybounds[2];
		bounds(offset, res, xbounds, ybounds);
        for(int j = xbounds[0]; j < xbounds[1]; j++) {
   		    double w1 = weight(offset(0) - j);
            for(int k = ybounds[0]; k < ybounds[1]; k++) {
    			double w = w1*weight(offset(1) - k);
				int index = j*res[1] + k;
                tmpPic = w * velStar[index];
                pic += tmpPic;

                tmpFlip = w * (velStar[index] - vel[index]);
                flip += tmpFlip;
            }
        }
        #ifndef NDEBUG
        if(pic.hasNaN()) {
            printf("\n\nPIC Vel has NaN\n");
            std::cout << pic << std::endl;
            exit(0);
        }
        if(flip.hasNaN()) {
            printf("FLIP Vel has NaN\n");
            std::cout << flip << std::endl;
            exit(0);
        }
        #endif
		p.v = (alpha * flip) + ((1 - alpha) * pic);
        //Mass proportional damping
        p.v *= 0.99999;
        #ifndef NDEBUG
        if(p.v.hasNaN()) {
            printf("Vel has NaN\n");
            std::cout << p.v << std::endl;
            exit(0);
        }
        #endif

        //Update Positions
        p.x += dt * p.v;
        #ifndef NDEBUG
        if(p.x.hasNaN()) {
            printf("Pos has NaN\n");
            std::cout << p.x << std::endl;
            exit(0);
        }
        #endif
        //Boundary collisions
        //Make sure there's a boundary of 2 cells for each side since the weighting
        //function will touch 2 cells out
        double lx = origin[0]+2*h;
        double ly = origin[1]+2*h;
        double ux = origin[0]+(res[0]-2)*h;
        double uy = origin[1]+(res[1]-2)*h;
        if(p.x(0) < lx) {
            p.x(0) = lx+EPS;
            p.v(0) *= -1;          //reverse velocity
        }
        if(p.x(1) < ly) {
            p.x(1) = ly+EPS;
            p.v(1) *= -1;
        }
        if(p.x(0) > ux) {
            p.x(0) = ux-EPS;
            p.v(0) *= -1;
        }
        if(p.x(1) > uy) {
            p.x(1) = uy-EPS;
            p.v(1) *= -1;
        }
    }
}


