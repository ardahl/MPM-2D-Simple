#include "world.hpp"
#include <cstdio>
#include <iostream>
#include <limits>
#include <cmath>
#include <Partio.h>
#include <Eigen/Geometry>
#include "json/json.h"
#include "range.hpp"
//For matrix sqrt(), used for polar decomposition
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;
using benlib::range;

#ifndef NDEBUG
std::ofstream debug;
#endif

#ifdef INFO
/// std::ofstream debug("debug.txt");
std::ofstream cond("cond.txt");
std::ofstream xdiff("xdiff.txt");
std::ofstream ydiff("ydiff.txt");
double rot1 = 0;
double rot2 = 0;
#endif


//Rotation = (speed*dt)/(mass of body) * pi/2

//APIC transfer of reference coordinates. Update the particle reference coordinates
//Apic for eulerian solids with elasto-plastic model

//calc.cs.umbc.edu. Nothing in home directory.
//ssh to cal[01-12], has access to data

//Steps:
//1. Are APIC transfers exact
//2. Measure difference between gradients and rotations
//3. Get sparse disc to work
//4. APIC transfer of reference

//APIC Transfers
//Lower epsilon
//x1. Check for getting zeros
//x2. Upping velocities, grid size see how error values scale 
//x3. Plot x and y difference for sign information
//4. Try different kernel and see what happens
//5. Look at other degeneracies
    //Delete particles with small error and rerun with just the outer circle
    //10 Particles
    //particles in a line
    //particles in a ring

//Look at grid values as well
//Try other cases
//Print weight values
//Simplest case where it fails
    //Particles in the same cell, different cells
//Plot grid values error
//Try different linear field, see if errors are similar
//Check matrix C
//Check if C is similar to gradient of Vel. Might be able to replace it in gradient update

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
        std::exit(1);
    }
    
    Vector3d c;
    std::string str, objType = "square";

	dt = root.get("dt", 1.0/30.0).asDouble();
	totalTime = root.get("totalTime", 5.0).asDouble();

    double lambda=0, mu=0;                      //Lame Constants for stress
    double compression; //critical compression (sec. 5 of stomahkin)
    double stretch; //critical stretch (sec. 5 of stomahkin)
    double massPropDamp, alpha;
	// read global material properties, these may be overridden per object
    auto lameIn = root["lame"];
    if (lameIn.size() != 2) {
        std::cout<< "bad lame, skipping" << std::endl;
    } 
    else {
        lambda = lameIn[0].asDouble();
        mu = lameIn[1].asDouble();
    }
	massPropDamp = root.get("massPropDamp",0.99999).asDouble();
	alpha = root.get("alpha",0.95).asDouble();


    auto stretchIn = root["stretch"];
    auto compressIn = root["compression"];
    if (stretchIn.isNull() || compressIn.isNull()){
        plasticEnabled = false;
        std::cout << "no plasticity" << std::endl;
        compression = 0.0;
        stretch = 0.0;
    } else {
        plasticEnabled = true;
        compression = compressIn.asDouble();
        stretch = stretchIn.asDouble();
        std::cout << "Compression: " << compression << std::endl;
        std::cout << "Stretch: " << stretch << std::endl;
    }
  
    double pmass = root["mass"].asDouble();
    auto colorIn = root["color"];
    if (colorIn.size() != 3) {
        c = Vector3d(1, 0, 0);
    } 
    else {
        c = Eigen::Vector3d(colorIn[0].asDouble(), colorIn[1].asDouble(), colorIn[2].asDouble());
    }


    auto objectsIn = root["objects"];
    for (auto i : range(objectsIn.size())) {
  	    objects.emplace_back();
        objType = objectsIn[i].get("type", "square").asString();
        objects[i].type = objType;
		if (objType == "square" || objType == "circle") {
		  auto locationIn = objectsIn[i]["location"];
		  if (locationIn.size() != 2) {
            std::cout<< "bad object location, skipping" << std::endl;
            continue;
		  }
		  objects[i].object(0) = locationIn[0].asDouble();
		  objects[i].object(1) = locationIn[1].asDouble();
		  auto sizeIn = objectsIn[i]["size"];
		  if (sizeIn.size() != 2) {
            std::cout<< "bad object size, skipping" << std::endl;
            objects[i].size[0] = 0;
            objects[i].size[1] = 0;
            continue;
		  }
		  objects[i].size[0] = sizeIn[0].asDouble();
		  objects[i].size[1] = sizeIn[1].asDouble();
		  auto resIn = objectsIn[i]["resolution"];
		  if (resIn.size() != 2) {
            std::cout<< "bad object resolution, skipping" << std::endl;
            objects[i].ores[0] = 1;
            objects[i].ores[0] = 1;
            continue;
		  }
		  objects[i].ores[0] = resIn[0].asInt();
		  objects[i].ores[1] = resIn[1].asInt();
		} else {
		  auto const pos = config.find_last_of('/');
		  std::string partfilename = config.substr(0, pos+1);
		  partfilename = partfilename + objectsIn[i].get("filename", "input.bgeo").asString();
		  std::cout<<"loading "<<partfilename<<std::endl;
		  readParticles(partfilename.c_str(), objects[i].particles);
		}
		MaterialProps &mp = objects[i].mp;
		auto lameIn = objectsIn[i]["lame"];
		if (lameIn.size() == 2) {
		  mp.lambda = lameIn[0].asDouble();
		  mp.mu = lameIn[1].asDouble();
		} else {
		  mp.lambda = lambda;
		  mp.mu = mu;
		}
		mp.massPropDamp = objectsIn[i].get("massPropDamp",massPropDamp).asDouble();
		mp.alpha = objectsIn[i].get("alpha",alpha).asDouble();
		mp.stretch = objectsIn[i].get("stretch", stretch).asDouble();
		mp.compression = objectsIn[i].get("compression", stretch).asDouble();
        
        mp.pmass = objectsIn[i].get("particleMass", pmass).asDouble();
        mp.pmass = objectsIn[i].get("mass", mp.pmass).asDouble();
        mp.mass = objectsIn[i].get("objectMass", -1.0).asDouble();
        
        auto colorIn = objectsIn[i]["color"];
        if (colorIn.size() == 3) {
            objects[i].color = Eigen::Vector3d(colorIn[0].asDouble(), colorIn[1].asDouble(), colorIn[2].asDouble());
        }
        else {
            objects[i].color = c;
        }
    }

    auto gridIn = root["grid"]; 
    {
 	    h = gridIn.get("h", 1.0).asDouble();
        auto originIn = gridIn["origin"];
        if (originIn.size() != 2) {
		  origin = Vector2d(0.0,0.0);
        } 
        else {
            origin = Vector2d(originIn[0].asDouble(),originIn[1].asDouble());
        }
		origin(0) += h/2.0; origin(1) += h/2.0;
        auto sizeIn = gridIn["size"];
        if (sizeIn.size() != 2) {
            std::cout<< "bad grid size, exiting" << std::endl;
			exit(-1);
        } 
        else {
            res[0] = sizeIn[0].asInt();
            res[1] = sizeIn[1].asInt();
        }
    }
  
    auto gravityIn = root["gravity"];
    if (gravityIn.isNull() || gravityIn.size() != 2) {
        std::cout<< "no gravity" << std::endl;
        gravityEnabled = false;
    }
    else {
        gravity(0)= gravityIn[0].asDouble();
        gravity(1)= gravityIn[1].asDouble();
        gravityEnabled = true;
    }

    auto rotationIn = root["rotate"];
    if (rotationIn.isNull()) {
        rotationEnabled = false;
        rotation = 0;
    } 
    else {
        rotationEnabled = true;
        rotation = rotationIn.asDouble();
    }
    
    origin(0) += h/2.0; origin(1) += h/2.0;
    

    for(size_t o = 0; o < objects.size(); o++) {
        Object& obj = objects[o];
        if(obj.type == "square") {
            center = obj.object + (Vector2d(obj.size[0],obj.size[1]) * 0.5);
            //Set up particles at each object vertex
            double diffx = obj.size[0] / (obj.ores[0]-1);
            double diffy = obj.size[1] / (obj.ores[1]-1);
            for(int i = 0; i < obj.ores[0]; i++) {
                for(int j = 0; j < obj.ores[1]; j++) {
                    Vector2d pos = obj.object + Vector2d(diffx*i, diffy*j);
                    Vector3d col = ((double)j/(obj.ores[1]-1))*obj.color;
                    Particle par(pos, Vector2d(0,0), col, obj.mp.pmass);
                    obj.particles.push_back(par);
                }
            }
            if(obj.mp.mass > 0) {
                double partMass = obj.mp.mass / obj.particles.size();
                for(size_t i = 0; i < obj.particles.size(); i++) {
                    obj.particles[i].m = partMass;
                }
            }
            else {
                obj.mp.mass = obj.mp.pmass * obj.particles.size();
            }
        }
        else if(obj.type == "circle") {
            center = obj.object;
            //non-randomly make a circle
            //technically this is an ellipse because I'm using the same data as the
            //square and just using the size vector as radius of the semi-major and semi-minor axes
            //just make a square and reject those outside the ellipse.
            //~78.5% of the resx*resy particles are accepted - Pi/4 * (l*w) particles 
            double diffx = 2*obj.size[0] / (obj.ores[0]-1);
            double diffy = 2*obj.size[1] / (obj.ores[1]-1);
            for(int i = 0; i < obj.ores[0]; i++) {
                for(int j = 0; j < obj.ores[1]; j++) {
                    Vector2d pos = center - Vector2d(obj.size[0], obj.size[1]) + Vector2d(diffx*i, diffy*j);
                    Vector3d col = ((double)j/(obj.ores[1]-1))*obj.color;
                    Vector2d ph = pos - obj.object;
                    if( ((ph(0)*ph(0))/(obj.size[0]*obj.size[0])) + ((ph(1)*ph(1))/(obj.size[1]*obj.size[1])) < 1+EPS) {
                        Particle par(pos, Vector2d(0,0), col, obj.mp.pmass);
                        obj.particles.push_back(par);
                    }
                }
            }
            if(obj.mp.mass > 0) {
                double partMass = obj.mp.mass / obj.particles.size();
                for(size_t i = 0; i < obj.particles.size(); i++) {
                    obj.particles[i].m = partMass;
                }
            }
            else {
                obj.mp.mass = obj.mp.pmass * obj.particles.size();
            }
        }
        else {
            obj.mp.pmass = obj.particles[0].m;
            obj.mp.mass = obj.mp.pmass * obj.particles.size();
        }
    }
  
    //Average position for center of mass
    for(size_t i = 0; i < objects.size(); i++) {
        Vector2d avePos = Vector2d::Zero();
        for(size_t j = 0; j < objects[i].particles.size(); j++) {
            avePos += objects[i].particles[j].x;
            objects[i].particles[j].c1 = objects[i].particles[j].color;
            objects[i].particles[j].c2 = Vector3d(0.0, 0.0, 1.0);
        }
        avePos /= objects[i].particles.size();
        objects[i].center = avePos;
        
        //initialize D
        objects[i].D = new Matrix2d[objects[i].particles.size()];
    }
  
    printf("Grid:\n");
    printf("Dimentions: %dx%d\n", res[0], res[1]);
    printf("Grid Spacing: %f\n", h);
    printf("X0: (%f, %f)\n", origin(0), origin(1));
    printf("X1: (%f, %f)\n", origin(0)+h*(res[0]-1), origin(1)+h*(res[1]-1));
    
    printf("\nConstants:\n");
    printf("Total Time: %f\n", totalTime);
    printf("dt: %f\n", dt);
    printf("Gravity: (%f, %f)\n", gravity(0), gravity(1));
    printf("Rotation: %f\n", rotation);
    
    for(size_t i = 0; i < objects.size(); i++) {
        printf("\nObject %d\n", (int)i);
        printf("Type: %s\n", objects[i].type.c_str());
        printf("Number of particles: %d\n", (int)objects[i].particles.size());
        printf("Particle Mass: %f\n", objects[i].particles[0].m);
        printf("Object Mass: %f\n", objects[i].mp.mass);
        printf("Lame Constants: %f, %f\n", objects[i].mp.lambda, objects[i].mp.mu);
        printf("Center of Mass: (%f, %f)\n", objects[i].center(0), objects[i].center(1));
        printf("Color: (%f, %f, %f)\n", objects[i].color(0), objects[i].color(1), objects[i].color(2));
    }
    printf("\n");

	// allocate grid values
    mass = new double[res[0]*res[1]];
    vel = new Vector2d[res[0]*res[1]];
    velStar = new Vector2d[res[0]*res[1]];
    frc = new Vector2d[res[0]*res[1]];
    
    /// #ifdef INFO
    /// polar.open("/nfs/scratch/adahl1/disc/polar.txt");
    /// kinetic.open("/nfs/scratch/adahl1/disc/kinetic.txt");
    /// #endif
}

inline double weight(double x) {
    double ax = std::abs(x);
    double x2 = x*x;
    double ax3 = x2*ax;
    if(ax < 1.0) {
        return 0.5*ax3 - x2 + (2.0/3.0);
    } else if(ax < 2.0) {
        return (-1.0/6.0)*ax3 + x2 - 2.0*ax + (4.0/3.0);
    }
    return 0;
}

inline double sign(double val) {
    return (0.0 < val) - (val < 0.0);
}

inline double gradweight1d(double x) {
    double ax = std::abs(x);
    double xax = x*ax;
    if(ax < 1.0) {
        return 1.5*xax - 2.0*x;
    } else if(ax < 2.0) {
        // x/ax is just the sign of x
        /// return (-0.5)*xax + 2.0*x - 2.0*(x/ax);
        return (-0.5)*xax + 2.0*x - 2*sign(x);
    }
    return 0;
}

inline Vector2d gradweight(const Vector2d &offset, double h) {
	return Vector2d(gradweight1d(offset(0))*weight(offset(1))/h, 
                    weight(offset(0))*gradweight1d(offset(1))/h);
}

//Keep the max and min? It should never reach past the ends since we're 
//putting an invisible wall around the last 2 cells in gridToParticles
inline void bounds(const Vector2d &offset, const int res[2], int *xbounds, int *ybounds) {
    /// xbounds[0] = std::max(0, ((int)std::ceil(BOUNDLOWER + offset(0))));
    /// xbounds[1] = std::min(res[0]-1, ((int)std::floor(BOUNDUPPER + offset(0)))+1);
    /// ybounds[0] = std::max(0, ((int)std::ceil(BOUNDLOWER + offset(1))));
    /// ybounds[1] = std::min(res[1]-1, ((int)std::floor(BOUNDUPPER + offset(1)))+1);
    xbounds[0] = ((int)(-2 + offset(0)))+1;
    xbounds[1] = ((int)( 2 + offset(0)))+1;
    ybounds[0] = ((int)(-2 + offset(1)))+1;
    ybounds[1] = ((int)( 2 + offset(1)))+1;
}

void World::init() {
    particleVolumesDensities();
    
    /******
     * Init Springs
     ******/
    int n = 10; //number of neighbors
    double ks = 300; //spring coeff
    double kd = 6; //damping coeff
    for(int o = 0; o < (int)objects.size(); o++) {
        std::vector<Particle> &parts = objects[o].particles;
        for(int i = 0; i < (int)parts.size(); i++) {
            Particle *p = &parts[i];
            //find n closest neighbors
            std::vector<Particle*> neighbors;
            std::vector<double> dists;
            for(int j = 0; j < (int)parts.size(); j++) {
                if(i == j) {
                    continue;
                }
                Particle *p2 = &parts[j];
                //Get Distance
                double dsq = (p->x - p2->x).squaredNorm();
                //Check if distance is smaller than anything existing
                for(int k = 0; k < (int)neighbors.size(); k++) {
                    if(dsq < dists[k]) {
                        neighbors.insert(neighbors.begin()+k, p2);
                        dists.insert(dists.begin()+k, dsq);
                    }
                    if((int)neighbors.size() > n) {
                        neighbors.pop_back();
                        dists.pop_back();
                    }
                }
            }
            //check if neighbors already have spring connected to this particle
            std::vector<Particle*> newNeigh;
            for(int j = 0; j < (int)neighbors.size(); j++) {
                std::vector<Spring> &springs = neighbors[j]->springs;
                if(springs.size() == 0) {
                    newNeigh.push_back(neighbors[j]);
                }
                bool exist = false;
                for(int k = 0; k < (int)springs.size(); k++) {
                    if(p->x == springs[k].p0->x || p->x == springs[k].p1->x) {
                        exist = true;
                    }
                }
                if(!exist) {
                    newNeigh.push_back(neighbors[j]);
                }
            }
            //otherwise, add it to the list
            for(int j = 0; j < (int)newNeigh.size(); j++) {
                Spring sp;
                sp.p0 = p;
                sp.p1 = newNeigh[j];
                sp.r = (p->x - newNeigh[j]->x).norm();
                sp.ks = ks;
                sp.kd = kd;
                p->springs.push_back(sp);
            }
        }
    }
    
    
    //Bootstrap on an APIC transfer to the grid and back
    #ifdef INFO
    for(int o = 0; o < (int)objects.size(); o++) {
        for(int i = 0; i < (int)objects[o].particles.size(); i++) {
            objects[o].particles[i].vold = objects[o].particles[i].v;
            objects[o].particles[i].xo = objects[o].particles[i].x;
        }
        //Do an APIC transfer to grid
        particlesToGrid();
        //Check grid values to see if grid position is recovered
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int index = i*res[1]+j;
                velStar[index] = vel[index];
            }
        }
        
        //Do an APIC transfer from grid
        gridToParticles();
        //Check particle values (should still be x,y)
        for(int i = 0; i < (int)objects[o].particles.size(); i++) {
            Vector2d diff = objects[o].particles[i].v - objects[o].particles[i].vold;
            /// double diffNorm = diff.norm();
            /// debug << "Old: (" << objects[o].particles[i].vold(0) << ", " << objects[o].particles[i].vold(1) << ")\n";
            /// debug << "New: (" << objects[o].particles[i].v(0) << ", " << objects[o].particles[i].v(1) << ")\n";
            /// debug << objects[o].particles[i].x(0) << " " << objects[o].particles[i].x(1) << " " << diffNorm << "\n";
            xdiff << objects[o].particles[i].x(0) << " " << objects[o].particles[i].x(1) << " " << diff(0) << "\n";
            ydiff << objects[o].particles[i].x(0) << " " << objects[o].particles[i].x(1) << " " << diff(1) << "\n";
            /// debug.flush();
            xdiff.flush();
            ydiff.flush();
        }
    }
    /// debug.close();
    xdiff.close();
    ydiff.close();
    /// std::exit(0);
    #endif
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
	for (unsigned int obj = 0; obj<objects.size(); obj++) {
	  std::vector<Particle> &particles = objects[obj].particles;
	  for(size_t i = 0; i < particles.size(); i++) {
        Particle &p = particles[i];
		p.rho = 0.0;
		Vector2d offset = (p.x - origin) / h;
		int xbounds[2], ybounds[2];
		bounds(offset, res, xbounds, ybounds);
        for(int j = xbounds[0]; j < xbounds[1]; j++) {
            double w1 = weight(offset(0) - j);
            for(int k = ybounds[0]; k < ybounds[1]; k++) {
                int index = j*res[1] + k;
                double w = w1 * weight(offset(1) - k);
                double r = mass[index] * w;
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
void World::step() {
    particlesToGrid();
    computeGridForces();
    updateGridVelocities();
    updateGradient();
    gridToParticles();
    #ifdef INFO
    if(stepNum % 200000 == 0) {
        /// double ang = 0.75*elapsedTime*M_PI_2;
        /// debug << "Frame " << stepNum / 50000 << " Angle: " << ang << ", ";
        /// double err = 0;
        /// Rotation2D<double> rad(ang);
        /// Matrix2d rotM = rad.toRotationMatrix();
        /// std::vector<Particle> &par = objects[0].particles;
        /// for(int i = 0; i < par.size(); i++) {
            /// Particle &p = par[i];
            /// Vector2d xc = p.xo - objects[0].center;
            /// err += ((rotM*xc+objects[0].center)-p.x).norm();
        /// }
        /// debug << err << "\n";
        //Check condition number
        //Loop through particles, find eigenvalues, divide larger by smaller
        //Output worst
        double maxCond = 0;
        for(int i = 0; i < (int)objects[0].particles.size(); i++) {
            Particle &p = objects[0].particles[i];
            //Eigenvalues
            /// EigenSolver<Matrix2d> es(p.gradientE*p.gradientP, false);
            /// double e1 = es.eigenvalues().real()(0);
            /// double e2 = es.eigenvalues().real()(1);
            /// std::cout << es.eigenvalues() << "\n" << e1 << "\n" << e2 << "\n";
            /// std::exit(0);
            Matrix2d grad = p.gradientE*p.gradientP;
            // [w x]
            // [y z]
            double w, x, y, z;
            w = grad(0,0);
            x = grad(0,1);
            y = grad(1,0);
            z = grad(1,1);
            //solve lambda^2 - (w+z)*lambda + (wz-xy) = 0
            double a, b, c;
            a = 1;
            b = -(w+z);
            c = w*z-x*y;
            double e1, e2;
            e1 = (-b + std::sqrt(b*b - 4*a*c)) / (2*a);
            e2 = (-b - std::sqrt(b*b - 4*a*c)) / (2*a);
            double cnum = 0;
            if(e1 > e2) {
                cnum = e1 / e2;
            }
            else {
                cnum = e2 / e1;
            }
            if(cnum > maxCond) {
                maxCond = cnum;
            }
        }
        cond << stepNum / 200000 << " " << maxCond << "\n";
        cond.flush();
    }
    #endif
    stepNum++;
	elapsedTime += dt;
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
    auto timer = prof.timeName("particlesToGrid");
    {Vector2d *v = vel; for (int i=0; i<res[0]*res[1]; i++, v++) (*v) = Vector2d(0.0,0.0);}
    {double *m = mass; for (int i=0; i<res[0]*res[1]; i++, m++) (*m) = 0.0;}
    
    Matrix2d tmpD = ((h*h)/3.0)*Matrix2d::Identity();
    Matrix2d tmpDinv = (3.0/(h*h))*Matrix2d::Identity();
    
    {auto timer2 = prof.timeName("ptg Object Loop");
	for (unsigned int obj = 0; obj<objects.size(); obj++) {
        std::vector<Particle> &particles = objects[obj].particles;
        Matrix2d *Dtensor = objects[obj].D;
        {auto timer3 = prof.timeName("ptg Particles Loop");
        for(size_t i = 0; i < particles.size(); i++) {
            Particle &p = particles[i];
            if(stepNum == 0) {
                Matrix2d C;
                C << 0, -0.75, 0.75, 0; //Rotational (rotation=1.0) Case
                /// C << 1, 0, 0, 1;    //Linear (x,y) Case
                p.B = C * tmpD;
            }
            Vector2d offset = (p.x - origin) / h;
            int xbounds[2], ybounds[2];
            bounds(offset, res, xbounds, ybounds);
            Vector2d xp = p.x;                                      //particle position
            #ifndef NDEBUG
            if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                if(i == 0) {
                    debug << "Particle " << i << "\nB:\n" << p.B << "\n";
                    debug << "X Bounds 1: " << xbounds[0] << ", " << xbounds[1] << "\n";
                    debug << "Y Bounds 1: " << ybounds[0] << ", " << ybounds[1] << "\n";
                }
            }
            #endif
            {auto timer4 = prof.timeName("ptg Bounded Loop 1");                     
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                double w1 = weight(offset(0) - j);
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    int index = j*res[1] + k;
                    double w = w1*weight(offset(1) - k);
                    mass[index] += w * p.m;
                    #ifndef NDEBUG
                    if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                        if(i == 0) {
                            debug << "Adding to mass at (" << j << "," << k << "): " << w*p.m << "\n";
                        }
                    }
                    #endif
                }
            }}
            {auto timer5 = prof.timeName("ptg Bounded Loop 2"); 
            #ifndef NDEBUG
            if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                if(i == 0) {
                    debug << "D:\n" << Dtensor[i] << "\n";
                }
            }
            #endif
            Matrix2d Dinv = tmpDinv;  
            Dtensor[i] = Dinv;
            Vector2d mv = p.m*p.v;
            Matrix2d mBD = p.m*p.B*Dinv;
            #ifndef NDEBUG
            if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                if(i == 0) {
                    debug << "Dinv:\n" << Dinv << "\n";
                    debug << "mv:\n" << mv << "\n";
                    debug << "mBD:\n" << mBD << "\n";
                }
            }
            #endif
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                double w1 = weight(offset(0) - j);
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    double w = w1*weight(offset(1) - k);
                    Vector2d xg = origin + h*Vector2d(j, k);
                    vel[j*res[1] + k] += w * (mv + mBD*(xg-xp).eval());
                    #ifndef NDEBUG
                    if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                        if(i == 0) {
                            debug << "Adding to vel at (" << j << "," << k << "):\n" << w*(mv+mBD*(xg-xp)) << "\n";
                        }
                    }
                    #endif
                }
            }}
        }}
        /// #ifdef INFO
        /// if(stepNum%5000 == 0) {
            /// double rotDet = 0;
            /// double angVel = 0;
            /// double inertia = 0, omega = 0;
            /// for(int i = 0; i < (int)particles.size(); i++) {
                /// Particle &p = particles[i];
                /// inertia += p.m * (p.x-objects[0].center).squaredNorm();
            /// }
            /// for(int i = 0; i < (int)particles.size(); i++) {
                /// Particle &p = particles[i];
                /// Matrix2d F = p.gradientE*p.gradientP;
                /// Matrix2d FT = F.transpose();
                /// //Decompose F = RS
                /// Matrix2d S = (FT*F).sqrt();
                /// /// Matrix2d R = F*S.inverse();
                /// //Want to look at forbeneous norm of symmetric part. So don't need rotation
                /// rotDet += (S-Matrix2d::Identity()).norm();
                /// omega += p.v.norm() / (objects[0].center-p.x).norm();
            /// }
            /// rotDet /= particles.size();
            /// omega /= particles.size();
            /// angVel = 0.5 * inertia * omega*omega;
            /// polar << elapsedTime << " " << rotDet << "\n";
            /// kinetic << elapsedTime << " " << angVel << "\n";
        /// }
        /// #endif
	}}
    #ifndef NDEBUG
    debug << "Mass:\n";
    for(int i = 0; i < res[0]; i++) {
		for(int j = 0; j < res[1]; j++) {
            debug << mass[i*res[1]+j] << " ";
        }
        debug << "\n";
    }
    debug << "\nMomentum:\n";
    for(int i = 0; i < res[0]; i++) {
		for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            debug << "(";
            if(vel[index](0) >= 0) {
                debug.precision(10);
            }
            else {
                debug.precision(9);
            }
            debug << vel[index](0) << ",";
            if(vel[index](1) >= 0) {
                debug.precision(10);
            }
            else {
                debug.precision(9);
            }
            debug << vel[index](1) << ") ";
        }
        debug << "\n";
    }
    #endif
    {auto timer6 = prof.timeName("ptg Vel Loop"); 
	for(int i = 0; i < res[0]; i++) {
		for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            if(mass[index] < EPS) {
                vel[index] = Vector2d(0.0, 0.0);
            }
            else {
                vel[index] /= mass[index];
            }
            #ifndef NDEBUG
            if(vel[index].hasNaN()) {
                printf("interpolated vel NaN at (%d, %d)\n", i, j);
                std::cout << vel[index] << std::endl;
                exit(0);
            }
            #endif
        }
	}}
    #ifndef NDEBUG
    debug << "\nVel:\n";
    for(int i = 0; i < res[0]; i++) {
		for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            debug << "(";
            if(vel[index](0) >= 0) {
                debug.precision(10);
            }
            else {
                debug.precision(9);
            }
            debug << vel[index](0) << ",";
            if(vel[index](1) >= 0) {
                debug.precision(10);
            }
            else {
                debug.precision(9);
            }
            debug << vel[index](1) << ") ";
        }
        debug << "\n";
    }
    #endif
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
    auto timer = prof.timeName("computeGridForces");
    {Vector2d *f = frc; for (int i=0; i<res[0]*res[1]; i++, f++) (*f) = Vector2d(0.0,0.0);}
	{for(int i = 0; i < (int)objects[0].particles.size(); i++) objects[0].particles[i].f = Vector2d::Zero();}
    
    {auto timer2 = prof.timeName("cgf Object Loop"); 
    for (unsigned int obj = 0; obj<objects.size(); obj++) {
        std::vector<Particle> &particles = objects[obj].particles;
        MaterialProps &mp = objects[obj].mp;
        {auto timer = prof.timeName("cgf Particles Loop"); 
        //*
        for(size_t i = 0; i < particles.size(); i++) {
            Particle &p = particles[i];
            for(int j = 0; i < p.springs.size(); j++) {
                Spring sp = p.springs[j];
                Vector2d I = sp.p0->x - sp.p1->x;
                Vector2d II = sp.p0->v - sp.p1->v;
                double inorm = I.norm();
                double idot = II.dot(I);
                Vector2d fa = -(sp.ks*(inorm-sp.r) + sp.kd*(idot/inorm)) * (I/inorm);
                sp.p0->f += fa;
                sp.p1->f -= fa;
            }
        }
        for(size_t i = 0; i < particles.size(); i++) {
            Particle &p = particles[i];
            Vector2d offset = (p.x - origin) / h;
            int xbounds[2], ybounds[2];
            bounds(offset, res, xbounds, ybounds);
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                double w1 = weight(offset(0) - j);
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    double w = w1*weight(offset(1) - k);
                    frc[j*res[1]+k] += w * p.f;
                }
            }
        }}
        //*/
        /*
        for(int i = 0; i < (int)particles.size(); i++) {
            Particle &p = particles[i];
            Matrix2d gradient = p.gradientE*p.gradientP; 
            /// Matrix2d gradient = p.gradientE;
            double J = gradient.determinant();
            if(J < 0) {
                /// printf("Negative Determinant: %f\n", J);
                /// std::cout << "Gradient:\n" << gradient << "\n";
                /// #ifdef INFO
                /// debug.close();
                /// xdiff.close();
                /// ydiff.close();
                /// #endif
                /// #ifndef NDEBUG
                /// debug.close();
                /// #endif
                /// exit(1);
                p.color = p.c2;
            }
            else {
                p.color = p.c1;
            }
            Matrix2d gradT = gradient.transpose();
            Matrix2d eps = 0.5 * (gradT * gradient - Matrix2d::Identity());
            double trace = eps.trace();
            Matrix2d stress = mp.lambda*trace*Matrix2d::Identity() + 2.0*mp.mu*eps;
            Matrix2d Ftmp = p.vol * J * stress;
            
            Vector2d offset = (p.x - origin) / h;
            int xbounds[2], ybounds[2];
            bounds(offset, res, xbounds, ybounds);
            #ifndef NDEBUG
            if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                if(i == 0) {
                    debug << "Gradient:\n" << gradient << "\n";
                    debug << "J: " << J << "\n";
                    debug << "Eps:\n" << eps << "\n";
                    debug << "Trace: " << trace << "\n";
                    debug << "Stress:\n" << stress << "\n";
                    debug << "Ftmp:\n" << Ftmp << "\n";
                    debug << "X Bounds 2: " << xbounds[0] << ", " << xbounds[1] << "\n";
                    debug << "Y Bounds 2: " << ybounds[0] << ", " << ybounds[1] << "\n";
                }
            }
            #endif
            {auto timer = prof.timeName("cgf Bound Loop"); 
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    Vector2d accumF = Ftmp * gradweight(Vector2d(offset(0)-j,offset(1)-k),h);
                    #ifndef NDEBUG
                    if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                        if(i == 0) {
                            debug << "Subbing to F:\n" << accumF << "\n";
                        }
                    }
                    #endif
                    frc[j*res[1] + k] -= accumF;
                    #ifndef NDEBUG
                    if(frc[j*res[1] + k].hasNaN()) {
                        printf("\nf NaN at (%d, %d)\n", j, k);
                        std::cout << "Force:\n" << frc[j*res[1] + k] << std::endl;
                        std::cout << "Volume: " << p.vol << std::endl;
                        std::cout << "Determinant: " << gradient.determinant() << std::endl;
                        std::cout << "Stress:\n" << stress << std::endl;
                        std::cout << "Gradient:\n" << gradweight(Vector2d(offset(0)-j, offset(1)-k),h) << std::endl;
                        exit(0);
                    }
                    #endif
                }
            }}
        }}
        */
	}}
}

/******************************
 * Update_Grid_Velocities
 *      for each grid node i do
 *          v^*_i = v_i^n + (dt * f_i)/m_i + dt * g
 *      end for
 *****************************/
void World::updateGridVelocities() {
    auto timer = prof.timeName("updateGridVelocities");
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1] + j;
            if(mass[index] < EPS) {
                velStar[index] = Vector2d(0, 0);
            }
            else {
                Vector2d extfrc(0.0,0.0);
                if(gravityEnabled) {
                    extfrc += mass[index]*gravity;
                }
                if(rotationEnabled) {
                    Vector2d d = origin+Vector2d(h*i,h*j)-center;
                    extfrc += mass[index]*rotation*Vector2d(-d(1), d(0));
                }
                velStar[index] = vel[index] + dt * (1.0/mass[index]) * (frc[index] + extfrc); //dt*g
            }
            #ifndef NDEBUG
            if(velStar[index].hasNaN()) {
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
void World::updateGradient() {
    auto timer = prof.timeName("updateGradient");
    for (unsigned int obj=0; obj<objects.size(); obj++) {
        std::vector<Particle> &particles = objects[obj].particles;
        MaterialProps &mp = objects[obj].mp;
        /// Matrix2d *D = objects[obj].D;

        {auto timer = prof.timeName("ug Particles Loop"); 
        for(size_t i = 0; i < particles.size(); i++) {
            Particle &p = particles[i];
            Matrix2d gradV = Matrix2d::Zero();
            Vector2d offset = (p.x - origin) / h;
            int xbounds[2], ybounds[2];
            bounds(offset, res, xbounds, ybounds);
            #ifndef NDEBUG
            if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                if(i == 0) {
                    debug << "X Bounds 3: " << xbounds[0] << ", " << xbounds[1] << "\n";
                    debug << "Y Bounds 3: " << ybounds[0] << ", " << ybounds[1] << "\n";
                }
            }
            #endif
            {auto timer = prof.timeName("ug Bound Loop"); 
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    int index = j*res[1]+k;
                    #ifndef NDEBUG
                    if(velStar[index].hasNaN()) {
                        printf("gradV velStar has NaN at (%d, %d)\n", j, k);
                        std::cout << velStar[index] << std::endl;
                        exit(0);
                    }
                    if(gradweight(Vector2d(offset(0)-j,offset(1)-k),h).transpose().hasNaN()) {
                        printf("gradV gradW has NaN at (%d, %d)\n", j, k);
                        std::cout << gradweight(Vector2d(offset(0)-j,offset(1)-k),h).transpose() << std::endl;
                        exit(0);
                    }
                    #endif
                    Matrix2d accumGrad = velStar[index] * gradweight(Vector2d(offset(0)-j,offset(1)-k),h).transpose();
                    gradV += accumGrad;
                    #ifndef NDEBUG
                    if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                        if(i == 0) {
                            debug << "Index: " << index << "\n";
                            debug << "VelStar:\n" << velStar[index] << "\n";
                            debug << "accumGrad:\n" << accumGrad << "\n";
                        }
                    }
                    #endif
                }
            }}
            #ifndef NDEBUG
            if(gradV.hasNaN()) {
                printf("gradV has NaN\n");
                std::cout << gradV << std::endl;
                exit(0);
            }
            #endif
            {auto timer = prof.timeName("ug Grad Update"); 
            /// Matrix2d C = p.B * D[i];
            /// Matrix2d fp = Matrix2d::Identity() + dt*C;
            Matrix2d fp = Matrix2d::Identity() + dt*gradV;
            
            Matrix2d tempGradE = fp*p.gradientE;
            #ifndef NDEBUG
            if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                if(i == 0) {
                    debug << "gradV:\n" << gradV << "\n";
                    debug << "C:\n" << C << "\n";
                    debug << "fp:\n" << fp << "\n";
                    debug << "Grad:\n" << tempGradE << "\n";
                }
            }
            #endif
            
            if (plasticEnabled){
                Matrix2d tempGrad = fp*p.gradientE*p.gradientP;
                
                JacobiSVD<Matrix2d> svd(tempGradE, ComputeFullU | ComputeFullV);
                
                Matrix2d svdU = svd.matrixU();
                Vector2d svdSV = svd.singularValues();
                Matrix2d svdV = svd.matrixV();
                
                Vector2d sVClamped;
                sVClamped << clamp(svdSV(0), 1-mp.compression, 1+mp.stretch), clamp(svdSV(1), 1-mp.compression, 1+mp.stretch);
                Matrix2d svdClamped = sVClamped.asDiagonal();
                
                p.gradientE = svdU * svdClamped * svdV.transpose();
                p.gradientP = svdV * svdClamped.inverse() * svdU.transpose() * tempGrad;
            } else {
                p.gradientE = tempGradE;
                p.gradientP = Matrix2d::Identity();
            }
            }
        }}
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
void World::gridToParticles() {
    auto timer = prof.timeName("gridToParticles");
    for (unsigned int obj=0; obj<objects.size(); obj++) {
        std::vector<Particle> &particles = objects[obj].particles;
        MaterialProps &mp = objects[obj].mp;
        {auto timer = prof.timeName("gtp Particles Loop"); 
        for(size_t i = 0; i < particles.size(); i++) {
            Particle &p = particles[i];
            //Update velocities
            p.B = Matrix2d::Zero();
            Vector2d apic = Vector2d::Zero();
            Vector2d offset = (p.x - origin) / h;
            int xbounds[2], ybounds[2];
            bounds(offset, res, xbounds, ybounds);
            #ifndef NDEBUG
            if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                if(i == 0) {
                    debug << "X Bounds 4: " << xbounds[0] << ", " << xbounds[1] << "\n";
                    debug << "Y Bounds 4: " << ybounds[0] << ", " << ybounds[1] << "\n";
                }
            }
            #endif
            {auto timer = prof.timeName("gtp Bound Loop"); 
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                double w1 = weight(offset(0) - j);
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    double w = w1*weight(offset(1) - k);
                    int index = j*res[1] + k;
                    Vector2d xg = origin + h*Vector2d(j, k);
                    Vector2d wvel = w * velStar[index];
                    apic += wvel;
                    p.B += wvel * (xg - p.x).eval().transpose().eval();
                    #ifndef NDEBUG
                    if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                        if(i == 0) {
                            debug << "Adding to apic:\n" << wvel << "\n";
                            debug << "Adding to B:\n" << wvel*(xg-p.x).transpose() << "\n";
                        }
                    }
                    #endif
                }
            }}
            {auto timer = prof.timeName("gtp Vel and Position Update"); 
            #ifndef NDEBUG
            if(apic.hasNaN()) {
                printf("\n\nAPIC Vel has NaN\n");
                std::cout << apic << std::endl;
                exit(0);
            }
            #endif
            p.v = apic;
            //Mass proportional damping
            p.v = mp.massPropDamp * p.v;
            #ifndef NDEBUG
            if(stepNum >= 50000*0 && stepNum < 50000*20 && stepNum%5000 == 0) {
                if(i == 0) {
                    debug << "APIC:\n" << apic << "\n";
                    debug << "Pos:\n" << p.x << "\n";
                    debug << "Adding:\n" << dt*apic << "\n";
                    debug << "New Pos:\n" << p.x+dt*p.v << "\n\n\n";
                }
            }
            #endif
            
            #ifndef NDEBUG
            if(p.v.hasNaN()) {
                printf("Vel has NaN\n");
                std::cout << p.v << std::endl;
                exit(0);
            }
            #endif

            //Update Positions
            p.x += dt * p.v;
            
            //Overriding particle velocities with expected values
            //ph = pos-center
            //expected value is rotation* (-ph(1), ph(0))
            /// Vector2d ph = p.x - objects[obj].center;
            /// p.v = 0.75 * Vector2d(-ph(1), ph(0));
            
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
            double ux = origin[0]+(res[0]-3)*h; //The last index is res-1 so it's -3
            double uy = origin[1]+(res[1]-3)*h; 
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
        }}
    }
}


void writeParticles(const char *fname, const std::vector<Particle> &particles) {
	Partio::ParticlesDataMutable *data = Partio::create();
	Partio::ParticleAttribute xattr;
	Partio::ParticleAttribute uattr;
	Partio::ParticleAttribute battr;
	Partio::ParticleAttribute geattr;
	Partio::ParticleAttribute gpattr;
	Partio::ParticleAttribute cattr;
	Partio::ParticleAttribute mattr;
	Partio::ParticleAttribute rattr;
	Partio::ParticleAttribute vattr;
	data->addParticles(particles.size());

	data->addAttribute("position", Partio::VECTOR, 2);
	data->addAttribute("velocity", Partio::VECTOR, 2);
	data->addAttribute("B", Partio::VECTOR, 4);
	data->addAttribute("gradientE", Partio::VECTOR, 4);
	data->addAttribute("gradientP", Partio::VECTOR, 4);
	data->addAttribute("color", Partio::VECTOR, 3);
	data->addAttribute("mass", Partio::FLOAT, 1);
	data->addAttribute("rho", Partio::FLOAT, 1);
	data->addAttribute("vol", Partio::FLOAT, 1);
	
	data->attributeInfo("position", xattr);
	data->attributeInfo("velocity", uattr);
	data->attributeInfo("B", battr);
	data->attributeInfo("gradientE", geattr);
	data->attributeInfo("gradientP", gpattr);
	data->attributeInfo("color", cattr);
	data->attributeInfo("mass", mattr);
	data->attributeInfo("rho", rattr);
	data->attributeInfo("vol", vattr);

	for (unsigned int i=0; i < particles.size(); i++) {
		const Particle &p = particles[i];
		float *x = data->dataWrite<float>(xattr, i);
		float *u = data->dataWrite<float>(uattr, i);
		float *b = data->dataWrite<float>(battr, i);
		float *ge = data->dataWrite<float>(geattr, i);
		float *gp = data->dataWrite<float>(gpattr, i);
		float *c = data->dataWrite<float>(cattr, i);
		float *m = data->dataWrite<float>(mattr, i);
		float *r = data->dataWrite<float>(rattr, i);
		float *v = data->dataWrite<float>(vattr, i);

		x[0] = p.x(0), x[1] = p.x(1);
        u[0] = p.v(0), u[1] = p.v(1);
		b[0] = p.B(0,0), b[1] = p.B(0,1), b[2] = p.B(1,0), b[3] = p.B(1,1);
		ge[0] = p.gradientE(0,0), ge[1] = p.gradientE(0,1), ge[2] = p.gradientE(1,0), ge[3] = p.gradientE(1,1);
		gp[0] = p.gradientP(0,0), gp[1] = p.gradientP(0,1), gp[2] = p.gradientP(1,0), gp[3] = p.gradientP(1,1);
		c[0] = p.color(0), c[1] = p.color(1), c[2] = p.color(2);
		m[0] = p.m;
		r[0] = p.rho;
		v[0] = p.vol;
	}

	Partio::write(fname, *data);
	data->release();
}


bool readParticles(const char *fname, std::vector<Particle> &particles) {
	Partio::ParticlesDataMutable *data = Partio::read(fname);
	if (data == 0) return 0;
	Partio::ParticleAttribute xattr;
	Partio::ParticleAttribute uattr;
	Partio::ParticleAttribute battr;
	Partio::ParticleAttribute geattr;
	Partio::ParticleAttribute gpattr;
	Partio::ParticleAttribute cattr;
	Partio::ParticleAttribute mattr;
	Partio::ParticleAttribute rattr;
	Partio::ParticleAttribute vattr;

	bool position = data->attributeInfo("position", xattr);
	bool velocity = data->attributeInfo("velocity", uattr);
	bool B = data->attributeInfo("B", battr);
	bool gradientE = data->attributeInfo("gradientE", geattr);
	bool gradientP = data->attributeInfo("gradientP", gpattr);
	bool color = data->attributeInfo("color", cattr);
	bool mass = data->attributeInfo("mass", mattr);
	bool rho = data->attributeInfo("rho", rattr);
	bool vol = data->attributeInfo("vol", vattr);

	particles.resize(data->numParticles());

	for (int i=0; i < data->numParticles(); i++) {
        Particle &p = particles[i];
        if (position) {
            float *x = data->dataWrite<float>(xattr, i);
            p.x[0] = x[0], p.x[1] = x[1];
        } else {
            p.x = Vector2d(0.0, 0.0);
        }
        if (velocity) {
            float *u = data->dataWrite<float>(uattr, i);
            p.v[0] = u[0], p.v[1] = u[1];
        } else {
            p.v = Vector2d(0.0, 0.0);
		}
        if (B) {
            float *b = data->dataWrite<float>(battr, i);
            p.B(0,0) = b[0], p.B(0,1) = b[1], p.B(1,0) = b[2], p.B(1,1) = b[3];
        } else {
            p.B = Matrix2d::Zero();
        }
        if (gradientE) {
            float *g = data->dataWrite<float>(geattr, i);
            p.gradientE(0,0) = g[0], p.gradientE(0,1) = g[1], p.gradientE(1,0) = g[2], p.gradientE(1,1) = g[3];
        } else {
            p.gradientE = Matrix2d::Identity();
        }
        if (gradientP) {
            float *g = data->dataWrite<float>(gpattr, i);
            p.gradientP(0,0) = g[0], p.gradientP(0,1) = g[1], p.gradientP(1,0) = g[2], p.gradientP(1,1) = g[3];
        } else {
            p.gradientP = Matrix2d::Identity();
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
	return true;
}

World::~World() {
    delete [] mass;
    delete [] vel;
    delete [] velStar;
    delete [] frc;
    prof.dump<std::chrono::duration<double>>(std::cout);
} 
