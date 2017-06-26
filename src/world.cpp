#include "world.hpp"
#include <cstdio>
#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>
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
std::ofstream xdiff("xdiff.txt");
std::ofstream ydiff("ydiff.txt");
#endif

//calc.cs.umbc.edu. Nothing in home directory.
//ssh to cal[01-12], has access to data

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
    prevVel = new Vector2d[res[0]*res[1]];
    frc = new Vector2d[res[0]*res[1]];
    mat = new Vector2d[res[0]*res[1]];
    weights = new double[res[0]*res[1]];
    stress = new Matrix2d[res[0]*res[1]]; 
    valid.assign(res[0]*res[1], 0);
    
    {Vector2d *v = vel; for (int i=0; i<res[0]*res[1]; i++, v++) (*v) = Vector2d::Zero();}
    {Vector2d *c = mat; for (int i=0; i<res[0]*res[1]; i++, c++) (*c) = Vector2d::Zero();}
    {Vector2d *vs = velStar; for (int i=0; i<res[0]*res[1]; i++, vs++) (*vs) = Vector2d::Zero();}
    {Vector2d *vo = prevVel; for (int i=0; i<res[0]*res[1]; i++, vo++) (*vo) = Vector2d::Zero();}
    {Vector2d *f = frc; for (int i=0; i<res[0]*res[1]; i++, f++) (*f) = Vector2d::Zero();}
    {Matrix2d *s = stress; for (int i=0; i<res[0]*res[1]; i++, s++) (*s) = Matrix2d::Zero();}
    {double *m = mass; for (int i=0; i<res[0]*res[1]; i++, m++) (*m) = 0.0;}
    {double *w = weights; for (int i=0; i<res[0]*res[1]; i++, w++) (*w) = 0.0;}
    
    #ifndef NDEBUG
    count = 0;
    sMax = 15;
    matTrans = new Vector2d*[sMax];
    for(int i = 0; i < sMax; i++) {
        matTrans[i] = new Vector2d[res[0]*res[1]];
    }
    #endif
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

#ifdef LAGRANGE
    /******
     * Init Springs
     ******/
    int n = 10; //number of neighbors
    double ks = 300; //spring coeff
    double kd = 6; //damping coeff
    for(int o = 0; o < (int)objects.size(); o++) {
        std::vector<Particle> &parts = objects[o].particles;
        std::vector<Spring> &springs = objects[o].springs;
        std::vector<int> counts(parts.size(), 0);
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
                if((int)neighbors.size() < n) {
                    int index;
                    for(index = 0; index < (int)neighbors.size(); index++) {
                        if(dsq < dists[index]) {
                            neighbors.insert(neighbors.begin()+index, p2);
                            dists.insert(dists.begin()+index, dsq);
                            break;
                        }
                    }
                    if(index == (int)neighbors.size()) {
                        neighbors.push_back(p2);
                        dists.push_back(dsq);
                    }
                    continue;
                }
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
                bool exist = false;
                for(int k = 0; k < (int)springs.size(); k++) {
                    if((p->x == springs[k].p0->x || p->x == springs[k].p1->x) && (neighbors[j]->x == springs[k].p0->x || neighbors[j]->x == springs[k].p1->x)) {
                        exist = true;
                        break;
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
                springs.push_back(sp);
            }
        }
        printf("Springs: %d\n", (int)springs.size());
        for(int i = 0; i < (int)springs.size(); i++) {
            Particle* p0 = springs[i].p0;
            Particle* p1 = springs[i].p1;
            for(int j = 0; j < (int)parts.size(); j++) {
                if(parts[j].x == p0->x || parts[j].x == p1->x) {
                    counts[j]++;
                }
            }
        }
    }
#endif //LAGRANGE
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
    
    //Init material coordinates
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int ind = i*res[1]+j;
            Vector2d xg = origin + h*Vector2d(i,j);
            mat[ind] = xg;
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
    /// #ifndef NDEBUG
    /// if(stepNum % 13333 == 0) {
        /// for(int i = 0; i < res[0]; i++) {
            /// for(int j = 0; j < res[1]; j++) {
                /// matTrans[count][i*res[1]+j] = mat[i*res[1]+j];
            /// }
        /// }
        /// count++;
        /// if(count == sMax) {
            /// for(int i = 0; i < res[0]; i++) {
                /// for(int j = 0; j < res[1]; j++) {
                    /// for(int k = 0; k < sMax; k++) {
                        /// debug << matTrans[k][i*res[1]+j](0) << " " << matTrans[k][i*res[1]+j](1) << " " << k << "\n";
                    /// }
                    /// debug << "\n\n";
                /// }
            /// }
            /// debug.close();
            /// std::exit(0);
        /// }
    /// }
    /// #endif
    particlesToGrid();
    computeGridForces();
    updateGridVelocities();
#ifndef MAT_TRANSFER
    updateGradient();
#else
    //Material coordinates are initialized to world coordinates
    //Extrapolate velocity field with fast march
        //March out ~5 times
    velExtrapolate();
    //Advect material coordinates
        //interpolate velocities at v[t] and v[t-dt] and average
        //step back by that amount and interpolate material coordinates around that point
    advect();
#endif
    gridToParticles();
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
    /// {Vector2d *c = mat; for (int i=0; i<res[0]*res[1]; i++, c++) (*c) = Vector2d(0.0,0.0);}
    {double *m = mass; for (int i=0; i<res[0]*res[1]; i++, m++) (*m) = 0.0;}
    {double *w = weights; for (int i=0; i<res[0]*res[1]; i++, w++) (*w) = 0.0;}
    
    Matrix2d tmpD = ((h*h)/3.0)*Matrix2d::Identity();
    Matrix2d tmpDinv = (3.0/(h*h))*Matrix2d::Identity();
    
	for (unsigned int obj = 0; obj<objects.size(); obj++) {
        std::vector<Particle> &particles = objects[obj].particles;
        Matrix2d *Dtensor = objects[obj].D;
        for(size_t i = 0; i < particles.size(); i++) {
            Particle &p = particles[i];
            if(stepNum == 0) {
                Matrix2d C;
                C << 0, -0.75, 0.75, 0; //Rotational (rotation=0.75) Case
                /// C << 0, 0, 0, 0;
                p.B = C * tmpD;
            }
            Vector2d offset = (p.x - origin) / h;
            int xbounds[2], ybounds[2];
            bounds(offset, res, xbounds, ybounds);
            Vector2d xp = p.x;                                      //particle position                   
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                double w1 = weight(offset(0) - j);
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    int index = j*res[1] + k;
                    double w = w1*weight(offset(1) - k);
                    mass[index] += w * p.m;
                    /// mat[index] += w * p.u;
                    /// weights[index] += w;
                }
            }
            Matrix2d Dinv = tmpDinv;  
            Dtensor[i] = Dinv;
            Vector2d mv = p.m*p.v;
            Matrix2d mBD = p.m*p.B*Dinv;
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                double w1 = weight(offset(0) - j);
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    double w = w1*weight(offset(1) - k);
                    Vector2d xg = origin + h*Vector2d(j, k);
                    vel[j*res[1] + k] += w * (mv + mBD*(xg-xp).eval());
                }
            }
        }
	}
	for(int i = 0; i < res[0]; i++) {
		for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            if(mass[index] < EPS) {
                vel[index] = Vector2d(0.0, 0.0);
                valid[index] = 0;
            }
            else {
                vel[index] /= mass[index];
                valid[index] = 1;
            }
            /// if(weights[index] < EPS) {
                /// mat[index] = Vector2d::Zero();
            /// }
            /// else {
                /// mat[index] /= weights[index];
            /// }
            #ifndef NDEBUG
            if(vel[index].hasNaN()) {
                printf("interpolated vel NaN at (%d, %d)\n", i, j);
                std::cout << vel[index] << std::endl;
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
    auto timer = prof.timeName("computeGridForces");
    {Vector2d *f = frc; for (int i=0; i<res[0]*res[1]; i++, f++) (*f) = Vector2d(0.0,0.0);}
#ifdef LAGRANGE
	{for(int i = 0; i < (int)objects[0].particles.size(); i++) objects[0].particles[i].f = Vector2d::Zero();}
#endif
    
#ifdef MAT_TRANSFER
    #ifndef NDEBUG
    Matrix2d *dgrad, *dF, *dFT, *dA, *dS, *dSt, *dstrain;
    double *gradDet = new double[res[0]*res[1]];
    double *cond = new double[res[0]*res[1]];
    dgrad = new Matrix2d[res[0]*res[1]];
    dF = new Matrix2d[res[0]*res[1]];
    dFT = new Matrix2d[res[0]*res[1]];
    dA = new Matrix2d[res[0]*res[1]];
    dS = new Matrix2d[res[0]*res[1]];
    dSt = new Matrix2d[res[0]*res[1]];
    dstrain = new Matrix2d[res[0]*res[1]];
    #endif
    //Do finite difference to make the deformation gradient at each grid point
    for(int obj = 0; obj < (int)objects.size(); obj++) {
        MaterialProps &mp = objects[obj].mp;
        /// #ifndef NDEBUG
        /// debug << "Mat:\n";
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int i = 0; i < res[0]; i++) {
                /// int index = i*res[1]+j;
                /// debug << "(";
                /// if(mat[index](0) < 0) {
                    /// debug << std::fixed << std::setprecision(5) << mat[index](0);
                /// }
                /// else {
                    /// debug << std::fixed << std::setprecision(6) << mat[index](0);
                /// }
                /// debug << ",";
                /// if(mat[index](1) < 0) {
                    /// debug << std::fixed << std::setprecision(5) << mat[index](1);
                /// }
                /// else {
                    /// debug << std::fixed << std::setprecision(6) << mat[index](1);
                /// }
                /// debug << ") ";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
        /// #endif
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int index = i*res[1]+j;
                if(!valid[index]) {
                    stress[index] = Matrix2d::Zero();
                    continue;
                }
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);
                Vector2d dx, dy;
                //First order central difference (x_n+1 - x_n-1) / 2h
                if(j == 0) {    //forward
                    dy = (mat[yp1] - mat[index]) / h;
                }
                else if(j == res[1]-1) {    //backward
                    dy = (mat[index] - mat[ym1]) / h;
                }
                else {  //center
                    dy = (mat[yp1] - mat[ym1]) / (2*h);
                }
                /// if(j == 0 || (!valid[ym1] && valid[yp1])) {    //forward
                    /// dy = (mat[yp1] - mat[index]) / h;
                /// }
                /// else if(j == res[1]-1 || (!valid[yp1] && valid[ym1])) {    //backward
                    /// dy = (mat[index] - mat[ym1]) / h;
                /// }
                /// else if(valid[yp1] && valid[ym1]) {  //center
                    /// dy = (mat[yp1] - mat[ym1]) / (2*h);
                /// }
                /// else {
                    /// dy = Vector2d(0, 1);
                /// }
                
                if(i == 0) {
                    dx = (mat[xp1] - mat[index]) / h;
                }
                else if(i == res[0]-1) {
                    dx = (mat[index] - mat[xm1]) / h;
                }
                else {
                    dx = (mat[xp1] - mat[xm1]) / (2*h);
                }
                /// if(i == 0 || (!valid[xm1] && valid[xp1])) {
                    /// dx = (mat[xp1] - mat[index]) / h;
                /// }
                /// else if(i == res[0]-1 || (!valid[xp1] && valid[xm1])) {
                    /// dx = (mat[index] - mat[xm1]) / h;
                /// }
                /// else if(valid[xp1] && valid[xm1]) {
                    /// dx = (mat[xp1] - mat[xm1]) / (2*h);
                /// }
                /// else {
                    /// dx = Vector2d(1, 0);
                /// }
                
                //Form inverse deformation gradient [dx dy]
                Matrix2d grad;
                grad.col(0) = dx;
                grad.col(1) = dy;
                #ifndef NDEBUG
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
                if(b*b - 4*a*c < 0) {
                    cnum = 0;
                }
                if(grad.hasNaN()) {
                    printf("\n%d, %d\n", i, j);
                    printf("dx: (%f, %f), (%f, %f)\n", dx(0), dx(1), grad.col(0)(0), grad.col(0)(1));
                    printf("dy: (%f, %f), (%f, %f)\n", dy(0), dy(1), grad.col(1)(0), grad.col(1)(1));
                    fflush(stdout);
                    std::exit(1);
                }
                gradDet[index] = grad.determinant();
                cond[index] = cnum;
                #endif  
                double det = grad.determinant();
                if(std::abs(det) < 1e-4) {
                    stress[index] = Matrix2d::Zero();
                    #ifndef NDEBUG
                    dgrad[index] = grad;
                    dF[index] = Matrix2d::Zero();
                    dFT[index] = Matrix2d::Zero();
                    dA[index] = Matrix2d::Zero();
                    dS[index] = Matrix2d::Zero();
                    dSt[index] = Matrix2d::Zero();
                    dstrain[index] = Matrix2d::Zero();
                    #endif
                }
                else {
                    Matrix2d F = grad.inverse();
                    #ifndef NDEBUG
                    if(F.hasNaN()) {
                        printf("%d, %d\n", i, j);
                        printf("Grad:\n");
                        std::cout << grad << "\n" << F << "\n";
                        fflush(stdout);
                        std::exit(1);
                    }
                    #endif
                    Matrix2d FT = F.transpose();
                    Matrix2d A = FT*F;
                    Matrix2d S = A.sqrt();
                    
                    /// EigenSolver<Matrix2d> es(A);
                    /// Matrix2d Dsqrt = es.eigenvalues().cwiseSqrt().asDiagonal();
                    /// Matrix2d S = es.eigenvectors() * Dsqrt * es.eigenvectors().inverse();
                    
                    /// double s = std::sqrt(A.determinant());
                    /// double t = std::sqrt(A.trace()+2*s);
                    /// Matrix2d S = A;
                    /// S(0,0) += s;
                    /// S(1,1) += s;
                    /// S /= t;
                    
                    #ifndef NDEBUG
                    if(S.hasNaN()) {
                        printf("S: %d, %d\n", i, j);
                        printf("[%f, %f]\n", S(0,0), S(0,1));
                        printf("[%f, %f]\n", S(1,0), S(1,1));
                        printf("F =[%f, %f]\n", F(0,0), F(0,1));
                        printf("   [%f, %f]\n", F(1,0), F(1,1));
                        printf("Ft=[%f, %f]\n", FT(0,0), FT(0,1));
                        printf("   [%f, %f]\n", FT(1,0), FT(1,1));
                        fflush(stdout);
                        std::exit(1);
                    }
                    #endif
                    Matrix2d St = S.transpose();
                    Matrix2d strain = 0.5 * (S + St) - Matrix2d::Identity();
                    #ifndef NDEBUG
                    if(strain.hasNaN()) {
                        printf("%d, %d\n", i, j);
                        printf("0.5 * ([%.5f,%.5f] + [%.5f,%.5f]) - I\n", S(0,0), S(0,1), St(0,0), St(0,1));
                        printf("      ([%.5f,%.5f]   [%.5f,%.5f])\n", S(1,0), S(1,1), St(1,0), St(1,1));
                        std::cout << strain << "\n";
                        fflush(stdout);
                        std::exit(1);
                    }
                    #endif
                    Matrix2d linstress = mp.lambda*strain.trace()*Matrix2d::Identity() + 2.0*mp.mu*strain;
                    #ifndef NDEBUG
                    if(linstress.hasNaN()) {
                        printf("%d, %d\n", i, j);
                        printf("lambda, mu: %f, %f\n", mp.lambda, mp.mu);
                        printf("trace: %f\n", strain.trace());
                        std::cout << strain << "\n";
                        std::cout << linstress << "\n";
                        fflush(stdout);
                        std::exit(1);
                    }
                    #endif
                    stress[index] = linstress;
                    
                    #ifndef NDEBUG
                    dgrad[index] = grad;
                    dF[index] = F;
                    dFT[index] = FT;
                    dA[index] = A;
                    dS[index] = S;
                    dSt[index] = St;
                    dstrain[index] = strain;
                    #endif
                }
            }
        }
        //Try smoothing stress field
        /// int numSmooth = 1; //just 1 pass for now
        /// Matrix2d *tmpstress = new Matrix2d[res[0]*res[1]];
        /// int xp1, xm1, yp1, ym1;
        /// for(int s = 0; s < numSmooth; s++) {
            /// for(int i = 0; i < res[0]; i++) {
                /// for(int j = 0; j < res[1]; j++) {
                    /// //Extend kernal past end
                    /// if(i == res[0]-1) {
                        /// xp1 = i;
                    /// }
                    /// else {
                        /// xp1 = i + 1;
                    /// }
                    /// if(i == 0) {
                        /// xm1 = i;
                    /// }
                    /// else {
                        /// xm1 = i - 1;
                    /// }
                    /// if(j == res[1]-1) {
                        /// yp1 = j;
                    /// }
                    /// else {
                        /// yp1 = j + 1;
                    /// }
                    /// if(j == 0) {
                        /// ym1 = j;
                    /// }
                    /// else {
                        /// ym1 = j - 1;
                    /// }
                    /// //Using kernal:
                    /// //[0  1  0]
                    /// //[1 -4  1]
                    /// //[0  1  0]
                    /// /// tmpstress[i*res[1]+j] = (stress[xm1*res[1]+j] + stress[xp1*res[1]+j] + stress[i*res[1]+ym1] + stress[i*res[1]+yp1] - 4*stress[i*res[1]+j]) / (h*h);
                    /// tmpstress[i*res[1]+j] = stress[xm1*res[1]+j] + stress[xp1*res[1]+j] + stress[i*res[1]+ym1] + stress[i*res[1]+yp1] - 4*stress[i*res[1]+j];
                    /// if(tmpstress[i*res[1]+j].hasNaN()) {
                        /// printf("Temp Stress NaN\n");
                        /// std::cout << tmpstress[i*res[1]+j] << "\n";
                        /// std::cout << stress[xm1*res[1]+j] << "\n";
                        /// std::cout << stress[xp1*res[1]+j] << "\n";
                        /// std::cout << stress[i*res[1]+ym1] << "\n";
                        /// std::cout << stress[i*res[1]+yp1] << "\n";
                        /// std::cout << 4*stress[i*res[1]+j] << "\n";
                        /// std::exit(1);
                    /// }
                /// }
            /// }
            /// for(int i = 0; i < res[0]*res[1]; i++) {
                /// stress[i] = tmpstress[i];
            /// }
        /// }
        /// delete[] tmpstress;
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int index = i*res[1]+j;
                if(!valid[index]) {
                    frc[index] = Vector2d::Zero();
                    continue;
                }
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);
                Vector2d fx, fy;
                /// if(j == 0) {    //forward
                    /// fy = (stress[yp1].row(1) - stress[index].row(1)) / h;
                    /// #ifndef NDEBUG
                    /// if(fy.hasNaN()) {
                        /// printf("forward fy: %d, %d\n", i, j);
                        /// std::cout << fy << "\n";
                        /// std::cout << stress[yp1] << "\n";
                        /// std::cout << stress[index] << "\n";
                        /// fflush(stdout);
                        /// std::exit(1);
                    /// }
                    /// #endif
                /// }
                /// else if(j == res[1]-1) {    //backward
                    /// fy = (stress[index].row(1) - stress[ym1].row(1)) / h;
                    /// #ifndef NDEBUG
                    /// if(fy.hasNaN()) {
                        /// printf("backward fy: %d, %d\n", i, j);
                        /// std::cout << fy << "\n";
                        /// std::cout << stress[index] << "\n";
                        /// std::cout << stress[ym1] << "\n";
                        /// fflush(stdout);
                        /// std::exit(1);
                    /// }
                    /// #endif
                /// }
                /// else {  //center
                    /// fy = (stress[yp1].row(1) - stress[ym1].row(1)) / (2*h);
                    /// #ifndef NDEBUG
                    /// if(fy.hasNaN()) {
                        /// printf("center fy: %d, %d\n", i, j);
                        /// std::cout << fy << "\n";
                        /// std::cout << stress[yp1] << "\n";
                        /// std::cout << stress[ym1] << "\n";
                        /// fflush(stdout);
                        /// std::exit(1);
                    /// }
                    /// #endif
                /// }
                if(j == 0 || (!valid[ym1] && valid[yp1])) {    //forward
                    fy = (stress[yp1].row(1) - stress[index].row(1)) / h;
                }
                else if(j == res[1]-1 || (!valid[yp1] && valid[ym1])) {    //backward
                    fy = (stress[index].row(1) - stress[ym1].row(1)) / h;
                }
                else if(valid[yp1] && valid[ym1]) {  //center
                    fy = (stress[yp1].row(1) - stress[ym1].row(1)) / (2*h);
                }
                else {
                    fy = Vector2d(0,0);
                }
                
                /// if(i == 0) {
                    /// fx = (stress[xp1].row(0) - stress[index].row(0)) / h;
                    /// #ifndef NDEBUG
                    /// if(fx.hasNaN()) {
                        /// printf("forward fx: %d, %d\n", i, j);
                        /// std::cout << fx << "\n";
                        /// std::cout << stress[xp1] << "\n";
                        /// std::cout << stress[index] << "\n";
                        /// fflush(stdout);
                        /// std::exit(1);
                    /// }
                    /// #endif
                /// }
                /// else if(i == res[0]-1) {
                    /// fx = (stress[index].row(0) - stress[xm1].row(0)) / h;
                    /// #ifndef NDEBUG
                    /// if(fx.hasNaN()) {
                        /// printf("backward fx: %d, %d\n", i, j);
                        /// std::cout << fx << "\n";
                        /// std::cout << stress[index] << "\n";
                        /// std::cout << stress[xm1] << "\n";
                        /// fflush(stdout);
                        /// std::exit(1);
                    /// }
                    /// #endif
                /// }
                /// else {
                    /// fx = (stress[xp1].row(0) - stress[xm1].row(0)) / (2*h);
                    /// #ifndef NDEBUG
                    /// if(fx.hasNaN()) {
                        /// printf("center fx: %d, %d\n", i, j);
                        /// std::cout << fx << "\n";
                        /// std::cout << stress[xp1] << "\n";
                        /// std::cout << stress[xm1] << "\n";
                        /// fflush(stdout);
                        /// std::exit(1);
                    /// }
                    /// #endif
                /// }
                if(i == 0 || (!valid[xm1] && valid[xp1])) {
                    fx = (stress[xp1].row(0) - stress[index].row(0)) / h;
                }
                else if(i == res[0]-1 || (!valid[xp1] && valid[xm1])) {
                    fx = (stress[index].row(0) - stress[xm1].row(0)) / h;
                }
                else if(valid[xp1] && valid[xm1]) {
                    fx = (stress[xp1].row(0) - stress[xm1].row(0)) / (2*h);
                }
                else {
                    fx = Vector2d(0,0);
                }
                
                Vector2d force(fx(0)+fy(0), fx(1)+fy(1));
                #ifndef NDEBUG
                if(force.hasNaN()) {
                    printf("%d, %d\n", i, j);
                    printf("(%f, %f) = (%f+%f, %f+%f)\n", force(0), force(0), fx(0), fx(1), fy(0), fy(1));
                    printf("fx: ((%f,%f) - (%f,%f)) / (2*%f)\n", stress[(i+1)*res[1]+j].col(0)(0), stress[(i+1)*res[1]+j].col(0)(1), stress[(i-1)*res[1]+j].col(0)(0), stress[(i-1)*res[1]+j].col(0)(1), h);
                    printf("fy: ((%f,%f) - (%f,%f)) / (2*%f)\n", stress[i*res[1]+(j+1)].col(1)(0), stress[i*res[1]+(j+1)].col(1)(1), stress[i*res[1]+(j-1)].col(1)(0), stress[i*res[1]+(j-1)].col(1)(1), h);
                    fflush(stdout);
                    std::exit(1);
                }
                #endif
                frc[i*res[1]+j] = force;
            }
        }
        #ifndef NDEBUG
        /// debug << "Grad:\n";
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int k = 0; k < 2; k++) {
                /// for(int i = 0; i < res[0]; i++) {
                    /// int ind = i*res[0]+j;
                    /// debug << "[";
                    /// if(dgrad[ind](k,0) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dgrad[ind](k,0);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dgrad[ind](k,0);
                    /// }
                    /// debug << ",";
                    /// if(dgrad[ind](k,1) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dgrad[ind](k,1);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dgrad[ind](k,1);
                    /// }
                    /// debug << "] ";
                /// }
                /// debug << "\n";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
        
        /// debug << "Determinant:\n";
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int i = 0; i < res[0]; i++) {
                /// int ind = i*res[0]+j;
                /// debug << gradDet[ind] << " ";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
        
        /// debug << "Condition:\n";
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int i = 0; i < res[0]; i++) {
                /// int ind = i*res[0]+j;
                /// debug << cond[ind] << " ";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
        
        /// debug << "F:\n";
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int k = 0; k < 2; k++) {
                /// for(int i = 0; i < res[0]; i++) {
                    /// int ind = i*res[0]+j;
                    /// debug << "[";
                    /// if(dF[ind](k,0) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dF[ind](k,0);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dF[ind](k,0);
                    /// }
                    /// debug << ",";
                    /// if(dF[ind](k,1) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dF[ind](k,1);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dF[ind](k,1);
                    /// }
                    /// debug << "] ";
                /// }
                /// debug << "\n";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
        
        /// debug << "FT:\n";
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int k = 0; k < 2; k++) {
                /// for(int i = 0; i < res[0]; i++) {
                    /// int ind = i*res[0]+j;
                    /// debug << "[";
                    /// if(dFT[ind](k,0) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dFT[ind](k,0);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dFT[ind](k,0);
                    /// }
                    /// debug << ",";
                    /// if(dFT[ind](k,1) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dFT[ind](k,1);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dFT[ind](k,1);
                    /// }
                    /// debug << "] ";
                /// }
                /// debug << "\n";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
        
        /// debug << "A:\n";
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int k = 0; k < 2; k++) {
                /// for(int i = 0; i < res[0]; i++) {
                    /// int ind = i*res[0]+j;
                    /// debug << "[";
                    /// if(dA[ind](k,0) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dA[ind](k,0);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dA[ind](k,0);
                    /// }
                    /// debug << ",";
                    /// if(dA[ind](k,1) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dA[ind](k,1);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dA[ind](k,1);
                    /// }
                    /// debug << "] ";
                /// }
                /// debug << "\n";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
        
        /// debug << "S:\n";
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int k = 0; k < 2; k++) {
                /// for(int i = 0; i < res[0]; i++) {
                    /// int ind = i*res[0]+j;
                    /// debug << "[";
                    /// if(dS[ind](k,0) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dS[ind](k,0);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dS[ind](k,0);
                    /// }
                    /// debug << ",";
                    /// if(dS[ind](k,1) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dS[ind](k,1);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dS[ind](k,1);
                    /// }
                    /// debug << "] ";
                /// }
                /// debug << "\n";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
        
        /// debug << "St:\n";
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int k = 0; k < 2; k++) {
                /// for(int i = 0; i < res[0]; i++) {
                    /// int ind = i*res[0]+j;
                    /// debug << "[";
                    /// if(dSt[ind](k,0) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dSt[ind](k,0);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dSt[ind](k,0);
                    /// }
                    /// debug << ",";
                    /// if(dSt[ind](k,1) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dSt[ind](k,1);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dSt[ind](k,1);
                    /// }
                    /// debug << "] ";
                /// }
                /// debug << "\n";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
        
        /// debug << "Strain:\n";
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int k = 0; k < 2; k++) {
                /// for(int i = 0; i < res[0]; i++) {
                    /// int ind = i*res[0]+j;
                    /// debug << "[";
                    /// if(dstrain[ind](k,0) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dstrain[ind](k,0);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dstrain[ind](k,0);
                    /// }
                    /// debug << ",";
                    /// if(dstrain[ind](k,1) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << dstrain[ind](k,1);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << dstrain[ind](k,1);
                    /// }
                    /// debug << "] ";
                /// }
                /// debug << "\n";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
        
        /// debug << "Stress:\n";
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int k = 0; k < 2; k++) {
                /// for(int i = 0; i < res[0]; i++) {
                    /// int ind = i*res[0]+j;
                    /// debug << "[";
                    /// if(stress[ind](k,0) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << stress[ind](k,0);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << stress[ind](k,0);
                    /// }
                    /// debug << ",";
                    /// if(stress[ind](k,1) < 0) {
                        /// debug << std::fixed << std::setprecision(5) << stress[ind](k,1);
                    /// }
                    /// else {
                        /// debug << std::fixed << std::setprecision(6) << stress[ind](k,1);
                    /// }
                    /// debug << "] ";
                /// }
                /// debug << "\n";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
        
        /// debug << "Force:\n";
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int i = 0; i < res[0]; i++) {
                /// int index = i*res[1]+j;
                /// debug << "(";
                /// if(frc[index](0) < 0) {
                    /// debug << std::fixed << std::setprecision(5) << frc[index](0);
                /// }
                /// else {
                    /// debug << std::fixed << std::setprecision(6) << frc[index](0);
                /// }
                /// debug << ",";
                /// if(frc[index](1) < 0) {
                    /// debug << std::fixed << std::setprecision(5) << frc[index](1);
                /// }
                /// else {
                    /// debug << std::fixed << std::setprecision(6) << frc[index](1);
                /// }
                /// debug << ") ";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n\n\n";
        delete[] dgrad;
        delete[] dF;
        delete[] dFT;
        delete[] dA;
        delete[] dS;
        delete[] dSt;
        delete[] dstrain;
        delete[] gradDet;
        delete[] cond;
        /// if(stepNum == 3) {
            /// std::exit(1);
        /// }
        #endif
    }
#else
    for (unsigned int obj = 0; obj<objects.size(); obj++) {
        std::vector<Particle> &particles = objects[obj].particles;
        MaterialProps &mp = objects[obj].mp;
    #ifdef LAGRANGE
        for(int j = 0; j < (int)objects[obj].springs.size(); j++) {
            Spring sp = objects[obj].springs[j];
            Vector2d I = sp.p0->x - sp.p1->x;
            Vector2d II = sp.p0->v - sp.p1->v;
            double inorm = I.norm();
            double idot = II.dot(I);
            Vector2d fa = -(sp.ks*(inorm/sp.r-1) + sp.kd*(idot/inorm)) * (I/inorm);
            sp.p0->f += fa;
            sp.p1->f -= fa;
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
        }
    #else
        for(int i = 0; i < (int)particles.size(); i++) {
            Particle &p = particles[i];
            Matrix2d gradient = p.gradientE*p.gradientP; 
            double J = gradient.determinant();
            if(J < 0) {
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
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    Vector2d accumF = Ftmp * gradweight(Vector2d(offset(0)-j,offset(1)-k),h);
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
            }
        }
    #endif
	}
#endif
}

/******************************
 * Update_Grid_Velocities
 *      for each grid node i do
 *          v^*_i = v_i^n + (dt * f_i)/m_i + dt * g
 *      end for
 *****************************/
void World::updateGridVelocities() {
    auto timer = prof.timeName("updateGridVelocities");
    /// #ifndef NDEBUG
    /// int it = 15399;
    /// if(stepNum > it) {
        /// debug << "force:\n";
        /// for(int i = 0; i < res[0]; i++) {
            /// for(int j = 0; j < res[1]; j++) {
                /// debug << "(";
                /// if(frc[i*res[1]+j](0) < 0) {
                    /// debug << std::fixed << std::setprecision(5) << frc[i*res[1]+j](0);
                /// } 
                /// else {
                    /// debug << std::fixed << std::setprecision(6) << frc[i*res[1]+j](0);
                /// } 
                /// debug << ",";
                /// if(frc[i*res[1]+j](1) < 0) {
                    /// debug << std::fixed << std::setprecision(5) << frc[i*res[1]+j](1);
                /// }
                /// else {
                    /// debug << std::fixed << std::setprecision(6) << frc[i*res[1]+j](1);
                /// }
                /// debug << ") ";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
    /// }
    /// #endif
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1] + j;
            if(stepNum > 0) {
                prevVel[index] = velStar[index];
                /// #ifndef NDEBUG
                /// if(stepNum > it && i == 31 && j == 19) {
                    /// debug << "previous: (" << prevVel[index](0) << ", " << prevVel[index](1) << ")\n";
                /// }
                /// #endif
            }
            /// if(mass[index] < EPS) {
            if(!valid[index]) {
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
                /// #ifndef NDEBUG
                /// if(stepNum > it && i == 31 && j == 19) {
                    /// debug << std::defaultfloat << "(" << vel[index](0) << "," << vel[index](1) << ") + " << dt << "*(1/" << mass[index] << ")*((" << frc[index](0) << "," << frc[index](1) << ")+(" << extfrc(0) << "," << extfrc(1) << "))\n";
                /// }
                /// #endif
            }
            if(stepNum == 0) {
                prevVel[index] = velStar[index];
            }
            /// #ifndef NDEBUG
            /// if(stepNum > it && i == 31 && j == 19) {
                /// debug << "vel: (" << velStar[index](0) << ", " << velStar[index](1) << ")\n";
            /// }
            /// #endif
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
                    if(gradweight(Vector2d(offset(0)-j,offset(1)-k),h).transpose().hasNaN()) {
                        printf("gradV gradW has NaN at (%d, %d)\n", j, k);
                        std::cout << gradweight(Vector2d(offset(0)-j,offset(1)-k),h).transpose() << std::endl;
                        exit(0);
                    }
                    #endif
                    Matrix2d accumGrad = velStar[index] * gradweight(Vector2d(offset(0)-j,offset(1)-k),h).transpose();
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
            /// Matrix2d C = p.B * D[i];
            /// Matrix2d fp = Matrix2d::Identity() + dt*C;
            Matrix2d fp = Matrix2d::Identity() + dt*gradV;
            
            Matrix2d tempGradE = fp*p.gradientE;
            
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
    }
}

void World::velExtrapolate() {
    /// #ifndef NDEBUG
    /// int it = 1442;
    /// if(stepNum > it) {
        /// debug << "Step " << stepNum << "\n";
        /// debug << "Vel:\n";
        /// for(int i = 0; i < res[0]; i++) {
            /// for(int j = 0; j < res[1]; j++) {
                /// debug << "(";
                /// if(velStar[i*res[1]+j](0) < 0) {
                    /// debug << std::fixed << std::setprecision(5) << velStar[i*res[1]+j](0);
                /// } 
                /// else {
                    /// debug << std::fixed << std::setprecision(6) << velStar[i*res[1]+j](0);
                /// } 
                /// debug << ",";
                /// if(velStar[i*res[1]+j](1) < 0) {
                    /// debug << std::fixed << std::setprecision(5) << velStar[i*res[1]+j](1);
                /// }
                /// else {
                    /// debug << std::fixed << std::setprecision(6) << velStar[i*res[1]+j](1);
                /// }
                /// debug << ") ";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
    /// }
    /// #endif
    //Applying simple extrapolation method from https://github.com/christopherbatty/Fluid3D/blob/master/fluidsim.cpp
    std::vector<char> oldValid;
    oldValid.resize(res[0]*res[1]);
    Vector2d *tmpVel = new Vector2d[res[0]*res[1]];
    //Perform this a couple of times
    for(int layer = 0; layer < 6; layer++) {
        oldValid = valid;
        /// #ifndef NDEBUG
        /// debug << "Iter: " << layer << "\n";
        /// #endif
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int ind = i*res[1]+j;
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);
                /// #ifndef NDEBUG
                /// debug << (int)oldValid[ind] << " ";
                /// #endif
                
                Vector2d sum(0,0);
                int count = 0;
                
                tmpVel[ind] = velStar[ind];
                //Only check around cells that are not already valid
                if(!oldValid[ind]) {
                    if(i+1 < res[0] && oldValid[xp1]) {
                        sum += velStar[xp1];
                        count++;
                    }
                    if(i-1 >= 0 && oldValid[xm1]) {
                        sum += velStar[xm1];
                        count++;
                    }
                    if(j+1 < res[1] && oldValid[yp1]) {
                        sum += velStar[yp1];
                        count++;
                    }
                    if(j-1 >= 0 && oldValid[ym1]) {
                        sum += velStar[ym1];
                        count++;
                    }
                    
                    //If neighboring cells were valid, average the values
                    if(count > 0) {
                        tmpVel[ind] = sum / count;
                        valid[ind] = 1;
                    }
                }
            }
            /// #ifndef NDEBUG
            /// debug << "\n";
            /// #endif
        }
        /// #ifndef NDEBUG
        /// debug << "\n";
        /// #endif
        //Replace grid with temp grid
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                velStar[i*res[1]+j] = tmpVel[i*res[1]+j];
            }
        }
    }
    /// #ifndef NDEBUG
    /// if(stepNum > it) {
        /// debug << "Extrapolated Vel:\n";
        /// for(int i = 0; i < res[0]; i++) {
            /// for(int j = 0; j < res[1]; j++) {
                /// debug << "(";
                /// if(velStar[i*res[1]+j](0) < 0) {
                    /// debug << std::fixed << std::setprecision(5) << velStar[i*res[1]+j](0);
                /// } 
                /// else {
                    /// debug << std::fixed << std::setprecision(6) << velStar[i*res[1]+j](0);
                /// } 
                /// debug << ",";
                /// if(velStar[i*res[1]+j](1) < 0) {
                    /// debug << std::fixed << std::setprecision(5) << velStar[i*res[1]+j](1);
                /// }
                /// else {
                    /// debug << std::fixed << std::setprecision(6) << velStar[i*res[1]+j](1);
                /// }
                /// debug << ") ";
            /// }
            /// debug << "\n";
        /// }
        /// debug << "\n";
        /// debug.close();
        /// std::exit(1);
    /// }
    /// #endif
    /// #ifndef NDEBUG
    /// debug.close();
    /// std::exit(1);
    /// #endif
    
    /// #ifndef NDEBUG
    /// int inc = 20000;
    /// if(stepNum % inc == 0) {
        /// std::string name = std::string("field-") + std::to_string(stepNum/inc);
        /// std::ofstream vf(name.c_str());
        /// double scale = 1;
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int i = 0; i < res[0]; i++) {
                /// Vector2d xg = origin+h*Vector2d(i,j);
                /// Vector2d v = scale*velStar[i*res[1]+j];
                /// Vector2d f = scale*frc[i*res[1]+j];
                /// vf << xg(0) << " " << xg(1) << " " << v(0) << " " << v(1) << "\n";
                /// vf << xg(0) << " " << xg(1) << " " << f(0) << " " << f(1) << "\n";
            /// }
        /// }
        /// vf.close();
        /// if(stepNum == 10*inc) {
            /// std::exit(0);
        /// }
    /// }
    /// #endif
    delete[] tmpVel;
}

void World::advect() {
    Vector2d *newMat = new Vector2d[res[0]*res[1]];
    //Advect material coordinates
        //interpolate velocities at v[t] and v[t-dt] and average
        //step back by that amount and interpolate material coordinates around that point
    /// #ifndef NDEBUG
    /// int it = 1400;
    /// #endif
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            if(!valid[index]) {
                newMat[index] = mat[index];
            }
            Vector2d xg = origin + h*Vector2d(i,j);
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "i, j: " << i << " " << j << "\n";
                /// debug << "xg: (" << xg(0) << ", " << xg(1) << ")\n";
                /// debug.flush();
            /// }
            /// #endif
            //v[t]
            Vector2d v1 = velStar[index];
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "v1: (" << v1(0) << ", " << v1(1) << "\n";
                /// debug.flush();
            /// }
            /// #endif
            //step position back by dt
            Vector2d ox = xg - dt*v1;
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "ox: (" << ox(0) << ", " << ox(1) << ")\n";
                /// debug.flush();
            /// }
            /// #endif
            //v[t-dt]
            Vector2d x = (ox - origin);
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "x: (" << x(0) << ", " << x(1) << ")\n";
                /// debug.flush();
            /// }
            /// #endif
            Vector2d ij = x / h;
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "ij: (" << ij(0) << ", " << ij(1) << ")\n";
                /// debug.flush();
            /// }
            /// #endif
            int i1 = (int)ij(0);
            int i2 = i1+1;
            int j1 = (int)ij(1);
            int j2 = j1+1;
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "i: [" << i1 << ", " << i2 << "]\n";
                /// debug << "j: [" << j1 << ", " << j2 << "]\n";
                /// debug.flush();
            /// }
            /// #endif
            Vector2d Q11 = prevVel[i1*res[1]+j1];
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "Q11: (" << Q11(0) << ", " << Q11(1) << ")\n";
                /// debug.flush();
            /// }
            /// #endif
            Vector2d Q21 = prevVel[i2*res[1]+j1];
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "Q21: (" << Q21(0) << ", " << Q21(1) << ")\n";
                /// debug.flush();
            /// }
            /// #endif
            Vector2d Q12 = prevVel[i1*res[1]+j2];
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "Q12: (" << Q12(0) << ", " << Q12(1) << ")\n";
                /// debug.flush();
            /// }
            /// #endif
            Vector2d Q22 = prevVel[i2*res[1]+j2];
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "Q22: (" << Q22(0) << ", " << Q22(1) << ")\n";
                /// debug.flush();
            /// }
            /// #endif
            double x1 = h*i1;
            double x2 = h*i2;
            double y1 = h*j1;
            double y2 = h*j2;
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "x: [" << x1 << ", " << x2 << "]\n";
                /// debug << "y: [" << y1 << ", " << y2 << "]\n";
                /// debug.flush();
            /// }
            /// #endif
            double xd1 = (x2-x(0))/(x2-x1);
            double xd2 = (x(0)-x1)/(x2-x1);
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "xd: [" << xd1 << ", " << xd2 << "]\n";
                /// debug.flush();
            /// }
            /// #endif
            Vector2d fy1 = xd1*Q11 + xd2*Q21;
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "fy1: (" << fy1(0) << ", " << fy1(1) << ")\n";
                /// debug.flush();
            /// }
            /// #endif
            Vector2d fy2 = xd1*Q12 + xd2*Q22;
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "fy2: (" << fy2(0) << ", " << fy2(1) << ")\n";
                /// debug.flush();
            /// }
            /// #endif
            Vector2d v2 = ((y2-x(1))/(y2-y1))*fy1 + ((x(1)-y1)/(y2-y1))*fy2;    
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "v2: (" << v2(0) << ", " << v2(1) << ")\n";
                /// debug.flush();
            /// }
            /// #endif    
            //New Vel
            Vector2d nv = (v1 + v2) / 2;
            /// #ifndef NDEBUG
            /// if(stepNum > it) {
                /// debug << "nv: (" << nv(0) << ", " << nv(1) << ")\n\n";
                /// debug.flush();
            /// }
            /// #endif
                
            //Step back and interpolate material coordinates
            ox = xg - dt*nv;
            x = (ox - origin);
            ij = x / h;
            i1 = (int)ij(0);
            i2 = i1+1;
            j1 = (int)ij(1);
            j2 = j1+1;
            Q11 = mat[i1*res[1]+j1];
            Q21 = mat[i2*res[1]+j1];
            Q12 = mat[i1*res[1]+j2];
            Q22 = mat[i2*res[1]+j2];
            x1 = h*i1;
            x2 = h*i2;
            y1 = h*j1;
            y2 = h*j2;
            xd1 = (x2-x(0))/(x2-x1);
            xd2 = (x(0)-x1)/(x2-x1);
            fy1 = xd1*Q11 + xd2*Q21;
            fy2 = xd1*Q12 + xd2*Q22;
            Vector2d nmat = ((y2-x(1))/(y2-y1))*fy1 + ((x(1)-y1)/(y2-y1))*fy2;   
            newMat[index] = nmat;   
        }
    }
    /// #ifndef NDEBUG
    /// if(stepNum > it) {
        /// debug << "\n\n";
        /// debug.close();
        /// std::exit(1);
    /// }
    /// #endif
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            mat[i*res[1]+j] = newMat[i*res[1]+j];
        }
    }
    delete[] newMat;
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
        for(size_t i = 0; i < particles.size(); i++) {
            Particle &p = particles[i];
            //Update velocities
            p.B = Matrix2d::Zero();
            Vector2d apic = Vector2d::Zero();
            Vector2d offset = (p.x - origin) / h;
            int xbounds[2], ybounds[2];
            bounds(offset, res, xbounds, ybounds);
            /// Matrix2d Bs, Be;
            /// Bs = Matrix2d::Zero();
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                double w1 = weight(offset(0) - j);
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    double w = w1*weight(offset(1) - k);
                    int index = j*res[1] + k;
                    Vector2d xg = origin + h*Vector2d(j, k);
                    Vector2d wvel = w * velStar[index];
                    apic += wvel;
                    p.B += wvel * (xg - p.x).transpose();
                    /// Bs += wvel * (xg - p.x).transpose();
                }
            }
            #ifndef NDEBUG
            if(apic.hasNaN()) {
                printf("\n\nAPIC Vel has NaN\n");
                std::cout << apic << std::endl;
                exit(0);
            }
            #endif
            p.v = apic;
            //For testing, do RK2 (trapezoidal) time integration. Use bilinear(2D) interpolation for values
            //Bilinear for starting value
            /// Vector2d x = (p.x - origin);
            /// Vector2d ij = x / h;
            /// int i1 = (int)ij(0);
            /// int i2 = i1+1;
            /// int j1 = (int)ij(1);
            /// int j2 = j1+1;
            /// Vector2d Q11 = velStar[i1*res[1]+j1];
            /// Vector2d Q21 = velStar[i2*res[1]+j1];
            /// Vector2d Q12 = velStar[i1*res[1]+j2];
            /// Vector2d Q22 = velStar[i2*res[1]+j2];
            /// double x1 = h*i1;
            /// double x2 = h*i2;
            /// double y1 = h*j1;
            /// double y2 = h*j2;
            /// double xd1 = (x2-x(0))/(x2-x1);
            /// double xd2 = (x(0)-x1)/(x2-x1);
            /// Vector2d fy1 = xd1*Q11 + xd2*Q21;
            /// Vector2d fy2 = xd1*Q12 + xd2*Q22;
            /// Vector2d vs = ((y2-x(1))/(y2-y1))*fy1 + ((x(1)-y1)/(y2-y1))*fy2;
            /// Vector2d vs = apic;
            //Do a temperary timestep for candidate position (trying backwards)
            /// Vector2d xn = p.x + dt * vs;
            /// Vector2d xn = p.x - dt * vs;
            //Bilinear for ending value
            /// x = (xn - origin);
            /// ij = x / h;
            /// i1 = (int)ij(0);
            /// i2 = i1+1;
            /// j1 = (int)ij(1);
            /// j2 = j1+1;
            /// Q11 = prevVel[i1*res[1]+j1];
            /// Q21 = prevVel[i2*res[1]+j1];
            /// Q12 = prevVel[i1*res[1]+j2];
            /// Q22 = prevVel[i2*res[1]+j2];
            /// x1 = h*i1;
            /// x2 = h*i2;
            /// y1 = h*j1;
            /// y2 = h*j2;
            /// xd1 = (x2-x(0))/(x2-x1);
            /// xd2 = (x(0)-x1)/(x2-x1);
            /// fy1 = xd1*Q11 + xd2*Q21;
            /// fy2 = xd1*Q12 + xd2*Q22;
            /// Vector2d ve = ((y2-x(1))/(y2-y1))*fy1 + ((x(1)-y1)/(y2-y1))*fy2;
            
            /// apic = Vector2d::Zero();
            /// offset = (xn - origin) / h;
            /// bounds(offset, res, xbounds, ybounds);
            /// Be = Matrix2d::Zero(); 
            /// for(int j = xbounds[0]; j < xbounds[1]; j++) {
                /// double w1 = weight(offset(0) - j);
                /// for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    /// double w = w1*weight(offset(1) - k);
                    /// int index = j*res[1] + k;
                    /// Vector2d xg = origin + h*Vector2d(j, k);
                    /// Vector2d wvel = w * prevVel[index];
                    /// apic += wvel;
                    /// Be += wvel * (xg - xn).transpose();
                /// }
            /// }
            /// Vector2d ve = apic;
            //Average and set velocity to that
            /// p.v = (vs + ve) / 2;
            /// p.B = (Bs + Be) / 2;
            //Mass proportional damping
            p.v = mp.massPropDamp * p.v;
            
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
    Partio::ParticleAttribute xoattr;
    data->addParticles(particles.size());

    data->addAttribute("rest", Partio::VECTOR, 2);
	data->addAttribute("position", Partio::VECTOR, 3);
	data->addAttribute("velocity", Partio::VECTOR, 2);
	data->addAttribute("B", Partio::VECTOR, 4);
	data->addAttribute("gradientE", Partio::VECTOR, 4);
	data->addAttribute("gradientP", Partio::VECTOR, 4);
	data->addAttribute("color", Partio::VECTOR, 3);
	data->addAttribute("mass", Partio::FLOAT, 1);
	data->addAttribute("rho", Partio::FLOAT, 1);
	data->addAttribute("vol", Partio::FLOAT, 1);
	
    data->attributeInfo("rest", xoattr);
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
        float *x0 = data->dataWrite<float>(xoattr, i);
		float *x = data->dataWrite<float>(xattr, i);
		float *u = data->dataWrite<float>(uattr, i);
		float *b = data->dataWrite<float>(battr, i);
		float *ge = data->dataWrite<float>(geattr, i);
		float *gp = data->dataWrite<float>(gpattr, i);
		float *c = data->dataWrite<float>(cattr, i);
		float *m = data->dataWrite<float>(mattr, i);
		float *r = data->dataWrite<float>(rattr, i);
		float *v = data->dataWrite<float>(vattr, i);

        x0[0] = p.u(0), x0[1] = p.u(1);
		x[0] = p.x(0), x[1] = p.x(1), x[2] = 0;
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
    Partio::ParticleAttribute xoattr;
	Partio::ParticleAttribute xattr;
	Partio::ParticleAttribute uattr;
	Partio::ParticleAttribute battr;
	Partio::ParticleAttribute geattr;
	Partio::ParticleAttribute gpattr;
	Partio::ParticleAttribute cattr;
	Partio::ParticleAttribute mattr;
	Partio::ParticleAttribute rattr;
	Partio::ParticleAttribute vattr;

    bool rest = data->attributeInfo("rest", xoattr);
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
            p.x(0) = x[0], p.x(1) = x[1];
        } else {
            p.x = Vector2d(0.0, 0.0);
        }
        if(rest) {
            float *x0 = data->dataWrite<float>(xoattr, i);
            p.u(0) = x0[0], p.u(1) = x0[1];
        } else {
            p.u = p.x;
        }
        if (velocity) {
            float *u = data->dataWrite<float>(uattr, i);
            p.v(0) = u[0], p.v(1) = u[1];
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
