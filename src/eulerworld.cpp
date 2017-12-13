#include "eulerworld.hpp"
#include <cstdio>
#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>
#include <Partio.h>
#include <Eigen/Geometry>
#include "json/json.h"
#include "range.hpp"

using namespace Eigen;
using benlib::range;

#ifndef NDEBUG
std::ofstream debug;
#endif

//TODO: Add option for start and end point of grid
//TODO: Go through and consolidate pointers, look into 2*X dimentional vectors instead

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

    int frames = root.get("frames", 30).asInt();
    auto dtIn = root["dt"];
    if (dtIn.isNull()) {
        //if dt not there, check if steps is
        auto stepsIn = root["steps"];
        if(stepsIn.isNull()) {
            //nothing there, just set a default
            dt = 1.0/30.0;
            steps = 1;
        }
        else {
            steps = stepsIn.asInt();
            dt = (1.0/(double)frames)/steps;
        }
    } 
    else {
        dt = dtIn.asDouble();
        steps = (int)std::round((1.0/(double)frames)/dt);
    }
    
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
        auto originIn = gridIn["origin"];
        if (originIn.size() != 2) {
		  origin = Vector2d(0.0,0.0);
        } 
        else {
            origin = Vector2d(originIn[0].asDouble(),originIn[1].asDouble());
        }
        auto sizeIn = gridIn["size"];
        if (sizeIn.size() != 2) {
            std::cout<< "bad grid size, exiting" << std::endl;
			exit(-1);
        } 
        else {
            res[0] = sizeIn[0].asInt();
            res[1] = sizeIn[1].asInt();
        }
        h = gridIn.get("h", 1.0).asDouble();
		Vector2d end(0.0,0.0);
		auto endIn = gridIn["end"];
		if(endIn.size() == 2) {
			end = Vector2d(endIn[0].asDouble(), endIn[1].asDouble());
			Vector2d dif = end-origin;
			h = dif(0)/(res[0]-1);
		}
        center = ((h*Vector2d(res[0]-1, res[1]-1)) / 2.0) + origin;
    }
  
    auto gravityIn = root["gravity"];
    if (gravityIn.isNull() || gravityIn.size() != 2) {
        std::cout<< "no gravity" << std::endl;
        gravity(0) = 0.0;
        gravity(1) = 0.0;
        gravityEnabled = false;
    }
    else {
        gravity(0)= gravityIn[0].asDouble();
        gravity(1)= gravityIn[1].asDouble();
        gravityEnabled = true;
    }

    auto rotationIn = root["rotate"];
    if (rotationIn.isNull()) {
        std::cout << "no rotation" << std::endl;
        rotationEnabled = false;
        rotation = 0;
    } 
    else {
        rotationEnabled = true;
        rotation = rotationIn.asDouble();
    }
    

    for(size_t o = 0; o < objects.size(); o++) {
        Object& obj = objects[o];
        if(obj.type == "square") {
            /// Vector2d objCenter = obj.object + (Vector2d(obj.size[0],obj.size[1]) * 0.5);
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
            Vector2d objCenter = obj.object;
            //non-randomly make a circle
            //technically this is an ellipse because I'm using the same data as the
            //square and just using the size vector as radius of the semi-major and semi-minor axes
            //just make a square and reject those outside the ellipse.
            //~78.5% of the resx*resy particles are accepted - Pi/4 * (l*w) particles 
            double diffx = 2*obj.size[0] / (obj.ores[0]-1);
            double diffy = 2*obj.size[1] / (obj.ores[1]-1);
            for(int i = 0; i < obj.ores[0]; i++) {
                for(int j = 0; j < obj.ores[1]; j++) {
                    Vector2d pos = objCenter - Vector2d(obj.size[0], obj.size[1]) + Vector2d(diffx*i, diffy*j);
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
    }
  
    printf("Grid:\n");
    printf("Dimentions: %dx%d\n", res[0], res[1]);
    printf("Grid Spacing: %f\n", h);
    printf("X0: (%f, %f)\n", origin(0), origin(1));
    printf("X1: (%f, %f)\n", origin(0)+h*(res[0]-1), origin(1)+h*(res[1]-1));
    printf("Center: (%f, %f)\n", center(0), center(1));
    
    printf("\nConstants:\n");
    printf("Total Time: %f\n", totalTime);
    printf("dt: %f\n", dt);
    printf("Steps per Frame: %d\n", (int)std::round(1.0/(30*dt)));
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
    tau = new Vector2d[res[0]*res[1]];
    mat = new Vector2d[res[0]*res[1]];
    dXx = new Vector2d[(res[0]-1)*res[1]];
    dXy = new Vector2d[res[0]*(res[1]-1)];
    stress = new Matrix2d[(res[0]-1)*(res[1]-1)];
    F = new Matrix2d[(res[0]-1)*(res[1]-1)]; 
    valid.assign(res[0]*res[1], 0);
    cvalid.assign(res[0]*res[1], 0);
    phi = new double[res[0]*res[1]];
    tmpVel = new Vector2d[res[0]*res[1]];
    #ifndef NDEBUG
    dD = new Matrix2d[res[0]*res[1]];
    #endif
    
    {double *m = mass; for (int i=0; i<res[0]*res[1]; i++, m++) (*m) = 0.0;}
    {Vector2d *v = vel; for (int i=0; i<res[0]*res[1]; i++, v++) (*v) = Vector2d::Zero();}
    {Vector2d *f = tau; for (int i=0; i<res[0]*res[1]; i++, f++) (*f) = Vector2d::Zero();}
    {Vector2d *c = mat; for (int i=0; i<res[0]*res[1]; i++, c++) (*c) = Vector2d::Zero();}
    {Vector2d *dx = dXx; for (int i=0; i<(res[0]-1)*res[1]; i++, dx++) (*dx) = Vector2d::Zero();}
    {Vector2d *dy = dXy; for (int i=0; i<res[0]*(res[1]-1); i++, dy++) (*dy) = Vector2d::Zero();}
    {Matrix2d *s = stress; for (int i=0; i<(res[0]-1)*(res[1]-1); i++, s++) (*s) = Matrix2d::Zero();}
    {Matrix2d *dg = F; for (int i=0; i<(res[0]-1)*(res[1]-1); i++, dg++) (*dg) = Matrix2d::Zero();}
    
    inc = steps;
    #ifndef NDEBUG
    inc = 1;
    count = 0;
    sMax = 11;
    matTrans = new Vector2d*[sMax];
    for(int i = 0; i < sMax; i++) {
        matTrans[i] = new Vector2d[res[0]*res[1]];
    }
    vdiff = VectorXd(2*res[0]*res[1]);
    vdiff2 = VectorXd(2*res[0]*res[1]);
    forces = new Vector2d[res[0]*res[1]];
    advFrc = new Vector2d[res[0]*res[1]];
    elastFrc = new Vector2d[res[0]*res[1]];
    extFrc = new Vector2d[res[0]*res[1]];
    #endif
}

inline double interpolate(Vector2d point, double* field, Vector2d origin, int res[2], double h) {
    /// wx = x - x1;
    /// wy = y - y1;
    /// v = (1-wx)*(1-wy)*v[i1][j1] + (wx)*(1-wy)*v[i1+1][j1] + (1-wx)*(wy)*v[i1][j1+1] + (wx)*(wy)*v[i1+1][j1+1];
    Vector2d x = (point - origin);
    Vector2d ij = x / h;    
    int i1 = (int)ij(0);
    int j1 = (int)ij(1);
    int i2 = i1+1;
    int j2 = j1+1;
    double x1 = h*i1;
    double y1 = h*j1;
    double wx = (x(0) - x1) / h;
    double wy = (x(1) - y1) / h;
    double v = 0;
    if(i1*res[1]+j1 >= 0 && i1*res[1]+j1 < res[0]*res[1]) {
        v += (1-wx)*(1-wy)*field[i1*res[1]+j1];
    }
    if(i2*res[1]+j1 >= 0 && i2*res[1]+j1 < res[0]*res[1]) {
        v += (wx)*(1-wy)*field[i2*res[1]+j1];
    }
    if(i1*res[1]+j2 >= 0 && i1*res[1]+j2 < res[0]*res[1]) {
        v += (1-wx)*(wy)*field[i1*res[1]+j2];
    }
    if(i2*res[1]+j2 >= 0 && i2*res[1]+j2 < res[0]*res[1]) {
        v += (wx)*(wy)*field[i2*res[1]+j2];
    }
    /// double v = (1-wx)*(1-wy)*field[i1*res[1]+j1] + 
               /// (wx)  *(1-wy)*field[i2*res[1]+j1] + 
               /// (1-wx)*(wy)  *field[i1*res[1]+j2] + 
               /// (wx)  *(wy)  *field[i2*res[1]+j2];
    return v;
}

inline Vector2d interpolate(Vector2d point, Vector2d* field, Vector2d origin, int res[2], double h) {
    /// wx = x - x1;
    /// wy = y - y1;
    /// v = (1-wx)*(1-wy)*v[i1][j1] + (wx)*(1-wy)*v[i1+1][j1] + (1-wx)*(wy)*v[i1][j1+1] + (wx)*(wy)*v[i1+1][j1+1];
    Vector2d x = (point - origin);
    Vector2d ij = x / h;    
    int i1 = (int)ij(0);
    int j1 = (int)ij(1);
    int i2 = i1+1;
    int j2 = j1+1;
    double x1 = h*i1;
    double y1 = h*j1;
    double wx = (x(0) - x1) / h;
    double wy = (x(1) - y1) / h;
    Vector2d v = Vector2d::Zero();
    if(i1*res[1]+j1 >= 0 && i1*res[1]+j1 < res[0]*res[1]) {
        v += (1-wx)*(1-wy)*field[i1*res[1]+j1];
    }
    if(i2*res[1]+j1 >= 0 && i2*res[1]+j1 < res[0]*res[1]) {
        v += (wx)*(1-wy)*field[i2*res[1]+j1];
    }
    if(i1*res[1]+j2 >= 0 && i1*res[1]+j2 < res[0]*res[1]) {
        v += (1-wx)*(wy)*field[i1*res[1]+j2];
    }
    if(i2*res[1]+j2 >= 0 && i2*res[1]+j2 < res[0]*res[1]) {
        v += (wx)*(wy)*field[i2*res[1]+j2];
    }
    /// Vector2d v = (1-wx)*(1-wy)*field[i1*res[1]+j1] + 
                 /// (wx)  *(1-wy)*field[i2*res[1]+j1] + 
                 /// (1-wx)*(wy)  *field[i1*res[1]+j2] + 
                 /// (wx)  *(wy)  *field[i2*res[1]+j2];
    
    return v;
}

inline Vector3d interpolate(Vector2d point, Vector3d* field, Vector3d bg, Vector2d origin, int res[2], double h) {
    /// wx = x - x1;
    /// wy = y - y1;
    /// v = (1-wx)*(1-wy)*v[i1][j1] + (wx)*(1-wy)*v[i1+1][j1] + (1-wx)*(wy)*v[i1][j1+1] + (wx)*(wy)*v[i1+1][j1+1];
    Vector2d x = (point - origin);
    Vector2d ij = x / h;    
    int i1 = (int)ij(0);
    int j1 = (int)ij(1);
    int i2 = i1+1;
    int j2 = j1+1;
    double x1 = h*i1;
    double y1 = h*j1;
    double wx = (x(0) - x1) / h;
    double wy = (x(1) - y1) / h;
    Vector3d v = Vector3d::Zero();
    if(i1*res[1]+j1 >= 0 && i1*res[1]+j1 < res[0]*res[1]) {
        v += (1-wx)*(1-wy)*field[i1*res[1]+j1];
    }
    else {
		v += (1-wx)*(1-wy)*bg;
	}
    if(i2*res[1]+j1 >= 0 && i2*res[1]+j1 < res[0]*res[1]) {
        v += (wx)*(1-wy)*field[i2*res[1]+j1];
    }
    else {
		v += (wx)*(1-wy)*bg;
	}
    if(i1*res[1]+j2 >= 0 && i1*res[1]+j2 < res[0]*res[1]) {
        v += (1-wx)*(wy)*field[i1*res[1]+j2];
    }
    else {
		v += (1-wx)*(wy)*bg;
	}
    if(i2*res[1]+j2 >= 0 && i2*res[1]+j2 < res[0]*res[1]) {
        v += (wx)*(wy)*field[i2*res[1]+j2];
    }
    else {
		v += (wx)*(wy)*bg;
	}
    /// Vector2d v = (1-wx)*(1-wy)*field[i1*res[1]+j1] + 
                 /// (wx)  *(1-wy)*field[i2*res[1]+j1] + 
                 /// (1-wx)*(wy)  *field[i1*res[1]+j2] + 
                 /// (wx)  *(wy)  *field[i2*res[1]+j2];
    
    return v;
}

void World::init() {
    color = new Vector3d[res[0]*res[1]];
    interpColor.resize(3*res[0]*res[1]);
    //Init material coordinates
    //(Temp for now) set up mass field
    /// Vector2d middle = origin + h*Vector2d((res[0]-1)/2, (res[1]-1)/2);
    double size[2] = {0.15, 0.15};	// ~1ft wide
    int count = 0;
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            Vector2d pos = origin + h*Vector2d(i, j);
            Vector2d ph = pos - center;
            // For Testing:
            // rotate material field and test D(X)
            Rotation2D<double> rot(-0.785398);
            Vector2d rpos = rot * ph;
            mat[index] = rpos + center;
            color[index] = Vector3d(135,206,255);
            vel[index] = 0.5*Vector2d(-ph(1), ph(0));
            //Create a square in the center
            /// if(i >= (res[0]/2)-4 && i <= (res[0]/2)+4 && j >= (res[1]/2)-4 && j <= (res[1]/2)+4) {
            if( ((ph(0)*ph(0))/(size[0]*size[0])) + ((ph(1)*ph(1))/(size[1]*size[1])) < 1+EPS) {	//circle
                /// mass[index] = 0.067;  //.067kg * 9x9 = 5.44kg ~ 12lbs
                count++;
                //Start it out rotating
                vel[index] = 0.5*Vector2d(-ph(1), ph(0));
                /// vel[index] = 0.5*Vector2d(1,0);
                Vector3d col = Vector3d(255, 0, 0);
                if(ph(0) < 0 && ph(1) < 0) {
                    col = Vector3d(0, 255, 0);
                }
                if(ph(0) >= 0 && ph(1) < 0) {
                    col = Vector3d(0, 0, 255);
                }
                if(ph(0) < 0 && ph(1) >= 0) {
                    col = Vector3d(255, 255, 0);
                }
                color[index] = col;
                valid[index] = 1;
            }
        }
    }
    double pmass = 9.07/count; // ~20lb / count
    for(int i = 0; i < res[0]; i++) {
		for(int j = 0; j < res[1]; j++) {
			int ind = i*res[1]+j;
			if(valid[ind]) {
				mass[ind] = pmass;
				//Add a passive particle at the cell
                #ifndef NDEBUG
                Particle p(mat[ind], Vector2d::Zero(), color[ind], pmass);
                particles.push_back(p);
                #endif
			}
		}
	}
    interpMass.resize(res[0]*res[1]);
    Fstar.resize(2*res[0]*res[1]);
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
    if(stepNum == 0) {
        #ifndef NDEBUG
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                matTrans[count][i*res[1]+j] = mat[i*res[1]+j];
            }
        }
        count++;
        #endif
        //write initial mass
        std::ofstream mout("./euler/mass0");
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                mout << mass[i*res[1]+j] << " ";
            }
            mout << "\n";
        }
        mout.close();
        //output vel for testing
        std::ofstream vout("./euler/vel0");
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                Vector2d p = origin+h*Vector2d(i,j);
                int ind = i*res[1]+j;
                vout << p(0) << " " << p(1) << " " << vel[ind](0) << " " << vel[ind](1) << " 0 0 255" << "\n";
            }
        }
        vout.close();
        
        std::ofstream colout("./euler/col0");
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                /// Vector3d c = color[i*res[1]+j];
                Vector2d coord = mat[i*res[1]+j];
                Vector3d c = interpolate(coord, color, Vector3d(135, 206, 255), origin, res, h);
                Vector2d pos = origin+h*Vector2d(i,j);
                colout << pos(0) << " " << pos(1) << " " << (int)std::round(c(0)) << " " << (int)std::round(c(1)) << " " << (int)std::round(c(2)) << "\n";
            }
        }
        colout.close();
    }
    getMass();
    cvalid = valid;
    #ifdef MAT_EXTRAP
    extrapolate_mat();
    valid = cvalid;
    #endif
    //Calculates deformation gradient and stress
    getDeformationGradient();
    computeStress();
    //Calculates the forces vector
    computeForce();
    //Solve M*v_t+1 = f*
    solveVelocity();
    velExtrapolate();
    //Advect material coordinates
    advect();
    
    stepNum++;
	elapsedTime += dt;
    
    if(stepNum % inc == 0) {
        //Write current mass
        std::ofstream mout("./euler/mass"+std::to_string(stepNum/inc));
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                mout << interpMass(i*res[1]+j) << " ";
            }
            mout << "\n";
        }
        mout.close();
        
        std::ofstream colout("./euler/col"+std::to_string(stepNum/inc));
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                Vector3i c = interpColor.segment(3*(i*res[1]+j), 3);
                Vector2d pos = origin+h*Vector2d(i,j);
                colout << pos(0) << " " << pos(1) << " " << c(0) << " " << c(1) << " " << c(2) << "\n";
            }
        }
        colout.close();
        
        std::ofstream vout("./euler/vel"+std::to_string(stepNum/inc));
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                Vector2d p = origin+h*Vector2d(i,j);
                int ind = i*res[1]+j;
                Vector3i col(255, 0, 0);
                if(cvalid[i*res[1]+j]) {
                    col = Vector3i(0, 0, 255);
                }
                vout << p(0) << " " << p(1) << " " << vel[ind](0) << " " << vel[ind](1) << " " << col(0) << " " << col(1) << " " << col(2) << "\n";
            }
        }
        vout.close();
        
        std::ofstream fout("./euler/frc"+std::to_string(stepNum/inc));
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                Vector2d p = origin+h*Vector2d(i,j);
                int ind = i*res[1]+j;
                Vector3i col(255, 0, 0);
                if(cvalid[i*res[1]+j]) {
                    col = Vector3i(0, 0, 255);
                }
                fout << p(0) << " " << p(1) << " " << forces[ind](0) << " " << forces[ind](1) << " " << col(0) << " " << col(1) << " " << col(2) << "\n";
            }
        }
        fout.close();
        
        std::ofstream afout("./euler/advfrc"+std::to_string(stepNum/inc));
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                Vector2d p = origin+h*Vector2d(i,j);
                int ind = i*res[1]+j;
                Vector3i col(255, 0, 0);
                if(cvalid[i*res[1]+j]) {
                    col = Vector3i(0, 0, 255);
                }
                afout << p(0) << " " << p(1) << " " << advFrc[ind](0) << " " << advFrc[ind](1) << " " << col(0) << " " << col(1) << " " << col(2) << "\n";
            }
        }
        afout.close();
        
        std::ofstream efout("./euler/elastfrc"+std::to_string(stepNum/inc));
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                Vector2d p = origin+h*Vector2d(i,j);
                int ind = i*res[1]+j;
                Vector3i col(255, 0, 0);
                if(cvalid[i*res[1]+j]) {
                    col = Vector3i(0, 0, 255);
                }
                efout << p(0) << " " << p(1) << " " << elastFrc[ind](0) << " " << elastFrc[ind](1) << " " << col(0) << " " << col(1) << " " << col(2) << "\n";
            }
        }
        efout.close();
        
        #ifndef NDEBUG
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                matTrans[count][i*res[1]+j] = mat[i*res[1]+j];
            }
        }
        count++;
        if(count == sMax) {
            std::ofstream mats("./euler/mats.txt");
            for(int i = 0; i < res[0]; i++) {
                for(int j = 0; j < res[1]; j++) {
                    if(valid[i*res[1]+j]) {
                        for(int k = 0; k < sMax; k++) {
                            mats << matTrans[k][i*res[1]+j](0) << " " << matTrans[k][i*res[1]+j](1) << " " << k << "\n";
                        }
                        mats << "\n\n";
                    }
                }
            }
            mats.close();
            
            std::ofstream par("./euler/particles");
            for(int i = 0; i < (int)particles.size(); i++) {
				Particle& p = particles[i];
				int j;
				for(j = 0; j < (int)p.hist.size(); j++) {
					par << p.hist[j](0) << " " << p.hist[j](1) << " " << j << "\n";
				}
				par << p.x(0) << " " << p.x(1) << " " << j << "\n\n\n";
			}
            par.close();
            std::exit(0);
        }
        #endif
    }
}

//Loop through grid nodes
//  Lookup material coordinate at that point
//  Interpolate initial mass field around that point. Form vector, then matrix
void World::getMass() {
    numNodesMass = 0;
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            //Get mat coords
            Vector2d coord = mat[index];
            //interpolate mass at that point
            double m = interpolate(coord, mass, origin, res, h);
            interpMass(index) = m;
            //Interpolate color for rendering
            Vector3d cd = interpolate(coord, color, Vector3d(135, 206, 255), origin, res, h);
            Vector3i c((int)std::round(cd(0)), (int)std::round(cd(1)), (int)std::round(cd(2)));
            interpColor.segment(3*index, 3) = c;
            if(m > EPS) {
                numNodesMass++;
                valid[index] = 1;
                /// if(stepNum == 0) {
					/// Vector2d pos = origin + h*Vector2d(i, j);
					/// Vector2d ph = pos - center;
					/// vel[index] = 0.5*Vector2d(-ph(1), ph(0));
				/// }
            }
            else {
                interpMass(index) = 0;
                valid[index] = 0;
				/// if(stepNum == 0) {
					/// vel[index] = Vector2d::Zero();
				/// }
            }
        }
    }
    
    #ifndef NDEBUG
    if(stepNum % inc == 0) {
        debug << "\n\nStep " << stepNum << "\n";
        debug << "Non-0 Mass: " << numNodesMass << "\n";
        debug << "\nX:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                debug << "(" << mat[i*res[1]+j](0) << "," << mat[i*res[1]+j](1) << ") ";
            }
            debug << "\n";
        }
        debug << "Mass:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                debug << mass[i*res[1]+j] << " ";
            }
            debug << "\n";
        }
        debug << "InterpMass:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                debug << interpMass(i*res[1]+j) << " ";
            }
            debug << "\n";
        }
        debug << "Valid:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                debug << (int)valid[i*res[1]+j] << " ";
            }
            debug << "\n";
        }
        debug << "Vel:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                debug << "(" << vel[i*res[1]+j](0) << "," << vel[i*res[1]+j](1) << ") ";
            }
            debug << "\n";
        }
    }
    /// debug.close();
    /// std::exit(0);
    #endif
}

#ifdef MAT_EXTRAP
//Extrapolate material coordinates from the center object
void World::extrapolate_mat() {
	//First create phi
    //As a hack for now, just copy the values so that the invalid ones are a very large value
    for(int i = 0; i < res[0]*res[1]; i++) {
        phi[i] = std::numeric_limits<double>::infinity();
        if(valid[i]) {
            phi[i] = -h/2.0;
        }
    }
    std::vector<char> tmpvalid = valid;
    //First sweep to extend distance field
    distSweep(phi, tmpvalid, 20, 1e-12);
    
	//Extrapolate out in a fast marching approach
	std::vector<std::tuple<double,int,int>> sortedphi;
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            if(!valid[i*res[1]+j]) {
                //Could switch but then sort requires custom function
                sortedphi.push_back(std::make_tuple(phi[i*res[1]+j], i, j));
            }
        }
    }
    std::sort(sortedphi.begin(), sortedphi.end());
    
    //Go through array, calculating new material for each element
    std::vector<std::tuple<double,int,int>>::iterator it;
	double xh = 0, yh = 0;
	Vector2d xm, ym;
    for(it = sortedphi.begin(); it != sortedphi.end(); it++) {
        int x = std::get<1>(*it);
        int y = std::get<2>(*it);
        int index = x*res[1]+y;
        int xp1 = (x+1)*res[1]+y;
        int xm1 = (x-1)*res[1]+y;
        int yp1 = x*res[1]+(y+1);
        int ym1 = x*res[1]+(y-1);
        
        //Look in x direction
        if(x < res[0]-1 && valid[xp1]) {	//right valid, -h
			xh = -h;
			xm = mat[xp1];
		}
		else if(x > 0 && valid[xm1]) {		//left valid, +h
			xh = h;
			xm = mat[xm1];
		}
		else {								//0
			xh = 0;
			xm = Vector2d::Zero();
		}
        
        //look in y direction
        if(y < res[1]-1 && valid[yp1]) {	//right valid, -h
			yh = -h;
			ym = mat[yp1];
		}
		else if(y > 0 && valid[ym1]) {		//left valid, +h
			yh = h;
			ym = mat[ym1];
		}
		else {								//0
			yh = 0;
			ym = Vector2d::Zero();
		}
		
		//If one direction is empty, we want to keep the value from the other direction
		if(xh == 0 && yh != 0) {
			xm = ym;
		}
		else if(yh == 0 && xh != 0) {
			ym = xm;
		}
		else if(xh == 0 && yh == 0){
			printf("Really shouldn't happen");
		}
		mat[index] = Vector2d(xm(0)+xh, ym(1)+yh);
		valid[index] = 1;
    }
}
#endif

//Loop through grid nodes from i,j to m-1,n-1
//  Do partial derivatives along i,j and i+1,j+1 boundaries. Store partial vector results on staggered grid
//Loop through cell centers
//  Interpolate partial results in x and y direction
//  Form deformation gradient F
void World::getDeformationGradient() {
    //Partials along X  (i+1,j)-(i,j)
    for(int i = 0; i < res[0]-1; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            int ind2 = (i+1)*res[1]+j;
            Vector2d dx = (mat[ind2] - mat[index]) / h;
            dXx[index] = dx;
        }
    }
    //Partials along Y  (i,j+1)-(i,j)
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]-1; j++) {
            int index = i*(res[1]-1)+j;
            int ind2 = i*(res[1]-1)+(j+1);
            Vector2d dy = (mat[ind2] - mat[index]) / h;
            dXy[index] = dy;
        }
    }
    //Cell Centers
    for(int i = 0; i < res[0]-1; i++) {
        for(int j = 0; j < res[1]-1; j++) {
            int index1 = i*res[1]+j;
            int index2 = i*(res[1]-1)+j;
            //Average dx values along (i,j) and (i,j+1)
            Vector2d dx = (dXx[index1] + dXx[i*res[1]+(j+1)]) / 2.0;
            //Average dy values along (i,j) and (i+1,j)
            Vector2d dy = (dXy[index2] + dXy[(i+1)*(res[1]-1)+j]) / 2.0;
            Matrix2d grad;
            grad.col(0) = dx;
            grad.col(1) = dy;
            if(std::abs(grad.determinant()) > 1e-4) {
                F[index2] = grad.inverse();
            }
            else {
                F[index2] = Matrix2d::Zero();                
            }
        }
    }
    
    #ifndef NDEBUG
    if(stepNum % inc == 0) {
        //Print gradients
        debug << "F: \n";
        for(int j = res[1]-2; j >= 0; j--) {
            for(int k = 0; k < 2; k++) {
                for(int i = 0; i < res[0]-1; i++) {
                    int ind = i*res[0]+j;
                    debug << "[";
                    if(F[ind](k,0) < 0) {
                        debug << F[ind](k,0);
                    }
                    else {
                        debug << F[ind](k,0);
                    }
                    debug << ",";
                    if(F[ind](k,1) < 0) {
                        debug << F[ind](k,1);
                    }
                    else {
                        debug << F[ind](k,1);
                    }
                    debug << "] ";
                }
                debug << "\n";
            }
            debug << "\n";
        }
        debug << "\n";
    }
    /// debug.close();
    /// std::exit(0);
    #endif
}

//Loop through cell centers
//  Calculate stress from F using stress-strain model
void World::computeStress() {
    //Only dealing with 1 object right now
    MaterialProps &mp = objects[0].mp;  
    for(int i = 0; i < res[0]-1; i++) {
        for(int j = 0; j < res[1]-1; j++) {
            int index = i*res[1]+j;
            Matrix2d Ft = F[index].transpose();
            Matrix2d strain = 0.5 * (Ft*F[index] - Matrix2d::Identity());
            Matrix2d linstress = mp.lambda*strain.trace()*Matrix2d::Identity() + 2.0*mp.mu*strain;
            stress[index] = linstress;
        }
    }
    
    #ifndef NDEBUG
    if(stepNum % inc == 0) {
        //Print gradients
        debug << "Stress: \n";
        for(int j = res[1]-2; j >= 0; j--) {
            for(int k = 0; k < 2; k++) {
                for(int i = 0; i < res[0]-1; i++) {
                    int ind = i*res[0]+j;
                    debug << "[";
                    if(stress[ind](k,0) < 0) {
                        debug << stress[ind](k,0);
                    }
                    else {
                        debug << stress[ind](k,0);
                    }
                    debug << ",";
                    if(stress[ind](k,1) < 0) {
                        debug << stress[ind](k,1);
                    }
                    else {
                        debug << stress[ind](k,1);
                    }
                    debug << "] ";
                }
                debug << "\n";
            }
            debug << "\n";
        }
        debug << "\n";
    }
    /// debug.close();
    /// std::exit(0);
    #endif
}

//Loop through grid nodes
//  Compute tau_c for each node direction
//  Compute tau from summing them
//  Compute D(v) using upwinding scheme
//  f* = Mv-dt*(M*D(v)*v - tau - f_ext)
//  add f* to F*. Each entry in F* corresponds to the same node in M
void World::computeForce() {
	#ifndef NDEBUG
	double maxErr = 0;
	#endif
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            /// if(!valid[index]) {
                /// tau[index] = Vector2d::Zero();
                /// Fstar.segment(2*index,2) = Vector2d::Zero();
                /// #ifndef NDEBUG
                /// dD[index] = Matrix2d::Zero();
                /// #endif
                /// continue;
            /// }
            //Compute tau by summing up all the tau_c for the 4 neighbors
            //-.5L(sum of n * stress)
            Vector2d ftau(0,0);
            Vector2d n1, n2;
            if(i != 0) { //left side
                if(j != 0) {    //bottom left, n=(-1,-1)
                    n1 = Vector2d(1, 0);
                    n2 = Vector2d(0, 1);
                    Matrix2d s = stress[(i-1)*res[1]+(j-1)];
                    ftau += (s*n1 + s*n2);
                }
                if(j != res[1]-1) { //upper left, n=(-1,1)
                    n1 = Vector2d(1, 0);
                    n2 = Vector2d(0, -1);
                    Matrix2d s = stress[(i-1)*res[1]+j];
                    ftau += (s*n1 + s*n2);
                }
            }
            if(i != res[0]-1) {    //right side
                if(j != 0) {    //bottom right, n=(1,-1)
                    n1 = Vector2d(-1, 0);
                    n2 = Vector2d(0, 1);
                    Matrix2d s = stress[i*res[1]+(j-1)];
                    ftau += (s*n1 + s*n2);
                }
                if(j != res[1]-1) { //upper right, n=(1,1)
                    n1 = Vector2d(-1, 0);
                    n2 = Vector2d(0, -1);
                    Matrix2d s = stress[i*res[1]+j];
                    ftau += (s*n1 + s*n2);
                }
            }
            ftau *= -0.5*h; 
            tau[index] = ftau;
            
            //Compute discrete jacobian
            Matrix2d D;
            if(!valid[index]) {
				D = Matrix2d::Zero();
			}
			else {
				D = upwindJac(vel, i, j);
			}
            #ifndef NDEBUG
            dD[index] = D;
            Matrix2d exactD;
            if(cvalid[index]) {
				exactD << 0.0, -0.5, 0.5, 0.0;
			}
			else {
				exactD << 0.0, 0.0, 0.0, 0.0;
			}
            double err = (D-exactD).norm();
            if(err > maxErr) {
				maxErr = err;
			}
            #endif
            //Compute force
            Vector2d fext(0.0,0.0);
            if(gravityEnabled) {
                fext += interpMass(index)*gravity;
            }
            if(rotationEnabled) {
                Vector2d d = origin+Vector2d(h*i,h*j)-center;
                fext += interpMass(index)*rotation*Vector2d(-d(1), d(0));
            }
            Matrix2d Mdiag = interpMass(index)*Matrix2d::Identity();
            /// D << 0, -0.5, 0.5, 0;
            Vector2d frc = Mdiag*D*vel[index] - tau[index] - fext;
            #ifndef NDEBUG
            forces[index] = frc;
            advFrc[index] = Mdiag*D*vel[index];
            elastFrc[index] = -tau[index];
            extFrc[index] = -fext;
            #endif
            //Fstar is actually momentum
            Vector2d fstar = Mdiag*vel[index] - dt*frc;
            Fstar.segment(2*index,2) = fstar;
        }
    }
    
    #ifndef NDEBUG
    if(stepNum % inc == 0) {
		std::ofstream eout("./euler/dv", std::ios_base::app | std::ios_base::out);
		eout << maxErr << std::endl;
		eout.close();
	}
    if(stepNum % inc == 0) {
        //Print tau, D, and forces
        debug << "Tau:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                debug << "(" << tau[i*res[1]+j](0) << "," << tau[i*res[1]+j](1) << ") ";
            }
            debug << "\n";
        }
        debug << "\n";
        debug << "Dv: \n";
        for(int j = res[1]-2; j >= 0; j--) {
            for(int k = 0; k < 2; k++) {
                for(int i = 0; i < res[0]-1; i++) {
                    int ind = i*res[0]+j;
                    if(i == 31 && j == 48) {
						printf("%f %f\n", dD[ind](k, 0), dD[ind](k, 1));
					}
                    debug << "[";
                    if(dD[ind](k,0) < 0) {
                        debug << dD[ind](k,0);
                    }
                    else {
                        debug << dD[ind](k,0);
                    }
                    debug << ",";
                    if(dD[ind](k,1) < 0) {
                        debug << dD[ind](k,1);
                    }
                    else {
                        debug << dD[ind](k,1);
                    }
                    debug << "] ";
                }
                debug << "\n";
            }
            debug << "\n";
        }
        debug << "\n";
        debug << "Advection Forces:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                int index = i*res[1]+j;
                debug << "(" << advFrc[index](0) << "," << advFrc[index](1) << ") ";
            }
            debug << "\n";
        }
        debug << "\n";
        debug << "Elastic Forces:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                int index = i*res[1]+j;
                debug << "(" << elastFrc[index](0) << "," << elastFrc[index](1) << ") ";
            }
            debug << "\n";
        }
        debug << "\n";
        debug << "External Forces:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                int index = i*res[1]+j;
                debug << "(" << extFrc[index](0) << "," << extFrc[index](1) << ") ";
            }
            debug << "\n";
        }
        debug << "\n";
        debug << "Forces:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                int index = i*res[1]+j;
                debug << "(" << forces[index](0) << "," << forces[index](1) << ") ";
            }
            debug << "\n";
        }
        debug << "\n";
        debug << "F*:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                int index = i*res[1]+j;
                debug << "(" << Fstar.segment(2*index,2)(0) << "," << Fstar.segment(2*index,2)(1) << ") ";
            }
            debug << "\n";
        }
        debug << "\n";
    }
    /// debug.close();
    /// std::exit(0);
    #endif
}

//Remove inactive nodes (M(i,j) = 0)
//Solve Mv=F*
void World::solveVelocity() {
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int ind = i*res[1]+j;
            if(interpMass(ind) > EPS) {
                Vector2d m = interpMass(ind)*Vector2d::Ones();
                Vector2d mi = m.cwiseInverse();
                Matrix2d Mi = mi.asDiagonal();
                Vector2d newVel = Mi*Fstar.segment(2*ind,2);
                #ifndef NDEBUG
                if(stepNum/steps > 9 && stepNum/steps < 16 && stepNum % inc == 0) {
					vdiff.segment(2*ind, 2) = newVel;
					if(stepNum/steps == 10) {
						vdiff2.segment(2*ind, 2) = newVel;
					}
				}
                #endif
                vel[ind] = newVel;
            }
        }
    }
    
    #ifndef NDEBUG
    if(stepNum/steps > 10 && stepNum/steps < 16 && stepNum % inc == 0) {
		std::ofstream vdiffout("./euler/vdiff"+std::to_string(stepNum/steps));
		/// vdiffout << "Frame " << stepNum/steps << ", Step " << stepNum%steps << "\n";
		for(int i = 0; i < res[0]; i++) {
			for(int j = 0; j < res[1]; j++) {
				int ind = i*res[1]+j;
				Vector2d v = vdiff.segment(2*ind, 2)-vdiff2.segment(2*ind, 2);
				/// vdiffout << "(" << v(0) << "," << v(1) << ") ";
				Vector2d p = origin+h*Vector2d(i,j);
                Vector3i col(255, 0, 0);
                if(cvalid[i*res[1]+j]) {
                    col = Vector3i(0, 0, 255);
                }
                vdiffout << std::setprecision(8) << std::scientific << p(0) << " " << p(1) << " " << v(0) << " " << v(1) << " " << col(0) << " " << col(1) << " " << col(2) << "\n";
				vdiff2 = vdiff;
			}
			/// vdiffout << "\n";
		}
		vdiffout.close();
	}
    if(stepNum % inc == 0) {
        //print new velocity
        debug << "New Vel:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                int index = i*res[1]+j;
                debug << "(" << vel[index](0) << "," << vel[index](1) << ") ";
            }
            debug << "\n";
        }
        debug << "\n";
    }
    /// debug.close();
    /// std::exit(0);
    #endif
}

//Extrapolate new v using fast sweeping method
//First step is to create the distance field. Starting with the inside distances
//defined sweep them out to get outside distances to the surface. Then sweep the 
//velocity out using the distance field in order to solve the PDE.
void World::velExtrapolate() {
	//Initialize level set with merging spheres
	//Reinitialization, solving the eikenol equation - fast first order
	//Set boundary of +/- h/2
	//sqrt(3)/2 *h for spheres
	
    //As a hack for now, just copy the values so that the invalid ones are a very large value
    for(int i = 0; i < res[0]*res[1]; i++) {
        phi[i] = std::numeric_limits<double>::infinity();
        if(valid[i]) {
            phi[i] = -h/2.0;
        }
    }
    std::vector<char> tmpvalid = valid;
    
    //First sweep to extend distance field
    distSweep(phi, tmpvalid, 20, 1e-12);
    //Then sweep the velocities over that field solving delu dot delphi = 0
    velSweep(vel, phi, 20, 1e-12);
    
    //TODO: Test only advecting where the object is
    /// valid = tmpvalid;
    //Test regular averaging
    /// sweepAve(vel, 5);
    
    //Try Exact field
    /// Vector2d middle = origin + h*Vector2d((res[0]-1)/2, (res[1]-1)/2);
    /// for(int i = 0; i < res[0]; i++) {
		/// for(int j = 0; j < res[1]; j++) {
			/// /// if(valid[i*res[1]+j]) {
				/// /// continue;
			/// /// }
			/// Vector2d pos = origin+h*Vector2d(i,j);
			/// Vector2d ph = pos - middle;
			/// Vector2d v(-ph(1), ph(0));
			/// vel[i*res[1]+j] = 0.5*v;
		/// }
	/// }
    
    #ifndef NDEBUG
    if(stepNum % inc == 0) {
        //print phi
        debug << "Phi field:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                int ind = i*res[1]+j;
                debug << phi[ind] << " ";
            }
            debug << "\n";
        }
        //print vel
        debug << "\nExtrapolated Vel:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                int index = i*res[1]+j;
                debug << "(" << vel[index](0) << "," << vel[index](1) << ") ";
            }
            debug << "\n";
        }
    }
    /// debug.close();
    /// std::exit(0);
    #endif
}

//Advect material coordinates with new velocity
//X = X - dt*D(v)*X
///They might have switched v and X at the end so it would be
///X = X - dt*D(X)*v
void World::advect() {
	//Passively move particles
	#ifndef NDEBUG
	double maxErr = 0;
	for(int i = 0; i < (int)particles.size(); i++) {
		Particle& p = particles[i];
		//interpolate velocity field at position
		Vector2d v = interpolate(p.x, vel, origin, res, h);
		//move (Forward Euler)
		if(stepNum % inc == 0) {
			p.hist.push_back(p.x);
		}
		p.x += dt*v;
	}
	#endif
    double alpha = 0;
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int ind = i*res[1]+j;
            if(!valid[ind]) {
                #ifndef NDEBUG
                dD[ind] = Matrix2d::Zero();
                #endif
                continue;
            }
            Matrix2d D;
            D = upwindJac(mat, i, j);
            #ifndef NDEBUG
            dD[ind] = D;
            Matrix2d exactD;
            if(cvalid[ind]) {
				Rotation2D<double> rot(-0.785398-dt*0.5*elapsedTime);
				exactD = rot.toRotationMatrix();
			}
			else {
				exactD = Matrix2d::Identity();
			}
            double err = (D-exactD).norm();
            if(err > maxErr) {
				maxErr = err;
			}
            #endif
            double a = dt * std::max(std::abs(vel[ind](0))/h, std::abs(vel[ind](1))/h);
            if(a > alpha) {
                alpha = a;
            }
            //For testing override D to be exact.
            //This should be equal to the deformation gradient
            /// D = Matrix2d::Identity();
            mat[ind] = mat[ind] - dt*D*vel[ind];
        }
    }
    
    if(alpha >= 1.0) {
        printf("\nAlpha: %f\n", alpha);
        printf("Frame %d, Step %d\n", stepNum/steps, stepNum%steps);
        #ifndef NDEBUG
        debug.close();
        #endif
        std::exit(0);
    }
    #ifndef NDEBUG
    if(stepNum % inc == 0) {
		std::ofstream eout("./euler/dx", std::ios_base::app | std::ios_base::out);
		eout << maxErr << std::endl;
		eout.close();
	}
    if(stepNum % inc == 0) {
        debug << "\nDX: \n";
        for(int j = res[1]-2; j >= 0; j--) {
            for(int k = 0; k < 2; k++) {
                for(int i = 0; i < res[0]-1; i++) {
                    int ind = i*res[0]+j;
                    if(i == 2 && j == 16) {
						printf("%f %f\n", dD[ind](k, 0), dD[ind](k, 1));
					}
                    debug << "[";
                    if(dD[ind](k,0) < 0) {
                        debug << std::fixed << std::setprecision(6) << dD[ind](k,0);
                    }
                    else {
                        debug << std::fixed << std::setprecision(6) << dD[ind](k,0);
                    }
                    debug << ",";
                    if(dD[ind](k,1) < 0) {
                        debug << std::fixed << std::setprecision(6) << dD[ind](k,1);
                    }
                    else {
                        debug << std::fixed << std::setprecision(6) << dD[ind](k,1);
                    }
                    debug << "] ";
                }
                debug << "\n";
            }
            debug << "\n";
        }
        debug << "\n";
        debug << "Material Advection:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                int index = i*res[1]+j;
                debug << "(" << mat[index](0) << "," << mat[index](1) << ") ";
            }
            debug << "\n";
        }
    }
    debug.close();
    std::exit(0);
    #endif
}

/************************
 * Helpers
 ************************/
Matrix2d World::upwindJac(Vector2d* field, int i, int j, bool boundary) {
    //Compute D(field) with upwinding scheme
    int index = i*res[1]+j;
    int xp1 = (i+1)*res[1]+j;
    int xm1 = (i-1)*res[1]+j;
    int yp1 = i*res[1]+(j+1);
    int ym1 = i*res[1]+(j-1);
    //DX
    int di = 2, dj = 16;
    //Dv
    /// int di = 31, dj = 48;
    if(i == di && j == dj) {
		printf("(i, j): (%f, %f), index: %d\n", field[index](0), field[index](1), (int)valid[index]);
		printf("(i+1, j): (%f, %f), xp1: %d\n", field[xp1](0), field[xp1](1), (int)valid[xp1]);
		printf("(i-1, j): (%f, %f), xm1: %d\n", field[xm1](0), field[xm1](1), (int)valid[xm1]);
		printf("(i, j+1): (%f, %f), yp1: %d\n", field[yp1](0), field[yp1](1), (int)valid[yp1]);
		printf("(i, j-1): (%f, %f), ym1: %d\n", field[ym1](0), field[ym1](1), (int)valid[ym1]);
	}
    double xx = 0, xy = 0, yx = 0, yy = 0;      //Upwind scheme so look at each value individually
    if(vel[index](0) >= 0) {    //X Velocity is positive, so moving from left to right, look left
		if(i == di && j == dj) {
			printf("Look X -1\n");
		}
        if((i != 0) && ((boundary && valid[xm1]) || !boundary)) {  //Looking backwards in x and valid value
			xx = (field[index](0) - field[xm1](0)) / h;
			if(i == di && j == dj) {
				printf("xx1: (%f - %f) / %f = %f\n", field[index](0), field[xm1](0), h, xx);
			}
        }
        else if((i != res[0]-1) && ((boundary && valid[xp1]) || !boundary)){  //Can't look backwards in x so do forward difference instead
			xx = (field[xp1](0) - field[index](0)) / h;
			if(i == di && j == dj) {
				printf("xx2: (%f - %f) / %f = %f\n", field[xp1](0), field[index](0), h, xx);
			}
        }
        
        if((j != 0) && ((boundary && valid[ym1]) || !boundary)) {  //Looking backwards in y and valid value
			xy = (field[index](0) - field[ym1](0)) / h;
			if(i == di && j == dj) {
				printf("xy1: (%f - %f) / %f = %f\n", field[index](0), field[ym1](0), h, xy);
			}
        }
        else if((j != res[1]-1) && ((boundary && valid[yp1]) || !boundary)) {  //Can't look backwards in y so do forward difference instead
			xy = (field[yp1](0) - field[index](0)) / h;
			if(i == di && j == dj) {
				printf("xy2: (%f - %f) / %f = %f\n", field[yp1](0), field[index](0), h, xy);
			}
        }
    }
    else {                          //X Velocity is negative, so forward difference
		if(i == di && j == dj) {
			printf("Look X +1\n");
		}
        if((i != res[0]-1) && ((boundary && valid[xp1]) || !boundary)) {
			xx = (field[xp1](0) - field[index](0)) / h;
			if(i == di && j == dj) {
				printf("xx3: (%f - %f) / %f = %f\n", field[xp1](0), field[index](0), h, xx);
			}
        }
        else if((i != 0) && ((boundary && valid[xm1]) || !boundary)) {
			xx = (field[index](0) - field[xm1](0)) / h;
			if(i == di && j == dj) {
				printf("xx4: (%f - %f) / %f = %f\n", field[index](0), field[xm1](0), h, xx);
			}
        }
        
        if((j != res[1]-1) && ((boundary && valid[yp1]) || !boundary)) {
			xy = (field[yp1](0) - field[index](0)) / h;
			if(i == di && j == dj) {
				printf("xy3: (%f - %f) / %f = %f\n", field[yp1](0), field[index](0), h, xy);
			}
        }
        else if((j != 0) && ((boundary && valid[ym1]) || !boundary)) {
			xy = (field[index](0) - field[ym1](0)) / h;
			if(i == di && j == dj) {
				printf("xy4: (%f - %f) / %f = %f\n", field[index](0), field[ym1](0), h, xy);
			}
        }
    }
    if(vel[index](1) >= 0) {    //Y Velocity is positive, so backward difference
		if(i == di && j == dj) {
			printf("Look Y -1\n");
		}
        if((i != 0) && ((boundary && valid[xm1]) || !boundary)) {
			yx = (field[index](1) - field[xm1](1)) / h;
			if(i == di && j == dj) {
				printf("yx1: (%f - %f) / %f = %f\n", field[index](1), field[xm1](1), h, yx);
			}
        }
        else if((i != res[0]-1) && ((boundary && valid[xp1]) || !boundary)) {
			yx = (field[xp1](1) - field[index](1)) / h;
			if(i == di && j == dj) {
				printf("yx2: (%f - %f) / %f = %f\n", field[xp1](1), field[index](1), h, yx);
			}
        }
        
        if((j != 0) && ((boundary && valid[ym1]) || !boundary)) {
			yy = (field[index](1) - field[ym1](1)) / h;
			if(i == di && j == dj) {
				printf("yy1: (%f - %f) / %f = %f\n", field[index](1), field[ym1](1), h, yy);
			}
        }
        else if((j != res[1]-1) && ((boundary && valid[yp1]) || !boundary)) {
			yy = (field[yp1](1) - field[index](1)) / h;
			if(i == di && j == dj) {
				printf("yy2: (%f - %f) / %f = %f\n", field[yp1](1), field[index](1), h, yy);
			}
        }
    }
    else {                          //Y Velocity is negative, so forward difference
		if(i == di && j == dj) {
			printf("Look Y +1\n");
		}
        if((i != res[0]-1) && ((boundary && valid[xp1]) || !boundary)) {
			yx = (field[xp1](1) - field[index](1)) / h;
			if(i == di && j == dj) {
				printf("yx3: (%f - %f) / %f = %f\n", field[xp1](1), field[index](1), h, yx);
			}
        }
        else if((i != 0) && ((boundary && valid[xm1]) || !boundary)) {
			yx = (field[index](1) - field[xm1](1)) / h;
			if(i == di && j == dj) {
				printf("yx4: (%f - %f) / %f = %f\n", field[index](1), field[xm1](1), h, yx);
			}
        }
        
        if((j != res[1]-1) && ((boundary && valid[yp1]) || !boundary)) {
			yy = (field[yp1](1) - field[index](1)) / h;
			if(i == di && j == dj) {
				printf("yy3: (%f - %f) / %f = %f\n", field[yp1](1), field[index](1), h, yy);
			}	
        }
        else if((j != 0) && ((boundary && valid[ym1]) || !boundary)) {
			yy = (field[index](1) - field[ym1](1)) / h;
			if(i == di && j == dj) {
				printf("yy4: (%f - %f) / %f = %f\n", field[index](1), field[ym1](1), h, yy);
			}
        }
    }
    Matrix2d D;
    D << xx, xy, yx, yy;
    return D;
}

void World::distSweep(double *field, std::vector<char> initial, int iters, double eps) {
    double err = 100;
    int it = 0;
    while(it < iters && err > eps) {
        int lowi=0, highi=0, inci=0;
        int lowj=0, highj=0, incj=0;
        switch(it % 4) {
            case 0:
                lowi = 0;
                lowj = 0;
                highi = res[0];
                highj = res[1];
                inci = 1;
                incj = 1;
                break;
            case 1:
                lowi = res[0]-1;
                lowj = 0;
                highi = -1;
                highj = res[1];
                inci = -1;
                incj = 1;
                break;
            case 2:
                lowi = res[0]-1;
                lowj = res[1]-1;
                highi = -1;
                highj = -1;
                inci = -1;
                incj = -1;
                break;
            case 3:
                lowi = 0;
                lowj = res[1]-1;
                highi = res[0];
                highj = -1;
                inci = 1;
                incj = -1;
                break;
        }
        err = 0;
        for(int i = lowi; i != highi; i += inci) {
            for(int j = lowj; j != highj; j += incj) {
                int index = i*res[1]+j;
                if(initial[index]) {
                    continue;
                }
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);
                
                //[(x-a)+]^2 + [(x-b)+]^2 = f^2 * h^2
                //a = u_xmin
                //b = u_ymin
                double a, b;
                if(i == 0) {
                    a = field[xp1];
                }
                else if(i == res[0]-1) {
                    a = field[xm1];
                }
                else {
                    a = std::min(field[xp1], field[xm1]);
                }
                if(j == 0) {
                    b = field[yp1];
                }
                else if(j == res[0]-1) {
                    b = field[ym1];
                }
                else {
                    b = std::min(field[yp1], field[ym1]);
                }
                
                double f = 1;
                double fh = f*h;            //Seeing what happens with f=1
                //Eq 2.4
                double ab = a-b;
                double xbar;
                if(std::abs(ab) >= fh) {
                    xbar = std::min(a, b) + fh;
                }
                else {
                    xbar = (a+b+std::sqrt(2*f*f*h*h-ab*ab)) / 2;
                }
                
                double prevVal = field[index];
                field[index] = std::min(field[index], xbar);
                
                //Keep the max change
                double currDif = prevVal - field[index];
                if(currDif > err) {
                    err = currDif;
                }
            }
        }
        it++;
    }
}

void World::velSweep(Vector2d *vel, double *field, int iters, double eps) {
    double err = 100;
    int it = 0;
    std::vector<char> newvalid = valid;
    while(it < iters && err > eps) {
        int lowi=0, highi=0, inci=0;
        int lowj=0, highj=0, incj=0;
        switch(it % 4) {
            case 0:
                lowi = 0;
                lowj = 0;
                highi = res[0];
                highj = res[1];
                inci = 1;
                incj = 1;
                break;
            case 1:
                lowi = res[0]-1;
                lowj = 0;
                highi = -1;
                highj = res[1];
                inci = -1;
                incj = 1;
                break;
            case 2:
                lowi = res[0]-1;
                lowj = res[1]-1;
                highi = -1;
                highj = -1;
                inci = -1;
                incj = -1;
                break;
            case 3:
                lowi = 0;
                lowj = res[1]-1;
                highi = res[0];
                highj = -1;
                inci = 1;
                incj = -1;
                break;
        }
        err = 0;
        for(int i = lowi; i != highi; i += inci) {
            for(int j = lowj; j != highj; j += incj) {
                int index = i*res[1]+j;
                if(valid[index]) {
                    continue;
                }
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);
                
                //(Fxmin*(phix-phixmin) + Fymin*(phix-phiymin)) / ((phix-phixmin) + (phix-phiymin))
                double a, b;
                Vector2d vi, vj;
                if(i == 0) {
                    a = field[xp1];
                    vi = vel[xp1];
                }
                else if(i == res[0]-1) {
                    a = field[xm1];
                    vi = vel[xm1];
                }
                else {
                    if(field[xp1] < field[xm1]) {
                        a = field[xp1];
                        vi = vel[xp1];
                    }
                    else {
                        a = field[xm1];
                        vi = vel[xm1];
                    }
                }
                if(j == 0) {
                    b = field[yp1];
                    vj = vel[yp1];
                }
                else if(j == res[0]-1) {
                    b = field[ym1];
                    vj = vel[ym1];
                }
                else {
                    if(field[yp1] < field[ym1]) {
                        b = field[yp1];
                        vj = vel[yp1];
                    }
                    else {
                        b = field[ym1];
                        vj = vel[ym1];
                    }
                }
                
                //If neither values are less than x_ij then the contribution is 0
                double phixi = 0, phixj = 0;
                if(field[index] > a) {
                    phixi = field[index]-a;
                }
                if(field[index] > b) {
                    phixj = field[index]-b;
                }
                //Solving eq. at bottom of pg 13 in Adalsteinsson and Sethian 1999
                Vector2d vtmp = (vi*phixi + vj*phixj) / (phixi + phixj);
                double dnorm = (vel[index]-vtmp).norm();
                if(dnorm > err) {
                    err = dnorm;
                }
                if(dnorm < eps) {
                    newvalid[index] = 1;
                }
                vel[index] = vtmp;
            }
        }
        it++;
    }
    valid = newvalid;
}

void World::sweepAve(Eigen::Vector2d *vel, int iters) {
    std::vector<char> oldValid;
    oldValid.resize(res[0]*res[1]);
    //Perform this a couple of times
    int it = 0;
    while(it < iters) {
        oldValid = valid;
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int ind = i*res[1]+j;
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);
                
                Vector2d sum(0,0);
                int count = 0;
                
                tmpVel[ind] = vel[ind];
                //Only check around cells that are not already valid
                if(!oldValid[ind]) {
                    if(i+1 < res[0] && oldValid[xp1]) {
                        sum += vel[xp1];
                        count++;
                    }
                    if(i-1 >= 0 && oldValid[xm1]) {
                        sum += vel[xm1];
                        count++;
                    }
                    if(j+1 < res[1] && oldValid[yp1]) {
                        sum += vel[yp1];
                        count++;
                    }
                    if(j-1 >= 0 && oldValid[ym1]) {
                        sum += vel[ym1];
                        count++;
                    }
                    if(i+1 < res[0] && j+1 < res[1] && oldValid[(i+1)*res[1]+(j+1)]) {
                        sum += vel[(i+1)*res[1]+(j+1)];
                        count++;
                    }
                    if(i+1 < res[0] && j-1 >=0 && oldValid[(i+1)*res[1]+(j-1)]) {
                        sum += vel[(i+1)*res[1]+(j-1)];
                        count++;
                    }
                    if(i-1 >= 0 && j+1 < res[1] && oldValid[(i-1)*res[1]+(j+1)]) {
                        sum += vel[(i-1)*res[1]+(j+1)];
                        count++;
                    }
                    if(i-1 >= 0 && j-1 >= 0 && oldValid[(i-1)*res[1]+(j-1)]) {
                        sum += vel[(i-1)*res[1]+(j-1)];
                        count++;
                    }
                    
                    //If neighboring cells were valid, average the values
                    if(count > 0) {
                        tmpVel[ind] = sum / count;
                        valid[ind] = 1;
                    }
                }
            }
        }
        //Replace grid with temp grid
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int ind = i*res[1]+j;
                vel[ind] = tmpVel[ind];
            }
        }
        it++;
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
    //Causes a segfault
    /// delete [] mass;
    /// delete [] vel;
    /// delete [] tau;
    /// delete [] mat;
    /// delete [] dXx;
    /// delete [] dXy;
    /// delete [] stress;
    /// delete [] F;
    prof.dump<std::chrono::duration<double>>(std::cout);
} 
