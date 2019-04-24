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

// Regular MPM or testing material transfer
#define MAT_TRANSFER
// Splotting material differences or material coordiantes
#define DIFF
// Whether to do APIC mat transfer
#define APIC_MAT

using namespace Eigen;
using benlib::range;

#ifndef NDEBUG
std::ofstream debug;
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

	// dt = root.get("dt", 1.0/30.0).asDouble();
    int frames = root.get("frames", 30).asInt();
    auto dtIn = root["dt"];
    if(dtIn.isNull()) {
        //if dt not there, check if steps is
        auto stepsIn = root["steps"];
        if(stepsIn.isNull()) {
            //nothing there, just set a default
            printf("No steps or dt\n");
            dt = 1.0 / frames;
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
    if(lameIn.size() != 2) {
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
    if(stretchIn.isNull() || compressIn.isNull()){
        plasticEnabled = false;
        std::cout << "no plasticity" << std::endl;
        compression = 0.0;
        stretch = 0.0;
    }
    else {
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
    for(auto i : range(objectsIn.size())) {
  	    objects.emplace_back();
        objType = objectsIn[i].get("type", "square").asString();
        objects[i].type = objType;
		if(objType == "square" || objType == "circle") {
		    auto locationIn = objectsIn[i]["location"];
		    if(locationIn.size() != 2) {
                std::cout<< "bad object location, skipping" << std::endl;
                continue;
		    }
		    objects[i].object(0) = locationIn[0].asDouble();
		    objects[i].object(1) = locationIn[1].asDouble();
		    auto sizeIn = objectsIn[i]["size"];
		    if(sizeIn.size() != 2) {
                std::cout<< "bad object size, skipping" << std::endl;
                objects[i].size[0] = 0;
                objects[i].size[1] = 0;
                continue;
		    }
		    objects[i].size[0] = sizeIn[0].asDouble();
		    objects[i].size[1] = sizeIn[1].asDouble();
		    auto resIn = objectsIn[i]["resolution"];
		    if(resIn.size() != 2) {
                std::cout<< "bad object resolution, skipping" << std::endl;
                objects[i].ores[0] = 1;
                objects[i].ores[0] = 1;
                continue;
		    }
		    objects[i].ores[0] = resIn[0].asInt();
		    objects[i].ores[1] = resIn[1].asInt();
		}
        else {
		    auto const pos = config.find_last_of('/');
		    std::string partfilename = config.substr(0, pos+1);
		    partfilename = partfilename + objectsIn[i].get("filename", "input.bgeo").asString();
		    std::cout<<"loading "<<partfilename<<std::endl;
		    readParticles(partfilename.c_str(), objects[i].particles);
		}
		MaterialProps &mp = objects[i].mp;
		auto lameIn = objectsIn[i]["lame"];
		if(lameIn.size() == 2) {
		    mp.lambda = lameIn[0].asDouble();
		    mp.mu = lameIn[1].asDouble();
		}
        else {
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
        if(colorIn.size() == 3) {
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

    for(size_t o = 0; o < objects.size(); o++) {
        Object& obj = objects[o];
        if(obj.type == "square") {
            center = obj.object + (VectorN(obj.size[0],obj.size[1]) * 0.5);
            //Set up particles at each object vertex
            double diffx = obj.size[0] / (obj.ores[0]-1);
            double diffy = obj.size[1] / (obj.ores[1]-1);
            for(int i = 0; i < obj.ores[0]; i++) {
                for(int j = 0; j < obj.ores[1]; j++) {
                    VectorN pos = obj.object + VectorN(diffx*i, diffy*j);
                    Vector3d col = ((double)j/(obj.ores[1]-1))*obj.color;
                    Particle par(pos, VectorN(0,0), col, obj.mp.pmass);
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
                    VectorN pos = center - VectorN(obj.size[0], obj.size[1]) + VectorN(diffx*i, diffy*j);
                    Vector3d col = ((double)j/(obj.ores[1]-1))*obj.color;
                    VectorN ph = pos - obj.object;
                    if( ((ph(0)*ph(0))/(obj.size[0]*obj.size[0])) + ((ph(1)*ph(1))/(obj.size[1]*obj.size[1])) < 1+EPS) {
                        Particle par(pos, VectorN(0,0), col, obj.mp.pmass);
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
        VectorN avePos = VectorN::Zero();
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
    mass.assign(res[0]*res[1], 0);
    vel.assign(res[0]*res[1], VectorN::Zero());
    velStar.assign(res[0]*res[1], VectorN::Zero());
    prevVel.assign(res[0]*res[1], VectorN::Zero());
    frc.assign(res[0]*res[1], VectorN::Zero());
    mat.assign(res[0]*res[1], VectorN::Zero());
    matdiff.assign(res[0]*res[1], VectorN::Zero());
    weights.assign(res[0]*res[1], 0);
    stress.assign(res[0]*res[1], MatrixN::Zero());
    valid.assign(res[0]*res[1], 0);

    offsets[0] = 0;
    offsets[1] = 1;
    offsets[2] = res[0];
    offsets[3] = res[0]+1;
    G.resize(QUAD, DIM);
    for(int i = 0; i < QUAD; i++) {
        G(i, 0) = i % 2 == 0 ? -1 : 1;
        G(i, 1) = i % 4 <= 1 ? -1 : 1;
        #if DIM == 3
        G(i, 2) = i < 4 ? -1 : 1;
        #endif
    }
    G *= 0.25 * (1.0/h);

    inc = steps;
    #ifndef NDEBUG
    // inc = 125000;
    inc = 2000;
    count = 0;
    sMax = 31;
    matTrans = new VectorN*[sMax];
    for(int i = 0; i < sMax; i++) {
        matTrans[i] = new VectorN[res[0]*res[1]];
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

inline VectorN gradweight(const VectorN &offset, double h) {
	return VectorN(gradweight1d(offset(0))*weight(offset(1))/h,
                    weight(offset(0))*gradweight1d(offset(1))/h);
}

//Keep the max and min? It should never reach past the ends since we're
//putting an invisible wall around the last 2 cells in gridToParticles
inline void bounds(const VectorN &offset, const int res[2], int *xbounds, int *ybounds) {
    /// xbounds[0] = std::max(0, ((int)std::ceil(BOUNDLOWER + offset(0))));
    /// xbounds[1] = std::min(res[0]-1, ((int)std::floor(BOUNDUPPER + offset(0)))+1);
    /// ybounds[0] = std::max(0, ((int)std::ceil(BOUNDLOWER + offset(1))));
    /// ybounds[1] = std::min(res[1]-1, ((int)std::floor(BOUNDUPPER + offset(1)))+1);
    xbounds[0] = ((int)(-2 + offset(0)))+1;
    xbounds[1] = ((int)( 2 + offset(0)))+1;
    ybounds[0] = ((int)(-2 + offset(1)))+1;
    ybounds[1] = ((int)( 2 + offset(1)))+1;
}

inline VectorN interpolate(VectorN point, const std::vector<VectorN>& field, VectorN origin, int res[2], double h) {
    /// wx = x - x1;
    /// wy = y - y1;
    /// v = (1-wx)*(1-wy)*v[i1][j1] + (wx)*(1-wy)*v[i1+1][j1] + (1-wx)*(wy)*v[i1][j1+1] + (wx)*(wy)*v[i1+1][j1+1];
    VectorN x = (point - origin);
    VectorN ij = x / h;
    int i1 = (int)ij(0);
    int j1 = (int)ij(1);
    int i2 = i1+1;
    int j2 = j1+1;
    double x1 = h*i1;
    double y1 = h*j1;
    double wx = (x(0) - x1) / h;
    double wy = (x(1) - y1) / h;
    VectorN v = (1-wx)*(1-wy)*field[i1*res[1]+j1] +
                 (wx)  *(1-wy)*field[i2*res[1]+j1] +
                 (1-wx)*(wy)  *field[i1*res[1]+j2] +
                 (wx)  *(wy)  *field[i2*res[1]+j2];

    return v;
}

inline void smoothField(std::vector<MatrixN>& field, int iters, int res[2], std::vector<char> valid) {
    std::vector<MatrixN> tmpfield(res[0]*res[1]);
    int xp1, xm1, yp1, ym1;
    for(int s = 0; s < iters; s++) {
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                //Extend kernal past end
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
                if(i == res[0]-1 || (i < res[0]-1 && !valid[(i+1)*res[1]+j])) {
                    xp1 = i;
                }
                else {
                    xp1 = i + 1;
                }
                if(i == 0 || (i > 0 && !valid[(i-1)*res[1]+j])) {
                    xm1 = i;
                }
                else {
                    xm1 = i - 1;
                }
                if(j == res[1]-1 || (i < res[1]-1 && !valid[i*res[1]+(j+1)])) {
                    yp1 = j;
                }
                else {
                    yp1 = j + 1;
                }
                if(j == 0 || (j > 0 && !valid[i*res[1]+(j-1)])) {
                    ym1 = j;
                }
                else {
                    ym1 = j - 1;
                }
                //Using kernal:
                //      [0 1 0]
                //1/8 * [1 4 1]
                //      [0 1 0]
                tmpfield[i*res[1]+j] = 0.125*(field[xm1*res[1]+j] + field[xp1*res[1]+j] + field[i*res[1]+ym1] + field[i*res[1]+yp1]) - 0.5*field[i*res[1]+j];
                if(tmpfield[i*res[1]+j].hasNaN()) {
                    printf("Smooth Temp Field NaN\n");
                    std::cout << tmpfield[i*res[1]+j] << "\n";
                    std::cout << field[xm1*res[1]+j] << "\n";
                    std::cout << field[xp1*res[1]+j] << "\n";
                    std::cout << field[i*res[1]+ym1] << "\n";
                    std::cout << field[i*res[1]+yp1] << "\n";
                    std::cout << 4*field[i*res[1]+j] << "\n";
                    std::exit(1);
                }
            }
        }
        for(int i = 0; i < res[0]*res[1]; i++) {
            field[i] = tmpfield[i];
        }
    }
}

inline void smoothField(std::vector<VectorN>& field, int iters, int res[2], std::vector<char> valid) {
    std::vector<VectorN> tmpfield(res[0]*res[1]);
    int xp1, xm1, yp1, ym1;
    for(int s = 0; s < iters; s++) {
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                //Extend kernal past end
                //Using an extension edge case
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
                if(i == res[0]-1 || (i < res[0]-1 && !valid[(i+1)*res[1]+j])) {
                    xp1 = i;
                }
                else {
                    xp1 = i + 1;
                }
                if(i == 0 || (i > 0 && !valid[(i-1)*res[1]+j])) {
                    xm1 = i;
                }
                else {
                    xm1 = i - 1;
                }
                if(j == res[1]-1 || (i < res[1]-1 && !valid[i*res[1]+(j+1)])) {
                    yp1 = j;
                }
                else {
                    yp1 = j + 1;
                }
                if(j == 0 || (j > 0 && !valid[i*res[1]+(j-1)])) {
                    ym1 = j;
                }
                else {
                    ym1 = j - 1;
                }
                //Using kernal:
                //      [0 1 0]
                //1/8 * [1 4 1]
                //      [0 1 0]
                tmpfield[i*res[1]+j] = 0.125*(field[xm1*res[1]+j] + field[xp1*res[1]+j] + field[i*res[1]+ym1] + field[i*res[1]+yp1]) + 0.5*field[i*res[1]+j];
                if(tmpfield[i*res[1]+j].hasNaN()) {
                    printf("Smooth Temp Field NaN\n");
                    std::cout << tmpfield[i*res[1]+j] << "\n";
                    std::cout << 0.125*field[xm1*res[1]+j] << "\n";
                    std::cout << 0.125*field[xp1*res[1]+j] << "\n";
                    std::cout << 0.125*field[i*res[1]+ym1] << "\n";
                    std::cout << 0.125*field[i*res[1]+yp1] << "\n";
                    std::cout << 0.5*field[i*res[1]+j] << "\n";
                    std::exit(1);
                }
            }
        }
        for(int i = 0; i < res[0]*res[1]; i++) {
            field[i] = tmpfield[i];
        }
    }
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
		VectorN offset = (p.x - origin) / h;
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
    // for(int i = 0; i < res[0]; i++) {
    //     for(int j = 0; j < res[1]; j++) {
    //         int ind = i*res[1]+j;
    //         VectorN xg = origin + h*VectorN(i,j);
    //         mat[ind] = xg;
    //     }
    // }
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
    #ifndef NDEBUG
    if(stepNum == 0) {
        std::ostringstream ss;
        ss << std::setw(2) << std::setfill('0') << 0 << "." << std::setw(6)<< 0;
        std::string pframe(ss.str());
        std::string parOut = "mat/particles-" + pframe + ".bgeo";
        writeParticles(parOut.c_str(), objects[0].particles);
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                matTrans[count][i*res[1]+j] = mat[i*res[1]+j];
            }
        }
        count++;
    }
    #endif
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
    // slAdvect();
    // eAdvect();
#endif
    gridToParticles();
    stepNum++;
	elapsedTime += dt;
    #ifndef NDEBUG
    if(stepNum % inc == 0) {
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                matTrans[count][i*res[1]+j] = mat[i*res[1]+j];
            }
        }
        count++;
        if(count == sMax) {
            std::ofstream mats("./mat/invmatcoords.txt");
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
            printf("\nExiting Debug\n");
            std::exit(0);
        }
        std::ostringstream ss;
        ss << std::setw(2) << std::setfill('0') << 0 << "." << std::setw(6)<< stepNum/inc;
        std::string pframe(ss.str());
        std::string parOut = "mat/particles-" + pframe + ".bgeo";
        writeParticles(parOut.c_str(), objects[0].particles);
    }
    #endif
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
    vel = std::vector<VectorN>(res[0]*res[1], VectorN::Zero());
    mass = std::vector<double>(res[0]*res[1], 0);
    mat = std::vector<VectorN>(res[0]*res[1], VectorN::Zero());
    matdiff = std::vector<VectorN>(res[0]*res[1], VectorN::Zero());

    MatrixN tmpD = ((h*h)/3.0)*MatrixN::Identity();
    MatrixN tmpDinv = (3.0/(h*h))*MatrixN::Identity();

	for (unsigned int obj = 0; obj<objects.size(); obj++) {
        std::vector<Particle> &particles = objects[obj].particles;
        for(size_t i = 0; i < particles.size(); i++) {
            Particle &p = particles[i];
            if(stepNum == 0) {
                MatrixN C;
                C << 0, -0.75, 0.75, 0; //Rotational (rotation=0.75) Case
                // C << 0, 0, 0, 0;
                p.B = C * tmpD;
            }
            VectorN offset = (p.x - origin) / h;
            int xbounds[2], ybounds[2];
            bounds(offset, res, xbounds, ybounds);
            VectorN xp = p.x;                                      //particle position
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                double w1 = weight(offset(0) - j);
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    int index = j*res[1] + k;
                    double w = w1*weight(offset(1) - k);
                    mass[index] += w * p.m;
                    #ifndef APIC_MAT
                    mat[index] += w * p.u;
                    matdiff[index] += w * (p.x - p.u);
                    #endif
                    weights[index] += w;
                }
            }
            VectorN mv = p.m*p.v;
            MatrixN mBD = p.m*p.B*tmpDinv;
            #ifdef APIC_MAT
            VectorN mu = p.m*p.u;
            MatrixN mBuD = p.m*p.Bu*tmpDinv;
            VectorN mdisp = p.m*(p.x-p.u);
            MatrixN mBdispD = p.m*p.Bdisp*tmpDinv;
            #endif
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                double w1 = weight(offset(0) - j);
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    double w = w1*weight(offset(1) - k);
                    VectorN xg = origin + h*VectorN(j, k);
                    vel[j*res[1] + k] += w * (mv + mBD*(xg-xp).eval());
                    #ifdef APIC_MAT
                    int index = j*res[1] + k;
                    mat[index] += w * (mu + mBuD*(xg-xp));
                    matdiff[index] += w * (mdisp + mBdispD*(xg-xp));
                    // mat[index] += w * (mu + mBD*(xg-xp));
                    // matdiff[index] += w * (mdisp + mBD*(xg-xp));
                    #endif
                }
            }
        }
	}
	for(int i = 0; i < res[0]; i++) {
		for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            if(mass[index] < EPS) {
                vel[index] = VectorN(0.0, 0.0);
                valid[index] = 0;
                mass[index] = 0;
                mat[index] = VectorN::Zero();
                matdiff[index] = VectorN::Zero();
            }
            else {
                vel[index] /= mass[index];
                valid[index] = 1;
                #ifndef APIC_MAT
                mat[index] /= weights[index];
                matdiff[index] /= weights[index];
                #else
                mat[index] /= mass[index];
                matdiff[index] /= mass[index];
                #endif
            }
            // if(weights[index] < EPS) {
            //     mat[index] = VectorN::Zero();
            //     matdiff[index] = VectorN::Zero();
            // }
            // else {
            //     mat[index] /= weights[index];
            //     matdiff[index] /= weights[index];
            // }
            #ifndef NDEBUG
            if(vel[index].hasNaN()) {
                printf("interpolated vel NaN at (%d, %d)\n", i, j);
                std::cout << vel[index] << std::endl;
                exit(0);
            }
            #endif
        }
	}
    #ifndef NDEBUG
    if(stepNum % inc == 0) {
        debug << "Material Coordinates:\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                int index = i*res[1]+j;
                debug << "(";
                if(matdiff[index](0) < 0) {
                    debug << std::fixed << std::setprecision(5) << matdiff[index](0);
                }
                else {
                    debug << std::fixed << std::setprecision(6) << matdiff[index](0);
                }
                debug << ",";
                if(matdiff[index](1) < 0) {
                    debug << std::fixed << std::setprecision(5) << matdiff[index](1);
                }
                else {
                    debug << std::fixed << std::setprecision(6) << matdiff[index](1);
                }
                debug << ") ";
            }
            debug << "\n";
        }
        std::ofstream materr(std::string("mat/materror-")+std::to_string(stepNum/inc));
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int ind = i*res[1]+j;
                if(!valid[ind]) {
                    continue;
                }
                VectorN xg = origin + h*VectorN(i, j);
                // materr << xg(0) << " " << xg(1) << " " << mat[ind](0)-xg(0) << " " << mat[ind](1)-xg(1) << "\n";
                #ifndef DIFF
                materr << xg(0) << " " << xg(1) << " 0\n";
                materr << mat[ind](0) << " " << mat[ind](1) << " 1\n\n\n";
                #else
                materr << xg(0) << " " << xg(1) << " 0\n";
                materr << xg(0)+mat[ind](0) << " " << xg(1)+mat[ind](1) << " 1\n\n\n";
                #endif
            }
        }
        // std::vector<Particle> &particles = objects[0].particles;
        // for(size_t i = 0; i < particles.size(); i++) {
        //     Particle &p = particles[i];
        //     materr << p.x(0) << " " << p.x(1) << " 0\n";
        //     materr << p.u(0) << " " << p.u(1) << " 1\n\n\n";
        // }
        materr.close();
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
    frc = std::vector<VectorN>(res[0]*res[1], VectorN::Zero());
#ifdef LAGRANGE
	{for(int i = 0; i < (int)objects[0].particles.size(); i++) objects[0].particles[i].f = VectorN::Zero();}
#endif

#ifdef MAT_TRANSFER
    #ifndef NDEBUG
    MatrixN *dF, *dFT, *dA, *dS, *dSt, *dstrain;
    #ifndef DIFF
    MatrixN *dgrad;
    dgrad = new MatrixN[res[0]*res[1]];
    #endif
    dF = new MatrixN[res[0]*res[1]];
    dFT = new MatrixN[res[0]*res[1]];
    dA = new MatrixN[res[0]*res[1]];
    dS = new MatrixN[res[0]*res[1]];
    dSt = new MatrixN[res[0]*res[1]];
    dstrain = new MatrixN[res[0]*res[1]];
    #endif
    #ifdef DIFF
    std::ofstream diffout;
    if(stepNum % inc == 0) {
        diffout.open(std::string("mat/F2-")+std::to_string(stepNum/inc));
    }
    #endif
    //Do finite difference to make the deformation gradient at each grid point
    for(int obj = 0; obj < (int)objects.size(); obj++) {
        MaterialProps &mp = objects[obj].mp;
        #ifndef DIFF
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int index = i*res[1]+j;
                if(!valid[index]) {
                    stress[index] = MatrixN::Zero();
                    #ifndef NDEBUG
                    dgrad[index] = MatrixN::Zero();
                    dF[index] = MatrixN::Zero();
                    dFT[index] = MatrixN::Zero();
                    dA[index] = MatrixN::Zero();
                    dS[index] = MatrixN::Zero();
                    dSt[index] = MatrixN::Zero();
                    dstrain[index] = MatrixN::Zero();
                    #endif
                    continue;
                }
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);
                VectorN dx, dy;
                //First order central difference (x_n+1 - x_n-1) / 2h
                if((j == 0 && valid[yp1]) || (j != 0 && j != res[1]-1 && !valid[ym1] && valid[yp1])) {    //forward
                    dy = (mat[yp1] - mat[index]) / h;
                }
                else if((j == res[1]-1 && valid[ym1]) || (j != 0 && j != res[1]-1 && !valid[yp1] && valid[ym1])) {    //backward
                    dy = (mat[index] - mat[ym1]) / h;
                }
                else if(j != 0 && j != res[1]-1 && valid[yp1] && valid[ym1]) {  //center
                    dy = (mat[yp1] - mat[ym1]) / (2*h);
                }
                else {
                    dy = VectorN(1, 0);
                }

                if((i == 0 && valid[xp1]) || (i != 0 && i != res[0]-1 && !valid[xm1] && valid[xp1])) {    //forward
                    dx = (mat[xp1] - mat[index]) / h;
                }
                else if((i == res[0]-1 && valid[xm1]) || (i != 0 && i != res[0]-1 && !valid[xp1] && valid[xm1])) {    //backward
                    dx = (mat[index] - mat[xm1]) / h;
                }
                else if(i != 0 && i != res[0]-1 && valid[xp1] && valid[xm1]) {  //center
                    dx = (mat[xp1] - mat[xm1]) / (2*h);
                }
                else {
                    dx = VectorN(0, 1);
                }

                //Form inverse deformation gradient [dx dy]
                MatrixN grad;
                grad.col(0) = dx;
                grad.col(1) = dy;
                #ifndef NDEBUG
                if(grad.hasNaN()) {
                    printf("\n%d, %d\n", i, j);
                    printf("dx: (%f, %f), (%f, %f)\n", dx(0), dx(1), grad.col(0)(0), grad.col(0)(1));
                    printf("dy: (%f, %f), (%f, %f)\n", dy(0), dy(1), grad.col(1)(0), grad.col(1)(1));
                    fflush(stdout);
                    std::exit(1);
                }
                #endif
                double det = grad.determinant();
                if(std::abs(det) < 1e-4) {
                    stress[index] = MatrixN::Zero();
                    #ifndef NDEBUG
                    dgrad[index] = grad;
                    dF[index] = MatrixN::Zero();
                    dFT[index] = MatrixN::Zero();
                    dA[index] = MatrixN::Zero();
                    dS[index] = MatrixN::Zero();
                    dSt[index] = MatrixN::Zero();
                    dstrain[index] = MatrixN::Zero();
                    #endif
                }
                else {
                    MatrixN F = grad.inverse();

                    #ifndef NDEBUG
                    if(F.hasNaN()) {
                        printf("%d, %d\n", i, j);
                        printf("Grad:\n");
                        std::cout << grad << "\n" << F << "\n";
                        fflush(stdout);
                        std::exit(1);
                    }
                    #endif
                    MatrixN FT = F.transpose();
                    MatrixN A = FT*F;
                    MatrixN S = A.sqrt();

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
                    MatrixN St = S.transpose();
                    MatrixN strain = 0.5 * (S + St) - MatrixN::Identity();
                    /// MatrixN strain = 0.5 * (FT*F - MatrixN::Identity());
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
                    MatrixN linstress = mp.lambda*strain.trace()*MatrixN::Identity() + 2.0*mp.mu*strain;
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
        #else
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int index = i*res[1]+j;
                if(!valid[index]) {
                    stress[index] = MatrixN::Zero();
                    #ifndef NDEBUG
                    dF[index] = MatrixN::Zero();
                    dFT[index] = MatrixN::Zero();
                    dA[index] = MatrixN::Zero();
                    dS[index] = MatrixN::Zero();
                    dSt[index] = MatrixN::Zero();
                    dstrain[index] = MatrixN::Zero();
                    #endif
                    continue;
                }
                //Taken from the eulerian stuff
                //Computing F from the difference between X and u
                MatrixN F;
                MatrixX coord(DIM, QUAD);
                for (int l = 0; l < QUAD; l++) {
                    coord.col(l) = matdiff[index + offsets[l]];
                }
                MatrixN Finv = MatrixN::Identity() - coord * G;
                if(abs(Finv.determinant()) > 1e-4) {
                    F = Finv.inverse();
                }
                else {
                    F = MatrixN::Identity();
                }
                #ifndef NDEBUG
                if(F.hasNaN()) {
                    printf("%d, %d\n", i, j);
                    printf("F:\n");
                    std::cout << F << "\n";
                    fflush(stdout);
                    std::exit(1);
                }
                if(stepNum % inc == 0) {
                    diffout << i << " " << j << "\n" << F << "\n";
                }
                #endif
                double det = F.determinant();
                if(std::abs(det) < 1e-4) {
                    stress[index] = MatrixN::Zero();
                    #ifndef NDEBUG
                    dF[index] = MatrixN::Zero();
                    dFT[index] = MatrixN::Zero();
                    dA[index] = MatrixN::Zero();
                    dS[index] = MatrixN::Zero();
                    dSt[index] = MatrixN::Zero();
                    dstrain[index] = MatrixN::Zero();
                    #endif
                }
                else {
                    MatrixN FT = F.transpose();
                    MatrixN A = FT*F;
                    MatrixN S = A.sqrt();

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
                    MatrixN St = S.transpose();
                    MatrixN strain = 0.5 * (S + St) - MatrixN::Identity();
                    /// MatrixN strain = 0.5 * (FT*F - MatrixN::Identity());
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
                    MatrixN linstress = mp.lambda*strain.trace()*MatrixN::Identity() + 2.0*mp.mu*strain;
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
        #ifndef NDEBUG
        diffout.close();
        #endif
        #endif
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int index = i*res[1]+j;
                if(!valid[index]) {
                    frc[index] = VectorN::Zero();
                    continue;
                }
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);
                VectorN fx, fy;

                if((j == 0 && valid[yp1]) || (j != 0 && j != res[1]-1 && !valid[ym1] && valid[yp1])) {    //forward
                    fy = (stress[yp1].row(1) - stress[index].row(1)) / h;
                }
                else if((j == res[1]-1 && valid[ym1]) || (j != 0 && j != res[1]-1 && !valid[yp1] && valid[ym1])) {    //backward
                    fy = (stress[index].row(1) - stress[ym1].row(1)) / h;
                }
                else if(j != 0 && j != res[1]-1 && valid[yp1] && valid[ym1]) {  //center
                    fy = (stress[yp1].row(1) - stress[ym1].row(1)) / (2*h);
                }
                else {
                    fy = VectorN(0,0);
                }


                if((i == 0 && valid[xp1]) || (i != 0 && i != res[0]-1 && !valid[xm1] && valid[xp1])) {    //forward
                    fx = (stress[xp1].row(0) - stress[index].row(0)) / h;
                }
                else if((i == res[0]-1 && valid[xm1]) || (i != 0 && i != res[0]-1 && !valid[xp1] && valid[xm1])) {    //backward
                    fx = (stress[index].row(0) - stress[xm1].row(0)) / h;
                }
                else if(i != 0 && i != res[0]-1 && valid[xp1] && valid[xm1]) {  //center
                    fx = (stress[xp1].row(0) - stress[xm1].row(0)) / (2*h);
                }
                else {
                    fx = VectorN(0,0);
                }

                VectorN force(fx(0)+fy(0), fx(1)+fy(1));
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
        if(stepNum % inc == 0) {
            debug << "Valid\n";
            for(int j = res[1]-1; j >= 0; j--) {
                for(int i = 0; i < res[0]; i++) {
                    debug << (int)valid[i*res[1]+j] << " ";
                }
                debug << "\n";
            }
            debug << "\n";
            debug << "Mass\n";
            for(int j = res[1]-1; j >= 0; j--) {
                for(int i = 0; i < res[0]; i++) {
                    debug << std::defaultfloat << mass[i*res[1]+j] << " ";
                }
                debug << "\n";
            }
            debug << "\n";
            #ifndef DIFF
            debug << "Grad:\n";
            for(int j = res[1]-1; j >= 0; j--) {
                for(int k = 0; k < 2; k++) {
                    for(int i = 0; i < res[0]; i++) {
                        int ind = i*res[0]+j;
                        debug << "[";
                        if(dgrad[ind](k,0) < 0) {
                            debug << std::fixed << std::setprecision(5) << dgrad[ind](k,0);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dgrad[ind](k,0);
                        }
                        debug << ",";
                        if(dgrad[ind](k,1) < 0) {
                            debug << std::fixed << std::setprecision(5) << dgrad[ind](k,1);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dgrad[ind](k,1);
                        }
                        debug << "] ";
                    }
                    debug << "\n";
                }
                debug << "\n";
            }
            debug << "\n";
            #endif
            debug << "F:\n";
            for(int j = res[1]-1; j >= 0; j--) {
                for(int k = 0; k < 2; k++) {
                    for(int i = 0; i < res[0]; i++) {
                        int ind = i*res[0]+j;
                        debug << "[";
                        if(dF[ind](k,0) < 0) {
                            debug << std::fixed << std::setprecision(5) << dF[ind](k,0);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dF[ind](k,0);
                        }
                        debug << ",";
                        if(dF[ind](k,1) < 0) {
                            debug << std::fixed << std::setprecision(5) << dF[ind](k,1);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dF[ind](k,1);
                        }
                        debug << "] ";
                    }
                    debug << "\n";
                }
                debug << "\n";
            }
            debug << "\n";

            debug << "FT:\n";
            for(int j = res[1]-1; j >= 0; j--) {
                for(int k = 0; k < 2; k++) {
                    for(int i = 0; i < res[0]; i++) {
                        int ind = i*res[0]+j;
                        debug << "[";
                        if(dFT[ind](k,0) < 0) {
                            debug << std::fixed << std::setprecision(5) << dFT[ind](k,0);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dFT[ind](k,0);
                        }
                        debug << ",";
                        if(dFT[ind](k,1) < 0) {
                            debug << std::fixed << std::setprecision(5) << dFT[ind](k,1);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dFT[ind](k,1);
                        }
                        debug << "] ";
                    }
                    debug << "\n";
                }
                debug << "\n";
            }
            debug << "\n";

            debug << "A:\n";
            for(int j = res[1]-1; j >= 0; j--) {
                for(int k = 0; k < 2; k++) {
                    for(int i = 0; i < res[0]; i++) {
                        int ind = i*res[0]+j;
                        debug << "[";
                        if(dA[ind](k,0) < 0) {
                            debug << std::fixed << std::setprecision(5) << dA[ind](k,0);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dA[ind](k,0);
                        }
                        debug << ",";
                        if(dA[ind](k,1) < 0) {
                            debug << std::fixed << std::setprecision(5) << dA[ind](k,1);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dA[ind](k,1);
                        }
                        debug << "] ";
                    }
                    debug << "\n";
                }
                debug << "\n";
            }
            debug << "\n";

            debug << "S:\n";
            for(int j = res[1]-1; j >= 0; j--) {
                for(int k = 0; k < 2; k++) {
                    for(int i = 0; i < res[0]; i++) {
                        int ind = i*res[0]+j;
                        debug << "[";
                        if(dS[ind](k,0) < 0) {
                            debug << std::fixed << std::setprecision(5) << dS[ind](k,0);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dS[ind](k,0);
                        }
                        debug << ",";
                        if(dS[ind](k,1) < 0) {
                            debug << std::fixed << std::setprecision(5) << dS[ind](k,1);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dS[ind](k,1);
                        }
                        debug << "] ";
                    }
                    debug << "\n";
                }
                debug << "\n";
            }
            debug << "\n";

            debug << "St:\n";
            for(int j = res[1]-1; j >= 0; j--) {
                for(int k = 0; k < 2; k++) {
                    for(int i = 0; i < res[0]; i++) {
                        int ind = i*res[0]+j;
                        debug << "[";
                        if(dSt[ind](k,0) < 0) {
                            debug << std::fixed << std::setprecision(5) << dSt[ind](k,0);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dSt[ind](k,0);
                        }
                        debug << ",";
                        if(dSt[ind](k,1) < 0) {
                            debug << std::fixed << std::setprecision(5) << dSt[ind](k,1);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dSt[ind](k,1);
                        }
                        debug << "] ";
                    }
                    debug << "\n";
                }
                debug << "\n";
            }
            debug << "\n";

            debug << "Strain:\n";
            for(int j = res[1]-1; j >= 0; j--) {
                for(int k = 0; k < 2; k++) {
                    for(int i = 0; i < res[0]; i++) {
                        int ind = i*res[0]+j;
                        debug << "[";
                        if(dstrain[ind](k,0) < 0) {
                            debug << std::fixed << std::setprecision(5) << dstrain[ind](k,0);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dstrain[ind](k,0);
                        }
                        debug << ",";
                        if(dstrain[ind](k,1) < 0) {
                            debug << std::fixed << std::setprecision(5) << dstrain[ind](k,1);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << dstrain[ind](k,1);
                        }
                        debug << "] ";
                    }
                    debug << "\n";
                }
                debug << "\n";
            }
            debug << "\n";

            debug << "Stress:\n";
            for(int j = res[1]-1; j >= 0; j--) {
                for(int k = 0; k < 2; k++) {
                    for(int i = 0; i < res[0]; i++) {
                        int ind = i*res[0]+j;
                        debug << "[";
                        if(stress[ind](k,0) < 0) {
                            debug << std::fixed << std::setprecision(5) << stress[ind](k,0);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << stress[ind](k,0);
                        }
                        debug << ",";
                        if(stress[ind](k,1) < 0) {
                            debug << std::fixed << std::setprecision(5) << stress[ind](k,1);
                        }
                        else {
                            debug << std::fixed << std::setprecision(6) << stress[ind](k,1);
                        }
                        debug << "] ";
                    }
                    debug << "\n";
                }
                debug << "\n";
            }
            debug << "\n";
            debug << "Force:\n";
            for(int j = res[1]-1; j >= 0; j--) {
                for(int i = 0; i < res[0]; i++) {
                    int index = i*res[1]+j;
                    debug << "(";
                    if(frc[index](0) < 0) {
                        debug << std::fixed << std::setprecision(5) << frc[index](0);
                    }
                    else {
                        debug << std::fixed << std::setprecision(6) << frc[index](0);
                    }
                    debug << ",";
                    if(frc[index](1) < 0) {
                        debug << std::fixed << std::setprecision(5) << frc[index](1);
                    }
                    else {
                        debug << std::fixed << std::setprecision(6) << frc[index](1);
                    }
                    debug << ") ";
                }
                debug << "\n";
            }
            debug << "\n\n\n";
            /// if(stepNum == 10*iter) {
                /// debug.close();
                /// std::exit(1);
            /// }
        }
        #ifndef DIFF
        delete[] dgrad;
        #endif
        delete[] dF;
        delete[] dFT;
        delete[] dA;
        delete[] dS;
        delete[] dSt;
        delete[] dstrain;
        #endif
    }
#else
    for (unsigned int obj = 0; obj<objects.size(); obj++) {
        std::vector<Particle> &particles = objects[obj].particles;
        MaterialProps &mp = objects[obj].mp;
    #ifdef LAGRANGE
        for(int j = 0; j < (int)objects[obj].springs.size(); j++) {
            Spring sp = objects[obj].springs[j];
            VectorN I = sp.p0->x - sp.p1->x;
            VectorN II = sp.p0->v - sp.p1->v;
            double inorm = I.norm();
            double idot = II.dot(I);
            VectorN fa = -(sp.ks*(inorm/sp.r-1) + sp.kd*(idot/inorm)) * (I/inorm);
            sp.p0->f += fa;
            sp.p1->f -= fa;
        }
        for(size_t i = 0; i < particles.size(); i++) {
            Particle &p = particles[i];
            VectorN offset = (p.x - origin) / h;
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
            MatrixN gradient = p.gradientE*p.gradientP;
            double J = gradient.determinant();
            if(J < 0) {
                p.color = p.c2;
            }
            else {
                p.color = p.c1;
            }
            MatrixN gradT = gradient.transpose();
            MatrixN eps = 0.5 * (gradT * gradient - MatrixN::Identity());
            double trace = eps.trace();
            MatrixN stress = mp.lambda*trace*MatrixN::Identity() + 2.0*mp.mu*eps;
            MatrixN Ftmp = p.vol * J * stress;

            VectorN offset = (p.x - origin) / h;
            int xbounds[2], ybounds[2];
            bounds(offset, res, xbounds, ybounds);
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    VectorN accumF = Ftmp * gradweight(VectorN(offset(0)-j,offset(1)-k),h);
                    frc[j*res[1] + k] -= accumF;
                    #ifndef NDEBUG
                    if(frc[j*res[1] + k].hasNaN()) {
                        printf("\nf NaN at (%d, %d)\n", j, k);
                        std::cout << "Force:\n" << frc[j*res[1] + k] << std::endl;
                        std::cout << "Volume: " << p.vol << std::endl;
                        std::cout << "Determinant: " << gradient.determinant() << std::endl;
                        std::cout << "Stress:\n" << stress << std::endl;
                        std::cout << "Gradient:\n" << gradweight(VectorN(offset(0)-j, offset(1)-k),h) << std::endl;
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
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1] + j;
            if(stepNum > 0) {
                prevVel[index] = velStar[index];
            }
            /// if(mass[index] < EPS) {
            if(!valid[index]) {
                velStar[index] = VectorN(0, 0);
            }
            else {
                VectorN extfrc(0.0,0.0);
                if(gravityEnabled) {
                    extfrc += mass[index]*gravity;
                }
                if(rotationEnabled) {
                    VectorN d = origin+VectorN(h*i,h*j)-center;
                    extfrc += mass[index]*rotation*VectorN(-d(1), d(0));
                }
                velStar[index] = vel[index] + dt * (1.0/mass[index]) * (frc[index] + extfrc); //dt*g
            }
            if(stepNum == 0) {
                prevVel[index] = velStar[index];
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

        for(size_t i = 0; i < particles.size(); i++) {
            Particle &p = particles[i];
            MatrixN gradV = MatrixN::Zero();
            VectorN offset = (p.x - origin) / h;
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
                    if(gradweight(VectorN(offset(0)-j,offset(1)-k),h).transpose().hasNaN()) {
                        printf("gradV gradW has NaN at (%d, %d)\n", j, k);
                        std::cout << gradweight(VectorN(offset(0)-j,offset(1)-k),h).transpose() << std::endl;
                        exit(0);
                    }
                    #endif
                    MatrixN accumGrad = velStar[index] * gradweight(VectorN(offset(0)-j,offset(1)-k),h).transpose();
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
            /// MatrixN C = p.B * Dinv;
            /// MatrixN fp = MatrixN::Identity() + dt*C;
            MatrixN fp = MatrixN::Identity() + dt*gradV;

            MatrixN tempGradE = fp*p.gradientE;

            if (plasticEnabled){
                MatrixN tempGrad = fp*p.gradientE*p.gradientP;

                JacobiSVD<MatrixN> svd(tempGradE, ComputeFullU | ComputeFullV);

                MatrixN svdU = svd.matrixU();
                VectorN svdSV = svd.singularValues();
                MatrixN svdV = svd.matrixV();

                VectorN sVClamped;
                sVClamped << clamp(svdSV(0), 1-mp.compression, 1+mp.stretch), clamp(svdSV(1), 1-mp.compression, 1+mp.stretch);
                MatrixN svdClamped = sVClamped.asDiagonal();

                p.gradientE = svdU * svdClamped * svdV.transpose();
                p.gradientP = svdV * svdClamped.inverse() * svdU.transpose() * tempGrad;
            } else {
                p.gradientE = tempGradE;
                p.gradientP = MatrixN::Identity();
            }
        }
    }
}

void World::fastSweep(std::vector<double>& field, std::vector<char> valid, double h, int res[2], int iters, double eps) {
    double err = 100;
    int it = 0;
    while(it < iters && err > eps) {
        int lowi, highi, inci;
        int lowj, highj, incj;
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
void World::velExtrapolateFS(std::vector<VectorN>& vel, const std::vector<double>& field, std::vector<char> &valid, int res[2], int iters, double eps) {
    double err = 100;
    int it = 0;
    std::vector<char> newvalid = valid;
    while(it < iters && err > eps) {
        int lowi, highi, inci;
        int lowj, highj, incj;
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
                VectorN vi, vj;
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
                VectorN vtmp = (vi*phixi + vj*phixj) / (phixi + phixj);
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

void World::velExtrapolate() {
    //Pre-extended Velocity Smooth
    /// smoothField(velStar, 15, res, valid);

    #ifndef NDEBUG
    if(stepNum % inc == 0) {
        val = valid;
        std::string name = std::string("mat/oldvelfield-") + std::to_string(stepNum/inc);
        std::ofstream vf(name.c_str());
        double scale = 1;
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                VectorN xg = origin+h*VectorN(i,j);
                VectorN v = scale*velStar[i*res[1]+j];
                vf << xg(0) << " " << xg(1) << " " << v(0) << " " << v(1);
                if(val[i*res[1]+j]) {
                    vf << " 255 0 0";
                }
                else {
                    vf << " 0 0 255";
                }
                vf << "\n";
            }
        }
        vf.close();
    }
    #endif
#ifndef FASTSWEEP
    //Applying simple extrapolation method from https://github.com/christopherbatty/Fluid3D/blob/master/fluidsim.cpp
    std::vector<char> oldValid;
    oldValid.resize(res[0]*res[1]);
    std::vector<VectorN> tmpVel(res[0]*res[1], VectorN::Zero());
    //Perform this a couple of times
    int it = 0, iters = 15;
    double err = 100, eps = 1e-8;
    while(it < iters && err > eps) {
        oldValid = valid;
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int ind = i*res[1]+j;
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);

                VectorN sum(0,0);
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
        }
        err = 0;
        //Replace grid with temp grid
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int ind = i*res[1]+j;
                double dnorm = (velStar[ind]-tmpVel[ind]).norm();
                if(dnorm > err) {
                    err = dnorm;
                }
                velStar[i*res[1]+j] = tmpVel[i*res[1]+j];
            }
        }
        it++;
    }
    //For testing, try exact field
    /// for(int i = 0; i < res[0]; i++) {
        /// for(int j = 0; j < res[1]; j++) {
            /// if(valid[i*res[1]+j]) {
                /// continue;
            /// }
            /// int ind = i*res[1]+j;
            /// VectorN xg = origin+h*VectorN(i,j);
            /// VectorN diff = xg - VectorN(0.0,0.5);
            /// VectorN rot(-diff(1), diff(0));
            /// velStar[ind] = 0.75*rot;
            /// valid[ind] = 1;
        /// }
    /// }
#else
    //Fast Sweeping Zhao 2005
        //In order to do this I'm pretty sure we'll need some special way to track the surface
    std::vector<double> phi(res[0]*res[1], 100);
    for(int i = 0; i < res[0]*res[1]; i++) {
        if(valid[i]) {
            phi[i] = -h/2.0;
        }
    }
    fastSweep(phi, valid, h, res, 6, 1e-9);
    velExtrapolateFS(velStar, phi, valid, res, 6, 1e-9);
#endif
    //Post-extended Velocity Smooth
    /// smoothField(velStar, 15, res, valid);

    #ifndef NDEBUG
    if(stepNum % inc == 0) {
        std::string name = std::string("mat/velfield-") + std::to_string(stepNum/inc);
        std::ofstream vf(name.c_str());
        name = std::string("mat/frcfield-") + std::to_string(stepNum/inc);
        std::ofstream ff(name.c_str());
        double scale = 1;
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                VectorN xg = origin+h*VectorN(i,j);
                VectorN v = scale*velStar[i*res[1]+j];
                vf << xg(0) << " " << xg(1) << " " << v(0) << " " << v(1);
                if(val[i*res[1]+j]) {
                    vf << " 255 0 0";
                }
                else {
                    vf << " 0 0 255";
                }
                vf << "\n";
                VectorN f = scale*frc[i*res[1]+j];
                ff << xg(0) << " " << xg(1) << " " << f(0) << " " << f(1);
                if(val[i*res[1]+j]) {
                    ff << " 255 0 0";
                }
                else {
                    ff << " 0 0 255";
                }
                ff << "\n";
            }
        }
        vf.close();
        ff.close();
        debug << "New Valid\n";
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                debug << (int)valid[i*res[1]+j] << " ";
            }
            debug << "\n";
        }
        debug << "\n";
        // if(stepNum == 10*inc) {
        //     std::exit(0);
        // }
    }
    #endif
}

//Semi-Lagrangian Advection
void World::slAdvect() {
    std::vector<VectorN> newMat(res[0]*res[1]);
    //Advect material coordinates
        //interpolate velocities at v[t] and v[t-dt] and average
        //step back by that amount and interpolate material coordinates around that point
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            if(!valid[index]) {
                newMat[index] = mat[index];
                continue;
            }
            VectorN xg = origin + h*VectorN(i,j);
            //v[t]
            VectorN v1 = velStar[index];
            //step position back by dt
            VectorN ox = xg - dt*v1;
            //v[t-dt]
            VectorN v2 = interpolate(ox, prevVel, origin, res, h);
            //New Vel
            VectorN nv = (v1 + v2) / 2;

            //Step back and interpolate material coordinates
            ox = xg - dt*nv;
            VectorN nmat = interpolate(ox, mat, origin, res, h);
            newMat[index] = nmat;
        }
    }
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            mat[i*res[1]+j] = newMat[i*res[1]+j];
        }
    }
}

//Eulerian Advection
void World::eAdvect() {
    //X = X - dt*D(v)*X
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            if(!valid[index]) {
                continue;
            }
            int xp1 = (i+1)*res[1]+j;
            int xm1 = (i-1)*res[1]+j;
            int yp1 = i*res[1]+(j+1);
            int ym1 = i*res[1]+(j-1);

            #ifndef NDEBUG
            if(stepNum % inc == 0) {
                debug << "D: " << i << ", " << j << "\n";
                if(j != res[1]-1) {
                    debug << "\t\t(" << velStar[yp1](0) << ", " << velStar[yp1](1) << ")\n";
                }
                if(i != 0) {
                    debug << "(" << velStar[xm1](0) << ", " << velStar[xm1](1) << ")  ";
                }
                debug << "(" << velStar[index](0) << ", " << velStar[index](1) << ")  ";
                if(i != res[0]-1) {
                    debug << "(" << velStar[xp1](0) << ", " << velStar[xp1](1) << ")";
                }
                debug << "\n";
                if(j != 0) {
                    debug << "\t\t(" << velStar[ym1](0) << ", " << velStar[ym1](1) << ")\n";
                }
            }
            #endif
            //Form D(v) - Jacobian of velocity with respect to position
            MatrixN D;
            D = upwindJac(velStar, velStar, i, j);

            #ifndef NDEBUG
            if(stepNum % inc == 0) {
                debug << D << "\n";
            }
            #endif

            #ifndef NDEBUG
            if(stepNum % inc == 0) {
                if(D.hasNaN()) {
                    printf("D has NaN: %d, %d\n", i, j);
                    std::cout << D << "\n\n";
                    std::exit(1);
                }
            }
            #endif
            VectorN newMat = mat[index] - dt * D * mat[index];
            #ifndef NDEBUG
            if(newMat.hasNaN()) {
                printf("NewMat has NaN: %d, %d\n", i, j);
                std::cout << mat[index] << "\n-\n" << dt << "\n*\n" << D << "\n*\n" << mat[index] << "\n";
                std::exit(1);
            }
            #endif
            mat[index] = newMat;
        }
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
        for(size_t i = 0; i < particles.size(); i++) {
            Particle &p = particles[i];
            //Update velocities
            p.B = MatrixN::Zero();
            #ifdef APIC_MAT
            p.Bu = MatrixN::Zero();
            p.Bdisp = MatrixN::Zero();
            #endif
            VectorN apic = VectorN::Zero();
            VectorN offset = (p.x - origin) / h;
            int xbounds[2], ybounds[2];
            bounds(offset, res, xbounds, ybounds);
            for(int j = xbounds[0]; j < xbounds[1]; j++) {
                double w1 = weight(offset(0) - j);
                for(int k = ybounds[0]; k < ybounds[1]; k++) {
                    double w = w1*weight(offset(1) - k);
                    int index = j*res[1] + k;
                    if(!valid[index]) {
                        // printf("Particle grabbing outside\n");
                        continue;
                    }
                    VectorN xg = origin + h*VectorN(j, k);
                    VectorN wvel = w * velStar[index];
                    apic += wvel;
                    p.B += wvel * (xg - p.x).transpose();
                    #ifdef APIC_MAT
                    p.Bu += w * mat[index] * (xg - p.x).transpose();
                    p.Bdisp += w * matdiff[index] * (xg - p.x).transpose();
                    #endif
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
            //Do a temperary timestep for candidate position
            //Bilinear for ending value
            //Average and set velocity to that

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

MatrixN World::upwindJac(const std::vector<VectorN>& field, const std::vector<VectorN>& velocity, int i, int j, bool boundary) {
    //Compute D(field) with upwinding scheme
    int index = i*res[1]+j;
    int xp1 = (i+1)*res[1]+j;
    int xm1 = (i-1)*res[1]+j;
    int yp1 = i*res[1]+(j+1);
    int ym1 = i*res[1]+(j-1);
    double xx = 0, xy = 0, yx = 0, yy = 0;      //Upwind scheme so look at each value individually
    if(velocity[index](0) >= 0) {    //X Velocity is positive, so moving from left to right, look left
        if((i != 0) && ((boundary && valid[xm1]) || !boundary)) {  //Looking backwards in x and valid value
			xx = (field[index](0) - field[xm1](0)) / h;
        }
        else if((i != res[0]-1) && ((boundary && valid[xp1]) || !boundary)){  //Can't look backwards in x so do forward difference instead
			xx = (field[xp1](0) - field[index](0)) / h;
        }

        if((j != 0) && ((boundary && valid[ym1]) || !boundary)) {  //Looking backwards in y and valid value
			xy = (field[index](0) - field[ym1](0)) / h;
        }
        else if((j != res[1]-1) && ((boundary && valid[yp1]) || !boundary)) {  //Can't look backwards in y so do forward difference instead
			xy = (field[yp1](0) - field[index](0)) / h;
        }
    }
    else {                          //X Velocity is negative, so forward difference
        if((i != res[0]-1) && ((boundary && valid[xp1]) || !boundary)) {
			xx = (field[xp1](0) - field[index](0)) / h;
        }
        else if((i != 0) && ((boundary && valid[xm1]) || !boundary)) {
			xx = (field[index](0) - field[xm1](0)) / h;
        }

        if((j != res[1]-1) && ((boundary && valid[yp1]) || !boundary)) {
			xy = (field[yp1](0) - field[index](0)) / h;
        }
        else if((j != 0) && ((boundary && valid[ym1]) || !boundary)) {
			xy = (field[index](0) - field[ym1](0)) / h;
        }
    }
    if(velocity[index](1) >= 0) {    //Y Velocity is positive, so backward difference
        if((i != 0) && ((boundary && valid[xm1]) || !boundary)) {
			yx = (field[index](1) - field[xm1](1)) / h;
        }
        else if((i != res[0]-1) && ((boundary && valid[xp1]) || !boundary)) {
			yx = (field[xp1](1) - field[index](1)) / h;
        }

        if((j != 0) && ((boundary && valid[ym1]) || !boundary)) {
			yy = (field[index](1) - field[ym1](1)) / h;
        }
        else if((j != res[1]-1) && ((boundary && valid[yp1]) || !boundary)) {
			yy = (field[yp1](1) - field[index](1)) / h;
        }
    }
    else {                          //Y Velocity is negative, so forward difference
        if((i != res[0]-1) && ((boundary && valid[xp1]) || !boundary)) {
			yx = (field[xp1](1) - field[index](1)) / h;
        }
        else if((i != 0) && ((boundary && valid[xm1]) || !boundary)) {
			yx = (field[index](1) - field[xm1](1)) / h;
        }

        if((j != res[1]-1) && ((boundary && valid[yp1]) || !boundary)) {
			yy = (field[yp1](1) - field[index](1)) / h;
        }
        else if((j != 0) && ((boundary && valid[ym1]) || !boundary)) {
			yy = (field[index](1) - field[ym1](1)) / h;
        }
    }
    Matrix2d D;
    D << xx, xy, yx, yy;
    return D;
}


void writeParticles(const char *fname, const std::vector<Particle> &particles) {
	Partio::ParticlesDataMutable *data = Partio::create();
	Partio::ParticleAttribute xattr;
	Partio::ParticleAttribute uattr;
	Partio::ParticleAttribute battr;
	Partio::ParticleAttribute geattr;
	Partio::ParticleAttribute gpattr;
	Partio::ParticleAttribute cattr;
    Partio::ParticleAttribute cattr2;
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
	data->addAttribute("Cd", Partio::VECTOR, 3);
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
	data->attributeInfo("Cd", cattr);
    data->attributeInfo("color", cattr2);
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
        float *c2 = data->dataWrite<float>(cattr2, i);
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
        c2[0] = p.color(0), c2[1] = p.color(1), c2[2] = p.color(2);
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
	bool color = data->attributeInfo("Cd", cattr);
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
            p.x = VectorN(0.0, 0.0);
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
            p.v = VectorN(0.0, 0.0);
		}
        if (B) {
            float *b = data->dataWrite<float>(battr, i);
            p.B(0,0) = b[0], p.B(0,1) = b[1], p.B(1,0) = b[2], p.B(1,1) = b[3];
        } else {
            p.B = MatrixN::Zero();
        }
        if (gradientE) {
            float *g = data->dataWrite<float>(geattr, i);
            p.gradientE(0,0) = g[0], p.gradientE(0,1) = g[1], p.gradientE(1,0) = g[2], p.gradientE(1,1) = g[3];
        } else {
            p.gradientE = MatrixN::Identity();
        }
        if (gradientP) {
            float *g = data->dataWrite<float>(gpattr, i);
            p.gradientP(0,0) = g[0], p.gradientP(0,1) = g[1], p.gradientP(1,0) = g[2], p.gradientP(1,1) = g[3];
        } else {
            p.gradientP = MatrixN::Identity();
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
    prof.dump<std::chrono::duration<double>>(std::cout);
}
