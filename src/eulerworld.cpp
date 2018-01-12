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
    double massPropDamp;
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
    auto rayleighIn = root["rayleigh"];
    if(rayleighIn.size() == 2) {
        rayleighAlpha = rayleighIn[0].asDouble();
        rayleighBeta = rayleighIn[1].asDouble();
    }
    else {
        rayleighAlpha = 0;
        rayleighBeta = 0;
    }
    double alpha;
    alpha = root.get("alpha",0.95).asDouble();

    double stretch, compression;
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

    totalActiveCells = (res[0]-2) * (res[1]-2);
    totalCells = res[0] * res[1];
    subsampleDh = h / 10.0;
    quarterCellSolidMass = 1;
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
    //Sizing vectors
    X = std::vector<VectorN>(totalCells);
    u = std::vector<VectorN>(totalCells);
    mass = std::vector<double>(totalCells);
    velocity = std::vector<VectorN>(totalCells);
    valid = std::vector<char>(totalCells, 0);
    vecWorkspace = std::vector<VectorN>(totalCells);
    phi = std::vector<double>(totalCells);
    detFInv = std::vector<double>(totalCells);
    gridToSolution = std::vector<int>(totalCells);
    gridToSolid = std::vector<int>(totalCells);
    tripletList = std::vector<ETriplet>(totalCells);
    restMass = std::vector<double>(totalCells, 0);
    nodeSolidVolumeFraction = std::vector<double>(totalCells);

    vStar = VectorX(totalCells);                           //updated velocity from forces solve
    mv = VectorX(totalCells);
    solidMv = VectorX(totalCells);
    unitForce = VectorX(totalCells);
    solidUnitForce = VectorX(totalCells);
    nzMassVec = VectorX(totalCells);
    nzMassVecInv = VectorX(totalCells);
    dampedMass = VectorX(totalCells);
    quarterVolumeFractions = VectorX(QUAD*totalCells);

    xmin = 0;
    ymin = 0;
    xmax = res[0];
    ymax = res[1];
    #if DIM == 3
    zmin = 0;
    zmax = res[2];
    tightMin = VectorI(0, 0, 0);
    tightMax = VectorI(res[0], res[1], res[2]);
    #else
    tightMin = VectorI(0, 0);
    tightMax = VectorI(res[0], res[1]);
    #endif

    inc = steps;
    #ifndef NDEBUG
    inc = 1;
    count = 0;
    sMax = 11;
    matTrans = new Vector2d*[sMax];
    for(int i = 0; i < sMax; i++) {
        matTrans[i] = new Vector2d[res[0]*res[1]];
    }
    #endif
}

inline Vector3d interpolate(VectorN point, std::vector<Vector3d> field, Vector3d bg, VectorN origin, int res[2], double h) {
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

template <class T>
inline T lerp(const T& value0, const T& value1, double f) {
    return (1 - f) * value0 + f * value1;
}

template <class T>
inline T bilerp(const T& v00, const T& v10, const T& v01, const T& v11, double fx, double fy) {
    return lerp(lerp(v00, v10, fx), lerp(v01, v11, fx), fy);
}

template <class T>
inline T trilerp(const T& v000, const T& v100, const T& v010, const T& v110,
                       const T& v001, const T& v101, const T& v011, const T& v111,
                       double fx, double fy, double fz) {
    return lerp(bilerp(v000, v100, v010, v110, fx, fy),
                bilerp(v001, v101, v011, v111, fx, fy), fz);
}

inline void get_barycentric(double x, int& i, double& f, int i_low, int i_high, double dh) {
    x /= dh;
    double s = std::floor(x);
    i = int(s);
    if(i < i_low) {
        i = i_low;
        f = 0.5;
    }
    else if(i > i_high) {
        i = i_high;
        f = 0.5;
    }
    else {
        f = x - s;
    }
}

template <class F>
inline F interpolate_value(const VectorN& point, std::vector<F> field, VectorN origin, int res[DIM], double dh) {
    #if DIM == 3
    int i, j, k;
    double fx, fy, fz;
    #else
    int i, j;
    double fx, fy;
    #endif

    get_barycentric(point(0) - origin(0), i, fx, 0, res[0] - 2, dh);
    get_barycentric(point(1) - origin(1), j, fy, 0, res[1] - 2, dh);
    #if DIM == 3
    get_barycentric(point(2) - origin(2), k, fz, 0, res[2] - 2, dh);
    #endif

    #if DIM == 3
    const F& v000 = field(i, j, k);
    const F& v001 = field(i, j, k + 1);
    const F& v010 = field(i, j + 1, k);
    const F& v011 = field(i, j + 1, k + 1);
    const F& v100 = field(i + 1, j, k);
    const F& v101 = field(i + 1, j, k + 1);
    const F& v110 = field(i + 1, j + 1, k);
    const F& v111 = field(i + 1, j + 1, k + 1);
    return trilerp(v000, v100, v010, v110, v001, v101, v011, v111, fx, fy, fz);
    #else
    const F& v00 = field[i*res[1]+j];
    const F& v01 = field[i*res[1]+(j+1)];
    const F& v10 = field[(i+1)*res[1]+j];
    const F& v11 = field[(i+1)*res[1]+(j+1)];
    return bilerp(v00, v10, v01, v11, fx, fy);
    #endif
}

void World::init() {
    color = std::vector<Vector3d>(res[0]*res[1], Vector3d(135,206,255));
    //Init material coordinates
    //(Temp for now) set up mass field
    /// Vector2d middle = origin + h*Vector2d((res[0]-1)/2, (res[1]-1)/2);
    double size[2] = {0.15, 0.15};	// ~1ft wide
    int count = 0;
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            VectorN pos = origin + h*VectorN(i, j);
            VectorN ph = pos - center;
            // For Testing:
            // rotate material field and test D(X)
            // Rotation2D<double> rot(-0.785398);
            // Vector2d rpos = rot * ph;
            // X[index] = rpos + center;
            X[index] = pos;
            // color[index] = Vector3d(135,206,255);
            // velocity[index] = 0.5*Vector2d(-ph(1), ph(0));
            velocity[index] = VectorN::Zero();
            //Create a square in the center
            /// if(i >= (res[0]/2)-4 && i <= (res[0]/2)+4 && j >= (res[1]/2)-4 && j <= (res[1]/2)+4) {
            if( ((ph(0)*ph(0))/(size[0]*size[0])) + ((ph(1)*ph(1))/(size[1]*size[1])) < 1+EPS) {	//circle
                /// mass[index] = 0.067;  //.067kg * 9x9 = 5.44kg ~ 12lbs
                count++;
                //Start it out rotating
                // velocity[index] = 0.5*Vector2d(-ph(1), ph(0));
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
    // double pmass = 18.14/count; // ~40lb / count
    double pmass = 1000;
    for(int i = 0; i < res[0]; i++) {
		for(int j = 0; j < res[1]; j++) {
			int ind = i*res[1]+j;
            VectorN pos = VectorN(i, j) * h + origin;
            u[ind] = pos - X[ind];
			if(valid[ind]) {
				restMass[ind] = pmass;
				//Add a passive particle at the cell
                #ifndef NDEBUG
                Particle p(X[ind], VectorN::Zero(), color[ind], pmass);
                particles.push_back(p);
                #endif
			}
		}
	}
    //Init materials
    computeSolidMass();
    computeSDF();

    gatherSolidMasses();
    gridToSolution.resize(totalCells, -1);
    int idx = 0;
    for (int x = 1; x < res[0] - 1; x++) {
        for (int y = 1; y < res[1] - 1; y++) {
            int index = x*res[1] + y;
            gridToSolution[index] = idx++;
        }
    }
}


/******************************
 * Algorithm Process from Eulerian Solid-Fluid Coupling
 ******************************/
void World::step() {
    if(stepNum == 0) {
        #ifndef NDEBUG
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                matTrans[count][i*res[1]+j] = X[i*res[1]+j];
            }
        }
        count++;
        #endif
        std::ofstream colout("./euler/col0");
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                VectorN coord = X[i*res[1]+j];
                Vector3d c = interpolate(coord, color, Vector3d(135, 206, 255), origin, res, h);
                VectorN pos = origin+h*Vector2d(i,j);
                colout << pos(0) << " " << pos(1) << " " << (int)std::round(c(0)) << " " << (int)std::round(c(1)) << " " << (int)std::round(c(2)) << "\n";
            }
        }
        colout.close();
        printf("Total Cells: %d\n", totalCells);
        printf("Total Active Cells: %d\n", totalActiveCells);
        printf("Total Solid Cells: %d\n", totalSolidCells);
        std::exit(0);
    }
    //Reset and resize variables
    mv.conservativeResize(totalActiveCells * DIM);
    unitForce.conservativeResize(totalActiveCells * DIM);
    unitForce.setZero();
    nzMassVec.conservativeResize(totalActiveCells * DIM);
    nzMassVecInv.conservativeResize(totalActiveCells * DIM);

    solidUnitForce.conservativeResize(totalSolidCells * DIM);
    solidUnitForce.setZero();
    solidMv.conservativeResize(totalSolidCells * DIM);
    solidMv.setZero();
    dampedMass.conservativeResize(totalSolidCells * DIM);
    dampedMass.setZero();
    //Compute F and Stiffness Matrix
    computeKf();
    //Compute mv
    computeMV();
    //Mass part of rayleigh damping
    double massDamp = 1 + dt * rayleighAlpha;
    //Compute body Forces
    for(int i = 1; i < res[0]-1; i++) {
        for(int j = 1; j < res[1]-1; j++) {
            int index = i*res[1]+j;
            int idx = gridToSolution[index];
            int solidIdx = gridToSolid[index];
            double m = mass[index];

            if(solidIdx < 0) {
                unitForce.segment<DIM>(idx*DIM) += m * gravity;
                //Not sure what the nz means
                nzMassVec[idx * DIM] = nzMassVec[idx * DIM + 1] = m;
                #if DIM == 3
                nzMassVec[idx * DIM + 2] = mass;
                #endif
            }
            else {
                dampedMass[solidIdx * DIM] = dampedMass[solidIdx * DIM + 1] = m * massDamp;
                #if DIM == 3
                dampedMass[solidIdx * DIM + 2] = m * massDamp;
                #endif
                solidUnitForce.segment<DIM>(solidIdx * DIM) += m * gravity;
            }
            nzMassVecInv[idx * DIM] = nzMassVecInv[idx * DIM + 1] = 1.0 / m;
            #if DIM == 3
            nzMassVecInv[idx * DIM + 2] = 1.0 / m;
            #endif
        }
    }
    //Construct global system matrix
    systemMat.resize(totalSolidCells*DIM, totalSolidCells*DIM);
    systemMat.setFromTriplets(tripletList.begin(), tripletList.end());
    systemMat *= dt*dt + rayleighBeta*dt;
    for(int i = 0; i < dampedMass.size(); i++) {
        systemMat.coeffRef(i, i) += dampedMass[i];
    }
    systemMat.makeCompressed();
    //Do PCR (conjugate residual) solve
    pcrImplicit();
    //finalize solution
    finalizeSolution();

    stepNum++;
	elapsedTime += dt;

    if(stepNum % inc == 0) {
        std::ofstream colout("./euler/col"+std::to_string(stepNum/inc));
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                VectorN coord = X[i*res[1]+j];
                Vector3d c = interpolate(coord, color, Vector3d(135, 206, 255), origin, res, h);
                VectorN pos = origin+h*VectorN(i,j);
                colout << pos(0) << " " << pos(1) << " " << (int)std::round(c(0)) << " " << (int)std::round(c(1)) << " " << (int)std::round(c(2)) << "\n";
            }
        }
        colout.close();

        #ifndef NDEBUG
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                matTrans[count][i*res[1]+j] = X[i*res[1]+j];
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
            std::exit(0);
        }
        #endif
    }
}

//Compute Stiffness Matrix and Material Forces
void World::computeKf() {
    int startSize = tripletList.size();

    VectorX localForce(solidUnitForce.size());
    localForce.setZero();
    std::vector<ETriplet> localMat;

    for(int i = tightMin[0]; i < tightMax[0]-1; i++) {
        for(int j = tightMin[1]; j < tightMax[1]-1; j++) {
            int index = i*res[1] + j;
            if(!valid[index]) {
                continue;
            }

            MatrixX nodalForcesFEM = MatrixX::Constant(DIM, QUAD, 0);
            MatrixX stiffness = MatrixX::Constant(DIM*QUAD, DIM*QUAD, 0);
            MatrixX coord(DIM, QUAD);
            for (int i = 0; i < QUAD; i++) {
                coord.col(i) = u[index + offsets[i]];
            }
            VectorN restPos(0, 0);
            MatrixX Xloc(DIM, QUAD);
            for (int i = 0; i < QUAD; i++) {
                Xloc.col(i) = X[index + offsets[i]];
                restPos += X[index + offsets[i]];
            }
            restPos /= QUAD;


            MatrixN F;
            double Jinv = 1.0 / computeF(coord, G, F);
            materialInit(F);
            MatrixN P = firstPiolaKirchhoff();
            MatrixN cauchyStress = Jinv * P * F.transpose();

            MatrixX pFpxhat = MatrixX::Constant(DIM*DIM, DIM*QUAD, 0);
            computePFPxhat(F, G, pFpxhat);
            MatrixX GF = G * F;
            #if DIM == 3
            nodalForcesFEM -= cauchyStress * G.transpose() * h * h * h;
            #else
            nodalForcesFEM -= cauchyStress * G.transpose() * h * h;
            #endif

            for (int k = 0; k < DIM*QUAD; k++) {
                int col = index + offsets[k / DIM];
                if(mass[col] == 0) {
                    continue;
                }
                MatrixN dF;
                #if DIM == 3
                dF(0, 0) = pFpxhat(0, k);
                dF(0, 1) = pFpxhat(1, k);
                dF(0, 2) = pFpxhat(2, k);
                dF(1, 0) = pFpxhat(3, k);
                dF(1, 1) = pFpxhat(4, k);
                dF(1, 2) = pFpxhat(5, k);
                dF(2, 0) = pFpxhat(6, k);
                dF(2, 1) = pFpxhat(7, k);
                dF(2, 2) = pFpxhat(8, k);
                #else
                dF(0, 0) = pFpxhat(0, k);
                dF(0, 1) = pFpxhat(1, k);
                dF(1, 0) = pFpxhat(2, k);
                dF(1, 1) = pFpxhat(3, k);
                #endif

                MatrixN deltaP = firstPKDifferential(dF);
                #if DIM == 3
                MatrixX deltaForce = deltaP * GF.transpose() * Jinv * h * h * h;
                #else
                MatrixX deltaForce = deltaP * GF.transpose() * Jinv * h * h;
                #endif
                for (int m = 0; m < QUAD; m++) {
                  stiffness(m * DIM, k) += deltaForce(0, m);
                  stiffness(m * DIM + 1, k) += deltaForce(1, m);
                  #if DIM == 3
                  stiffness(m * DIM + 2, k) += deltaForce(2, m);
                  #endif
                }
            }

            for (int k = 0; k < QUAD; k++) {
                int idx = gridToSolid[index + offsets[k]];
                if (idx >= 0) {
                    localForce.segment<DIM>(idx * DIM) += nodalForcesFEM.col(k);
                }
            }
            for (int k = 0; k < DIM*QUAD; k++) {
                int col = index + offsets[k / DIM];
                if(mass[col] < EPS) {
                    continue;
                }
                if(gridToSolid[col] < 0) {
                    continue;
                }
                int newCol = gridToSolid[col] * DIM + k % DIM;
                for(int l = 0; l < DIM*QUAD; l++) {
                    int row = index + offsets[l / DIM];
                    if(mass[row] < EPS) {
                        continue;
                    }
                    if(gridToSolid[row] < 0) {
                        continue;
                    }
                    int newRow = gridToSolid[row] * DIM + l % DIM;
                    localMat.push_back(ETriplet(newRow, newCol, stiffness(l, k)));
                }
            }
        }
    }
    tripletList.resize(localMat.size() + startSize);
    std::vector<ETriplet>::iterator dest = tripletList.begin() + startSize;
    std::copy(localMat.begin(), localMat.end(), dest);
    solidUnitForce += localForce;
}

//Compute F and return the determinant
double World::computeF(const MatrixX& coord, const MatrixX& pNpF, MatrixN& F) {
    MatrixN Finv = MatrixN::Identity() - coord * pNpF;
    double Jinv = Finv.determinant();
    if(abs(Jinv) > 1e-4) {
        F = Finv.inverse();
    }
    else {
        F = MatrixN::Identity();
        Jinv = 1.0;
    }
    return 1.0 / Jinv;
}

//Initalize Material. Decompose F for use in stress tensor
//Keep R, S, L
void World::materialInit(const MatrixN& F) {
    svd(F, U, Fhat, V);
    R = U * V.transpose();
    S = V * Fhat * V.transpose();

    MatrixN I = MatrixN::Identity();
    MatrixN Ld = (mu * I) + (lambda * (Fhat - I).trace() - 2 * mu) * ((Fhat.trace() * I) - Fhat).inverse();

    Ld(0, 0) = Ld(0, 0) > 0 ? Ld(0, 0) : 0;
    Ld(1, 1) = Ld(1, 1) > 0 ? Ld(1, 1) : 0;
    #if DIM == 3
    Ld(2, 2) = Ld(2, 2) > 0 ? Ld(2, 2) : 0;
    #endif
    L = V * Ld * V.transpose();
}

//Compute the Singular Value Decomposition of F into U * Fhat * VT
//https://scicomp.stackexchange.com/questions/8899/robust-algorithm-for-2x2-svd
void World::svd(const MatrixN& F, MatrixN& U, MatrixN& Fhat, MatrixN& V) {
    // compute the SVD using a 2x2 specific algorithm
    double A = (F(0, 0) + F(1, 1)) / 2.0;
    double B = (F(0, 0) - F(1, 1)) / 2.0;
    double C = (F(1, 0) + F(0, 1)) / 2.0;
    double D = (F(1, 0) - F(0, 1)) / 2.0;

    double Q = std::sqrt(A*A + D*D);
    double R = std::sqrt(B*B + C*C);

    double sx = Q + R;
    double sy = Q - R;

    double a1 = std::atan2(C, B);
    double a2 = std::atan2(D, A);

    double theta = (a2 - a1) / 2.0;
    double phi = (a2 + a1) / 2.0;

    //Sign of sy. 1 if >= 0, -1 otherwise
    int S = (sy >= 0.0) ? 1 : -1;
    //V is rotation by theta
    V(0, 0) = V(1, 1) = std::cos(theta);
    V(1, 0) = std::sin(theta);
    V(0, 1) = -V(1, 0);

    //U is a rotation by phi
    U(0, 0) = std::cos(phi);
    U(1, 1) = S * U(0, 0);
    U(1, 0) = std::sin(phi);
    U(0, 1) = -S * U(1, 0);

    Fhat.setZero();
    Fhat(0, 0) = sx;
    Fhat(1, 1) = std::abs(sy);


    //More standardized way
    // MatrixN Fnormal3 = F.transpose() * F;
    // Eigen::JacobiSVD<MatrixN> eigenSystem(Fnormal3, Eigen::ComputeFullV);
    // VectorN eigenvalues = eigenSystem.singularValues();
    // V = eigenSystem.matrixV();
    //
    // VectorN oEig = eigenvalues;
    // MatrixN oV = V;
    //
    // // populate the diagonal
    // Fhat.setZero();
    // Fhat(0, 0) = std::sqrt(eigenvalues(0));
    // Fhat(1, 1) = std::sqrt(eigenvalues(1));
    // #if DIM == 3
    // Fhat(2, 2) = std::sqrt(eigenvalues(2));
    // #endif
    //
    // if(std::isnan(Fhat(0, 0))) Fhat(0, 0) = 0.0;
    // if(std::isnan(Fhat(1, 1))) Fhat(1, 1) = 0.0;
    // #if DIM == 3
    // if(std::isnan(Fhat(2, 2))) Fhat(2, 2) = 0.0;
    // #endif
    //
    // // if V is a reflection, multiply a column by -1
    // // to ensure it is a rotation
    // double detV = V.determinant();
    // if (detV < 0.0) {
    //     V(0, 0) *= -1.0;
    //     V(1, 0) *= -1.0;
    //     #if DIM == 3
    //     V(2, 0) *= -1.0;
    //     #endif
    // }
    //
    // // compute U
    // U.setZero();
    // U(0, 0) = (Fhat(0, 0) > 0.0f) ? 1.0f / Fhat(0, 0) : 0.0f;
    // U(1, 1) = (Fhat(1, 1) > 0.0f) ? 1.0f / Fhat(1, 1) : 0.0f;
    // #if DIM == 3
    // U(2, 2) = (Fhat(2, 2) > 0.0f) ? 1.0f / Fhat(2, 2) : 0.0f;
    // #endif
    // U = F * V * U;
    // orthogonalizeU(U, Fhat);
    //
    // // correct any reflections in U to ensure it is a rotation
    // if (F.determinant() < 0.0) removeUReflection(U, Fhat);
}

//First Piola Kirchhoff Stress Tensor
MatrixN World::firstPiolaKirchhoff() {
    MatrixN E = S - MatrixN::Identity();
    MatrixN tmp = 2.0 * mu * E;
    double tr = lambda * E.trace();
    tmp(0, 0) += tr;
    tmp(1, 1) += tr;
    #if DIM == 3
    tmp(2, 2) += tr;
    #endif
    return R * tmp;
}

MatrixN World::firstPKDifferential(const MatrixN& dF) {
    MatrixN dFhat = R.transpose() * dF;
    MatrixN dFsym = 0.5 * (dFhat + dFhat.transpose());
    #if DIM == 3
    MatrixN dFskew = 0.5 * (dFhat - dFhat.transpose());
    #endif

    MatrixN dPsym = (2 * mu * dFsym) + (lambda * dFsym.trace() * MatrixN::Identity());

    #if DIM == 3
    VectorN f(-dFskew(1, 2) + dFskew(2, 1), -dFskew(2, 0) + dFskew(0, 2), -dFskew(0, 1) + dFskew(1, 0));
    #else
    // VectorN f(-dFskew(1, 0) + dFskew(0, 1), dFskew(1, 0) - dFskew(0, 1));
    #endif

    double tr = lambda * dFhat.trace();
    dPsym(0, 0) += tr;
    dPsym(1, 1) += tr;
    #if DIM == 3
    dPsym(2, 2) += tr;
    #endif

    #if DIM == 3
    MatrixN dPskew = crossMatrix(L*f);
    #else
    MatrixN dPskew = MatrixN::Zero();
    #endif
    //z axis is not doing anything so the skew would be 0.
    //any part of dFskew(?,2) and dFskew(2,?) is going to be 0. f ends up with (0,0,x)
    //L is 0 in 3rd dimention, so L*f is 0. cross product matrix is thus a matrix of 0
    //Adding 0 doesn't change anything
    MatrixN deltaP = dPsym + dPskew;
    return R * deltaP;
}

#if DIM == 3
MatrixN World::crossMatrix(const VectorN& vec) {
    MatrixN final;
    final(0, 0) = 0.0;
    final(1, 0) = vec[2];
    final(2, 0) = -vec[1];

    final(0, 1) = -vec[2];
    final(1, 1) = 0.0;
    final(2, 1) = vec[0];

    final(0, 2) = vec[1];
    final(1, 2) = -vec[0];
    final(2, 2) = 0.0;

    return final;
 }
 #endif

void World::computePFPxhat(const MatrixN& F, const MatrixX& pNpx, MatrixX& pFpxhat) {
    MatrixX GF = pNpx * F;
    for(int i = 0; i < QUAD; i++) {
        for(int j = 0; j < DIM; j++) {
            pFpxhat(j * DIM, i * DIM + j) = GF(i, 0);
            pFpxhat(j * DIM + 1, i * DIM + j) = GF(i, 1);
            #if DIM == 3
            pFpxhat(j * DIM + 2, i * DIM + j) = GF(i, 2);
            #endif
        }
    }
}

void World::computeMV() {
    for(int i = xmin; i < xmax-1; i++) {
        for(int j = ymin; j < ymax-1; j++) {
            int index = i * res[1] + j;
            int idx = gridToSolid[index];
            if (idx >= 0) {
              solidMv.segment<DIM>(idx * DIM) += mass[index] * velocity[index];
            }
        }
    }
}

//Solve eq. 8 for the velocity from the forces
void World::pcrImplicit() {
    systemMatDiag = systemMat.diagonal();
    solidUnitForce *= dt;
    solidUnitForce += solidMv;
    unitForce *= dt;
    unitForce += mv;

    // Eqn (8) in the paper
    VectorX solidVelocity(solidUnitForce.size());
    solvePCR(systemMat, solidUnitForce, solidVelocity, systemMatDiag);

    vStar.conservativeResize(totalActiveCells * DIM);
    for(int index = 0; index < totalCells; index++) {
        int idx = gridToSolution[index];
        if (idx < 0) {
            continue;
        }
        int solidIdx = gridToSolid[index];
        if(solidIdx < 0) {
            vStar.segment<DIM>(idx*DIM) = unitForce.segment<DIM>(idx*DIM) * nzMassVecInv[idx*DIM];
        }
        else {
            vStar.segment<DIM>(idx*DIM) = solidVelocity.segment<DIM>(solidIdx*DIM);
        }
    }
}

//Preconditioned Conjugate Residual Method
//Preconditioned with Jacobi (M)
bool World::solvePCR(const SparseMat& A, const VectorX& b, VectorX& x, const VectorX& M) {
    double eps = 1e-3;
    int maxIteration = 25;
    x.conservativeResize(b.size());
    x.setZero();

    VectorX q, s, Ap, tmp, residual, direction;
    q.conservativeResize(b.size());
    q.setZero();
    residual = b;

    s.conservativeResize(b.size());
    s.setZero();
    direction.conservativeResize(b.size());
    direction.setZero();
    tmp.conservativeResize(b.size());
    tmp.setZero();

    ArrayN Ma = M.array();
    direction = (residual.array() / Ma).matrix();
    residual = direction;

    s = A * residual;

    Ap = s;

    double deltaNew = residual.dot(s);

    double bNorm = b.norm();

    bool converged = false;
    for (int i = 0; i < maxIteration && !converged; i++) {
        q = (Ap.array() / Ma).matrix();

        double alpha = Ap.dot(q);
        if(fabs(alpha) > 0.0) {
            alpha = deltaNew / alpha;
        }

        x += alpha * direction;

        residual -= alpha * q;

        s = A * residual;

        double deltaOld = deltaNew;

        deltaNew = residual.dot(s);

        double beta = deltaNew / deltaOld;

        direction *= beta;
        direction += residual;

        Ap *= beta;
        Ap += s;

        if(i > 0 && i % 5 == 0) {
            tmp = A * x;
            tmp -= b;
            double tmpNorm = tmp.norm();
            if(tmpNorm < eps * bNorm) {
                converged = true;
                break;
            }
        }
    }
    return converged;
}

void World::finalizeSolution() {
    //Copy new velocity to grid
    for(int i = 0; i < totalCells; i++) {
        int idx = gridToSolution[i];
        if(idx >= 0) {
            velocity[idx] = vStar.segment<DIM>(idx * DIM);
        }
    }
    //extend velocity field
    for (int i = 0; i < totalCells; i++) {
        u[i] += velocity[i] * dt;
    }
    advectField(u, vecWorkspace);
    //update solid displacements
    // Adjust boundary displacement
    for(int i = xmin; i < xmax; i++) {
        u[i*res[1]+ymin] = 2 * u[i*res[1]+(ymin+1)] - u[i*res[1]+(ymin+2)];
        u[i*res[1]+(ymax-1)] = 2 * u[i*res[1]+(ymax-2)] - u[i*res[1]+(ymax-3)];
    }
    for(int j = ymin; j < ymax; j++) {
        u[xmin*res[1]+j] = 2 * u[(xmin+1)*res[1]+j] - u[(xmin+2)*res[1]+j];
        u[(xmax-1)*res[1]+j] = 2 * u[(xmax-2)*res[1]+j] - u[(xmax-3)*res[1]+j];
    }

    // Recover material position field
    int index = 0;
    for(int x = 0; x < res[0]; x++) {
        for(int y = 0; y < res[1]; y++, index++) {
            VectorN pos(x * h, y * h);
            pos += origin - u[index];
            X[index] = pos;
        }
    }
    vecWorkspace = velocity;
    //compute solid masses
    computeSolidMass(true);
    //advect surface mesh

    //compute SDF
    computeSDF();
    //gather solid masses
    gatherSolidMasses();
    //boundary conditions
    for(int x = 0; x < res[0]; x++) {
        velocity[x*res[1]+0].setZero();
        velocity[x*res[1]+(res[1]-1)].setZero();
    }
    for(int y = 0; y < res[1]; y++) {
        velocity[0*res[1]+y].setZero();
        velocity[(res[0]-1)*res[1]+y].setZero();
    }
}

void World::computeSDF() {
    double large_distance = res[0] + res[1] + 2;

    xmax = ymax = 0;
    xmin = res[0];
    ymin = res[1];

    phi = std::vector<double>(res[0]*res[1], large_distance);

    int index = 0;
    for(int x = 0; x < res[0]; x++) {
        for(int y = 0; y < res[1]; y++, index++) {
            if(mass[index] > EPS) {
                valid[index] = 1;
                phi[index] = -0.5;
                if(x < xmin) {
                    xmin = x;
                }
                if(x > xmax) {
                    xmax = x;
                }
                if(y < ymin) {
                    ymin = y;
                }
                if(y > ymax) {
                    ymax = y;
                }
            }
            else {
                valid[index] = 0;
            }
        }
    }

    //Number of cells to extend the boarder
    int tightBorder = 3;
    int extrap = 20;

    #if DIM == 3
    tightMin = VectorI(std::max(xmin - tightBorder, 0), std::max(ymin - tightBorder, 0), std::max(zmin - tightBorder, 0));
    tightMax = VectorI(std::min(xmax + tightBorder, res[0]), std::min(ymax + tightBorder, res[1]), std::min(zmax + tightBorder, zRes));
    #else
    tightMin = VectorI(std::max(xmin - tightBorder, 0), std::max(ymin - tightBorder, 0));
    tightMax = VectorI(std::min(xmax + tightBorder, res[0]), std::min(ymax + tightBorder, res[1]));
    #endif

    xmin = std::max(xmin - extrap, 0);
    ymin = std::max(ymin - extrap, 0);
    #if DIM == 3
    zmin = std::max(zmin - extrap, 0);
    #endif

    xmax = std::min(xmax + extrap, res[0]);
    ymax = std::min(ymax + extrap, res[1]);
    #if DIM == 3
    zmax = std::min(zmax + extrap, res[2]);
    #endif

    //Sweep in 4 directions twice, so 8 iterations
    distSweep(phi, valid, 8, 1e-6);
}

void World::advectField(std::vector<VectorN>& field, std::vector<VectorN>& workspace) {
    for(int x = 0; x < res[0]; x++) {
        for (int y = 0; y < res[1]; y++) {
            int index = x * res[1] + y;
            VectorN node(x, y);
            node = node * h + origin;
            VectorN newPos = node - velocity[index] * dt;
            workspace[index] = interpolate_value(newPos, field, origin, res, h);
        }
    }
    field.swap(workspace);
}

void World::computeSolidMass(bool usesdf) {
    mass = std::vector<double>(totalCells, 0);
    valid = std::vector<char>(totalCells, 0);
    detFInv = std::vector<double>(totalCells, 0);
    quarterVolumeFractions.setZero();

    double subsampleVolume = subsampleDh * subsampleDh;
    int massSampleRes = h / subsampleDh;
    int halfMassSampleRes = massSampleRes >> 1;
    double fraction = 1.0 / (massSampleRes * massSampleRes);

    for(int x = tightMin[0]; x < tightMax[0]-1; x++) {
        for(int y = tightMin[1]; y < tightMax[1]-1; y++) {
            int index = x * res[1] + y;
            MatrixX coord(DIM, QUAD);
            for(int i = 0; i < QUAD; i++) {
                coord.col(i) = u[index + offsets[i]];
            }
            MatrixN F;
            detFInv[index] = std::abs(1.0 / computeF(coord, G, F));
        }
    }

    for(int x = tightMin[0]; x < tightMax[0]-1; x++) {
        for(int y = tightMin[1]; y < tightMax[1]-1; y++) {
            int index = x * res[1] + y;
            VectorN node(x*h+origin(0), y*h+origin(1));
            if(usesdf && phi[index] > 2.0) {
                continue;
            }
            for(int p = 0; p < massSampleRes; p++) {
                for(int q = 0; q < massSampleRes; q++) {
                    VectorN pos = node + VectorN((p+0.5)*subsampleDh, (q+0.5)*subsampleDh);
                    VectorN rPos = interpolate_value(pos, X, origin, res, h);
                    double d = interpolate_value(rPos, restMass, origin, res, h) * subsampleVolume;
                    if(d < 1e-6) {
                        continue;
                    }
                    int offset = p < halfMassSampleRes ? 0 : 1;
                    offset += q < halfMassSampleRes ? 0 : 2;
                    quarterVolumeFractions[index*QUAD+offset] += fraction;
                }
            }
            double sum = 0;
            for(int i = 0; i < QUAD; i++) {
                sum += quarterVolumeFractions[index*QUAD+i];
            }
            double threshold = 0.5;
            if(sum > threshold) {
                for(int i = 0; i < QUAD; i++) {
                    quarterVolumeFractions[index*QUAD+i] = 1;
                }
                valid[index] = 1;
            }
            else {
                if(sum > EPS) {
                    valid[index] = 2;
                }
                for(int i = 0; i < QUAD; i++) {
                    quarterVolumeFractions[index*QUAD+i] = 0;
                }
            }
        }
    }

    for(int i = 0; i < QUAD; i++) {
        for(int x = tightMin[0]; x < tightMax[0]-1; x++) {
            for(int y = tightMin[1]; y < tightMax[1]-1; y++) {
                int index = x * res[1] + y;
                if(!valid[index]) {
                    continue;
                }
                mass[index+offsets[i]] += quarterVolumeFractions[index*QUAD+i] * detFInv[index] * quarterCellSolidMass;
            }
        }
    }
}

void World::gatherSolidMasses() {
    // nodeSolidVolumeFraction.clear();
    //
    // for (int i = 0; i < QUAD; i++) {
    //     for (int x = 1; x < res[0] - 1; x++) {
    //         for (int y = 1; y < res[1] - 1; y++) {
    //             int index = x*res[1]+y;
    //             nodeSolidVolumeFraction[index+offsets[i]] += quarterVolumeFractions[index*QUAD+i];
    //         }
    //     }
    // }
    // for(int i = 0; i < (int)nodeSolidVolumeFraction.size(); i++) {
    //     nodeSolidVolumeFraction[i] *= 0.125;
    // }

    gridToSolid.resize(totalCells, -1);
    int idx = 0;
    for(int x = 1; x < res[0] - 1; x++) {
        for (int y = 1; y < res[1] - 1; y++) {
            int index = x*res[1]+y;
            if(mass[index] > 0) {
                gridToSolid[index] = idx++;
            }
            else {
                gridToSolid[index] = -1;
            }
        }
    }
    totalSolidCells = idx;
}

/************************
 * Helpers
 ************************/

void World::distSweep(std::vector<double>& field, const std::vector<char>& initial, int iters, double eps) {
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
