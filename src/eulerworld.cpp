#include "eulerworld.hpp"
#include <cstdio>
#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>
#include <Partio.h>
#include <Eigen/Geometry>
#include <random>
#include <unordered_map>
#include <utility>
#include "json/json.h"
#include "range.hpp"

using namespace Eigen;
using benlib::range;

#ifndef NDEBUG
std::ofstream debug;
#endif

World::World(std::string config) {
    Eigen::initParallel();
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

    testVel = root.get("testVel", -1).asInt();

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
    // subsampleDh = h / 10.0;
    subsampleDh = 0.00125;
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
    gridToSolution = std::vector<int>(totalCells, -1);
    gridToSolid = std::vector<int>(totalCells);
    tripletList = std::vector<ETriplet>(totalCells);
    restMass = std::vector<double>(totalCells, 0);
    nodeSolidVolumeFraction = std::vector<double>(totalCells);

    vStar = VectorX(DIM*totalCells);                           //updated velocity from forces solve
    solidMv = VectorX(DIM*totalCells);
    solidUnitForce = VectorX(DIM*totalCells);
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
    inc = 100;
    // inc = 1;
    count = 0;
    sMax = 31;
    matTrans = new Vector2d*[sMax];
    for(int i = 0; i < sMax; i++) {
        matTrans[i] = new Vector2d[res[0]*res[1]];
    }
    #endif
}

World::~World() {
    std::cout << "\n\nTimes:\n";
    prof.dump<std::chrono::duration<double>>(std::cout);
    std::cout << "Percentage:\n";
    prof.dumpPercentages(std::cout);
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
            ///Rotation
            // double scale = 1.5;
            double scale = 5;
            if(testVel == 0) {
                if(i == 0 && j == 0) {
                    printf("Vel: Rotation\n");
                }
                velocity[index] = scale*Vector2d(-ph(1), ph(0));
            }
            else if(testVel == 1) { ///Linear
                if(i == 0 && j == 0) {
                    printf("Vel: Linear\n");
                }
                velocity[index] = scale*VectorN(1, 0);
            }
            else { ///None
                if(i == 0 && j == 0) {
                    printf("Vel: None\n");
                }
                velocity[index] = VectorN::Zero();
            }
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
    // gridToSolution.resize(totalCells, -1);
    gridToSolution = std::vector<int>(totalCells, -1);
    int idx = 0;
    for (int x = 1; x < res[0] - 1; x++) {
        for (int y = 1; y < res[1] - 1; y++) {
            int index = x*res[1] + y;
            gridToSolution[index] = idx++;
        }
    }

    //Init material particles
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, std::nextafter(h, std::numeric_limits<double>::max()));
    int particlesPerCell = 32;
    double density = quarterCellSolidMass / (h * h * 0.25);
    double volume = (h * h) / particlesPerCell;

    for(int i = tightMin[0]; i < tightMax[0]; i++) {
        for(int j = tightMin[1]; j < tightMax[1]; j++) {
            int ind = i*res[1] + j;
            if(!valid[ind]) {
                continue;
            }
            VectorN grid = origin + h*VectorN(i, j);
            for(int k = 0; k < particlesPerCell; k++) {
                VectorN jitter = grid + VectorN(dis(gen), dis(gen));
                // double d = interpolate_value(jitter, restMass, origin, res, h);
                double d = density * volume;
                VectorN v = interpolate_value(jitter, velocity, origin, res, h);
                Particle p(jitter, v, color[ind], d);
                matParticles.push_back(p);
            }
        }
    }
    vecWorkspace = velocity;
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
        #ifdef CLUSTER
        std::ofstream colout("/nfs/scratch/adahl1/euler/col0");
        #else
        std::ofstream colout("./euler/col0");
        #endif
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                VectorN coord = X[i*res[1]+j];
                Vector3d c = interpolate(coord, color, Vector3d(135, 206, 255), origin, res, h);
                VectorN pos = origin+h*Vector2d(i,j);
                colout << pos(0) << " " << pos(1) << " " << (int)std::round(c(0)) << " " << (int)std::round(c(1)) << " " << (int)std::round(c(2)) << "\n";
            }
        }
        colout.close();
        #ifdef CLUSTER
        std::ofstream vout("/nfs/scratch/adahl1/euler/vel0");
        #else
        std::ofstream vout("./euler/vel0");
        #endif
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                Vector2d p = origin+h*Vector2d(i,j);
                int ind = i*res[1]+j;
                vout << p(0) << " " << p(1) << " " << velocity[ind](0) << " " << velocity[ind](1) << " 0 0 255" << "\n";
            }
        }
        vout.close();
        printf("Total Cells: %d\n", totalCells);
        printf("Total Active Cells: %d\n", totalActiveCells);
        printf("Total Solid Cells: %d\n", totalSolidCells);
    }
    //Reset and resize variables
    solidUnitForce.conservativeResize(totalSolidCells * DIM);
    solidUnitForce.setZero();
    solidMv.conservativeResize(totalSolidCells * DIM);
    solidMv.setZero();
    dampedMass.conservativeResize(totalSolidCells * DIM);
    dampedMass.setZero();
    tripletList.clear();
    //Compute F and Stiffness Matrix
    {
    auto timer = prof.timeName("ComputeKf");
    computeKf();
    }
    //Compute mv
    {
    auto timer = prof.timeName("computeMV");
    computeMV();
    }
    // for(int i = 0; i < (int)mass.size(); i++) {
    //     printf("%f, %d, %d\n", mass[i], gridToSolid[i], gridToSolution[i]);
    // }
    // std::exit(0);
    //Mass part of rayleigh damping
    {
    auto timer = prof.timeName("Body Forces");
    double massDamp = 1 + dt * rayleighAlpha;
    //Compute body Forces
    #pragma omp parallel for schedule(static) collapse(2)
    for(int i = 1; i < res[0]-1; i++) {
        for(int j = 1; j < res[1]-1; j++) {
            int index = i*res[1]+j;
            // int idx = gridToSolution[index];
            int solidIdx = gridToSolid[index];
            double m = mass[index];

            if(solidIdx < 0) {
                // unitForce.segment<DIM>(idx*DIM) += m * gravity;
            }
            else {
                dampedMass[solidIdx * DIM] = dampedMass[solidIdx * DIM + 1] = m * massDamp;
                #if DIM == 3
                dampedMass[solidIdx * DIM + 2] = m * massDamp;
                #endif
                if(rotationEnabled) {
                    VectorN d = (origin+h*Vector2d(i,j)) - center;
                    solidUnitForce.segment<DIM>(solidIdx * DIM) += m * rotation * VectorN(-d(1), d(0));
                }
                if(gravityEnabled) {
                    solidUnitForce.segment<DIM>(solidIdx * DIM) += m * gravity;
                }
            }
        }
    }
    std::ofstream fout("./euler/frc"+std::to_string(stepNum/inc));
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            Vector2d p = origin+h*Vector2d(i,j);
            int ind = i*res[1]+j;
            Vector3i col(255, 0, 0);
            if(gridToSolid[ind] >= 0) {
                int idx = gridToSolid[ind];
                fout << p(0) << " " << p(1) << " " << solidUnitForce.segment<DIM>(idx*DIM)(0) << " " << solidUnitForce.segment<DIM>(idx*DIM)(1) << " " << col(0) << " " << col(1) << " " << col(2) << "\n";
            }
        }
    }
    fout.close();
    //Construct global system matrix
    systemMat.resize(totalSolidCells*DIM, totalSolidCells*DIM);
    systemMat.setFromTriplets(tripletList.begin(), tripletList.end());
    systemMat *= dt*dt + rayleighBeta*dt;
    for(int i = 0; i < dampedMass.size(); i++) {
        systemMat.coeffRef(i, i) += dampedMass[i];
    }
    systemMat.makeCompressed();
    }
    //Do PCR (conjugate residual) solve
    {
    auto timer = prof.timeName("PCR");
    pcrImplicit();
    }
    //finalize solution
    {
    auto timer = prof.timeName("Finalize");
    finalizeSolution();
    }
    stepNum++;
	elapsedTime += dt;

    if(stepNum % inc == 0) {
        #ifdef CLUSTER
        std::ofstream colout("/nfs/scratch/adahl1/euler/col"+std::to_string(stepNum/inc));
        #else
        std::ofstream colout("./euler/col"+std::to_string(stepNum/inc));
        #endif
        for(int j = res[1]-1; j >= 0; j--) {
            for(int i = 0; i < res[0]; i++) {
                VectorN coord = X[i*res[1]+j];
                Vector3d c = interpolate(coord, color, Vector3d(135, 206, 255), origin, res, h);
                VectorN pos = origin+h*VectorN(i,j);
                colout << pos(0) << " " << pos(1) << " " << (int)std::round(c(0)) << " " << (int)std::round(c(1)) << " " << (int)std::round(c(2)) << "\n";
            }
        }
        colout.close();
        #ifdef CLUSTER
        std::ofstream vout("/nfs/scratch/adahl1/euler/vel"+std::to_string(stepNum/inc));
        #else
        std::ofstream vout("./euler/vel"+std::to_string(stepNum/inc));
        #endif
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                Vector2d p = origin+h*Vector2d(i,j);
                int ind = i*res[1]+j;
                Vector3i col(255, 0, 0);
                // if(cvalid[i*res[1]+j]) {
                //     col = Vector3i(0, 0, 255);
                // }
                vout << p(0) << " " << p(1) << " " << velocity[ind](0) << " " << velocity[ind](1) << " " << col(0) << " " << col(1) << " " << col(2) << "\n";
            }
        }
        vout.close();

        #ifndef NDEBUG
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                matTrans[count][i*res[1]+j] = X[i*res[1]+j];
            }
        }
        count++;
        if(count == sMax) {
            #ifdef CLUSTER
            std::ofstream mats("/nfs/scratch/adahl1/euler/mats.txt");
            #else
            std::ofstream mats("./euler/mats.txt");
            #endif
            for(int i = 0; i < res[0]; i++) {
                for(int j = 0; j < res[1]; j++) {
                    // if(valid[i*res[1]+j]) {
                        for(int k = 0; k < sMax; k++) {
                            mats << matTrans[k][i*res[1]+j](0) << " " << matTrans[k][i*res[1]+j](1) << " " << k << "\n";
                        }
                        mats << "\n\n";
                    // }
                }
            }
            mats.close();
            std::cout << "\n\nTimes:\n";
            prof.dump<std::chrono::duration<double>>(std::cout);
            std::cout << "\n\nPercentage:\n";
            prof.dumpPercentages(std::cout);
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

    double a1, a2;
    //Safeguard against atan2(0, 0) = NaN
    if(std::abs(C) < 1e-12 && std::abs(B) < 1e-12) {
        a1 = 0;
    }
    else {
        a1 = std::atan2(C, B);
    }
    if(std::abs(D) < 1e-12 && std::abs(A) < 1e-12) {
        a2 = 0;
    }
    else {
        a2 = std::atan2(D, A);
    }

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
    #pragma omp parallel for schedule(static)
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

    // Eqn (8) in the paper
    VectorX solidVelocity(solidUnitForce.size());
    solidVelocity.setZero();
    solvePCR(systemMat, solidUnitForce, solidVelocity, systemMatDiag);

    vStar.conservativeResize(totalActiveCells * DIM);
    #pragma omp parallel for schedule(static)
    for(int index = 0; index < totalCells; index++) {
        int idx = gridToSolution[index];
        if (idx < 0) {
            continue;
        }
        int solidIdx = gridToSolid[index];
        if(solidIdx < 0) {
            vStar.segment<DIM>(idx*DIM) = VectorN::Zero();
        }
        else {
            vStar.segment<DIM>(idx*DIM) = solidVelocity.segment<DIM>(solidIdx*DIM);
            // printf("S, %d: (%f, %f) = (%f, %f)\n", idx, vStar.segment<DIM>(idx*DIM)(0), vStar.segment<DIM>(idx*DIM)(1), solidVelocity.segment<DIM>(solidIdx*DIM)(0), solidVelocity.segment<DIM>(solidIdx*DIM)(1));
        }
    }
}

//Preconditioned Conjugate Residual Method
//Preconditioned with Jacobi (M)
bool World::solvePCR(const SparseMat& A, const VectorX& b, VectorX& x, const VectorX& M) {
    double eps = 1e-3;
    int maxIteration = 25;

    residual = b;

    ArrayX Ma = M.array();
    direction = (residual.array() / Ma).matrix();
    residual = direction;

    // s = A * residual;
    parallelMatVec(A, residual, s);
    Ap = s;

    double deltaNew = residual.dot(s);
    double bNorm = b.norm();

    bool converged = false;
    for (int i = 0; i < maxIteration && !converged; i++) {
        q = (Ap.array() / Ma).matrix();

        double alpha = Ap.dot(q);
        if(fabs(alpha) > 1e-9) {
            alpha = deltaNew / alpha;
        }

        x += alpha * direction;
        residual -= alpha * q;
        // s = A * residual;
        parallelMatVec(A, residual, s);
        double deltaOld = deltaNew;

        deltaNew = residual.dot(s);
        double beta = 0;
        if(deltaOld > 1e-9) {
            beta = deltaNew / deltaOld;
        }
        direction *= beta;
        direction += residual;
        Ap *= beta;
        Ap += s;

        // if(i > 0 && i % 5 == 0) {
            // tmp = A * x;
            parallelMatVec(A, x, tmp);
            tmp -= b;
            double tmpNorm = tmp.norm();
            if(tmpNorm < eps * bNorm) {
                converged = true;
                break;
            }
        // }
    }
    return converged;
}

void World::parallelMatVec(const SparseMat& A, const VectorX& b, VectorX& x) {
    x.conservativeResize(A.rows());
    x.setZero();
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < A.rows(); ++k) {
        for (SparseMat::InnerIterator it(A, k); it; ++it) {
            x[k] += it.value() * b[it.col()];
        }
    }
}

void World::finalizeSolution() {
    //Copy new velocity to grid
    {
    auto timer = prof.timeName("\tVelocity Copy");
    #pragma omp parallel for
    for(int i = 0; i < totalCells; i++) {
        int idx = gridToSolution[i];
        if(idx >= 0) {
            velocity[i] = vStar.segment<DIM>(idx * DIM);
        }
    }
    }
    //grid to particles and particles to grid
    {
    auto timer = prof.timeName("\tParticles and Grid");
    #ifndef APICVEL
    gridToParticles();
    particlesToGrid();
    #else
    apicg2p();
    apicp2g();
    #endif
    }
    //extend velocity field
    {
    auto timer = prof.timeName("\tExtend Velocity Field");
    std::vector<double> phifield(res[0]*res[1], res[0]+res[1]+2);
    for(int x = 0; x < res[0]; x++) {
        for(int y = 0; y < res[1]; y++) {
            int ind = x*res[1]+y;
            if(valid[ind]) {
                phifield[ind] = -h/2.0;
            }
        }
    }
    distSweep(phifield, valid, 8, 1e-8);
    velSweep(velocity, phifield, 8, 1e-8);
    }
    //advect displacement field
    {
    auto timer = prof.timeName("\tDisplacement Advection");
    double max_alpha = 0;
    #pragma omp parallel for reduction(max : max_alpha)
    for (int i = 0; i < totalCells; i++) {
        u[i] += velocity[i] * dt;
        double alpha = dt * std::max(std::abs(velocity[i](0))/h, std::abs(velocity[i](1))/h);
        if(alpha > max_alpha) {
            max_alpha = alpha;
        }
    }
    if(max_alpha >= 1.0) {
        printf("\nCFL Violated with alpha of %f\n", max_alpha);
        std::exit(0);
    }
    advectField(u, vecWorkspace);
    }
    //update solid displacements
    // Adjust boundary displacement
    {
    auto timer = prof.timeName("\tBoundary Material");
    for(int i = xmin; i < xmax; i++) {
        u[i*res[1]+ymin] = 2 * u[i*res[1]+(ymin+1)] - u[i*res[1]+(ymin+2)];
        u[i*res[1]+(ymax-1)] = 2 * u[i*res[1]+(ymax-2)] - u[i*res[1]+(ymax-3)];
    }
    for(int j = ymin; j < ymax; j++) {
        u[xmin*res[1]+j] = 2 * u[(xmin+1)*res[1]+j] - u[(xmin+2)*res[1]+j];
        u[(xmax-1)*res[1]+j] = 2 * u[(xmax-2)*res[1]+j] - u[(xmax-3)*res[1]+j];
    }
    }
    // Recover material position field
    {
    auto timer = prof.timeName("\tMaterial Recovery");
    int index = 0;
    for(int x = 0; x < res[0]; x++) {
        for(int y = 0; y < res[1]; y++, index++) {
            VectorN pos(x * h, y * h);
            pos += origin - u[index];
            X[index] = pos;
        }
    }
    vecWorkspace = velocity;
    }
    //compute solid masses
    {
    auto timer = prof.timeName("\tCompute Solid Mass");
    computeSolidMass(true);
    }
    //advect surface mesh

    //compute SDF
    {
    auto timer = prof.timeName("\tCompute SDF");
    computeSDF();
    }
    //gather solid masses
    {
    auto timer = prof.timeName("\tGather Solid Masses");
    gatherSolidMasses();
    }
    //boundary conditions
    {
    auto timer = prof.timeName("\tBoundary Velocity");
    for(int x = 0; x < res[0]; x++) {
        velocity[x*res[1]+0].setZero();
        velocity[x*res[1]+(res[1]-1)].setZero();
    }
    for(int y = 0; y < res[1]; y++) {
        velocity[0*res[1]+y].setZero();
        velocity[(res[0]-1)*res[1]+y].setZero();
    }
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
    int tightBorder = 5;
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
    #pragma omp parallel for collapse(2)
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

    #pragma omp parallel for schedule(static) collapse(2)
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

    #pragma omp parallel for schedule(static) collapse(2)
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
        #pragma omp parallel for schedule(static) collapse(2)
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
    gridToSolid = std::vector<int>(totalCells, -1);
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

void World::gridToParticles() {
    VectorN wallNormals[4] = {VectorN(1, 0), VectorN(0, 1),
                              VectorN(-1, 0), VectorN(0, -1)};
    double boundaries[4] = {origin(0),
                            origin(1),
                            origin(0) + h * (res[0] - 1),
                            origin(1) + h * (res[1] - 1)};
    double alpha = 0.99;
    #pragma omp parallel for schedule(static)
    for(int c = 0; c < (int)matParticles.size(); c++) {
        Particle& p = matParticles[c];
        double x = (p.x(0) - origin(0)) / h;
        double y = (p.x(1) - origin(1)) / h;
        int i = std::floor(x) - 1;
        int j = std::floor(y) - 1;
        VectorN pic(0, 0);
        VectorN flip(0, 0);
        double weights = 0;

        for(int r = i; r < i+4; r++) {
            for(int t = j; t < j+4; t++) {
                if(r < 0 || t < 0 || r > res[0]-1 || t > res[1]-1) {
                    continue;
                }
                int ind = r*res[1]+t;
                if(mass[ind] == 0) {
                    continue;
                }

                VectorN diff = (1.0/h) * (p.x - origin);
                diff -= VectorN(r, t);
                VectorN n = VectorN::Zero();
                for(int a = 0; a < DIM; a++) {
                    if(std::abs(diff(a)) < 1) {
                        n(a) = 0.5 * std::abs(std::pow(diff(a), 3)) - diff(a)*diff(a) + 2.0/3.0;
                    }
                    else if(std::abs(diff(a)) < 2) {
                        n(a) = -1.0/6.0 * std::abs(std::pow(diff(a), 3)) + diff(a)*diff(a) - 2*std::abs(diff(a)) + 4.0/3.0;
                    }
                }
                double w = n(0) * n(1);
                weights += w;
                pic += velocity[ind] * w;
                flip += (velocity[ind] - vecWorkspace[ind]) * w;
            }
        }

        flip += p.v;
        p.v = (1 - alpha) * pic + alpha * flip;
    }
    double wallThickness = 0.015;
    double springConstant = 10000;
    double friction = 0.5;
    #pragma omp parallel for schedule(dynamic)
    for(int c = 0; c < (int)matParticles.size(); c++) {
        Particle& p = matParticles[c];
        VectorN futurePos = p.x;
        double fx = (futurePos(0) - origin(0)) / h;
        double fy = (futurePos(1) - origin(1)) / h;
        int i = std::floor(fx);
        int j = std::floor(fy);
        i = std::min(std::max(i, 0), res[0] - 2);
        j = std::min(std::max(j, 0), res[1] - 2);
        // handle collision with walls
        for(int b = 0; b < 4; b++) {
            double dist = std::abs(p.x(b % 2) - boundaries[b]);
            if(dist < wallThickness) {
                double overlap = wallThickness - dist;
                int xIndex = i;
                if(b == 0) {
                    xIndex = 0;
                }
                else if(b == 2) {
                    xIndex = res[0] - 1;
                }
                int yIndex = j;
                if(b == 1) {
                    yIndex = 0;
                }
                else if(b == 3) {
                    yIndex = res[1] - 1;
                }
                VectorN v_other = velocity[xIndex*res[1]+yIndex];
                VectorN vrel = p.v - v_other;
                double vn = vrel.dot(wallNormals[b]);
                if(vn <= 0) {
                    double impulse = -std::min(dt*springConstant*overlap, p.m*(0.1*overlap/dt-vn));
                    double delta_vn = -impulse / p.m;
                    VectorN vt = vrel - vn * wallNormals[b];
                    VectorN new_vt = std::max(0.1-friction*delta_vn/vt.norm(), 0.0) * vt;
                    p.v = (vn + delta_vn) * wallNormals[b] + new_vt + v_other;
                    break;
                }
            }
        }
        p.x += p.v * dt;
    }
}

void World::particlesToGrid() {
    std::unordered_map<int, std::pair<VectorN, double> > indexToWeight;
    #pragma omp parallel
    {
    std::unordered_map<int, std::pair<VectorN, double> > localMap;
    #pragma omp parallel for schedule(static)
    for(int c = 0; c < (int)matParticles.size(); c++) {
        Particle& p = matParticles[c];
        double x = (p.x(0) - origin(0)) / h;
        double y = (p.x(1) - origin(1)) / h;
        int i = std::floor(x) - 1;
        int j = std::floor(y) - 1;

        for(int s = i; s < i+4; s++) {
            for(int t = j; t < j+4; t++) {
                if(s < 0 || t < 0 || s > res[0]-1 || t > res[1]-1) {
                    continue;
                }
                VectorN diff = (1.0/h) * (p.x - origin);
                diff -= VectorN(s, t);
                VectorN n(0, 0);
                for(int k = 0; k < 2; k++) {
                    if(std::abs(diff[k]) < 1) {
                        n[k] = 0.5*std::abs(std::pow(diff[k], 3)) - diff[k]*diff[k] + 2.0/3.0;
                    }
                    else if(std::abs(diff[k]) < 2) {
                        n[k] = -1.0/6.0*std::abs(std::pow(diff[k], 3)) + diff[k]*diff[k] - 2*std::abs(diff[k]) + 4.0/3.0;
                    }
                }
                double w = n[0] * n[1];
                int index = s*res[1]+t;
                if(localMap.find(index) == localMap.end()) {
                    localMap[index].first = w * p.v;
                    localMap[index].second = w;
                }
                else {
                    localMap[index].first += w * p.v;
                    localMap[index].second += w;
                }
            }
        }
    }
    #pragma omp barrier
    #pragma omp critical
    {
        for(auto iToW : localMap) {
            if(indexToWeight.find(iToW.first) == indexToWeight.end()) {
                indexToWeight[iToW.first] = iToW.second;
            }
            else {
                indexToWeight[iToW.first].first += iToW.second.first;
                indexToWeight[iToW.first].second += iToW.second.second;
            }
        }
    }
    }
    mass.clear();
    static int cutoff = 3;
    for(auto iToW : indexToWeight) {
        if(iToW.second.second < cutoff) {
            continue;
        }
        velocity[iToW.first] = iToW.second.first / iToW.second.second;
        mass[iToW.first] = 1;
    }
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

inline void bounds(const VectorN &offset, const int res[2], int *xbounds, int *ybounds) {
    xbounds[0] = ((int)(-2 + offset(0)))+1;
    xbounds[1] = ((int)( 2 + offset(0)))+1;
    ybounds[0] = ((int)(-2 + offset(1)))+1;
    ybounds[1] = ((int)( 2 + offset(1)))+1;
}

void World::apicg2p() {
    for(size_t i = 0; i < matParticles.size(); i++) {
        Particle &p = matParticles[i];
        //Update velocities
        p.B = MatrixN::Zero();
        VectorN apic = VectorN::Zero();
        VectorN offset = (p.x - origin) / h;
        int xbounds[2], ybounds[2];
        bounds(offset, res, xbounds, ybounds);
        for(int j = xbounds[0]; j < xbounds[1]; j++) {
            double w1 = weight(offset(0) - j);
            for(int k = ybounds[0]; k < ybounds[1]; k++) {
                double w = w1*weight(offset(1) - k);
                int index = j*res[1] + k;
                // if(!valid[index]) {
                //     printf("Particle grabbing outside\n");
                // }
                VectorN xg = origin + h*VectorN(j, k);
                VectorN wvel = w * velocity[index];
                apic += wvel;
                p.B += wvel * (xg - p.x).transpose();
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
        // p.v = mp.massPropDamp * p.v;

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

void World::apicp2g() {
    MatrixN tmpD = ((h*h)/3.0)*MatrixN::Identity();
    MatrixN tmpDinv = (3.0/(h*h))*MatrixN::Identity();
    mass.clear();

    for(size_t i = 0; i < matParticles.size(); i++) {
        Particle &p = matParticles[i];
        if(stepNum == 0) {
            MatrixN C;
            C << 0, -5.0, 5.0, 0; //Rotational (rotation=0.75) Case
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
            }
        }
        VectorN mv = p.m*p.v;
        MatrixN mBD = p.m*p.B*tmpDinv;
        for(int j = xbounds[0]; j < xbounds[1]; j++) {
            double w1 = weight(offset(0) - j);
            for(int k = ybounds[0]; k < ybounds[1]; k++) {
                double w = w1*weight(offset(1) - k);
                VectorN xg = origin + h*VectorN(j, k);
                velocity[j*res[1] + k] += w * (mv + mBD*(xg-xp).eval());
            }
        }
    }
	for(int i = 0; i < res[0]; i++) {
		for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            if(mass[index] < EPS) {
                velocity[index] = VectorN(0.0, 0.0);
                // valid[index] = 0;
                mass[index] = 0;
            }
            else {
                velocity[index] /= mass[index];
                // valid[index] = 1;
            }
            #ifndef NDEBUG
            if(velocity[index].hasNaN()) {
                printf("interpolated vel NaN at (%d, %d)\n", i, j);
                std::cout << velocity[index] << std::endl;
                exit(0);
            }
            #endif
        }
	}
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

void World::velSweep(std::vector<VectorN>& vel, const std::vector<double>& field, int iters, double eps) {
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
	Partio::ParticleAttribute cattr;
	Partio::ParticleAttribute mattr;
	Partio::ParticleAttribute rattr;
	Partio::ParticleAttribute vattr;
    data->addParticles(particles.size());

	data->addAttribute("position", Partio::VECTOR, 3);
	data->addAttribute("velocity", Partio::VECTOR, 2);
	data->addAttribute("B", Partio::VECTOR, 4);
	data->addAttribute("color", Partio::VECTOR, 3);
	data->addAttribute("mass", Partio::FLOAT, 1);
	data->addAttribute("rho", Partio::FLOAT, 1);
	data->addAttribute("vol", Partio::FLOAT, 1);

	data->attributeInfo("position", xattr);
	data->attributeInfo("velocity", uattr);
	data->attributeInfo("B", battr);
	data->attributeInfo("color", cattr);
	data->attributeInfo("mass", mattr);
	data->attributeInfo("rho", rattr);
	data->attributeInfo("vol", vattr);

	for (unsigned int i=0; i < particles.size(); i++) {
		const Particle &p = particles[i];
		float *x = data->dataWrite<float>(xattr, i);
		float *u = data->dataWrite<float>(uattr, i);
		float *b = data->dataWrite<float>(battr, i);
		float *c = data->dataWrite<float>(cattr, i);
		float *m = data->dataWrite<float>(mattr, i);
		float *r = data->dataWrite<float>(rattr, i);
		float *v = data->dataWrite<float>(vattr, i);

		x[0] = p.x(0), x[1] = p.x(1), x[2] = 0;
        u[0] = p.v(0), u[1] = p.v(1);
		b[0] = p.B(0,0), b[1] = p.B(0,1), b[2] = p.B(1,0), b[3] = p.B(1,1);
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
	Partio::ParticleAttribute cattr;
	Partio::ParticleAttribute mattr;
	Partio::ParticleAttribute rattr;
	Partio::ParticleAttribute vattr;

	bool position = data->attributeInfo("position", xattr);
	bool velocity = data->attributeInfo("velocity", uattr);
	bool B = data->attributeInfo("B", battr);
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
