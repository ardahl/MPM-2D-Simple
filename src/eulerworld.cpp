#include "eulerworld.hpp"
#include <cstdio>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <Partio.h>
// #include <Eigen/Geometry>
#include <unordered_map>
#include <utility>
#include <tuple>
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
	// read global material properties, these may be overridden per object
    auto lameIn = root["lame"];
    if (lameIn.size() != 2) {
        std::cout<< "bad lame, skipping" << std::endl;
    }
    else {
        lambda = lameIn[0].asDouble();
        mu = lameIn[1].asDouble();
    }
    auto rayleighIn = root["rayleigh"];
    if(rayleighIn.size() == 2) {
        rayleighAlpha = rayleighIn[0].asDouble();
        rayleighBeta = rayleighIn[1].asDouble();
    }
    else {
        rayleighAlpha = 0;
        rayleighBeta = 0;
    }

    auto colorIn = root["color"];
    if (colorIn.size() != 3) {
        c = Vector3d(1, 0, 0);
    }
    else {
        c = Eigen::Vector3d(colorIn[0].asDouble(), colorIn[1].asDouble(), colorIn[2].asDouble());
    }

    auto gridIn = root["grid"];
    {
        auto originIn = gridIn["origin"];
        if(originIn.size() != 2) {
		    origin = VectorN(0.0,0.0);
        }
        else {
            origin = VectorN(originIn[0].asDouble(),originIn[1].asDouble());
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
        dh = gridIn.get("dh", 1.0).asDouble();
        dhInv = 1.0 / dh;
		VectorN end(0.0,0.0);
		auto endIn = gridIn["end"];
		if(endIn.size() == 2) {
			end = VectorN(endIn[0].asDouble(), endIn[1].asDouble());
			VectorN dif = end-origin;
			dh = dif(0)/(res[0]-1);
            dhInv = 1.0 / dh;
		}
        center = ((dh*VectorN(res[0]-1, res[1]-1)) / 2.0) + origin;
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

    auto objectsIn = root["objects"];
    for(auto i : range(objectsIn.size())) {
        //Object Properties
        //lambda, mu
        //density
        //plastic yield
        //pmass
        //size
        //center
        //test vel
        double olambda, omu;
		auto lameIn = objectsIn[i]["lame"];
		if(lameIn.size() == 2) {
		    olambda = lameIn[0].asDouble();
		    omu = lameIn[1].asDouble();
		}
        else {
		    olambda = lambda;
		    omu = mu;
		}
        double density = objectsIn[i].get("density", 1000).asDouble();
        double yield = objectsIn[i].get("plasticYield", -1).asDouble();
        double pmass = objectsIn[i].get("pmass", 1000).asDouble();
        double size[2] = {0.15, 0.15};
        auto sizeIn = objectsIn[i]["size"];
        if (sizeIn.size() == 2) {
            size[0] = sizeIn[0].asDouble();
            size[1] = sizeIn[1].asDouble();
        }
        VectorN center;
        auto centerIn = objectsIn[i]["center"];
        if(centerIn.size() != 2) {
		    center = Vector2d(0.0,0.0);
        }
        else {
            center = Vector2d(centerIn[0].asDouble(),centerIn[1].asDouble());
        }
        std::vector<double> restMass(res[0]*res[1], 0);
        for(int i = 0; i < res[0]; i++) {
    		for(int j = 0; j < res[1]; j++) {
    			int index = i*res[1]+j;
                VectorN pos = VectorN(i, j) * dh + origin;
                VectorN ph = pos - center;
                if( ((ph(0)*ph(0))/(size[0]*size[0])) + ((ph(1)*ph(1))/(size[1]*size[1])) < 1+EPS) {
                    restMass[index] = pmass;
                }
    		}
    	}
		Object* solid = new Object(res, dh, origin, restMass);
        solid->center = center;
        solid->size[0] = size[0];
        solid->size[1] = size[1];
        solid->setMaterial(olambda, omu);
        solid->setTestVel(testVel);
        solid->setMassDensity(density);
        if(yield > 0) {
            solid->setPlasticYield(yield);
        }
        solid->setInitialState();
        solid->setInitialVelocity();
        solids.push_back(solid);
    }

    printf("Grid:\n");
    printf("Dimentions: %dx%d\n", res[0], res[1]);
    printf("Grid Spacing: %f\n", dh);
    printf("X0: (%f, %f)\n", origin(0), origin(1));
    printf("X1: (%f, %f)\n", origin(0)+dh*(res[0]-1), origin(1)+dh*(res[1]-1));
    printf("Center: (%f, %f)\n", center(0), center(1));

    printf("\nConstants:\n");
    printf("Total Time: %f\n", totalTime);
    printf("dt: %f\n", dt);
    printf("Steps per Frame: %d\n", (int)std::round(1.0/(30*dt)));
    printf("Gravity: (%f, %f)\n", gravity(0), gravity(1));
    printf("Rotation: %f\n\n", rotation);

    for(int i = 0; i < (int)solids.size(); i++) {
        printf("Solid %d\n", i);
        printf("\tCenter: (%f, %f)\n", solids[i]->center(0), solids[i]->center(1));
        printf("\tLame Constants: %f, %f\n", solids[i]->solidMaterial->lambda, solids[i]->solidMaterial->mu);
        if(solids[i]->isPlastic) {
            printf("\tPlastic Yield: %f\n", solids[i]->plasticYield);
        }
        else {
            printf("\tNo Plasticity\n");
        }
    }
    printf("\n");

    fluidDensity = 1;
    fluidViscocity = 0;
    bouyancy = 0;
    intensity = std::vector<double>(res[0]*res[1], 0);

    totalActiveCells = (res[0]-2) * (res[1]-2);
    totalCells = res[0] * res[1];
    offsets[0] = 0;
    offsets[1] = 1;
    offsets[2] = res[0];
    offsets[3] = res[0]+1;
    //Sizing vectors
    solidMass = std::vector<double>(totalCells);
    workspace = std::vector<double>(totalCells);
    nodeSolidVolumeFraction = std::vector<double>(totalCells);
    quarterSolidVolumeFractions = VectorX(QUAD*totalCells);
    velocity = std::vector<VectorN>(totalCells);
    vecWorkspace = std::vector<VectorN>(totalCells);

    inc = steps;
    #ifndef NDEBUG
    inc = 1;
    // inc = 10;
    sMax = 11;
    // sMax = 61;
    count = 0;
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

MatrixN interpolate_value(const std::vector<MatrixN>& data, const VectorN& point, int& index, VectorN origin, int res[DIM], double dh) {
    #if DIM == 3
    int i, j, k;
    double fx, fy, fz;
    #else
    int i, j;
    double fx, fy;
    #endif
    get_barycentric(point(0) - origin(0), i, fx, 0, res[0] - 1, dh);
    get_barycentric(point(1) - origin(1), j, fy, 0, res[1] - 1, dh);
    #if DIM == 3
    get_barycentric(point(2) - origin(2), k, fz, 0, res[2] - 1, dh);
    #endif
    index = i*res[1]+j;

    return data[index];
}

void World::init() {
    for(int i = 0; i < (int)solids.size(); i++) {
        solids[i]->init();
    }
    gatherSolidMasses();
    for(int i = 0; i < (int)solids.size(); i++) {
        solids[i]->recoverVelocityField(velocity);
    }

    gridToSolution.resize(totalCells, -1);
    int idx = 0;
    for(int x = 1; x < res[0] - 1; x++) {
        for (int y = 1; y < res[1] - 1; y++) {
            int index = x*res[1]+y;
            gridToSolution[index] = idx++;
        }
    }
    computeDivergenceConstraintJacobian();
}

void World::computeDivergenceConstraintJacobian() {
    std::vector<Tripletd> jacobiTripletList;
    int cnt = 0;
    for(int x = 0; x < res[0]-1; x++) {
        for(int y = 0; y < res[1]-1; y++) {
            int index = x*res[1]+y;
            for(int i = 0; i < QUAD; i++) {
                int idx = gridToSolution[index+offsets[i]];
                int signX = i % 2 == 0 ? -1 : 1;
                int signY = i % 4 <= 1 ? -1 : 1;
                #if DIM == 3
                int signZ = i < 4 ? -1 : 1;
                #endif
                if (idx >= 0) {
                  jacobiTripletList.push_back(Tripletd(cnt, idx*DIM, signX));
                  jacobiTripletList.push_back(Tripletd(cnt, idx*DIM+1, signY));
                  #if DIM == 3
                  jacobiTripletList.push_back(Tripletd(cnt, idx*DIM+2, signZ));
                  #endif
                }
            }
            cnt++;
        }
    }
    divergenceConstraintJacobian.resize(cnt, totalActiveCells*DIM);
    divergenceConstraintJacobian.setFromTriplets(jacobiTripletList.begin(), jacobiTripletList.end());
    divergenceConstraintJacobian.makeCompressed();


    JJTDiagEntries.clear();
    for(int k = 0; k < divergenceConstraintJacobian.rows(); ++k) {
        std::vector<std::pair<int, double> > entries;
        for(SparseMat::InnerIterator it(divergenceConstraintJacobian, k); it; ++it) {
            entries.push_back(std::make_pair(it.col(), it.value()*it.value()));
        }
        JJTDiagEntries.push_back(entries);
    }
    JT = divergenceConstraintJacobian.transpose();
    JT.makeCompressed();
}

/******************************
 * Algorithm Process from Eulerian Solid-Fluid Coupling
 ******************************/
void World::step() {
    if(stepNum == 0) {
        #ifndef NDEBUG
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                matTrans[count][i*res[1]+j] = solids[0]->X[i*res[1]+j];
            }
        }
        count++;
        #endif
        #ifdef CLUSTER
        writeParticles("/nfs/scratch/adahl1/euler/part-00.00000.bgeo", solids[0]->vertices);
        #else
        writeParticles("./euler/part-00.00000.bgeo", solids[0]->vertices);
        #endif
        #ifdef CLUSTER
        std::ofstream vout("/nfs/scratch/adahl1/euler/vel0");
        #else
        std::ofstream vout("./euler/vel0");
        #endif
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                Vector2d p = origin+dh*Vector2d(i,j);
                int ind = i*res[1]+j;
                vout << p(0) << " " << p(1) << " " << velocity[ind](0) << " " << velocity[ind](1) << " 0 0 255\n";
            }
        }
        vout.close();
        printf("Total Cells: %d\n", totalCells);
        printf("Total Active Cells: %d\n", totalActiveCells);
        printf("Total Solid Cells: %d\n", totalSolidCells);
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
    dampedSolidMass.conservativeResize(totalSolidCells * DIM);
    dampedSolidMass.setZero();
    std::vector<Tripletd> tripletList;
    //Compute F and Stiffness Matrix
    {
    auto timer = prof.timeName("Compute Kf and MV");
    for(int i = 0; i < (int)solids.size(); i++) {
        solids[i]->computeKf(solidUnitForce, tripletList, gridToSolid);
        solids[i]->computeMV(solidMv, gridToSolid);
    }
    }
    //Skipping vorticity stuff since it's always 0
    //Mass part of rayleigh damping
    {
    auto timer = prof.timeName("Body Forces");
    double massDamp = 1 + dt * rayleighAlpha;
    //Compute body Forces
    #pragma omp parallel for schedule(static) collapse(2)
    for(int i = 1; i < res[0]-1; i++) {
        for(int j = 1; j < res[1]-1; j++) {
            int index = i*res[1]+j;
            int idx = gridToSolution[index];
            int solidIdx = gridToSolid[index];

            double fluidMass = std::max(1.0 - nodeSolidVolumeFraction[index], 0.0) * fluidDensity * dh * dh;
            double mass = fluidMass + solidMass[index];

            if(solidIdx < 0) {
                unitForce.segment<DIM>(idx*DIM) += mass * gravity;
                // buoyancy
                unitForce[idx*DIM+1] += bouyancy * fluidMass * intensity[index];
                nzMassVec[idx*DIM] = nzMassVec[idx*DIM+1] = mass;
                #if DIM == 3
                nzMassVec[idx*DIM+2] = mass;
                #endif
                mv.segment<DIM>(idx*DIM) = fluidMass * velocity[index];
            }
            else {
                dampedSolidMass[solidIdx * DIM] = dampedSolidMass[solidIdx * DIM + 1] = solidMass[index] * massDamp + fluidMass;
                #if DIM == 3
                dampedSolidMass[solidIdx * DIM + 2] = mass * massDamp + fluidMass;
                #endif
                solidMv.segment<DIM>(solidIdx*DIM) += fluidMass * velocity[index];
                if(rotationEnabled) {
                    VectorN d = (origin+dh*Vector2d(i,j)) - center;
                    solidUnitForce.segment<DIM>(solidIdx * DIM) += mass * rotation * VectorN(-d(1), d(0));
                }
                if(gravityEnabled) {
                    solidUnitForce.segment<DIM>(solidIdx * DIM) += mass * gravity;
                }
                // buoyancy
                solidUnitForce[solidIdx*DIM+1] += bouyancy * fluidMass * intensity[index];
            }
            nzMassVecInv[idx*DIM] = nzMassVecInv[idx*DIM+1] = 1.0 / mass;
            #if DIM == 3
            nzMassVecInv[idx*DIM+2] = 1.0 / mass;
            #endif
        }
    }
    //Construct global system matrix
    systemMatrix.resize(totalSolidCells*DIM, totalSolidCells*DIM);
    systemMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    systemMatrix *= dt*dt + rayleighBeta*dt;
    for(int i = 0; i < dampedSolidMass.size(); i++) {
        systemMatrix.coeffRef(i, i) += dampedSolidMass[i];
    }
    systemMatrix.makeCompressed();
    }
    printf("Triplet: %d\n", (int)tripletList.size());
    printf("Damped Solid: %d\n", (int)dampedSolidMass.size());
    printf("System Matrix: (%d, %d)\n", (int)systemMatrix.rows(), (int)systemMatrix.cols());

    {
    auto timer = prof.timeName("Divergence RHS");
    computeDivergenceConstraintRHS();
    }
    //Do PCR (conjugate residual) solve
    {
        std::ofstream vb("./euler/v0");
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int index = i*res[1]+j;
                VectorN xg = dh*VectorN(i,j) + origin;
                VectorN v = velocity[index];
                vb << xg(0) << " " << xg(1) << " " << v(0) << " " << v(1) << "\n";
            }
        }
        vb.close();
    auto timer = prof.timeName("PCR");
    pcrImplicit();
    }
    std::ofstream fout("./euler/frc"+std::to_string(stepNum/inc));
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            Vector2d p = origin+dh*Vector2d(i,j);
            int ind = i*res[1]+j;
            Vector3i col(255, 0, 0);
            if(gridToSolid[ind] >= 0) {
                int idx = gridToSolid[ind];
                fout << p(0) << " " << p(1) << " " << solidUnitForce.segment<DIM>(idx*DIM)(0) << " " << solidUnitForce.segment<DIM>(idx*DIM)(1) << " " << col(0) << " " << col(1) << " " << col(2) << "\n";
            }
        }
    }
    fout.close();
    //finalize solution
    {
    auto timer = prof.timeName("Finalize");
    finalizeSolution();
    }
    stepNum++;
	elapsedTime += dt;

    if(stepNum % inc == 0) {
        std::ostringstream ss;
        ss << std::setw(2) << std::setfill('0') << 0 << "." << std::setw(5) << stepNum/inc;
        #ifdef CLUSTER
        std::string outLoc = "/nfs/scratch/adahl1/euler/part-"+ss.str()+".bgeo";
        #else
        std::string outLoc = "./euler/part-"+ss.str()+".bgeo";
        #endif
        writeParticles(outLoc.c_str(), solids[0]->vertices);
        #ifdef CLUSTER
        std::ofstream vout("/nfs/scratch/adahl1/euler/vel"+std::to_string(stepNum/inc));
        #else
        std::ofstream vout("./euler/vel"+std::to_string(stepNum/inc));
        #endif
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                Vector2d p = origin+dh*Vector2d(i,j);
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
                matTrans[count][i*res[1]+j] = solids[0]->X[i*res[1]+j];
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

void World::computeDivergenceConstraintRHS() {
    divergenceConstraintsRHS.conservativeResize(divergenceConstraintJacobian.rows());

    // #pragma omp parallel for schedule(static)
    int cnt = 0;
    for(int x = 0; x < res[0]-1; x++) {
        for(int y = 0; y < res[1]-1; y++) {
            int index = x*res[1]+y;
            divergenceConstraintsRHS[cnt] = 0;
            for(int i = 0; i < QUAD; i++) {
                int idx = gridToSolution[index+offsets[i]];
                int signX = i % 2 == 0 ? -1 : 1;
                int signY = i % 4 <= 1 ? -1 : 1;
                #if DIM == 3
                int signZ = i < 4 ? -1 : 1;
                #endif
                if(idx < 0) {
                    divergenceConstraintsRHS[cnt] -= signX * velocity[index+offsets[i]](0);
                    divergenceConstraintsRHS[cnt] -= signY * velocity[index+offsets[i]](1);
                    #if DIM == 3
                    divergenceConstraintsRHS[cnt] -= signZ * velocity[index+offsets[i]](2);
                    #endif
                }
            }
            cnt++;
        }
    }
}

//Solve eq. 8 for the velocity from the forces
void World::pcrImplicit() {
    systemMatDiag = systemMatrix.diagonal();
    solidUnitForce *= dt;
    solidUnitForce += solidMv;
    unitForce *= dt;
    unitForce += mv;

    std::ofstream sf("./euler/sf1");
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            if(gridToSolid[index] < 0) {
                continue;
            }
            VectorN xg = dh*VectorN(i,j) + origin;
            VectorN v = solidUnitForce.segment<DIM>(gridToSolid[index]*DIM);
            sf << xg(0) << " " << xg(1) << " " << v(0) << " " << v(1) << "\n";
        }
    }
    sf.close();

    // Eqn (8) in the paper
    solvePCR(systemMatrix, solidUnitForce, solidVelocityAtDt0, systemMatDiag);
    std::ofstream vpcr0("./euler/vpcr0");
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            if(gridToSolid[index] < 0) {
                continue;
            }
            VectorN xg = dh*VectorN(i,j) + origin;
            VectorN v = solidVelocityAtDt0.segment<DIM>(gridToSolid[index]*DIM);
            vpcr0 << xg(0) << " " << xg(1) << " " << v(0) << " " << v(1) << "\n";
        }
    }
    vpcr0.close();
    velocityAtDt0.conservativeResize(totalActiveCells * DIM);
    #pragma omp parallel for schedule(static)
    for(int index = 0; index < totalCells; index++) {
        int idx = gridToSolution[index];
        if (idx < 0) {
            continue;
        }
        int solidIdx = gridToSolid[index];
        if(solidIdx < 0) {
            velocityAtDt0.segment<DIM>(idx*DIM) = unitForce.segment<DIM>(idx*DIM) * nzMassVecInv[idx*DIM];
        }
        else {
            velocityAtDt0.segment<DIM>(idx*DIM) = solidVelocityAtDt0.segment<DIM>(solidIdx*DIM);
        }
    }
    std::ofstream vpcr1("./euler/vpcr1");
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            if(gridToSolution[index] < 0) {
                continue;
            }
            VectorN xg = dh*VectorN(i,j) + origin;
            VectorN v = velocityAtDt0.segment<DIM>(gridToSolution[index]*DIM);
            vpcr1 << xg(0) << " " << xg(1) << " " << v(0) << " " << v(1) << "\n";
        }
    }
    vpcr1.close();

    JMInvJTDiag.conservativeResize(divergenceConstraintJacobian.rows());
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < (int)JJTDiagEntries.size(); i++) {
        double sum = 0;
        for(std::pair<int, double>& entry : JJTDiagEntries[i]) {
            sum += nzMassVecInv[entry.first] * entry.second;
        }
        JMInvJTDiag[i] = sum;
    }
    VectorX b1 = divergenceConstraintJacobian*velocityAtDt0-divergenceConstraintsRHS;

    VectorX lambda1;
    //Div Free solver
    solvePCRDivFree(b1, lambda1, JMInvJTDiag);
    b1 = -divergenceConstraintJacobian.transpose() * lambda1;
    VectorX solidB1(totalSolidCells * DIM);
    #pragma omp parallel for schedule(static)
    for(int x = 0; x < totalCells; x++) {
        int idx = gridToSolution[x];
        int solidIdx = gridToSolid[x];
        if(solidIdx >= 0) {
            solidB1.segment<DIM>(solidIdx * DIM) = b1.segment<DIM>(idx * DIM);
        }
    }
    solvePCR(systemMatrix, solidB1, solidPressureCorrection, systemMatDiag);
    pressureCorrection.conservativeResize(totalActiveCells * DIM);
    #pragma omp parallel for schedule(static)
    for(int index = 0; index < totalCells; index++) {
        int idx = gridToSolution[index];
        if(idx < 0) {
            continue;
        }
        int solidIdx = gridToSolid[index];
        if(solidIdx < 0) {
            pressureCorrection.segment<DIM>(idx * DIM) = b1.segment<DIM>(idx * DIM) * nzMassVecInv[idx * DIM];
        }
        else {
            pressureCorrection.segment<DIM>(idx * DIM) = solidPressureCorrection.segment<DIM>(solidIdx * DIM);
        }
    }
    std::ofstream vpcr2("./euler/vpcr2");
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            if(gridToSolution[index] < 0) {
                continue;
            }
            VectorN xg = dh*VectorN(i,j) + origin;
            VectorN v = pressureCorrection.segment<DIM>(gridToSolution[index]*DIM);
            vpcr2 << xg(0) << " " << xg(1) << " " << v(0) << " " << v(1) << "\n";
        }
    }
    vpcr2.close();

    velocityAtDt1 = velocityAtDt0 + pressureCorrection;
    std::ofstream vpcr3("./euler/vpcr3");
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            if(gridToSolution[index] < 0) {
                continue;
            }
            VectorN xg = dh*VectorN(i,j) + origin;
            VectorN v = velocityAtDt1.segment<DIM>(gridToSolution[index]*DIM);
            vpcr3 << xg(0) << " " << xg(1) << " " << v(0) << " " << v(1) << "\n";
        }
    }
    vpcr3.close();
    // velocity clear;
    velocity = std::vector<VectorN>(totalCells, VectorN::Zero());
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < totalCells; i++) {
        int idx = gridToSolution[i];
        if(idx >= 0) {
            velocity[i] = velocityAtDt1.segment<DIM>(idx * DIM);
        }
    }
}

//Preconditioned Conjugate Residual Method
//Preconditioned with Jacobi
bool World::solvePCR(const SparseMat& A, const VectorX& b, VectorX& x, const VectorX& precond) {
    // double eps = 1e-8;
    double eps = 0.02;
    // int maxIteration = 10;
    int maxIteration = 500;

    x.conservativeResize(b.size());
    x.setZero();
    q.conservativeResize(b.size());
    q.setZero();
    residual = b;

    ArrayX Ma = precond.array();
    direction = (residual.array() / Ma).matrix();
    residual = direction;

    // s = A * residual;
    parallelMatVec(A, residual, s);
    Ap = s;

    double deltaNew = residual.dot(s);
    double bNorm = b.norm();

    bool converged = false;
    for(int i = 0; i < maxIteration && !converged; i++) {
        q = (Ap.array() / Ma).matrix();

        double alpha = Ap.dot(q);
        if(std::fabs(alpha) > 0.0) {
            alpha = deltaNew / alpha;
        }

        x += alpha * direction;
        residual -= alpha * q;
        // s = A * residual;
        parallelMatVec(A, residual, s);
        double deltaOld = deltaNew;

        deltaNew = residual.dot(s);
        // double beta = 0;
        // if(deltaOld > 1e-9) {
        //     beta = deltaNew / deltaOld;
        // }
        double beta = deltaNew / deltaOld;
        direction *= beta;
        direction += residual;
        Ap *= beta;
        Ap += s;

        if(i > 0 && i % 5 == 0) {
            // tmp = A * x;
            parallelMatVec(A, x, tmp);
            tmp -= b;
            double tmpNorm = tmp.norm();
            if(tmpNorm < eps * bNorm) {
                converged = true;
                printf("Number of iterations to converge: %d\n", i);
                break;
            }
        }
    }
    return converged;
}

void World::getMatVecMult(const VectorX& input, VectorX& output) {
    VectorX midResult(JT.rows());
    #pragma omp parallel for schedule(static)
    for(int k = 0; k < JT.rows(); ++k) {
        midResult[k] = 0;
        for(SparseMat::InnerIterator it(JT, k); it; ++it) {
            midResult[k] += it.value() * input[it.col()];
        }
    }
    output.conservativeResize(divergenceConstraintJacobian.rows());
    #pragma omp parallel for schedule(static)
    for(int k = 0; k < divergenceConstraintJacobian.rows(); ++k) {
        output[k] = 0;
        for(SparseMat::InnerIterator it(divergenceConstraintJacobian, k); it; ++it) {
            output[k] += it.value() * midResult[it.col()] * nzMassVecInv[it.col()];
        }
    }
}

bool World::solvePCRDivFree(const VectorX& b, VectorX& x, const VectorX& precond) {
    // double eps = 1e-8;
    double eps = 0.02;
    // int maxIteration = 10;
    int maxIteration = 500;

    if(x.size() != b.size()) {
        x.conservativeResize(b.size());
        x.setZero();
    }

    getMatVecMult(x, q);
    residual = b - q;

    ArrayX Ma = precond.array();
    direction = (residual.array() / Ma).matrix();
    residual = direction;

    getMatVecMult(residual, s);
    Ap = s;

    double deltaNew = residual.dot(s);

    double bNorm = b.norm();

    bool converged = false;
    for(int i = 0; i < maxIteration && !converged; i++) {
        q = (Ap.array() / Ma).matrix();

        double alpha = Ap.dot(q);
        if(std::fabs(alpha) > 0.0) {
            alpha = deltaNew / alpha;
        }
        x += alpha * direction;
        residual -= alpha * q;

        getMatVecMult(residual, s);
        double deltaOld = deltaNew;

        deltaNew = residual.dot(s);
        double beta = deltaNew / deltaOld;

        direction *= beta;
        direction += residual;
        Ap *= beta;
        Ap += s;

        if(i > 0 && i % 5 == 0) {
            getMatVecMult(x, tmp);
            tmp -= b;
            double tmpNorm = tmp.norm();
            if(tmpNorm < eps * bNorm) {
                converged = true;
                printf("Number of iterations to converge: %d\n", i);
                break;
            }
        }
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
    auto timer = prof.timeName("\tFinalize Solids");
    for(int i = 0; i < (int)solids.size(); i++)
        solids[i]->finalizeSolution(velocity, dt);
    }
    //grid to particles and particles to grid
    {
    auto timer = prof.timeName("\tParticles and Grid");
    for(int i = 0; i < (int)solids.size(); i++) {
        solids[i]->gridToParticles();
    }
    for(int i = 0; i < (int)solids.size(); i++) {
        solids[i]->particlesToGrid();
    }
    }
    //extend velocity field
    {
    auto timer = prof.timeName("\tExtend Velocity Field");
    for (int i = 0; i < (int)solids.size(); i++) {
        solids[i]->extendVectorField();
    }
    }
    {
    auto timer = prof.timeName("\tSolid Displacements");
    for(int i = 0; i < (int)solids.size(); i++) {
        solids[i]->finalizeSolution();
    }
    }
    {
    auto timer = prof.timeName("\tCompute Solid Mass");
    for(int i = 0; i < (int)solids.size(); i++) {
        solids[i]->computeSolidMass(true);
    }
    }
    {
    auto timer = prof.timeName("\tAdvect Surface Particles");
    for(int i = 0; i < (int)solids.size(); i++) {
        solids[i]->advectSurface();
    }
    }
    {
    auto timer = prof.timeName("\tSolid SDF and Particle Init");
    for(int i = 0; i < (int)solids.size(); i++) {
        solids[i]->computeSDF();
        solids[i]->initMaterialParticles();
    }
    }
    {
    auto timer = prof.timeName("\tGather Solid Masses");
    gatherSolidMasses();
    }
    {
    auto timer = prof.timeName("\tAdvecting Fields");
    advectFluidVelocityField(vecWorkspace, dt);
    advectField(intensity, workspace, dt);
    }

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < (int)solids.size(); i++) {
        solids[i]->recoverVelocityField(velocity);
    }

    // Set Dirichlet Boundary condition
    for(int x = 0; x < res[0]; x++) {
        velocity[x*res[1]].setZero();
        velocity[x*res[1]+res[1]-1].setZero();
    }
    for(int y = 0; y < res[1]; y++) {
        velocity[y].setZero();
        velocity[(res[0]-1)*res[1]+y].setZero();
    }
}

void World::gatherSolidMasses() {
    // solidMass clear
    solidMass = std::vector<double>(totalCells, 0);
    quarterSolidVolumeFractions.setZero();
    for(Object* ob : solids) {
        for(int i = 0; i < (int)solidMass.size(); i++) {
            solidMass[i] += ob->mass[i];
        }
        quarterSolidVolumeFractions += ob->quarterVolumeFractions;
    }
    // nodeSolidVolumeFraction clear
    nodeSolidVolumeFraction = std::vector<double>(totalCells, 0);
    for(int i = 0; i < QUAD; i++) {
        #pragma omp parallel for schedule(static) collapse(2)
        for(int x = 1; x < res[0]-1; x++) {
            for(int y = 1; y < res[1]-1; y++) {
                int index = x*res[1]+y;
                nodeSolidVolumeFraction[index+offsets[i]] += quarterSolidVolumeFractions[index*QUAD+i];
            }
        }
    }
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < (int)nodeSolidVolumeFraction.size(); i++) {
        nodeSolidVolumeFraction[i] *= 0.125;
    }

    gridToSolid = std::vector<int>(totalCells, -1);
    int idx = 0;
    for(int x = 1; x < res[0] - 1; x++) {
        for (int y = 1; y < res[1] - 1; y++) {
            int index = x*res[1]+y;
            if(solidMass[index] > 0) {
                gridToSolid[index] = idx++;
            }
            else {
                gridToSolid[index] = -1;
            }
        }
    }
    totalSolidCells = idx;
}

void World::advectFluidVelocityField(std::vector<VectorN>& workspace, double dt) {
    #pragma omp parallel for schedule(static)
    for(int x = 1; x < res[0]-1; x++) {
        for(int y = 1; y < res[1]-1; y++) {
            int index = x*res[1]+y;
            if(solidMass[index] > 0) {
                continue;
            }
            VectorN node(x*dh+origin(0), y*dh+origin(1));
            VectorN v = interpolate_value(node-0.5*velocity[index]*dt, velocity, origin, res, dh);
            VectorN newPos = node - v*dt;

            double mass = interpolate_value(newPos, solidMass, origin, res, dh);
            double newdt = dt/2.0;
            while(mass > 0 && newdt > dt/16.0) {
                newPos += v*newdt;
                newdt /= 2;
                mass = interpolate_value(newPos, solidMass, origin, res, dh);
            }
            workspace[index] = interpolate_value(newPos, velocity, origin, res, dh);
        }
    }
    velocity.swap(workspace);
}

void World::advectField(std::vector<double>& field, std::vector<double>& workspace, double dt) {
    // workspace clear
    workspace = std::vector<double>(totalCells, 0);
    #pragma omp parallel for schedule(static) collapse(2)
    for(int x = 0; x < res[0]; x++) {
        for(int y = 0; y < res[1]; y++) {
            int index = x*res[1]+y;
            if(solidMass[index] > 0) {
                continue;
            }
            VectorN node(x*dh+origin(0), y*dh+origin(1));
            VectorN v = velocity[index];
            VectorN newPos = node - v*dt;

            double mass = interpolate_value(newPos, solidMass, origin, res, dh);
            double newdt = dt / 2.0;
            while(mass > 0 && newdt > dt/16.0) {
                newPos += v*newdt;
                newdt /= 2;
                mass = interpolate_value(newPos, solidMass, origin, res, dh);
            }
            workspace[index] = interpolate_value(newPos, field, origin, res, dh);
        }
    }
    field.swap(workspace);
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

void World::apicAdvect() {
    //Grid to particles
    double alpha = 0.99;
    #pragma omp parallel for schedule(static)
    for(size_t c = 0; c < solids[0]->matParticles.size(); c++) {
        Particle &p = solids[0]->matParticles[c];
        // p.vo = p.v;
        //Update velocities
        p.B = MatrixN::Zero();
        // p.u = VectorN::Zero();
        VectorN apic = VectorN::Zero();
        VectorN flip = VectorN::Zero();
        double x = (p.x(0) - origin(0)) / dh;
        double y = (p.x(1) - origin(1)) / dh;
        int i = std::floor(x) - 1;
        int j = std::floor(y) - 1;
        for(int r = i; r < i+4; r++) {
            for(int t = j; t < j+4; t++) {
                if(r < 0 || t < 0 || r > res[0]-1 || t > res[1]-1) {
                    continue;
                }
                int ind = r*res[1]+t;
                if(solids[0]->mass[ind] == 0) {
                    continue;
                }
                VectorN diff = dhInv * (p.x - origin);
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

                VectorN xg = origin + dh*VectorN(r, t);
                VectorN wu = w * solids[0]->u[ind];
                apic += wu;
                flip += (solids[0]->u[ind] - solids[0]->matWorkspace[ind]) * w;
                p.B += wu * (xg - p.x).transpose();

                // VectorN wvel = w * velocity[ind];
                // apic += wvel;
                // flip += (velocity[ind] - vecWorkspace[ind]) * w;
                // p.B += wvel * (xg - p.x).transpose();
                // p.u += u[ind];
            }
        }
        #ifndef NDEBUG
        if(apic.hasNaN()) {
            printf("\n\nAPIC Vel has NaN\n");
            std::cout << apic << std::endl;
            exit(0);
        }
        #endif
        p.u = (1 - alpha) * apic + alpha * flip;
        // flip += p.v;
        // p.v = (1 - alpha) * apic + alpha * flip;
        //Mass proportional damping
        // p.v = mp.massPropDamp * p.v;
    }

    #pragma omp parallel for schedule(dynamic)
    for(int c = 0; c < (int)solids[0]->matParticles.size(); c++) {
        Particle& p = solids[0]->matParticles[c];
        #ifndef NDEBUG
        if(p.v.hasNaN()) {
            printf("Vel has NaN\n");
            std::cout << p.v << std::endl;
            exit(0);
        }
        #endif
        //Update Positions
        // p.xo = p.x;
        p.x += dt * p.v;

        #ifndef NDEBUG
        if(p.x.hasNaN()) {
            printf("Pos has NaN\n");
            std::cout << p.x << std::endl;
            exit(0);
        }
        #endif
    }

    //Particles to Grid
    MatrixN tmpDinv = (3.0/(dh*dh))*MatrixN::Identity();
    std::unordered_map<int, std::tuple<VectorN, double, double> > indexToWeight;
    #pragma omp parallel
    {
    std::unordered_map<int, std::tuple<VectorN, double, double> > localMap;
    #pragma omp parallel for schedule(static)
    for(size_t c = 0; c < solids[0]->matParticles.size(); c++) {
        Particle &p = solids[0]->matParticles[c];
        double x = (p.x(0) - origin(0)) / dh;
        double y = (p.x(1) - origin(1)) / dh;
        int i = std::floor(x) - 1;
        int j = std::floor(y) - 1;
        VectorN mv = p.m*p.u;
        MatrixN mBD = p.m*p.B*tmpDinv;
        for(int s = i; s < i+4; s++) {
            for(int t = j; t < j+4; t++) {
                if(s < 0 || t < 0 || s > res[0]-1 || t > res[1]-1) {
                    continue;
                }
                VectorN diff = dhInv * (p.x - origin);
                diff -= VectorN(s, t);
                VectorN n(0, 0);
                for(int k = 0; k < DIM; k++) {
                    if(std::abs(diff(k)) < 1) {
                        n(k) = 0.5*std::abs(std::pow(diff(k), 3)) - diff(k)*diff(k) + 2.0/3.0;
                    }
                    else if(std::abs(diff(k)) < 2) {
                        n(k) = -1.0/6.0*std::abs(std::pow(diff(k), 3)) + diff(k)*diff(k) - 2*std::abs(diff(k)) + 4.0/3.0;
                    }
                }
                double w = n(0) * n(1);
                int index = s*res[1]+t;
                VectorN xg = origin + dh*VectorN(s, t);
                VectorN update = w * (mv + mBD*(xg-p.x));
                if(localMap.find(index) == localMap.end()) {
                    std::get<0>(localMap[index]) = update;
                    std::get<1>(localMap[index]) = w*p.m;
                    std::get<2>(localMap[index]) = w;
                }
                else {
                    std::get<0>(localMap[index]) += update;
                    std::get<1>(localMap[index]) += w*p.m;
                    std::get<2>(localMap[index]) += w;
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
                std::get<0>(indexToWeight[iToW.first]) += std::get<0>(iToW.second);
                std::get<1>(indexToWeight[iToW.first]) += std::get<1>(iToW.second);
                std::get<2>(indexToWeight[iToW.first]) += std::get<2>(iToW.second);
            }
        }
    }
    }
    // solids[0]->mass.clear();
    solids[0]->mass = std::vector<double>(totalCells, 0);
    int cutoff = 3;
    for(auto iToW : indexToWeight) {
        if(std::get<2>(iToW.second) < cutoff) {
            continue;
        }
        solids[0]->u[iToW.first] = std::get<0>(iToW.second) / std::get<1>(iToW.second);
        solids[0]->mass[iToW.first] = 1;
    }
}

void World::apicg2p() {
    // std::ofstream advectout("./euler/advect");
    // double alpha = 0.99;
    #pragma omp parallel for schedule(static)
    for(size_t c = 0; c < solids[0]->matParticles.size(); c++) {
        Particle &p = solids[0]->matParticles[c];
        p.vo = p.v;
        //Update velocities
        p.B = MatrixN::Zero();
        VectorN apic = VectorN::Zero();
        double x = (p.x(0) - origin(0)) / dh;
        double y = (p.x(1) - origin(1)) / dh;
        int i = std::floor(x) - 1;
        int j = std::floor(y) - 1;
        for(int r = i; r < i+4; r++) {
            for(int t = j; t < j+4; t++) {
                if(r < 0 || t < 0 || r > res[0]-1 || t > res[1]-1) {
                    continue;
                }
                int ind = r*res[1]+t;
                if(solids[0]->mass[ind] == 0) {
                    continue;
                }
                VectorN diff = dhInv * (p.x - origin);
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

                VectorN xg = origin + dh*VectorN(r, t);
                VectorN wvel = w * solids[0]->velocity[ind];
                // VectorN cleft = w * (xg - p.x);
                // VectorN cright = p.m * solids[0]->velocity[ind];
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
        //Mass proportional damping
        // p.v = mp.massPropDamp * p.v;
    }

    VectorN wallNormals[4] = {VectorN(1, 0), VectorN(0, 1),
                              VectorN(-1, 0), VectorN(0, -1)};
    double boundaries[4] = {origin(0),
                            origin(1),
                            origin(0) + dh * (res[0] - 1),
                            origin(1) + dh * (res[1] - 1)};
    double wallThickness = 0.015;
    double springConstant = 10000;
    double friction = 0.5;
    #pragma omp parallel for schedule(dynamic)
    for(int c = 0; c < (int)solids[0]->matParticles.size(); c++) {
        Particle& p = solids[0]->matParticles[c];
        VectorN futurePos = p.x;
        double fx = (futurePos(0) - origin(0)) / dh;
        double fy = (futurePos(1) - origin(1)) / dh;
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
                VectorN v_other = solids[0]->velocity[xIndex*res[1]+yIndex];
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

        #ifndef NDEBUG
        if(p.v.hasNaN()) {
            printf("Vel has NaN\n");
            std::cout << p.v << std::endl;
            exit(0);
        }
        #endif

        //Update Positions
        p.xo = p.x;
        p.x += dt * p.v;

        #ifndef NDEBUG
        if(p.x.hasNaN()) {
            printf("Pos has NaN\n");
            std::cout << p.x << std::endl;
            exit(0);
        }
        #endif
    }
}

void World::apicp2g() {
    MatrixN tmpDinv = (3.0/(dh*dh))*MatrixN::Identity();
    std::unordered_map<int, std::tuple<VectorN, double, double> > indexToWeight;
    #pragma omp parallel
    {
    std::unordered_map<int, std::tuple<VectorN, double, double> > localMap;
    #pragma omp parallel for schedule(static)
    for(size_t c = 0; c < solids[0]->matParticles.size(); c++) {
        Particle &p = solids[0]->matParticles[c];
        double x = (p.x(0) - origin(0)) / dh;
        double y = (p.x(1) - origin(1)) / dh;
        int i = std::floor(x) - 1;
        int j = std::floor(y) - 1;
        MatrixN Dinv = MatrixN::Zero();
        for(int s = i; s < i+4; s++) {
            for(int t = j; t < j+4; t++) {
                if(s < 0 || t < 0 || s > res[0]-1 || t > res[1]-1) {
                    continue;
                }
                VectorN diff = dhInv * (p.x - origin);
                diff -= VectorN(s, t);
                VectorN n(0, 0);
                for(int k = 0; k < DIM; k++) {
                    if(std::abs(diff(k)) < 1) {
                        n(k) = 0.5*std::abs(std::pow(diff(k), 3)) - diff(k)*diff(k) + 2.0/3.0;
                    }
                    else if(std::abs(diff(k)) < 2) {
                        n(k) = -1.0/6.0*std::abs(std::pow(diff(k), 3)) + diff(k)*diff(k) - 2*std::abs(diff(k)) + 4.0/3.0;
                    }
                }
                double w = n(0) * n(1);
                VectorN xg = origin + dh*VectorN(s, t);
                VectorN xixp = xg - p.x;
                Dinv += w * xixp * xixp.transpose();
            }
        }
        VectorN mv = p.m*p.v;
        // if(Dinv.determinant() < 1e-12) {
        //     printf("Singular D: %f\n", Dinv.determinant());
        // }
        // MatrixN mBD = p.m*p.B*Dinv.inverse();
        MatrixN mBD = p.m*p.B*tmpDinv;
        for(int s = i; s < i+4; s++) {
            for(int t = j; t < j+4; t++) {
                if(s < 0 || t < 0 || s > res[0]-1 || t > res[1]-1) {
                    continue;
                }
                VectorN diff = dhInv * (p.x - origin);
                diff -= VectorN(s, t);
                VectorN n(0, 0);
                for(int k = 0; k < DIM; k++) {
                    if(std::abs(diff(k)) < 1) {
                        n(k) = 0.5*std::abs(std::pow(diff(k), 3)) - diff(k)*diff(k) + 2.0/3.0;
                    }
                    else if(std::abs(diff(k)) < 2) {
                        n(k) = -1.0/6.0*std::abs(std::pow(diff(k), 3)) + diff(k)*diff(k) - 2*std::abs(diff(k)) + 4.0/3.0;
                    }
                }
                double w = n(0) * n(1);
                int index = s*res[1]+t;
                VectorN xg = origin + dh*VectorN(s, t);
                VectorN update = w * (mv + mBD*(xg-p.x));
                if(localMap.find(index) == localMap.end()) {
                    std::get<0>(localMap[index]) = update;
                    std::get<1>(localMap[index]) = w*p.m;
                    std::get<2>(localMap[index]) = w;
                }
                else {
                    std::get<0>(localMap[index]) += update;
                    std::get<1>(localMap[index]) += w*p.m;
                    std::get<2>(localMap[index]) += w;
                }
            }
        }
        p.x = p.xo;
    }
    #pragma omp barrier
    #pragma omp critical
    {
        for(auto iToW : localMap) {
            if(indexToWeight.find(iToW.first) == indexToWeight.end()) {
                indexToWeight[iToW.first] = iToW.second;
            }
            else {
                std::get<0>(indexToWeight[iToW.first]) += std::get<0>(iToW.second);
                std::get<1>(indexToWeight[iToW.first]) += std::get<1>(iToW.second);
                std::get<2>(indexToWeight[iToW.first]) += std::get<2>(iToW.second);
            }
        }
    }
    }
    // solids[0]->mass.clear();
    solids[0]->mass = std::vector<double>(totalCells, 0);
    int cutoff = 3;
    for(auto iToW : indexToWeight) {
        if(std::get<2>(iToW.second) < cutoff) {
            continue;
        }
        solids[0]->velocity[iToW.first] = std::get<0>(iToW.second) / std::get<1>(iToW.second);
        solids[0]->mass[iToW.first] = 1;
    }
}

/************************
 * Helpers
 ************************/

void writeParticles(const char *fname, const std::vector<Particle> &particles) {
	Partio::ParticlesDataMutable *data = Partio::create();
	Partio::ParticleAttribute xattr;
	Partio::ParticleAttribute vattr;
    Partio::ParticleAttribute uattr;
	Partio::ParticleAttribute battr;
	Partio::ParticleAttribute cattr;
	Partio::ParticleAttribute mattr;
	Partio::ParticleAttribute rattr;
	Partio::ParticleAttribute voattr;
    data->addParticles(particles.size());

	data->addAttribute("position", Partio::VECTOR, 3);
	data->addAttribute("velocity", Partio::VECTOR, 2);
    data->addAttribute("X", Partio::VECTOR, 2);
	data->addAttribute("B", Partio::VECTOR, 4);
	data->addAttribute("Cd", Partio::VECTOR, 3);
	data->addAttribute("mass", Partio::FLOAT, 1);
	data->addAttribute("rho", Partio::FLOAT, 1);
	data->addAttribute("vol", Partio::FLOAT, 1);

	data->attributeInfo("position", xattr);
	data->attributeInfo("velocity", vattr);
    data->attributeInfo("X", uattr);
	data->attributeInfo("B", battr);
	data->attributeInfo("Cd", cattr);
	data->attributeInfo("mass", mattr);
	data->attributeInfo("rho", rattr);
	data->attributeInfo("vol", voattr);

	for (unsigned int i=0; i < particles.size(); i++) {
		const Particle &p = particles[i];
		float *x = data->dataWrite<float>(xattr, i);
        float *v = data->dataWrite<float>(vattr, i);
		float *u = data->dataWrite<float>(uattr, i);
		float *b = data->dataWrite<float>(battr, i);
		float *c = data->dataWrite<float>(cattr, i);
		float *m = data->dataWrite<float>(mattr, i);
		float *r = data->dataWrite<float>(rattr, i);
		float *vo = data->dataWrite<float>(voattr, i);

		x[0] = p.x(0), x[1] = p.x(1), x[2] = 0;
        v[0] = p.v(0), v[1] = p.v(1);
        u[0] = p.u(0), u[1] = p.u(1);
		b[0] = p.B(0,0), b[1] = p.B(0,1), b[2] = p.B(1,0), b[3] = p.B(1,1);
		c[0] = p.color(0), c[1] = p.color(1), c[2] = p.color(2);
		m[0] = p.m;
		r[0] = p.rho;
		vo[0] = p.vol;
	}

	Partio::write(fname, *data);
	data->release();
}


bool readParticles(const char *fname, std::vector<Particle> &particles) {
	Partio::ParticlesDataMutable *data = Partio::read(fname);
	if (data == 0) return 0;
	Partio::ParticleAttribute xattr;
    Partio::ParticleAttribute vattr;
	Partio::ParticleAttribute uattr;
	Partio::ParticleAttribute battr;
	Partio::ParticleAttribute cattr;
	Partio::ParticleAttribute mattr;
	Partio::ParticleAttribute rattr;
	Partio::ParticleAttribute voattr;

	bool position = data->attributeInfo("position", xattr);
	bool velocity = data->attributeInfo("velocity", vattr);
    bool X = data->attributeInfo("X", uattr);
	bool B = data->attributeInfo("B", battr);
	bool color = data->attributeInfo("Cd", cattr);
	bool mass = data->attributeInfo("mass", mattr);
	bool rho = data->attributeInfo("rho", rattr);
	bool vol = data->attributeInfo("vol", voattr);

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
            float *v = data->dataWrite<float>(vattr, i);
            p.v(0) = v[0], p.v(1) = v[1];
        } else {
            p.v = Vector2d(0.0, 0.0);
		}
        if (X) {
            float *u = data->dataWrite<float>(uattr, i);
            p.u(0) = u[0], p.u(1) = u[1];
        } else {
            p.u = p.x;
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
            float *vo = data->dataWrite<float>(voattr, i);
            p.vol = vo[0];
        } else {
            p.vol = 1.0;
        }
	}
	return true;
}
