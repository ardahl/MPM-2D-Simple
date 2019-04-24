#include "eulerworld.hpp"
#include <cstdio>
#include <iostream>
#include <limits>
#include <cmath>
#include <random>
#include <unordered_map>
#include <utility>

Object::Object(int resolution[DIM], double dh, VectorN origin, const std::vector<double>& restMass) :
    solidMaterial(NULL), dh(dh), origin(origin), restMass(restMass), quarterCellSolidMass(1), isPlastic(false)
{
    res[0] = resolution[0];
    res[1] = resolution[1];
    verticesInit = true;

    totalCells = res[0]*res[1];
    dhInv = 1.0 / dh;

    /***
    ***/
    subsampleDh = 0.00125;
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
    G *= 0.25 * dhInv;
    xmin = 0;
    ymin = 0;
    xmax = res[0];
    ymax = res[1];
    tightMin = VectorI(0, 0);
    tightMax = VectorI(res[0], res[1]);

    X = std::vector<VectorN>(totalCells);
    u = std::vector<VectorN>(totalCells);
    mass = std::vector<double>(totalCells);
    velocity = std::vector<VectorN>(totalCells);
    valid = std::vector<char>((res[0]-1)*(res[1]-1), 0);
    phi = std::vector<double>(totalCells);
    vecWorkspace = std::vector<VectorN>(totalCells);
    workspace = std::vector<double>(totalCells);
    detFInv = std::vector<double>(totalCells);

    quarterVolumeFractions = VectorX::Zero(QUAD*totalCells);

    scale = 1.0;
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
F Object::interpolate_value(const VectorN& point, std::vector<F> field, VectorN origin, int res[DIM], double dh) {
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

MatrixN Object::interpolate_value(const std::vector<MatrixN>& data, const VectorN& point, int& index, VectorN origin, int res[DIM], double dh) {
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

void Object::setInitialState() {
    int index = 0;
    for(int x = 0; x < res[0]; x++) {
        for (int y = 0; y < res[1]; y++, index++) {
            VectorN pos(x*dh, y*dh);
            pos += origin;
            u[index] = VectorN::Zero();
            X[index] = pos;
        }
    }
}

void Object::setInitialVelocity() {
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            VectorN pos = origin + dh*VectorN(i, j);
            VectorN ph = pos - center;
            if(testVel == 0) {
                if(i == 0 && j == 0) {
                    printf("Vel: Rotation\n");
                }
                velocity[index] = scale*VectorN(-ph(1), ph(0));
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
            if(i == 0 || j == 0 || i == res[0]-1 || j == res[1]-1) {
                velocity[index] = VectorN::Zero();
            }
        }
    }
}

void Object::setMaterial(double lamda, double mu) {
    solidMaterial = new Corotation(lamda, mu);
    #if defined(_OPENMP)
    int totalCores = omp_get_max_threads();
    #else
    int totalCores = 1;
    #endif
    materialCopies.resize(totalCores);
    for(int i = 0; i < totalCores; i++) {
        materialCopies[i] = new Corotation(lambda, mu);
    }
}

void Object::init() {
    computeSolidMass();
    computeSDF();
    initMaterialParticles();
    advectSurface();
}

void Object::computeSolidMass(bool usesdf) {
    mass = std::vector<double>(totalCells, 0);
    valid = std::vector<char>(totalCells, 0);
    detFInv = std::vector<double>(totalCells, 0);
    quarterVolumeFractions.setZero();

    double subsampleVolume = subsampleDh * subsampleDh;
    int massSampleRes = dh / subsampleDh;
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
            VectorN node(x*dh+origin(0), y*dh+origin(1));
            if(usesdf && phi[index] > 2.0) {
                continue;
            }
            for(int p = 0; p < massSampleRes; p++) {
                for(int q = 0; q < massSampleRes; q++) {
                    VectorN pos = node + VectorN((p+0.5)*subsampleDh, (q+0.5)*subsampleDh);
                    VectorN rPos = interpolate_value(pos, X, origin, res, dh);
                    double d = interpolate_value(rPos, restMass, origin, res, dh) * subsampleVolume;
                    if(d < 1e-6) {
                        continue;
                    }
                    // printf("d: %f, (%d, %d)\n", d, x, y);
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
                // printf("Sum: %f, (%d, %d)\n", sum, x, y);
                for(int i = 0; i < QUAD; i++) {
                    quarterVolumeFractions[index*QUAD+i] = 1;
                }
                valid[index] = 1;
            }
            else {
                if(sum > 0) {
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
    double sum = 0;
    for(int i = 0; i < (int)mass.size(); i++) {
        sum += mass[i];
    }
    printf("Mass Sum: %f\n", sum);
}

void Object::computeSDF() {
    double large_distance = res[0] + res[1] + 2;

    xmax = ymax = 0;
    xmin = res[0];
    ymin = res[1];

    phi = std::vector<double>(res[0]*res[1], large_distance);

    int index = 0;
    for(int x = 0; x < res[0]; x++) {
        for(int y = 0; y < res[1]; y++, index++) {
            if(mass[index] > 0) {
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
    distSweep(phi, valid, 8, 1e-8);
}

void Object::initMaterialParticles() {
    matParticles.clear();
    //Init material particles
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, std::nextafter(dh, std::numeric_limits<double>::max()));
    int particlesPerCell = 32;
    double density = quarterCellSolidMass / (dh * dh * 0.25);
    double volume = (dh * dh) / particlesPerCell;
    if(verticesInit) {
        printf("Particle Mass: %f\n", density*volume);
    }

    MatrixN tmpD = ((dh*dh)/3.0)*MatrixN::Identity();
    for(int i = tightMin[0]; i < tightMax[0]; i++) {
        for(int j = tightMin[1]; j < tightMax[1]; j++) {
            int ind = i*res[1] + j;
            if(!valid[ind]) {
                continue;
            }
            VectorN grid = origin + dh*VectorN(i, j);
            for(int k = 0; k < particlesPerCell; k++) {
                VectorN jitter = grid + VectorN(dis(gen), dis(gen));
                // double d = interpolate_value(jitter, restMass, origin, res, h);
                double d = density * volume;
                VectorN v = interpolate_value(jitter, velocity, origin, res, dh);
                Particle p(jitter, v, Eigen::Vector3d(0, 255, 0), d);
            	MatrixN C;
            	if(testVel == 0) {
		             C << 0, -scale, scale, 0; //Rotational (rotation=0.75) Case
            	}
		        else { //Linear and no velocity cases both have 0 gradient
		             C << 0, 0, 0, 0;
		        }
		        p.B = C * tmpD;
                p.u = VectorN::Zero();
		        matParticles.push_back(p);
            }
            //Vertex Particles for rendering, don't need as many
            if(verticesInit) {
                for(int k = 0; k < 4; k++) {
                    VectorN jitter = grid + VectorN(dis(gen), dis(gen));
                    VectorN ph = jitter - center;
                    Eigen::Vector3d col = Eigen::Vector3d::Zero();
                    if( ((ph(0)*ph(0))/(size[0]*size[0])) + ((ph(1)*ph(1))/(size[1]*size[1])) < 1+EPS) {
                        col = Eigen::Vector3d(255, 0, 0);
                        if(ph(0) < 0 && ph(1) < 0) {
                            col = Eigen::Vector3d(0, 255, 0);
                        }
                        if(ph(0) >= 0 && ph(1) < 0) {
                            col = Eigen::Vector3d(0, 0, 255);
                        }
                        if(ph(0) < 0 && ph(1) >= 0) {
                            col = Eigen::Vector3d(255, 255, 0);
                        }
                    }
                    Particle p(jitter, VectorN::Zero(), col, 0);
                    p.u = jitter;
                    p.xo = jitter;
                    p.vo = VectorN::Zero();
                    vertices.push_back(p);
                }
            }
        }
    }
    if(verticesInit) {
        printf("Number of Particles: %d\n", (int)matParticles.size());
    }
    verticesInit = false;
}

void Object::recoverVelocityField(std::vector<VectorN>& v) {
    for(int x = 0; x < res[0]; x++) {
        for(int y = 0; y < res[1]; y++) {
            int index = x*res[1]+y;
            if(mass[index] > 0) {
                v[index] = velocity[index];
            }
        }
    }
}

void Object::finalizeSolution(std::vector<VectorN>& vel, double step) {
    velocity = std::vector<VectorN>(totalCells, VectorN::Zero());
    #pragma omp parallel for schedule(static)
    for(int x = xmin; x < xmax; x++) {
        for(int y = ymin; y < ymax; y++) {
            int index = x*res[1]+y;
            velocity[index] = vel[index];
        }
    }
    dt = step;
}

void Object::computeKf(VectorX& force, std::vector<Tripletd>& tripletList, const std::vector<int>& indexMap) {
    double flowRate = 0.5;
    std::vector<int> matOffsets(materialCopies.size()+1, 0);
    int startSize = tripletList.size();

    #pragma omp parallel
    {
    #if defined(_OPENMP)
    int id = omp_get_thread_num();
    #else
    int id = 0;
    #endif
    VectorX localForce(force.size());
    localForce.setZero();
    std::vector<Tripletd> localMat;
    std::unordered_map<int, MatrixN> updatedFps;
    #pragma omp parallel for schedule(static) collapse(2)
    for(int x = tightMin[0]; x < tightMax[0]-1; x++) {
        for(int y = tightMin[1]; y < tightMax[1]-1; y++) {
            int index = x*res[1]+y;
            if(valid[index]) {
                continue;
            }
            MatrixX nodalForcesFEM = MatrixX::Constant(DIM, QUAD, 0);
            MatrixX stiffness = MatrixX::Constant(DIM*QUAD, DIM*QUAD, 0);

            MatrixX coord(DIM, QUAD);
            for(int i = 0; i < QUAD; i++) {
                coord.col(i) = u[index + offsets[i]];
            }
            VectorN restPos(0, 0);
            MatrixX Xtmp(DIM, QUAD);
            for(int i = 0; i < QUAD; i++) {
                Xtmp.col(i) = X[index + offsets[i]];
                restPos += X[index + offsets[i]];
            }
            restPos /= 8;

            Corotation* material = materialCopies[id];
            MatrixN F;
            double Jinv = 1.0 / computeF(coord, G, F);
            material->init(F);
            MatrixN P = material->firstPiolaKirchhoff();
            MatrixN cauchyStress = Jinv * P * F.transpose();
            if(isPlastic) {
                if(cauchyStress.norm() > plasticYield) {
                    int materialGridIndex = 0;
                    MatrixN FpInv = interpolate_value(FpInvVec, restPos, materialGridIndex, origin, res, dh);
                    F *= FpInv;
                    material->init(F);
                    P = material->firstPiolaKirchhoff();
                    cauchyStress = Jinv * P * F.transpose();
                    double cNorm = cauchyStress.norm();
                    if(cNorm > plasticYield) {
                        double gamma = flowRate * (cNorm - plasticYield) / cNorm;
                        MatrixN& Fhat = material->Fhat;
                        MatrixN& V = material->V;

                        double detFhatCubicRoot = std::cbrt(1.0 / std::max(std::abs(Jinv), 1e-3));
                        Fhat(0, 0) = pow(Fhat(0, 0) / detFhatCubicRoot, -gamma);
                        Fhat(1, 1) = pow(Fhat(1, 1) / detFhatCubicRoot, -gamma);
                        #if DIM == 3
                        Fhat(2, 2) = pow(Fhat(2, 2) / detFhatCubicRoot, -gamma);
                        #endif

                        MatrixN deltaFpInv = V * Fhat * V.transpose();
                        FpInv = deltaFpInv * FpInv;
                        if((FpInv.array() != FpInv.array()).any()) {
                            std::cout << "Bad update!!! FpInv has nan!!!" << std::endl;
                        }
                        else {
                            updatedFps[materialGridIndex] = FpInv;
                        }
                    }
                }
            }

            MatrixX pFpxhat = MatrixX::Constant(DIM*DIM, DIM*QUAD, 0);
            computePFPxhat(F, G, pFpxhat);
            MatrixX GF = G * F;
            #if DIM == 3
            nodalForcesFEM -= cauchyStress * G.transpose() * dh * dh * dh;
            #else
            nodalForcesFEM -= cauchyStress * G.transpose() * dh * dh;
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

                MatrixN deltaP = material->firstPKDifferential(dF);
                #if DIM == 3
                MatrixX deltaForce = deltaP * GF.transpose() * Jinv * dh * dh * dh;
                #else
                MatrixX deltaForce = deltaP * GF.transpose() * Jinv * dh * dh;
                #endif
                for (int m = 0; m < QUAD; m++) {
                    stiffness(m * DIM, k) += deltaForce(0, m);
                    stiffness(m * DIM + 1, k) += deltaForce(1, m);
                    #if DIM == 3
                    stiffness(m * DIM + 2, k) += deltaForce(2, m);
                    #endif
                }
            }

            for(int i = 0; i < QUAD; i++) {
                int idx = indexMap[index + offsets[i]];
                if(idx >= 0) {
                    localForce.segment<DIM>(idx * DIM) += nodalForcesFEM.col(i);
                }
            }
            for (int i = 0; i < DIM*QUAD; i++) {
                int col = index + offsets[i / DIM];
                if(mass[col] == 0) {
                    continue;
                }
                if(indexMap[col] < 0) {
                    continue;
                }
                int newCol = indexMap[col] * DIM + i % DIM;
                for(int k = 0; k < DIM*QUAD; k++) {
                    int row = index + offsets[k / DIM];
                    if(mass[row] == 0) {
                        continue;
                    }
                    if(indexMap[row] < 0) {
                        continue;
                    }
                    int newRow = indexMap[row] * DIM + k % DIM;
                    localMat.push_back(Tripletd(newRow, newCol, stiffness(k, i)));
                }
            }
        }
    }
    matOffsets[id + 1] = localMat.size();
    #pragma omp barrier
    #pragma omp single
    {
        for(int t = 1; t < (int)matOffsets.size(); t++) {
            matOffsets[t] += matOffsets[t - 1];
        }
        tripletList.resize(matOffsets.back() + startSize);
    }
    std::vector<Tripletd>::iterator dest = tripletList.begin() + startSize;
    std::advance(dest, matOffsets[id]);
    std::copy(localMat.begin(), localMat.end(), dest);
    #pragma omp critical
    {
        force += localForce;
        for(const std::pair<int, MatrixN>& ele : updatedFps)
            FpInvVec[ele.first] = ele.second;
    }
    }
}

//Compute F and return the determinant
double Object::computeF(const MatrixX& coord, const MatrixX& pNpF, MatrixN& F) {
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

void Object::computePFPxhat(const MatrixN& F, const MatrixX& pNpx, MatrixX& pFpxhat) {
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

void Object::computeMV(VectorX mv, const std::vector<int>& indexMap) {
    #pragma omp parallel for schedule(static)
    for (int x = xmin; x < xmax - 1; x++) {
        for (int y = ymin; y < ymax - 1; y++) {
            int index = x*res[1]+y;
            int idx = indexMap[index];
            if(idx >= 0) {
                mv.segment<DIM>(idx * DIM) += mass[index] * velocity[index];
            }
        }
    }
}

void Object::gridToParticles() {
    std::ofstream vb("./euler/v1");
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            VectorN xg = dh*VectorN(i,j) + origin;
            VectorN v = velocity[index];
            vb << xg(0) << " " << xg(1) << " " << v(0) << " " << v(1) << "\n";
        }
    }
    vb.close();

    VectorN wallNormals[4] = {VectorN(1, 0), VectorN(0, 1),
                              VectorN(-1, 0), VectorN(0, -1)};
    double boundaries[4] = {origin(0),
                            origin(1),
                            origin(0) + dh * (res[0] - 1),
                            origin(1) + dh * (res[1] - 1)};
    double alpha = 0.99;
    #pragma omp parallel for schedule(static)
    for(int c = 0; c < (int)matParticles.size(); c++) {
        Particle& p = matParticles[c];
        p.vo = p.v;
        double x = (p.x(0) - origin(0)) / dh;
        double y = (p.x(1) - origin(1)) / dh;
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
        double fx = (futurePos(0) - origin(0)) / dh;
        double fy = (futurePos(1) - origin(1)) / dh;
        int i = std::floor(fx);
        int j = std::floor(fy);
        i = std::min(std::max(i, 0), res[0] - 2);
        j = std::min(std::max(j, 0), res[1] - 2);
        // handle collision with walls
        for(int b = 0; b < 4; b++) {
            double dist = std::abs(p.x(b % DIM) - boundaries[b]);
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

    std::ofstream vm("./euler/v2");
    for(int i = 0; i < (int)matParticles.size(); i++) {
        VectorN xg = matParticles[i].x;
        VectorN v = matParticles[i].v;
        vm << xg(0) << " " << xg(1) << " " << v(0) << " " << v(1) << "\n";
    }
    vm.close();
}

void Object::particlesToGrid() {
    std::unordered_map<int, std::pair<VectorN, double> > indexToWeight;
    #pragma omp parallel
    {
    std::unordered_map<int, std::pair<VectorN, double> > localMap;
    #pragma omp parallel for schedule(static)
    for(int c = 0; c < (int)matParticles.size(); c++) {
        Particle& p = matParticles[c];
        double x = (p.x(0) - origin(0)) / dh;
        double y = (p.x(1) - origin(1)) / dh;
        int i = std::floor(x) - 1;
        int j = std::floor(y) - 1;

        for(int s = i; s < i+4; s++) {
            for(int t = j; t < j+4; t++) {
                if(s < 0 || t < 0 || s > res[0]-1 || t > res[1]-1) {
                    continue;
                }
                VectorN diff = dhInv * (p.x - origin);
                diff -= VectorN(s, t);
                VectorN n(0, 0);
                for(int k = 0; k < DIM; k++) {
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
    // mass clear
    mass = std::vector<double>(totalCells, 0);
    int cutoff = 3;
    for(auto iToW : indexToWeight) {
        if(iToW.second.second < cutoff) {
            continue;
        }
        velocity[iToW.first] = iToW.second.first / iToW.second.second;
        mass[iToW.first] = 1;
    }
    std::ofstream va("./euler/v3");
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            VectorN xg = dh*VectorN(i,j) + origin;
            VectorN v = velocity[index];
            va << xg(0) << " " << xg(1) << " " << v(0) << " " << v(1) << "\n";
        }
    }
    va.close();
    std::exit(0);
}

void Object::extendVectorField() {
    std::vector<double> phifield(res[0]*res[1], res[0]+res[1]+2);
    for(int x = 0; x < res[0]; x++) {
        for(int y = 0; y < res[1]; y++) {
            int ind = x*res[1]+y;
            if(valid[ind]) {
                phifield[ind] = -dh/2.0;
            }
        }
    }
    distSweep(phifield, valid, 8, 1e-8);
    velSweep(velocity, phifield, 8, 1e-8);
}

void Object::finalizeSolution() {
    for(int i = 0; i < totalCells; i++) {
        u[i] += velocity[i] * dt;
    }
    advectField(u, vecWorkspace, dt);

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
            VectorN pos(x*dh, y*dh);
            pos += origin - u[index];
            X[index] = pos;
        }
    }
    vecWorkspace = velocity;
}

void Object::advectField(std::vector<VectorN>& field, std::vector<VectorN>& workspace, double dt) {
    for(int x = 0; x < res[0]; x++) {
        for(int y = 0; y < res[1]; y++) {
            int index = x*res[1]+y;
            VectorN node(x, y);
            node = node * dh + origin;
            VectorN newPos = node - velocity[index] * dt;
            workspace[index] = interpolate_value(newPos, field, origin, res, dh);
        }
    }
    field.swap(workspace);
}

void Object::advectSurface() {
    double epsilon = 1e-4;
    #pragma omp parallel for schedule(static)
    for(unsigned int x = 0; x < vertices.size(); x++) {
        Particle& p = vertices[x];
        VectorN& vertex = p.x;
        VectorN v = interpolate_value(vertex, velocity, origin, res, dh);
        vertex += v * dt;
        double residual = 1;
        int count = 0;
        while(residual > epsilon && count++ < 10) {
            VectorN ui = interpolate_value(vertex, u, origin, res, dh);
            VectorN prest = vertex - ui;
            VectorN diff = prest - p.u;
            residual = diff.norm();

            int i = (vertex[0] - origin[0]) / dh;
            if(i < 0) {
                i = 0;
            }
            if(i > res[0] - 2) {
                i = res[1] - 2;
            }
            int j = (vertex[1] - origin[1]) / dh;
            if(j < 0) {
                j = 0;
            }
            if(j > res[1] - 2) {
                j = res[1] - 2;
            }

            int index = i*res[1]+j;
            MatrixX coord(DIM, QUAD);
            for (int i = 0; i < QUAD; i++) {
                coord.col(i) = u[index + offsets[i]];
            }
            MatrixN F;
            computeF(coord, G, F);
            vertex -= F * diff;
        }
        vertex[0] = std::max(origin[0], vertex[0]);
        vertex[1] = std::max(origin[1], vertex[1]);
        #if DIM == 3
        vertex[2] = std::max(origin[2], vertex[2]);
        #endif

        vertex[0] = std::min(origin[0] + dh * (res[0] - 1), vertex[0]);
        vertex[1] = std::min(origin[1] + dh * (res[1] - 1), vertex[1]);
        #if DIM == 3
        vertex[2] = std::min(origin[2] + dh * (res[2] - 1), vertex[2]);
        #endif
    }
}

void Object::distSweep(std::vector<double>& field, const std::vector<char>& initial, int iters, double eps) {
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
                double fh = f*dh;            //Seeing what happens with f=1
                //Eq 2.4
                double ab = a-b;
                double xbar;
                if(std::abs(ab) >= fh) {
                    xbar = std::min(a, b) + fh;
                }
                else {
                    xbar = (a+b+std::sqrt(2*f*f*dh*dh-ab*ab)) / 2;
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

void Object::velSweep(std::vector<VectorN>& vel, const std::vector<double>& field, int iters, double eps) {
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
