#ifndef EULERIAN_HPP_
#define EULERIAN_HPP_

#include "defines.hpp"
#include <vector>
#include <algorithm>
#include <fstream>
#include "profiler.hpp"

class Particle {
public:
    //Temp
    MatrixN B;         //B matrix from APIC paper

    VectorN x, v, xo, vo, u;      //postition, velocity
    Eigen::Vector3d color;
    double m;                  //mass
    double rho;                //density
    double vol;                //volume

    Particle(VectorN x, VectorN v, Eigen::Vector3d color, double m):
	  B(Eigen::Matrix2d::Zero()), x(x), v(v), color(color), m(m), rho(0.0), vol(0.0) {}
    Particle() {}
};

class Corotation {
public:
    Corotation(double lambda, double mu);
    void init(const MatrixN& F);
    void svd(const MatrixN& F, MatrixN& U, MatrixN& Fhat, MatrixN& V);
    MatrixN firstPiolaKirchhoff();
    MatrixN firstPKDifferential(const MatrixN& dF);
    #if DIM == 3
    MatrixN crossMatrix(const VectorN& vec);
    #endif
    Corotation* copy();

    double lambda, mu;
    //From SVD of F
    MatrixN R, S, L, U, Fhat, V;
};

class Object {
public:
    std::vector<VectorN> X;
    std::vector<VectorN> u;
    std::vector<VectorN> vecWorkspace, matWorkspace;
    int offsets[QUAD];

    Corotation* solidMaterial;
    std::vector<Corotation*> materialCopies;
    double lambda, mu;

    int res[DIM];
    double size[DIM];
    double dh, dhInv;
    VectorN origin;
    std::vector<double> restMass;
    std::vector<MatrixN> FpInvVec;
    VectorI tightMax;
    VectorI tightMin;
    int xmin, xmax;
    int ymin, ymax;
    std::vector<int> activeCells;
    std::vector<int> gridToActiveCell;
    SparseMat systemMatrix;

    double subsampleDh;
    double quarterCellSolidMass;
    double dt;
    int totalCells;
    std::vector<char> valid;
    std::vector<double> phi;
    std::vector<double> mass;

    VectorX quarterVolumeFractions;
    std::vector<double> workspace;
    std::vector<VectorN> velocity;
    VectorX mv;
    VectorX unitForce;
    bool isPlastic;
    double plasticYield;

    VectorX velocityAtDt0;
    VectorX velocityAtDt1;
    VectorX gradient;
    VectorX constraintObjective;
    double objective;

    std::vector<double> detFInv;
    MatrixX G;
    std::vector<Particle> matParticles;
    std::vector<Particle> vertices;
    int testVel;
    double scale;
    VectorN center;
    bool verticesInit;

    Object(int resolution[DIM], double dh, VectorN origin, const std::vector<double>& restMass);
    void setTestVel(int v) { testVel = v; }
    void setInitialState();
    void setInitialVelocity();
    void setMaterial(double lambda, double mu);
    void setPlasticYield(double v) {
        isPlastic = true;
        plasticYield = v;
    }
    void setMassDensity(double density) {
        quarterCellSolidMass = density * dh * dh * 0.25;
    }

    void init();
    void computeSolidMass(bool usesdf=false);
    void computeSDF();
    void initMaterialParticles();
    void recoverVelocityField(std::vector<VectorN>& v);
    void finalizeSolution(std::vector<VectorN>& vel, double step);
    void extendVectorField();
    void distSweep(std::vector<double>& field, const std::vector<char>& initial, int iters, double eps);
    void velSweep(std::vector<Eigen::Vector2d>& vel, const std::vector<double>& field, int iters, double eps);
    void finalizeSolution();
    void advectField(std::vector<VectorN>& field, std::vector<VectorN>& workspace, double dt);
    void advectSurface();

    void computeKf(VectorX& force, std::vector<Tripletd>& tripletList, const std::vector<int>& indexMap);
    double computeF(const MatrixX& coord, const MatrixX& pNpF, MatrixN& F);
    void computePFPxhat(const MatrixN& F, const MatrixX& pNpx, MatrixX& pFpxhat);
    void computeMV(VectorX mv, const std::vector<int>& indexMap);

    //Base version from paper
    void gridToParticles();
    void particlesToGrid();

    template <class F>
    F interpolate_value(const VectorN& point, std::vector<F> field, VectorN origin, int res[DIM], double dh);
    MatrixN interpolate_value(const std::vector<MatrixN>& data, const VectorN& point, int& index, VectorN origin, int res[DIM], double dh);
};

class World {
public:
    int stepNum, steps;
    double elapsedTime, dt, totalTime;
    std::string filename;
    VectorN origin;                 //lower left node position
    int res[DIM];                             //Grid node dimensions
    double dh, dhInv;                               //Grid node spacing
    int totalCells;
    int totalSolidCells;
    int totalActiveCells;
    int offsets[QUAD];

    std::vector<Object*> solids;

    std::vector<VectorN> velocity;
    std::vector<VectorN> vecWorkspace;
    VectorX mv, solidMv;
    VectorX unitForce, solidUnitForce;
    double rayleighAlpha, rayleighBeta;
    VectorX dampedSolidMass;
    VectorX nzMassVec, nzMassVecInv;
    std::vector<double> solidMass;
    std::vector<double> nodeSolidVolumeFraction;
    VectorX quarterSolidVolumeFractions;
    std::vector<int> gridToSolution;
    std::vector<int> gridToSolid;
    VectorX velocityAtDt0;
    VectorX velocityAtDt1;
    VectorX pressureCorrection;
    VectorX solidVelocityAtDt0;
    VectorX solidVelocityAtDt1;
    VectorX solidPressureCorrection;
    SparseMat systemMatrix;
    VectorX systemMatDiag;
    SparseMat JT;
    SparseMat JMInvJT;
    VectorX JMInvJTDiag;
    std::vector<std::vector<std::pair<int, double> > > JJTDiagEntries;
    double fluidDensity;
    double fluidViscocity;
    SparseMat divergenceConstraintJacobian;
    VectorX divergenceConstraintsRHS;
    std::vector<double> intensity;
    std::vector<double> workspace;
    double bouyancy;

    //PCR
    VectorX q, s, Ap, tmp, residual, direction;

    // external forces
    VectorN gravity;
    double rotation;
    bool rotationEnabled, gravityEnabled;
    VectorN center;

    #ifndef NDEBUG
    VectorN **matTrans;
    int count, sMax;
    #endif
    int inc;

    int testVel;
    World(std::string config);
    ~World();
    //Functions
    void init();
    //Perform a step of length dt
    void step();
    void computeDivergenceConstraintRHS();
    void pcrImplicit();
    void finalizeSolution();

    void gatherSolidMasses();
    void computeDivergenceConstraintJacobian();
    bool solvePCR(const SparseMat& A, const VectorX& b, VectorX& x, const VectorX& M);
    bool solvePCRDivFree(const VectorX& b, VectorX& x, const VectorX& precond);
    void parallelMatVec(const SparseMat& A, const VectorX& b, VectorX& x);
    void getMatVecMult(const VectorX& input, VectorX& output);
    void advectFluidVelocityField(std::vector<VectorN>& workspace, double dt);
    void advectField(std::vector<double>& field, std::vector<double>& workspace, double dt);

    //Apic version for testing
    void apicg2p();
    void apicp2g();
    void apicAdvect();

    //Helpers
    void sweepAve(Eigen::Vector2d *vel, int iters);

    benlib::Profiler prof;
};

inline double clamp(double x, double low, double high)
{
    return std::min(high, std::max(x, low));
}

void writeParticles(const char *fname, const std::vector<Particle> &particles);
bool readParticles(const char *fname, std::vector<Particle> &particles);

#endif
