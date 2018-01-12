#ifndef EULERIAN_HPP_
#define EULERIAN_HPP_

#include "defines.hpp"
#include <vector>
#include <algorithm>
#include <fstream>
#include "profiler.hpp"
#define MAT_EXTRAP
class Particle;

class Particle {
public:
    //Temp
    Eigen::Matrix2d B;         //B matrix from APIC paper

    Eigen::Vector2d u;        //Initial position (rest state)
    Eigen::Vector2d x, v;      //postition, velocity
    #ifndef NDEBUG
    std::vector<Eigen::Vector2d> hist;
    #endif
    Eigen::Vector3d color;
    Eigen::Vector3d c1, c2;
    Eigen::Matrix2d gradientE; //elastic portion of deformation gradient
    Eigen::Matrix2d gradientP; //plastic portion of deformation gradient
    double m;                  //mass
    double rho;                //density
    double vol;                //volume

    Particle(Eigen::Vector2d x, Eigen::Vector2d v, Eigen::Vector3d color, double m):
	  B(Eigen::Matrix2d::Zero()), u(x), x(x), v(v), color(color), gradientE(Eigen::Matrix2d::Identity()), gradientP(Eigen::Matrix2d::Identity()), m(m), rho(0.0), vol(0.0) {}
    Particle() {}
};

struct MaterialProps {
    double lambda, mu;                  //Lame Constants for stress
    double compression;                 //critical compression (sec. 5 of stomahkin)
    double stretch;                     //critical stretch (sec. 5 of stomahkin)
    double massPropDamp, mass, pmass, alpha;
};

struct Object {
    MaterialProps mp;
    Eigen::Vector3d color;
    std::string type;
    double size[2];
    int ores[2];
    Eigen::Vector2d object, center;
    std::vector<Particle> particles;
};


class World {
public:
    int stepNum, steps;
    double elapsedTime, dt, totalTime;
    std::string filename;
    Eigen::Vector2d origin;                 //lower left node position
    int res[DIM];                             //Grid node dimensions
    double h;                               //Grid node spacing
    VectorI tightMax;
    VectorI tightMin;
    int xmin, xmax;
    int ymin, ymax;
    std::vector<Object> objects;

    std::vector<double> restMass;
    int offsets[QUAD];
    MatrixX G;
    //Structures used for calculations
    int totalCells;
    int totalSolidCells;
    int totalActiveCells;
    std::vector<VectorN> X;                 //material coordinates
    std::vector<VectorN> u;                 //displacements of each coordinate to their current location
    std::vector<double> mass;
    std::vector<VectorN> velocity;
    std::vector<char> valid;
    std::vector<VectorN> vecWorkspace;
    std::vector<double> phi;
    std::vector<double> detFInv;
    double subsampleDh;
    double quarterCellSolidMass;
    VectorX vStar;                           //updated velocity from forces solve
    VectorX mv;
    VectorX solidMv;
    VectorX unitForce;
    VectorX solidUnitForce;
    VectorX nzMassVec;
    VectorX nzMassVecInv;
    VectorX dampedMass;
    std::vector<int> gridToSolution;
    std::vector<int> gridToSolid;
    std::vector<ETriplet> tripletList;
    double mu;
    double lambda;
    double rayleighAlpha;
    double rayleighBeta;
    VectorX quarterVolumeFractions;
    std::vector<double> nodeSolidVolumeFraction;
    //From SVD of F
    //These are all used at the same time and nowhere else, no need for an array
    MatrixN R;
    MatrixN S;
    MatrixN L;
    MatrixN U;
    MatrixN Fhat;
    MatrixN V;
    SparseMat systemMat;
    VectorX systemMatDiag;

    //Color for Rendering
    std::vector<Eigen::Vector3d> color;

    // external forces
    VectorN gravity;
    double rotation;
    bool rotationEnabled, gravityEnabled, plasticEnabled;
    VectorN center;

    #ifndef NDEBUG
    VectorN **matTrans;
    int count, sMax;
    std::vector<Particle> particles;
    #endif
    int inc;

    World(std::string config);
    ~World();
    //Functions
    void init();                            //Do any configurations
    //Perform a step of length dt
    void step();
    void computeKf();
    void computeMV();
    void pcrImplicit();
    void finalizeSolution();

    double computeF(const MatrixX& coord, const MatrixX& pNpF, MatrixN& F);
    void materialInit(const MatrixN& F);
    void svd(const MatrixN& F, MatrixN& U, MatrixN& Fhat, MatrixN& V);
    MatrixN firstPiolaKirchhoff();
    MatrixN firstPKDifferential(const MatrixN& dF);
    void computePFPxhat(const MatrixN& F, const MatrixX& pNpx, MatrixX& pFpxhat);
    bool solvePCR(const SparseMat& A, const VectorX& b, VectorX& x, const VectorX& M);
    void computeSDF();
    void advectField(std::vector<VectorN>& field, std::vector<VectorN>& workspace);
    void computeSolidMass(bool usesdf=false);
    void gatherSolidMasses();
    MatrixN crossMatrix(const VectorN& vec);

    //Helpers
    void distSweep(std::vector<double>& field, const std::vector<char>& initial, int iters, double eps);
    void velSweep(Eigen::Vector2d *vel, double *field, int iters, double eps);
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
