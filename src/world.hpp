#ifndef WORLD_HPP_
#define WORLD_HPP_

#include <Eigen/Dense>
#include "defines.hpp"
#include <vector>
#include <algorithm>
#include <fstream>
#include "profiler.hpp"

class Particle;

class Spring {
public:
    Particle *p0, *p1;
    double r;           //rest length
    double ks, kd;      //spring const, damping const
};

class Particle {
public:
    //Temp
    MatrixN B;         //B matrix from APIC paper

    VectorN u;        //Initial position (rest state)
    VectorN x, v;      //postition, velocity
    Eigen::Vector3d color;
    Eigen::Vector3d c1, c2;
    MatrixN gradientE; //elastic portion of deformation gradient
    MatrixN gradientP; //plastic portion of deformation gradient
    VectorN f;
    double m;                  //mass
    double rho;                //density
    double vol;                //volume

    Particle(VectorN x, VectorN v, Eigen::Vector3d color, double m):
	  B(MatrixN::Zero()), u(x), x(x), v(v), color(color), gradientE(MatrixN::Identity()), gradientP(MatrixN::Identity()), f(VectorN::Zero()), m(m), rho(0.0), vol(0.0) {}
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
  VectorN object, center;
  std::vector<Particle> particles;      //B*Dinv is an approximation of gradV, so maybe store it for later
  std::vector<Spring> springs;
};


class World {
public:
    int stepNum, steps;
    double elapsedTime, dt, totalTime;
    std::string filename;
    VectorN origin;                 //lower left position
    int res[2];                             //Grid dimensions
    double h;                               //Grid spacing
    //Structures used for calculations
    std::vector<double> mass;
    // double *mass;
    //No pressure projection, so no staggered grid needed
    std::vector<VectorN> vel, prevVel, velStar, frc;
    std::vector<VectorN> mat, matdiff;
    std::vector<MatrixN> stress;
    std::vector<double> weights;
    MatrixX G;
    int offsets[4];
    // VectorN *vel, *velStar, *frc;   //previous velocity, new velocity, grid forces
    // VectorN *prevVel;
    //Material coordinates
    // VectorN *mat;
    // MatrixN *stress;
    // double *weights;

    // particle object
    std::vector<Object> objects;

    // external forces
    VectorN gravity;
    double rotation;
    bool rotationEnabled, gravityEnabled, plasticEnabled;
    VectorN center;

    World(std::string config);
    ~World();
    //Functions
    void init();                            //Do any configurations, also call Compute_Particle_Volumes_And_Densities
    //Perform a step of length dt
    void step();
    void particleVolumesDensities();        //Compute_Particle_Volumes_And_Densities
    void particlesToGrid();                 //Rasterize_Particle_Data_To_Grid
    void computeGridForces();               //Compute_Grid_Forces
    void updateGridVelocities();            //Update_Grid_Velocities
    void updateGradient();                  //Update_Deformation_Gradient
    void gridToParticles();                 //Update_Particle_Velocities and Update_Particle_Positions

    MatrixN upwindJac(const std::vector<VectorN>& field, const std::vector<VectorN>& velocity, int i, int j, bool boundary=true);

    //Semi-lagrangian advection
    std::vector<char> valid;
    void velExtrapolate();
    void slAdvect();
    //Eulerian Advection
    void eAdvect();
    void fastSweep(std::vector<double>& field, std::vector<char> valid, double h, int res[2], int iters, double eps);
    void velExtrapolateFS(std::vector<VectorN>& vel, const std::vector<double>& field, std::vector<char> &valid, int res[2], int iters, double eps);

    //optimization stuff (not part of the actual algorithm)
    //Mostly for testing purposes right now
    #ifndef NDEBUG
    VectorN **matTrans;
    int count, sMax, inc;
    std::vector<char> val;
    #endif

  benlib::Profiler prof;
};

inline double clamp(double x, double low, double high)
{
    return std::min(high, std::max(x, low));
}

void writeParticles(const char *fname, const std::vector<Particle> &particles);
bool readParticles(const char *fname, std::vector<Particle> &particles);

#endif
