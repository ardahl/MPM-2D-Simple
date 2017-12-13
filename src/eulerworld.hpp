#ifndef EULERIAN_HPP_
#define EULERIAN_HPP_

#include <Eigen/Dense>
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
    int res[2];                             //Grid node dimensions
    double h;                               //Grid node spacing
    //Structures used for calculations
    double *mass, *phi;
    //dXx and dXy are staggered by 1. Stress is staggered inside by 1
    Eigen::Vector2d *vel, *tau, *mat, *dXx, *dXy;   //grid velocity, traction forces, material coordinates, derivatives of X
    Eigen::Matrix2d *stress, *F;                    //Stress and Deformation Gradient
    //Meta info for mass
    int numNodesMass;               //Number of nodes with non 0 mass
    Eigen::VectorXd interpMass;
    Eigen::VectorXd Fstar;          //Force on right side of system to solve
    std::vector<char> valid, cvalid;
    //Color for Rendering
    Eigen::Vector3d *color;
    Eigen::VectorXi interpColor;
    
    // particle object
    std::vector<Object> objects;
    
    // external forces
    Eigen::Vector2d gravity;
    double rotation;
    bool rotationEnabled, gravityEnabled, plasticEnabled;
    Eigen::Vector2d center;

    World(std::string config);
    ~World();
    //Functions
    void init();                            //Do any configurations
    //Perform a step of length dt
    void step();
    
    void getMass();                         //Form diagonal mass matrix from interpolation
    void getDeformationGradient();          //Calculate deformation gradient for each cell center
    void computeStress();                   //Compute stress from deformation gradients
    void computeForce();                    //Compute tau and assemble force vector
    void solveVelocity();                   //Solve Mv=f* for new grid node velocities
    void velExtrapolate();                  //Extrapolate new velocity to inactive cells
    void advect();                          //Advect material coordinates
	#ifdef MAT_EXTRAP
	void extrapolate_mat();
	#endif

    //Helpers
    Eigen::Matrix2d upwindJac(Eigen::Vector2d* field, int i, int j, bool boundary=true);       //Forms jacobian using upwinding scheme on the velocity field at position (i,j)
    void distSweep(double *field, std::vector<char> initial, int iters, double eps);
    void velSweep(Eigen::Vector2d *vel, double *field, int iters, double eps);
    void sweepAve(Eigen::Vector2d *vel, int iters);
    
    #ifndef NDEBUG
    Eigen::Matrix2d *dD;
    Eigen::Vector2d **matTrans;
    int count, sMax;
    Eigen::VectorXd vdiff, vdiff2;
    std::vector<Particle> particles;
    Eigen::Vector2d *forces, *advFrc, *elastFrc, *extFrc;
    #endif
    int inc;
    Eigen::Vector2d *tmpVel;

    benlib::Profiler prof;
};

inline double clamp(double x, double low, double high)
{
    return std::min(high, std::max(x, low));
}

void writeParticles(const char *fname, const std::vector<Particle> &particles);
bool readParticles(const char *fname, std::vector<Particle> &particles);

#endif
