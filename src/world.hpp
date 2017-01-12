#ifndef MATERIAL_HPP_
#define MATERIAL_HPP_

#include <Eigen/Dense>
#include "defines.hpp"
#include <vector>
#include "profiler.hpp"

class Force;

class Particle {
public:
    #ifndef NDEBUG
    Eigen::Vector2d interpPos;
    Eigen::Vector2d x0, x1;
    #endif
    //Temp
    Eigen::Matrix2d stress;
    
    Eigen::Vector2d x, v;      //postition, velocity
    Eigen::Vector3d color;
    double m;                  //mass
    Eigen::Matrix2d gradient;  //deformation gradient
    double rho;         //density
    double vol;         //volume
    Particle(Eigen::Vector2d x, Eigen::Vector2d v, Eigen::Vector3d color, double m): 
	  x(x), v(v), color(color), m(m), gradient(Eigen::Matrix2d::Identity()), rho(0), vol(0) {
          #ifndef NDEBUG
          x0 = x;
          #endif
      }
    Particle() {}
};

class World {
public:
    int stepNum;
    double elapsedTime, dt, totalTime;
    std::string filename;
    std::vector<Particle> particles;
    Eigen::Vector2d origin;                 //lower left and upper right positions
    int res[2];                             //Grid dimensions
    double h;                               //Grid spacing
    double lambda, mu;                      //Lame Constants for stress
    double massPropDamp;
    //Structures used for calculations
    double *mass;
    //No pressure projection, so no staggered grid needed
    Eigen::Vector2d *vel, *velStar, *frc;   //previous velocity, new velocity, grid forces

    // external forces
    Eigen::Vector2d gravity;
    double rotation;
    bool rotationEnabled, gravityEnabled;
    Eigen::Vector2d center;

    //Debugging stuff
    #ifndef NDEBUG
    Eigen::Vector2d *worldPos;
    double angle;
    double m;
    #endif

    World(std::string config);
  ~World();
     //Functions
    void init();                            //Do any configurations, also call Compute_Particle_Volumes_And_Densities
    std::vector<Particle> getParticles() { return particles; }
    //Perform a step of length dt
    void step();
    void particleVolumesDensities();        //Compute_Particle_Volumes_And_Densities
    void particlesToGrid();                 //Rasterize_Particle_Data_To_Grid
    void computeGridForces();               //Compute_Grid_Forces
    void updateGridVelocities();            //Update_Grid_Velocities
    void updateGradient();                  //Update_Deformation_Gradient
    void gridToParticles();                 //Update_Particle_Velocities and Update_Particle_Positions

  benlib::Profiler prof;
};

void writeParticles(const char *fname, const std::vector<Particle> &particles);
bool readParticles(const char *fname, std::vector<Particle> &particles);

#endif
