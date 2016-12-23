#ifndef MATERIAL_HPP_
#define MATERIAL_HPP_

#include <Eigen/Dense>
#include "defines.hpp"
#include <vector>

class Force;

class Particle {
public:
    #ifndef NDEBUG
    Eigen::Vector2d interpPos;
    #endif
    //Temp
    Eigen::Matrix2d stress;
    
    Eigen::Vector2d x, v;      //postition, velocity
    Eigen::Vector3d color;
    double m;           //mass
    Eigen::Matrix2d gradient;  //deformation gradient
    double rho;         //density
    double vol;         //volume
    Particle(Eigen::Vector2d x, Eigen::Vector2d v, Eigen::Vector3d color, double m): 
	  x(x), v(v), color(color), m(m), gradient(Eigen::Matrix2d::Identity()), rho(0), vol(0) {}
    Particle() {}
};

class World {
public:
    int stepNum;
    double elapsedTime;
    std::string filename;
    std::vector<Particle> particles;
    Eigen::Vector2d origin;                 //lower left and upper right positions
    int res[2];                             //Grid dimensions
    double h;                               //Grid spacing
    double lambda, mu;                      //Lame Constants for stress
    //Structures used for calculations
    double *mass;
    //No pressure projection, so no staggered grid needed
    Eigen::Vector2d *vel, *velStar, *frc;   //previous velocity, new velocity, grid forces

    // external forces
    Eigen::Vector2d gravity;
    double rotation;
    bool rotationEnabled, gravityEnabled;
    Eigen::Vector2d center;


    #ifndef NDEBUG
    Eigen::Vector2d *worldPos;
    #endif

    World(std::string config);
     //Functions
    void init();                            //Do any configurations, also call Compute_Particle_Volumes_And_Densities
    void getMesh();
    //Perform a step of length dt
    void step(double dt);
    void particleVolumesDensities();        //Compute_Particle_Volumes_And_Densities
    void particlesToGrid();                 //Rasterize_Particle_Data_To_Grid
    void computeGridForces();               //Compute_Grid_Forces
    void updateGridVelocities(double dt);   //Update_Grid_Velocities
    /// void velocitySolve();
    void updateGradient(double dt);         //Update_Deformation_Gradient
    void gridToParticles(double dt);        //Update_Particle_Velocities and Update_Particle_Positions
};

#endif
