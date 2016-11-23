#ifndef MATERIAL_HPP_
#define MATERIAL_HPP_

#include "grid.hpp"

class Force;

class Particle {
public:
    Eigen::Vector2d x, v;      //postition, velocity
    Eigen::Vector3d color;     
    double m;           //mass
    Eigen::Matrix2d gradient;  //deformation gradient
    double rho;         //density
    double vol;         //volume
    Particle(Eigen::Vector2d x, Eigen::Vector2d v, Eigen::Vector3d color, double m): x(x), v(v), color(color), m(m), gradient(Eigen::Matrix2d::Identity()), rho(0), vol(0) {}
};

class Material {
public:
    std::vector<Particle*> particles;
    std::vector<Force*> forces;
    Eigen::Vector2d x0;
    int m, n;                               //Grid dimensions
    double h;                               //Grid spacing
    double lambda, mu;                      //Lame Constants for stress
    //Structures used for calculations
    Grid<double> mass;
    //No pressure projection, so no staggered grid needed
    Grid<Eigen::Vector2d> vel, velStar, f;         //previous velocity, new velocity, grid forces
    
    Material(std::string config);
    Eigen::Vector2d getExtForces(double dt, int i, int j);    //External forces for grid cell (i, j)
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

class Force {
public:
    virtual Eigen::Vector2d addForces(Material *mat, double dt, int i, int j) = 0;
};

class Gravity : public Force {
public:
    Eigen::Vector2d g;
    bool enabled;
    Gravity(Eigen::Vector2d g): g(g), enabled(true) {}
    Eigen::Vector2d addForces(Material *mat, double dt, int i, int j);
};

#endif
