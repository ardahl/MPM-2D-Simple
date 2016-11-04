#ifndef MATERIAL_HPP_
#define MATERIAL_HPP_

#include "grid.hpp"
/// #include "trimesh2/include/TriMesh.h"

class Force;

class Particle {
public:
    Vector2d x, v;      //postition, velocity
    Vector3d color;     
    double m;           //mass
    Matrix2d gradient;  //deformation gradient
    double rho;         //density
    double vol;         //volume
    Particle(Vector2d x, Vector2d v, Vector3d color, double m): x(x), v(v), color(color), m(m), gradient(Matrix2d::Identity()), rho(0), vol(0) {}
};

class Material {
public:
    /// trimesh::TriMesh* m;
    std::vector<Particle*> particles;
    /// std::vector<Force*> forces;
    int m, n;                               //Grid dimensions
    double h;                               //Grid spacing
    double lambda, mu;                      //Lame Constants for stress
    //Structures used for calculations
    Grid<double> mass;
    //No pressure projection, so no staggered grid needed
    Grid<Vector2d> vel, velStar, f;         //previous velocity, new velocity, grid forces
    
    Material(std::string config);
    //Functions
    void init();                            //Do any configurations, also call Compute_Particle_Volumes_And_Densities
    void getMesh();
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
    virtual void addForces(Material *fluid, double dt) = 0;
};

#endif
