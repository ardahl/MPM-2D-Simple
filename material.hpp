#ifndef MATERIAL_HPP_
#define MATERIAL_HPP_

#include "grid.hpp"

class Force;

class Particle {
public:
    Vector2d x, v;      //postition, velocity
    Vector3d color;     
    double m;           //mass
    Matrix3d gradient, gradE, gradP;  //deformation gradient
    Particle(Vector2d x, Vector2d v, Vector3d color, double m): x(x), v(v), color(color), m(m), gradient(Matrix2d::Identity()), gradE(Matrix2d::Identity()), gradP(Matrix2d::Identity()) {}
};

class Material {
public:
    std::vector<Particle*> particles;
    std::vector<Force*> forces;
    //Structures used for calculations
    StaggeredGrid vel, velOld;
    VectorXd velStar;
    Grid<double> mass;                          //Staggered or regular?
    int m, n;
    
    Material();
    //Functions
    void init();
    void getMesh();
    void step(double dt);
    void particlesToGrid();
    void gridForces();
    void updateGridVelocity(double dt);
    void velocitySolve();
    void updateGradient(double dt);
    void gridToParticles(double dt);
};

class Force {
public:
    virtual void addForces(Material *fluid, double dt) = 0;
};

#endif
