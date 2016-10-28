#ifndef MATERIAL_HPP_
#define MATERIAL_HPP_

#include "grid.hpp"

class Force;

class Particle {
public:
    Vector2d x, v;      //postition, velocity
    Vector3d color;     
    double m;           //mass
    Matrix3d gradient;  //deformation gradient
    Particle(Vector2d x, Vector2d v, Vector3d color, double m): x(x), v(v), color(color), m(m) {}
};

class Material {
public:
    int m, n;
    std::vector<Particle*> particles;
    std::vector<Force*> forces;
    StaggeredGrid vel, velOld;
    Material();
};

class Force {
public:
    virtual void addForces(Material *fluid, double dt) = 0;
};

#endif
