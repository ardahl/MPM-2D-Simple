#include "material.hpp"

void Material::step(double dt) {
    //particlesToGrid
    velOld = vel;
    //gridForces
    //updateGridVelocity
    //velocitySolve
    //updateGradient
    gridToParticles();
}

void Material::particlesToGrid() {
    mass.assign(0);
    vel.assign(Vector2d(0,0));
    velWeight.assign(Vector2d(1e-3,1e-3)); // avoid division by zero
    // TODO: for each particle, add its velocity to its nearby grid cells.
    // Use velWeight to determine the weighted average denominator
    for(int i = 0; i < particles.size(); i++) {
		Particle* p = particles[i];
        mass.addInterpolated(p->x, p->m);
		vel.addInterpolated(p->x, p->v * p->m);
		/// velWeight.addInterpolated(p->x, Vector2d(1, 1));
	}
	for(int i = 0; i < vel.u.m; i++) {
		for(int j = 0; j < vel.u.n; j++) {
			/// vel.u(i, j) /= velWeight.u.get(i, j);
            vel.u(i, j) /= mass(i, j);
		}
	}
	for(int i = 0; i < vel.v.m; i++) {
		for(int j = 0; j < vel.v.n; j++) {
			/// vel.v(i, j) /= velWeight.v.get(i, j);
            vel.v(i, j) /= mass(i, j);
		}
	}
}

void Material::gridForces() {
    
}

void Material::updateGridVelocity() {
    
}

void Material::velocitySolve() {
    
}

void Material::updateGradient(double dt) {
    for(int i = 0; i < particles.size(); i++) {
        Particle* p = particles[i];
        Matrix3d delV = Matrix3d::Identity();
        Matrix3d Fp = (Matrix3d::Identity() + dt * delV) * p->gradient;
    }
}

void Material::gridToParticles() {
    // TODO: update particle velocities using vel and velOld grids
    double alpha = 0.95;
    for(int i = 0; i < particles.size(); i++) {
		Particle* p = particles[i];
		Vector2d pic = vel.interpolate(p->x);
		Vector2d flip = p->v + (vel.interpolate(p->x) - velOld.interpolate(p->x));
		p->v = (alpha * flip) + ((1 - alpha) * pic);
	}
}
