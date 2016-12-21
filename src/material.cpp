#include "material.hpp"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>

using namespace Eigen;

#ifndef NDEBUG
std::ofstream debug;
int count = 0;
#endif

//TODO: Different colored particles
//TODO: Make circle

//TODO: Check if forces to cells are what you should be getting 
//applying it to particles
//Plot 2norm of (F-R)
//Rotation = (speed*dt)/(mass of body) * pi/2

Material::Material(std::string config) {
    stepNum = 0;
    std::ifstream scn;
    scn.open(config);
    if(scn.fail()) {
        printf("File could not be opened: %s\n", config.c_str());
        exit(0);
    }
    std::string objFile;
    double pmass, l, w, rotspeed;
    Vector3d c;
    Vector2d object, grav;
    int resx, resy;
    std::string str;
    x0 = Vector2d(0, 0);
    Force* frc = NULL;
    while(scn >> str) {
        //Ignore empty lines
        if(str == "") {
            continue;
        }
        //Ignore Comments
        else if(str[0] == '#') {
            scn.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }
        else if(str == "object") {
            scn >> object(0) >> object(1) >> l >> w >> resx >> resy;
        }
        else if(str == "grid") {
            scn >> x0(0) >> x0(1) >> m >> n;
        }
        else if(str == "mass") {
            scn >> pmass;
        }
        else if(str == "color") {
            scn >> c(0) >> c(1) >> c(2);
        }
        else if(str == "h") {
            scn >> h;
        }
        else if(str == "lame") {
            scn >> lambda >> mu;
        }
        else if(str == "g") {
            scn >> grav(0) >> grav(1);
            frc = new Gravity(grav);
            forces.push_back(frc);
        }
        else if(str == "rotate") {
            scn >> rotspeed;
            center = object+Vector2d(l/2.0, w/2.0);
            frc = new Rotate(center, rotspeed);
            forces.push_back(frc);
        }
    }
    scn.close();
    Vector2d xGrid = x0 + Vector2d(h/2.0,h/2.0);
    x0 = xGrid;
    x1 = x0 + h*Vector2d(m, n);
    //Set up particles at each object vertex
    double diffx = l / (resx-1);
    double diffy = w / (resy-1);
    for(int i = 0; i < resx; i++) {
        for(int j = 0; j < resy; j++) {
            Vector2d pos = object + Vector2d(diffx*i, diffy*j);
            Particle* par = new Particle(pos, Vector2d(0,0), c, pmass);
            particles.push_back(par);
        }
    }
    //Average position for center of mass
    Vector2d avePos = Vector2d::Zero();
    for(size_t i = 0; i < particles.size(); i++) {
        avePos += particles[i]->x;
    }
    avePos /= particles.size();

    printf("Number of particles: %d\n", (int)particles.size());
    printf("Dimentions: %dx%d\n", m, n);
    printf("Grid Spacing: %f\n", h);
    printf("Lame Constants: %f, %f\n", lambda, mu);
    printf("Gravity: (%f, %f)\n", grav(0), grav(1));
    printf("Object Spacing: (%f, %f)\n", diffx, diffy);
    printf("Center of Mass: (%f, %f), (%f, %f)\n", center(0), center(1), avePos(0), avePos(1));
    printf("X0: (%f, %f)\n", x0(0), x0(1));
    printf("X1: (%f, %f)\n", x1(0), x1(1));
    mass = Grid<double>(m, n, xGrid, h);
    vel = velStar = f = Grid<Vector2d>(m, n, xGrid, h);
    
    //Test interpolation
    /// #ifndef NDEBUG
    /// //Test grid interpolation
    /// //Set grid values to their world position
    /// worldPos = Grid<Vector2d>(m, n, xGrid, h);
    /// for(int i = 0; i < m; i++) {
        /// for(int j = 0; j < n; j++) {
            /// worldPos(i, j) = worldPos.x0 + Vector2d(i*h, j*h);
            /// debug << "(" << worldPos(i, j)(0) << ", " << worldPos(i, j)(1) << ")  ";
        /// }
        /// debug << "\n";
    /// }
    /// //Try and recover particle positions
    /// debug << "Interpolated Particle Positions\n";
    /// for(int i = 0; i < particles.size(); i++) {
        /// Particle* p = particles[i];
        /// p->interpPos = Vector2d(0.0, 0.0);
        /// // debug << "\nWeights:\n";
        /// for(int j = worldPos.lower(p->x, 0); j < worldPos.upper(p->x, 0); j++) {
            /// for(int k = worldPos.lower(p->x, 1); k < worldPos.upper(p->x, 1); k++) {
                /// p->interpPos += worldPos.weight(p->x, j, k) * worldPos(j, k);
            /// }
        /// }
        /// Vector2d difference = p->x - p->interpPos;
        /// if(!difference.isZero(1e8)) {
            /// debug << "Particle " << i << ": (" << p->x(0) << ", " << p->x(1) << ") -> (" << p->interpPos(0) << ", " << p->interpPos(1) << ")\n";
        /// }
    /// }
    /// #endif
    
    //TODO: Check gradient
}

void Material::init() {
    particleVolumesDensities();
}

/******************************
 * Compute_Particle_Volumes_And_Densities
 *      Raserise_Particle_Data_To_Grid
 *      for each particle p do
 *          rho_p = 0
 *          for each grid node i s.t. w_ip > 0
 *              rho_p += w_ip * m_i
 *          end for
 *          rho_p /= h^3
 *          V_p = m_p / rho_p
 *      end for
 *****************************/
void Material::particleVolumesDensities() {
    particlesToGrid();
    for(size_t i = 0; i < particles.size(); i++) {
        Particle* p = particles[i];
        for(int j = mass.lower(p->x, 0); j < mass.upper(p->x, 0); j++) {
            for(int k = mass.lower(p->x, 1); k < mass.upper(p->x, 1); k++) {
                double r = mass(j, k) * mass.weight(p->x, j, k);
                p->rho += r;
            }
        }
        p->rho /= (h*h);            //3D: h*h*h
        #ifndef NDEBUG
        if(std::isnan(p->rho) || std::isinf(p->rho)) {
            printf("Paricle %d rho has NaN: %f\n", (int)i, p->rho);
            exit(0);
        }
        #endif
        p->vol = p->m / p->rho;
        #ifndef NDEBUG
        if(std::isnan(p->vol) || std::isinf(p->vol)) {
            printf("Paricle %d volume has NaN: %f\nMass: %f\nDensity: %f\n", (int)i, p->vol, p->m, p->rho);
            exit(0);
        }
        #endif
    }
}

/******************************
 * Algorithm Process from Stomakhin et al. 2013 and Yue et al. 2015
 *
 * Compute_Particle_Volumes_And_Densities
 * while true
 *      Rasterize_Particle_Data_to_Grid
 *      Compute_Grid_Forces
 *      Update_Grid_Velocities
 *      Update_Deformation_Gradient
 *      Update_Particle_Velocities
 *      Update_Particle_Positions
 ******************************/
void Material::step(double dt) {
    particlesToGrid();
    computeGridForces();
    updateGridVelocities(dt);
    updateGradient(dt);
    gridToParticles(dt);
    stepNum++;
}

/******************************
 * Rasterize_Particle_Data_To_Grid
 *      for each grid node i do
 *          m_i = 0
 *          v_i = 0
 *      end for
 *      for each particle p do
 *          for each grid node i s.t. w_ip > 0 do
 *              m_i += w_ip * m_p
 *              v_i += w_ip * m_p * v_p
 *          end for
 *      end for
 *      for each grid cell i do
 *          v_i /= m_i
 *      end for
 *
 *****************************/
void Material::particlesToGrid() {
    mass.assign(0.0);
    vel.assign(Vector2d(0,0));
    for(size_t i = 0; i < particles.size(); i++) {
		Particle* p = particles[i];
        for(int j = mass.lower(p->x, 0); j < mass.upper(p->x, 0); j++) {
            for(int k = mass.lower(p->x, 1); k < mass.upper(p->x, 1); k++) {
                mass(j, k) += mass.weight(p->x, j, k) * p->m;
                vel(j, k) += vel.weight(p->x, j, k) * p->m * p->v;
            }
        }
	}
    /// #pragma omp parallel for collapse(2)
	for(int i = 0; i < vel.m; i++) {
		for(int j = 0; j < vel.n; j++) {
            if(std::abs(mass(i, j)) < EPS) {
                vel(i, j) = Vector2d(0.0, 0.0);
            }
            else {
                vel(i, j) /= mass(i, j);
            }
            #ifndef NDEBUG
            if(vel(i, j).hasNaN()) {
                printf("interpolated vel NaN at (%d, %d)\n", i, j);
                std::cout << vel(i, j) << std::endl;
                exit(0);
            }
            #endif
		}
	}
}

/******************************
 * Compute_Grid_Forces
 *      for each grid node i do
 *          f_i = 0
 *      end for
 *      for each particle p do
 *          J_p = det(F_p)
 *          epsilon_p = (1/2)*(F^T * F - I)
 *          stress_p = lambda * trace(epsilon_p) * I + 2 * mu * epsilon_p
 *          for each grid cell i s.t. w_ip > 0
 *              f_i -= V_p * J_p * stress_p * grad(w_ip)
 *          end for
 *      end for
 *****************************/
void Material::computeGridForces() {
    // TODO: Check Force direction
    f.assign(Vector2d(0, 0));
    for(size_t i = 0; i < particles.size(); i++) {
        Particle* p = particles[i];
        double J = p->gradient.determinant();
        Matrix2d eps = 0.5 * (p->gradient.transpose() * p->gradient - Matrix2d::Identity());
        double trace = eps.trace();
        Matrix2d stress = lambda*trace*Matrix2d::Identity() + 2.0*mu*eps;
        for(int j = f.lower(p->x, 0); j < f.upper(p->x, 0); j++) {
            for(int k = f.lower(p->x, 1); k < f.upper(p->x, 1); k++) {
                Vector2d accumF = p->vol * J * stress * f.gradWeight(p->x, j, k);
                f(j, k) -= accumF;
                #ifndef NDEBUG
                if(f(j, k).hasNaN()) {
                    printf("\nf NaN at (%d, %d)\n", j, k);
                    std::cout << "Force:\n" << f(j, k) << std::endl;
                    std::cout << "Volume: " << p->vol << std::endl;
                    std::cout << "Determinant: " << p->gradient.determinant() << std::endl;
                    std::cout << "Stress:\n" << p->stress << std::endl;
                    std::cout << "Gradient:\n" << f.gradWeight(p->x, j, k) << std::endl;
                    exit(0);
                }
                #endif
            }
        }
    }
}

/******************************
 * Update_Grid_Velocities
 *      for each grid node i do
 *          v^*_i = v_i^n + (dt * f_i)/m_i + dt * g
 *      end for
 *****************************/
void Material::updateGridVelocities(double dt) {
    /// #pragma omp parallel for collapse(2)
    for(int i = 0; i < velStar.m; i++) {
        for(int j = 0; j < velStar.n; j++) {
            if(std::abs(mass(i, j)) < EPS) {
                velStar(i, j) = Vector2d(0, 0);
            }
            else {
                velStar(i, j) = vel(i, j) + dt * (1.0/mass(i, j)) * f(i, j) + getExtForces(dt, i, j); //dt*g
            }
            #ifndef NDEBUG
            if(velStar(i, j).hasNaN()) {
                printf("velStar NaN at (%d, %d)\n", i, j);
                std::cout << velStar(i, j) << std::endl;
                exit(0);
            }
            #endif
        }
    }
}

/******************************
 * Update_Deformation_Gradient
 *      for each particle p do
 *          grad(v_p) = 0
 *          for each grid node i s.t. w_ip > 0
 *              grad(v_p) += v^*_i * (grad(w_ip))^T
 *          end for
 *          f_p,n+1 = I + dt*grad(v_p)
 *          F_p,n+1 *= f_p,n+1
 *      end for
 *****************************/
void Material::updateGradient(double dt) {
    /// #pragma omp parallel for 
    for(size_t i = 0; i < particles.size(); i++) {
        Particle* p = particles[i];
        Matrix2d gradV = Matrix2d::Zero();
        for(int j = velStar.lower(p->x, 0); j < velStar.upper(p->x, 0); j++) {
            for(int k = velStar.lower(p->x, 1); k < velStar.upper(p->x, 1); k++) {
                #ifndef NDEBUG
                if(velStar(j, k).hasNaN()) {
                    printf("gradV velStar has NaN at (%d, %d)\n", j, k);
                    std::cout << velStar(j, k) << std::endl;
                    exit(0);
                }
                if(velStar.gradWeight(p->x, j, k).transpose().hasNaN()) {
                    printf("gradV gradW has NaN at (%d, %d)\n", j, k);
                    std::cout << velStar.gradWeight(p->x, j, k).transpose() << std::endl;
                    exit(0);
                }
                #endif
                Matrix2d accumGrad = velStar(j, k) * velStar.gradWeight(p->x, j, k).transpose();
                gradV += accumGrad;
            }
        }
        #ifndef NDEBUG
        if(gradV.hasNaN()) {
            printf("gradV has NaN\n");
            std::cout << gradV << std::endl;
            exit(0);
        }
        #endif
        Matrix2d fp = Matrix2d::Identity() + dt*gradV;
        p->gradient = fp*p->gradient;
    }
}

/******************************
 * Update_Particle_Velocities
 *      for each particle p do
 *          v_pic = 0
 *          v_flip = v_p
 *          for all grid cells i s.t. w_ip > 0
 *              v_pic += w_ip * v^*_i
 *              v_flip += w_ip * (v^*_i - v_i)
 *          end for
 *          v_p = (1-alpha)*v_pic + alpha*v_flip
 *      end for
 *
 * Update_Particle_Positions
 *      for each particle p do
 *          x_p += dt * v_p
 *      end for
 *****************************/
void Material::gridToParticles(double dt) {
    double alpha = 0.95;
    /// #pragma omp parallel for
    for(size_t i = 0; i < particles.size(); i++) {
		Particle* p = particles[i];
		//Update velocities
        Vector2d pic = Vector2d::Zero();
        Vector2d flip = p->v;
        Vector2d tmpPic, tmpFlip;
        for(int j = velStar.lower(p->x, 0); j < velStar.upper(p->x, 0); j++) {
            for(int k = velStar.lower(p->x, 1); k < velStar.upper(p->x, 1); k++) {
                tmpPic = velStar.weight(p->x, j, k) * velStar(j, k);
                pic += tmpPic;

                tmpFlip = velStar.weight(p->x, j, k) * (velStar(j, k) - vel(j, k));
                flip += tmpFlip;
            }
        }
        #ifndef NDEBUG
        if(pic.hasNaN()) {
            printf("\n\nPIC Vel has NaN\n");
            std::cout << pic << std::endl;
            exit(0);
        }
        if(flip.hasNaN()) {
            printf("FLIP Vel has NaN\n");
            std::cout << flip << std::endl;
            exit(0);
        }
        #endif
		p->v = (alpha * flip) + ((1 - alpha) * pic);
        //Mass proportional damping
        p->v *= 0.99999;
        #ifndef NDEBUG
        if(p->v.hasNaN()) {
            printf("Vel has NaN\n");
            std::cout << p->v << std::endl;
            exit(0);
        }
        #endif

        //Update Positions
        p->x += dt * p->v;
        #ifndef NDEBUG
        if(p->x.hasNaN()) {
            printf("Pos has NaN\n");
            std::cout << p->x << std::endl;
            exit(0);
        }
        #endif
        //Boundary collisions
        //Make sure there's a boundary of 2 cells for each side since the weighting
        //function will touch 2 cells out
        double lx = x0(0)+2*h;
        double ly = x0(1)+2*h;
        double ux = x0(0)+(m-2)*h;
        double uy = x0(1)+(n-2)*h;
        if(p->x(0) < lx) {
            p->x(0) = lx+EPS;
            p->v(0) *= -1;          //reverse velocity
        }
        if(p->x(1) < ly) {
            p->x(1) = ly+EPS;
            p->v(1) *= -1;
        }
        if(p->x(0) > ux) {
            p->x(0) = ux-EPS;
            p->v(0) *= -1;
        }
        if(p->x(1) > uy) {
            p->x(1) = uy-EPS;
            p->v(1) *= -1;
        }
    }
}


Vector2d Material::getExtForces(double dt, int i, int j) {
    Vector2d force(0, 0);
    for(size_t k = 0; k < forces.size(); k++) {
        force += forces[k]->addForces(this, dt, i, j);
    }
    return force;
}

Vector2d Gravity::addForces(Material *mat, double dt, int i, int j) {
    if(!enabled) {
        return Vector2d(0, 0);
    }
    return dt*g;
}

Vector2d Rotate::addForces(Material *mat, double dt, int i, int j) {
    double x = (mat->x0(0)+mat->h*i)-center(0);
    double y = (mat->x0(1)+mat->h*j)-center(1);
    return speed*dt*Vector2d(-y, x);
}