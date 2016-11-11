#include "material.hpp"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>

Material::Material(std::string config) {
    std::ifstream scn;
    scn.open(config);
    if(scn.fail()) {
        printf("File could not be opened: %s\n", config.c_str());
        exit(0);
    }
    std::string objFile;
    double pmass, l, w;
    Vector3d c;
    Vector2d object, grav;
    int resx, resy;
    std::string str;
    x0 = Vector2d(0, 0);
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
            /// scn >> objFile;
            /// m = trimesh::TriMesh::read(objFile);
            scn >> object(0) >> object(1) >> l >> w >> resx >> resy;
        }
        else if(str == "dim") {
            scn >> m >> n;
        }
        else if(str == "mass") {
            scn >> pmass;
        }
        else if(str == "color") {
            scn >> c(0) >> c(1) >> c(2);
        }
        /// else if(str == "h") {
            /// scn >> h;
        /// }
        else if(str == "lame") {
            scn >> lambda >> mu;
        }
        else if(str == "g") {
            scn >> grav(0) >> grav(1);
        }
    }
    scn.close();
    h = 1.0/m;      //Grid bounds always to from 0 to 1 so we can set this explicitly
    //Set up particles at each object vertex
    for(int i = 0; i <= resx; i++) {
        for(int j = 0; j <= resy; j++) {
            Vector2d pos = object + Vector2d(l*((double)i/resx), w*((double)j/resy));
            Particle* par = new Particle(pos, Vector2d(0,0), c, pmass);
            particles.push_back(par);
        }
    }
    printf("Number of particles: %d\n", (int)particles.size());
    printf("Dimentions: %dx%d\n", m, n);
    printf("Spacing: %f\n", h);
    printf("Lame Constants: %f, %f\n", lambda, mu);
    printf("Gravity: (%f, %f)\n", grav(0), grav(1));
    Vector2d xGrid = x0 + Vector2d(h/2.0,h/2.0);
    mass = Grid<double>(m, n, xGrid, h);
    vel = velStar = f = Grid<Vector2d>(m, n, xGrid, h);
    Force* g = new Gravity(grav);
    forces.push_back(g);
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
    for(int i = 0; i < particles.size(); i++) {
        Particle* p = particles[i];
        for(int j = mass.lower(p->x(0), 0); j < mass.upper(p->x(0), 0); j++) {
            for(int k = mass.lower(p->x(1), 1); k < mass.upper(p->x(1), 1); k++) {
                p->rho += mass.interpolate(p->x);
            }
        }
        p->rho /= (h*h*h);
        if(std::isnan(p->rho)) {
            printf("Paricle %d rho has NaN: %f\n", i, p->rho);
            exit(0);
        }
        p->vol = p->m / p->rho;
        if(std::isnan(p->vol)) {
            printf("Paricle %d volume has NaN: %f\n", i, p->vol);
            exit(0);
        }
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
    mass.assign(1e-6); //prevent divide by 0
    vel.assign(Vector2d(0,0));
    for(int i = 0; i < particles.size(); i++) {
		Particle* p = particles[i];
        mass.addInterpolated(p->x, p->m);
		vel.addInterpolated(p->x, p->v * p->m);
	}
	for(int i = 0; i < vel.m; i++) {
		for(int j = 0; j < vel.n; j++) {
            vel(i, j) /= mass(i, j);
            if(vel(i, j).hasNaN()) {
                printf("interpolated vel NaN at (%d, %d)\n", i, j);
                std::cout << vel(i, j) << std::endl;
                exit(0);
            }
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
    f.assign(Vector2d(0, 0));
    for(int i = 0; i < particles.size(); i++) {
        Particle* p = particles[i];
        double J = p->gradient.determinant();
        Matrix2d eps = 0.5 * (p->gradient.transpose() * p->gradient - Matrix2d::Identity());
        Matrix2d stress = lambda*eps.trace()*Matrix2d::Identity() + 2.0*mu*eps;
        for(int j = f.lower(p->x(0), 0); j < f.upper(p->x(0), 0); j++) {
            for(int k = f.lower(p->x(1), 1); k < f.upper(p->x(1), 1); k++) {
                f(j, k) -= p->vol * J * stress * f.gradWeight(p->x, j, k);
                if(f(j, k).hasNaN()) {
                    printf("f NaN at (%d, %d)\n", j, k);
                    std::cout << f(j, k) << std::endl;
                    exit(0);
                }
            }
        }
    }
}

/******************************
 * Update_Grid_Velocities
 *      for each grid node i do
 *          v^*_i = (dt * f_i)/m_i + dt * g
 *      end for
 *****************************/
void Material::updateGridVelocities(double dt) {
    velStar.assign(Vector2d(0,0));
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < velStar.m; i++) {
        for(int j = 0; j < velStar.n; j++) {
            velStar(i, j) = dt * (1.0/mass(i, j)) * f(i, j);// + getExtForces(dt, i, j); //dt*g
            if(velStar(i, j).hasNaN()) {
                printf("velStar NaN at (%d, %d)\n", i, j);
                std::cout << velStar(i, j) << std::endl;
                exit(0);
            }
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
    #pragma omp parallel for
    for(int i = 0; i < particles.size(); i++) {
        Particle* p = particles[i];
        Matrix2d gradV = Matrix2d::Identity();
        for(int j = velStar.lower(p->x(0), 0); j < velStar.upper(p->x(0), 0); j++) {
            for(int k = velStar.lower(p->x(1), 1); k < velStar.upper(p->x(1), 1); k++) {
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
                gradV += velStar(j, k) * velStar.gradWeight(p->x, j, k).transpose();
            }
        }
        if(gradV.hasNaN()) {
            printf("gradV has NaN\n");
            std::cout << gradV << std::endl;
            exit(0);
        }
        Matrix2d fp = Matrix2d::Identity() + dt*gradV;
        p->gradient *= fp;
    }
}

/******************************
 * Update_Particle_Velocities
 *      for each particle p do
 *          v_pic = 0
 *          v_flip = 0
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
    #pragma omp parallel for
    for(int i = 0; i < particles.size(); i++) {
		Particle* p = particles[i];
		//Update velocities
        Vector2d pic = velStar.interpolate(p->x);
		Vector2d flip = p->v + (velStar.interpolate(p->x) - vel.interpolate(p->x));
        if(pic.hasNaN()) {
            printf("PIC Vel has NaN\n");
            std::cout << pic << std::endl;
            exit(0);
        }
        if(flip.hasNaN()) {
            printf("FLIP Vel has NaN\n");
            std::cout << flip << std::endl;
            exit(0);
        }
		p->v = (alpha * flip) + ((1 - alpha) * pic);
        if(p->v.hasNaN()) {
            printf("Vel has NaN\n");
            std::cout << p->v << std::endl;
            exit(0);
        }
        
        //Update Positions
        p->x += dt * p->v;
        if(p->x.hasNaN()) {
            printf("Pos has NaN\n");
            std::cout << p->x << std::endl;
            exit(0);
        }
    }
}


Vector2d Material::getExtForces(double dt, int i, int j) {
    Vector2d force(0, 0);
    for(int i = 0; i < forces.size(); i++) {
        force += forces[i]->addForces(this, dt, i, j);
    }
    return force;
}

Vector2d Gravity::addForces(Material *mat, double dt, int i, int j) {
    if(!enabled) {
        return Vector2d(0, 0);
    }
    return dt * g;
}
