#include "material.hpp"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>

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
    Vector2d x0;
    int resx, resy;
    std::string str;
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
            scn >> x0(0) >> x0(1) >> l >> w >> resx >> resy;
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
        else if(str == "h") {
            scn >> h;
        }
        else if(str == "lame") {
            scn >> lambda >> mu;
        }
    }
    scn.close();
    
    //Set up particles at each object vertex
    for(int i = 0; i <= resx; i++) {
        for(int j = 0; j <= resy; j++) {
            Vector2d pos = x0 + Vector2d(l*(i/resx), w*(j/resy));
            Particle* par = new Particle(pos, Vector2d(0,0), c, pmass);
            particles.push_back(par);
        }
    }
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
        for(int j = mass.lower(p->x(0)); j < mass.upper(p->x(0)); j++) {
            for(int k = mass.lower(p->x(1)); k < mass.upper(p->x(1)); k++) {
                p->rho += mass.interpolate(p->x);
            }
        }
        p->rho /= (h*h*h);
        p->vol = p->m / p->rho;
    }
}

/******************************
 * Algorithm Process from Stomakhin et al. 2013 and Yue et al. 2015
 * 
 * Compute_Particle_Volumes_And_Densities
 * while not ended
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
    mass.assign(0);
    vel.assign(Vector2d(0,0));
    for(int i = 0; i < particles.size(); i++) {
		Particle* p = particles[i];
        mass.addInterpolated(p->x, p->m);
		vel.addInterpolated(p->x, p->v * p->m);
	}
	for(int i = 0; i < vel.m; i++) {
		for(int j = 0; j < vel.n; j++) {
            vel(i, j) /= mass(i, j);
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
        for(int j = f.lower(p->x(0)); j < f.upper(p->x(0)); j++) {
            for(int k = f.lower(p->x(1)); k < f.upper(p->x(1)); k++) {
                f(i, j) -= p->vol * J * stress * f.gradWeight(p->x, i, j);
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
    for(int i = 0; i < vel.m; i++) {
        for(int j = 0; j < vel.n; j++) {
            velStar(i, j) = dt * (1.0/mass(i, j)) * f(i, j);    ///+ dt * g;
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
    for(int i = 0; i < particles.size(); i++) {
        Particle* p = particles[i];
        Matrix2d gradV = Matrix2d::Identity();
        for(int j = velStar.lower(p->x(0)); j < velStar.upper(p->x(0)); j++) {
            for(int k = velStar.lower(p->x(1)); k < velStar.upper(p->x(1)); k++) {
                gradV += velStar(i, j) * velStar.gradWeight(p->x, i, j).transpose();
            }
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
    for(int i = 0; i < particles.size(); i++) {
		Particle* p = particles[i];
		//Update velocities
        Vector2d pic = velStar.interpolate(p->x);
		Vector2d flip = p->v + (velStar.interpolate(p->x) - vel.interpolate(p->x));
		p->v = (alpha * flip) + ((1 - alpha) * pic);
        
        //Update Positions
        p->x += dt * p->v;
    }
}
