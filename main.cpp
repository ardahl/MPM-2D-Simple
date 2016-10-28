#include "material.hpp"

int main(int argc, char** argv) {
    
    /******************************
     * Algorithm Process
     * 
     * while not ended
     *      Rasterize particle data to the grid
     *      Compute Grid Forces
     *      Update grid velocities
     *      Solve semi-implicit integration system
     *      Update Deformation Gradient
     *      Update particle velocities
     *      Update particle positions
     * 
     * 
     * 
     * Rasterize particle data to the grid
     *      Transfer mass       m^n_i = Ep m_p * w^n_ip
     *      Transfer velocities v^n_i = Ep (v^n_p * m_p * w^n_ip) / m^n_i
     * 
     * Compute grid forces
     *      
     * 
     * Update grid velocities
     *      Update to v^*_i = v^n_i + dt * m^-1_i * f^n_i
     * 
     * Solve semi-implicit integration system
     * 
     * 
     * Update deformation gradient
     * 
     * 
     * Update particle velocities
     *      v^n+1_p = (1-a)*v_pic + a*v_flip
     * 
     * Update particle positions
     *      x^n+1_p = x^n_p + dt * v^n+1_p
     * 
     ******************************/
    
    return 0;
}
