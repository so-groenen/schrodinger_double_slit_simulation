#ifndef INTERFERRO_HPP
#define INTERFERRO_HPP

#include<iostream>
#include <eigen3/Eigen/SparseLU>



inline void zero_out(Eigen::VectorXcf& vec, size_t iy, size_t jx, size_t Ny)
{
    size_t k = (jx-1)*(Ny-2) + (iy-1);
    vec(k) = 0;
}


struct Interferometer
{
    size_t N_x{};
    size_t N_y{};
    size_t thickness{};
    size_t width{};
    size_t opening{};
    size_t height{};
    Interferometer(size_t Nx, size_t Ny)
        :   N_x(Nx), N_y(Ny)
    {
    }
    void set_param(size_t thickness_, size_t opening_, size_t width_, size_t height_)
    {
        thickness  = thickness_;
        width      = width_;
        opening    = opening_;
        height     = height_;
    }
    void activate_interaction(Eigen::VectorXcf& psi) const
    {
        size_t midX   = N_x/2;
        size_t midY   = N_y/2;
        size_t top    = 0; 
        size_t bottom = (N_y-1);
        for (size_t dx = 1; dx < (thickness+1); dx++)
        {
            for (size_t dy = 1; dy < (height+1); dy++)
            {        
                zero_out(psi, top    + (dy-1), midX - thickness/2  + (dx-1), N_y);     // from "top    of screen" to top    slit opening
                zero_out(psi, bottom - (dy-1), midX - thickness/2  + (dx-1), N_y);     // from "bottom of screen" to bottom slit opening
            }    
            for (size_t dy = 1; dy < (width+1); dy++)
            {
                zero_out(psi, midY - (width/2) + (dy-1), midX - thickness/2  + (dx-1), N_y);
            }
        }
    }
};

#endif