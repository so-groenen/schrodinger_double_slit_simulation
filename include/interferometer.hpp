#ifndef INTERFERRO_HPP
#define INTERFERRO_HPP

#include <iostream>
#include <eigen3/Eigen/SparseLU>


struct Interferometer
{
    size_t N_x{};
    size_t N_y{};
    size_t thickness{};
    size_t width{};
    size_t opening{};
    size_t height{};
    Interferometer(size_t Nx, size_t Ny);
    void set_param(size_t thickness_, size_t opening_, size_t width_, size_t height_);
    void activate_interaction(Eigen::VectorXcf& psi) const;
private:
    static void zero_out(Eigen::VectorXcf& vec, size_t iy, size_t jx, size_t Ny);
};

#endif