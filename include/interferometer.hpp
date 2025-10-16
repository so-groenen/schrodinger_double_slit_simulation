#ifndef INTERFERRO_HPP
#define INTERFERRO_HPP

#include <iostream>
#include "raylib.h"
#include "Eigen/SparseLU"
struct Point
{
    size_t x{};
    size_t y{};
};

struct Interferometer
{
    size_t m_Nx{};
    size_t m_Ny{};
    size_t m_thickness{};
    size_t m_width{};
    size_t m_height{};
    size_t m_top{};
    size_t m_bottom{};
    Point m_mid_point{};
    Point m_start{};
    Point m_slit_pt{};

    Interferometer(size_t Nx, size_t Ny);
    void set_param(size_t thickness, size_t width, size_t height);
    void activate_interaction(Eigen::VectorXcf& psi) const;
    void draw(const RenderTexture2D& tile, Point start, size_t y, size_t x) const;
private:
    static void zero_out(Eigen::VectorXcf& vec, size_t iy, size_t jx, size_t Ny);
};

#endif