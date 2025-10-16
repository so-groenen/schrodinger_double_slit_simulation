#include "interferometer.hpp"


Interferometer::Interferometer(size_t Nx, size_t Ny)
    : m_Nx{Nx}, m_Ny{Ny}, m_top{0}, m_bottom{Ny-1}, m_mid_point{Nx/2,Ny/2}
{
}
void Interferometer::set_param(size_t thickness_, size_t width_, size_t height_)
{
    m_thickness  = thickness_;
    m_width      = width_;
    m_height     = height_;

    size_t slit_x {m_mid_point.x - m_thickness/2};
    size_t slit_y {m_mid_point.y - m_width/2};
    m_slit_pt =   {slit_x,slit_y};
}
void Interferometer::activate_interaction(Eigen::VectorXcf& psi) const
{
    for (size_t dx = 1; dx < (m_thickness+1); dx++)
    {
        for (size_t dy = 1; dy < (m_height+1); dy++)
        {        
            zero_out(psi, m_top    + (dy-1), m_mid_point.x - m_thickness/2  + (dx-1), m_Ny);     // from "top    of screen" to top    slit opening
            zero_out(psi, m_bottom - (dy-1), m_mid_point.x - m_thickness/2  + (dx-1), m_Ny);     // from "bottom of screen" to bottom slit opening
        }    
        for (size_t dy = 1; dy < (m_width+1); dy++)
        {
            zero_out(psi, m_mid_point.y - (m_width/2) + (dy-1), m_mid_point.x - m_thickness/2  + (dx-1), m_Ny);
        }
    }
}
void Interferometer::draw(const RenderTexture2D& tile, Point start, size_t y, size_t x) const
{
    int dy = tile.texture.height;
    int dx = tile.texture.width;
 
    // Draw the part inside between the two slits
    if( (y < m_width ) && (x < m_thickness) )
    {
        DrawTexture(tile.texture, start.x + (m_slit_pt.x + x)*dx, start.y + (m_slit_pt.y + y)*dy, RAYWHITE);
    }

    // Draw the upper & lower parts
    if( (y < m_height ) && (x < m_thickness))
    {                    
        DrawTexture(tile.texture, start.x + (m_slit_pt.x + x)*dx, start.y +  (m_top    + y)*dy, RAYWHITE);
        DrawTexture(tile.texture, start.x + (m_slit_pt.x + x)*dx, start.y +  (m_bottom - y)*dy, RAYWHITE);
    }
}
void Interferometer::zero_out(Eigen::VectorXcf& vec, size_t iy, size_t jx, size_t Ny)
{
    size_t k = (jx-1)*(Ny-2) + (iy-1);
    vec(k) = 0;
}

