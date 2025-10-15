#include <iostream>
#include <print>
#include "raylib.h"
#include "eigen3/Eigen/SparseLU"
#include "interferometer.hpp"
#include "wavefunction_builder.hpp"
#include "helper_functions.hpp"
#include "crank_nicolson_manager.hpp"
#include "schrodinger_equation.hpp"


SchodingerEquation::SchodingerEquation(Vector2 L, Vector2 dr, Vector2 init_pos)
    : m_initial_pos{init_pos}, Lx{L.x}, Ly{L.y}, Nx{get_num_elements(0, L.x, dr.x)}, Ny{get_num_elements(0, L.y, dr.y)}
{
    float dx = dr.x;
    float dy = dr.y;
    float dt = (dx*dx)/4.f;
    m_rx = - dt / ( 2.f*imaginary_unit*(dx*dx));
    m_ry = - dt / ( 2.f*imaginary_unit*(dy*dy));
    m_a0 = (1.0f + 2.0f*m_rx + 2.0f*m_ry);
    m_b0 = (1.0f - 2.0f*m_rx - 2.0f*m_ry);

    init_wave_function();
    init_solver();
}
float SchodingerEquation::get_wf_modulus(size_t k) const
{
    std::complex<float> z = m_psi(k);
    return sqrt(   pow(z.real(),2)  + pow(z.imag(),2)   );
}
void SchodingerEquation::evolve()
{
    m_psi_temp.noalias() = m_sparse_M*m_psi;
    m_psi = m_solver.solve(m_psi_temp);
}
void SchodingerEquation::interact(const Interferometer& double_slit)
{
    double_slit.activate_interaction(m_psi);
}
float SchodingerEquation::get_max_amplitude() const
{
    float vmax{-1};
    float value{};
    for (size_t k = 0; k < (size_t)m_psi.size(); k++)
    {
        value = get_wf_modulus(k);
        if (value > vmax)
        {
            vmax = value;
        }
    }
    return vmax;
}
void SchodingerEquation::init_wave_function()
{
    float x0 = m_initial_pos.x;
    float y0 = m_initial_pos.y;
    try 
    {
        WaveFunctionBuilder wf_builder{Lx, Ly};
        wf_builder.set_initial_pos(x0, y0);
        wf_builder.set_deviation(0.2);
        auto psi_temp = wf_builder.build_wavefunction(Ny, Nx);
        m_psi         = std::move(psi_temp);
        std::println("Wavefunction allocated: size: {}byes.",  get_size(m_psi));

    }
    catch(const std::exception& e)
    {
        std::println("Wavefunction allocation failure: {}", e.what());
        exit(EXIT_FAILURE);
    }
}
void SchodingerEquation::init_solver()
{
    try
    {
        CrankNicolsonManager cn_manager{Nx, Ny};
        cn_manager.set_diagonal_elements(m_a0, m_b0);
        cn_manager.set_off_diag_elements(m_rx, m_ry);
        auto[sparse_A, sparse_M_temp] = cn_manager.get_sparse_matrices();

        std::println("Sparse matrix allocated: {}bytes.", get_size(sparse_A) );  
        m_solver.compute(sparse_A);                                   // Schrodinger solver is initilized here.        
        m_sparse_M = std::move(sparse_M_temp);
    }
    catch(const std::exception& e)
    {
        std::println("Sparse-Matrix allocation failure: {}", e.what());
        exit(EXIT_FAILURE);
    }
}
 