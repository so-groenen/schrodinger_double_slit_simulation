// #include <iostream>
// #include <print>
// #include "raylib.h"
// #include "Eigen/SparseLU"
#include "schrodinger_equation_builder.hpp"
// #include "schrodinger_equation.hpp"
// #include "crank_nicolson_builder.hpp"
// #include "schrodinger_equation_builder.hpp"

constexpr float SIGMA_DEVIATION = 0.2f;
constexpr std::complex<float> imaginary_unit{0, 1}; 

SchodingerEquationBuilder::SchodingerEquationBuilder(const Vector2& L, const Vector2& dr, const Vector2& init_pos,
                                                    std::unique_ptr<IMatrixBuilder> matrix_builder, 
                                                    std::unique_ptr<IWaveFunctionBuilder> wf_builder)

    : m_initial_pos{init_pos}, m_Lx{L.x}, m_Ly{L.y}, m_Nx{get_num_elements(0, L.x, dr.x)}, m_Ny{get_num_elements(0, L.y, dr.y)}
      ,m_sparse_mat_buidler{std::move(matrix_builder)}, m_wf_builder{std::move(wf_builder)}
{
    float dx = dr.x;
    float dy = dr.y;
    float dt = (dx*dx)/4.f;
    m_rx = - dt / ( 2.f*imaginary_unit*(dx*dx));
    m_ry = - dt / ( 2.f*imaginary_unit*(dy*dy));
    m_a0 = (1.0f + 2.0f*m_rx + 2.0f*m_ry);
    m_b0 = (1.0f - 2.0f*m_rx - 2.0f*m_ry);
}

void SchodingerEquationBuilder::init_wave_function()
{
    float x0 = m_initial_pos.x;
    float y0 = m_initial_pos.y;
    try 
    {
        m_wf_builder->set_system_size(m_Lx, m_Ly);
        m_wf_builder->set_initial_pos(x0, y0);
        m_wf_builder->set_deviation(SIGMA_DEVIATION);
        auto psi_temp = m_wf_builder->build_wavefunction(m_Ny, m_Nx);
        m_psi         = std::move(psi_temp);
        std::println("Wavefunction allocated: size: {}byes.",  get_size(m_psi));

    }
    catch(const std::exception& e)
    {
        std::println("Wavefunction allocation failure: {}", e.what());
        exit(EXIT_FAILURE);
    }
}
auto SchodingerEquationBuilder::build_equation() -> SchodingerEquation
{
    init_wave_function();
    init_sparse_matrices();
    return SchodingerEquation(m_Nx, m_Ny, (m_sparse_A), std::move(m_sparse_M), std::move(m_psi) );
}
void SchodingerEquationBuilder::init_sparse_matrices()
{
    try
    {
        m_sparse_mat_buidler->set_num_elements(m_Nx, m_Ny);
        m_sparse_mat_buidler->set_diagonal_elements(m_a0, m_b0);
        m_sparse_mat_buidler->set_off_diag_elements(m_rx, m_ry);
        auto[sparse_A, sparse_M] = m_sparse_mat_buidler->get_sparse_matrices();

        std::println("Sparse matrix allocated: {}bytes.", get_size(sparse_A) );  
        m_sparse_A = std::move(sparse_A);                            
        m_sparse_M = std::move(sparse_M);
    }
    catch(const std::exception& e)
    {
        std::println("Sparse-Matrix allocation failure: {}", e.what());
        exit(EXIT_FAILURE);
    }
}