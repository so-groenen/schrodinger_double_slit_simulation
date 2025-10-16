#ifndef EQUATION_BUILDER_HPP
#define EQUATION_BUILDER_HPP

#include <iostream>
#include <print>
#include "raylib.h"
// #include "eigen3/Eigen/SparseLU"
#include "Eigen/SparseLU"
#include "interferometer.hpp"
#include "interface_wavefunction_builder.hpp"
#include "schrodinger_equation.hpp"
#include "helper_functions.hpp"
#include "interface_matrix_builder.hpp"

#include <memory>
 

class SchodingerEquationBuilder
{      
    std::complex<float> m_rx{};
    std::complex<float> m_ry{};
    std::complex<float> m_a0{};
    std::complex<float> m_b0{};
    Vector2 m_initial_pos{};
    std::unique_ptr<IMatrixBuilder> m_sparse_mat_buidler{};
    std::unique_ptr<IWaveFunctionBuilder> m_wf_builder{};
    SparseMatrix m_sparse_A{};
    SparseMatrix m_sparse_M{};
    Eigen::VectorXcf m_psi{};
    float m_Lx{};
    float m_Ly{};
    size_t m_Nx{};
    size_t m_Ny{};
public:
    explicit SchodingerEquationBuilder(const Vector2& L, const Vector2& dr, const Vector2& init_pos,
                                        std::unique_ptr<IMatrixBuilder> matrix_builder,
                                        std::unique_ptr<IWaveFunctionBuilder> wf_builder);
    auto build_equation() -> SchodingerEquation;
private:
    void init_wave_function();
    void init_sparse_matrices();
};

#endif