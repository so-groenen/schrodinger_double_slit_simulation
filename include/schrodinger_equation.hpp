#ifndef SCHRODINGER_HPP
#define SCHRODINGER_HPP

#include <iostream>
#include <print>
#include "raylib.h"
#include "eigen3/Eigen/SparseLU"
#include "interferometer.hpp"
#include "wavefunction_builder.hpp"
#include "helper_functions.hpp"
#include "crank_nicolson_manager.hpp"


using SparseMatrix = Eigen::SparseMatrix<std::complex<float>>;

class SchodingerEquation
{      
    std::complex<float> m_rx{};
    std::complex<float> m_ry{};
    std::complex<float> m_a0{};
    std::complex<float> m_b0{};
    Eigen::SparseLU<SparseMatrix> m_solver{};
    Vector2 m_initial_pos{};
    Eigen::VectorXcf m_psi_temp{};    
public:
    float Lx{};
    float Ly{};
    size_t Nx{};
    size_t Ny{};
    SparseMatrix m_sparse_M;
    Eigen::VectorXcf m_psi{};
public:
    explicit SchodingerEquation(Vector2 L, Vector2 dr, Vector2 init_pos);
    float get_wf_modulus(size_t k) const;
    void evolve();
    void interact(const Interferometer& double_slit);
    float get_max_amplitude() const;
private:
    void init_wave_function();
    void init_solver();
};

#endif