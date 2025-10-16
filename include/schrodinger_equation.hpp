#ifndef SCHRODINGER_HPP
#define SCHRODINGER_HPP

#include <iostream>
#include <print>
#include "Eigen/SparseLU"
#include "interferometer.hpp"
#include "helper_functions.hpp"
#include <memory>

using SparseMatrix = Eigen::SparseMatrix<std::complex<float>>;

class SchodingerEquation
{      
    Eigen::SparseLU<SparseMatrix> m_solver{};
    SparseMatrix m_sparse_M{};
    Eigen::VectorXcf m_psi{};
    Eigen::VectorXcf m_psi_backup{};
    Eigen::VectorXcf m_psi_temp{};    
public:
    explicit SchodingerEquation(size_t Nx_, size_t Ny_, const SparseMatrix& sparse_A , SparseMatrix&& sparse_M, Eigen::VectorXcf&& initial_wf);
    size_t Nx{};
    size_t Ny{};
    float get_wf_modulus(size_t k) const;
    void evolve();
    void interact(const Interferometer& double_slit);
    float get_max_amplitude() const;
    void reset();
};

 

#endif