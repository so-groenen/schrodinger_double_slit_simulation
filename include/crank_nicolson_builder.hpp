#ifndef CN_MANAGER_HPP
#define CN_MANAGER_HPP

#include <iostream>
#include "Eigen/SparseLU"
#include "interface_matrix_builder.hpp"
using SparseMatrix = Eigen::SparseMatrix<std::complex<float>>;

 
class CrankNicolsonBuilder : public IMatrixBuilder
{
    std::complex<float> m_a0{};
    std::complex<float> m_b0{};
    std::complex<float> m_rx{};
    std::complex<float> m_ry{};
    size_t N_center{};
    size_t m_Nx{};
    size_t m_Ny{};
public:
    explicit CrankNicolsonBuilder();
    void set_num_elements(size_t Nx, size_t Ny) override;
    void set_diagonal_elements(std::complex<float> a0, std::complex<float> b0) override; // define
    void set_off_diag_elements(std::complex<float> rx, std::complex<float> ry) override;
    auto get_sparse_matrices() const -> std::tuple<SparseMatrix,SparseMatrix>  override;
private:
    void init_dense_matrices(Eigen::MatrixXcf& A, Eigen::MatrixXcf& M) const;
    void set_dense_matrices(Eigen::MatrixXcf& A, Eigen::MatrixXcf& M) const;
};

#endif