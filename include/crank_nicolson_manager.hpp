#ifndef CN_MANAGER_HPP
#define CN_MANAGER_HPP

#include <iostream>
#include <eigen3/Eigen/SparseLU>

using SparseMatrix = Eigen::SparseMatrix<std::complex<float>>;

class CrankNicolsonManager
{
    std::complex<float> m_a0{};
    std::complex<float> m_b0{};
    std::complex<float> m_rx{};
    std::complex<float> m_ry{};
    size_t N_center{};
    size_t N_x{};
    size_t N_y{};
public:
    CrankNicolsonManager(size_t Nx, size_t Ny);
    void set_diagonal_elements(std::complex<float> a0, std::complex<float> b0); // define
    void set_off_diag_elements(std::complex<float> rx, std::complex<float> ry);
    auto get_sparse_matrices() const -> std::tuple<SparseMatrix,SparseMatrix>;
private:
    void init_temp_matrices(Eigen::MatrixXcf& A, Eigen::MatrixXcf& M) const;
    void set_temp_matrices(Eigen::MatrixXcf& A, Eigen::MatrixXcf& M) const;
};

#endif