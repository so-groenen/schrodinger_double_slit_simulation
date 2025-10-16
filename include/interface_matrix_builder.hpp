#ifndef IMATRIX_BUILDER_HPP
#define IMATRIX_BUILDER_HPP

#include <iostream>
#include "Eigen/SparseLU"

using SparseMatrix = Eigen::SparseMatrix<std::complex<float>>;

class IMatrixBuilder
{
public:
    virtual void set_num_elements(size_t Nx, size_t Ny) = 0;
    virtual void set_diagonal_elements(std::complex<float> a0, std::complex<float> b0) = 0;
    virtual void set_off_diag_elements(std::complex<float> rx, std::complex<float> ry) = 0;
    virtual auto get_sparse_matrices() const -> std::tuple<SparseMatrix,SparseMatrix> = 0;
    virtual ~IMatrixBuilder() = default;
};

 

#endif